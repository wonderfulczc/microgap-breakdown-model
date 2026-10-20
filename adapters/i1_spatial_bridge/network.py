"""Axisymmetric cell-face FV conduction; symmetric series-resistance faces, no floor."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve


def harmonic(a,b):
    a,b=np.broadcast_arrays(np.asarray(a,dtype=float),np.asarray(b,dtype=float))
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or np.any(a<0) or np.any(b<0):
        raise ValueError('INVALID_CONDUCTIVITY')
    result=np.zeros_like(a);positive=(a>0)&(b>0)
    low=np.minimum(a[positive],b[positive]);high=np.maximum(a[positive],b[positive])
    result[positive]=2*low/(1+low/high)
    return result


@dataclass
class Network:
    nodes: np.ndarray
    shape: tuple
    left: np.ndarray
    right: np.ndarray
    conductance: np.ndarray
    HV: np.ndarray
    RETURN: np.ndarray
    labels: np.ndarray
    component_HV: np.ndarray
    component_RETURN: np.ndarray


def build(sigma,mask,dr,dz):
    sigma=np.asarray(sigma,dtype=float);mask=np.asarray(mask)
    if sigma.shape!=mask.shape or sigma.ndim!=2 or dr<=0 or dz<=0:raise ValueError('GRID_MISMATCH')
    if np.any(~np.isfinite(sigma)) or np.any(sigma<0):raise ValueError('INVALID_CONDUCTIVITY')
    if np.any(~np.isin(mask,[0,1,2])):raise ValueError('UNKNOWN_MATERIAL')
    nz,nr=sigma.shape;positive=(mask==0)&(sigma>0)
    nodes=np.flatnonzero(positive.ravel());idx=np.full(mask.shape,-1,dtype=int);idx.ravel()[nodes]=np.arange(len(nodes))
    left=[];right=[];weight=[];hv=np.zeros(len(nodes));ret=hv.copy()
    # Each gas-gas face occurs once; gas-conductor uses gas sigma and half-cell distance.
    for axis in [0,1]:
        a=(slice(None,-1),slice(None)) if axis==0 else (slice(None),slice(None,-1))
        b=(slice(1,None),slice(None)) if axis==0 else (slice(None),slice(1,None))
        ma,mb=mask[a],mask[b];sa,sb=sigma[a],sigma[b];ia,ib=idx[a],idx[b]
        area=np.broadcast_to(np.pi*((np.arange(nr)+1)**2-np.arange(nr)**2)*dr**2,(nz-1,nr)) if axis==0 else np.broadcast_to(2*np.pi*np.arange(1,nr)*dr*dz,(nz,nr-1))
        distance=dz if axis==0 else dr
        both=(ia>=0)&(ib>=0)
        left.extend(ia[both]);right.extend(ib[both]);weight.extend(harmonic(sa[both],sb[both])*area[both]/distance)
        for gas_side,metal_side,gas_sigma,gas_ids in [(ma,mb,sa,ia),(mb,ma,sb,ib)]:
            for role,collection in [(1,hv),(2,ret)]:
                contact=(gas_ids>=0)&(gas_side==0)&(metal_side==role)
                np.add.at(collection,gas_ids[contact],gas_sigma[contact]*area[contact]/(.5*distance))
    left=np.asarray(left,dtype=int);right=np.asarray(right,dtype=int);weight=np.asarray(weight,dtype=float)
    adjacency=coo_matrix((np.ones(2*len(left)),(np.r_[left,right],np.r_[right,left])),shape=(len(nodes),len(nodes))).tocsr()
    count,labels=connected_components(adjacency,directed=False)
    h=np.zeros(count,dtype=bool);r=h.copy()
    h[labels[hv>0]]=True;r[labels[ret>0]]=True
    return Network(nodes,mask.shape,left,right,weight,hv,ret,labels,h,r)


def energy(net,u):
    return float(np.sum(net.conductance*(u[net.left]-u[net.right])**2,dtype=np.longdouble)
                 +np.sum(net.HV*(1-u)**2,dtype=np.longdouble)+np.sum(net.RETURN*u*u,dtype=np.longdouble))


def solve(net,trial=None):
    n=len(net.nodes)
    bridge=net.component_HV&net.component_RETURN
    active=bridge[net.labels]
    u=np.zeros(n)
    u[net.component_HV[net.labels]&~net.component_RETURN[net.labels]]=1
    ids=np.flatnonzero(active)
    status='OPEN_CIRCUIT';matrix_nnz=0
    if len(ids):
        diag=np.bincount(np.r_[net.left,net.right],weights=np.r_[net.conductance,net.conductance],minlength=n)+net.HV+net.RETURN
        matrix=coo_matrix((np.r_[diag,-net.conductance,-net.conductance],
                    (np.r_[np.arange(n),net.left,net.right],np.r_[np.arange(n),net.right,net.left])),shape=(n,n)).tocsr()[ids][:,ids]
        matrix_nnz=matrix.nnz
        # Exact diagonal equilibration only: no edge change, no added diagonal regularizer.
        scale=1/np.sqrt(diag[ids]);D=diags(scale)
        u[ids]=scale*spsolve(D@matrix@D,scale*net.HV[ids])
        status='SOLVED_BRIDGED_NUMERICAL_SUBGRAPH'
    flux=net.conductance*(u[net.left]-u[net.right])
    balance=np.bincount(np.r_[net.left,net.right],weights=np.r_[flux,-flux],minlength=n)+net.HV*(u-1)+net.RETURN*u
    ih=float(np.sum(net.HV*(1-u),dtype=np.longdouble));ir=float(np.sum(net.RETURN*u,dtype=np.longdouble))
    dissipated=energy(net,u)
    norm=dissipated if dissipated>0 else None
    result=dict(status=status,positive_sigma_nodes=n,gas_gas_links=len(net.left),
        connected_components=len(net.component_HV),bridging_components=int(bridge.sum()),solved_unknowns=len(ids),matrix_nnz=matrix_nnz,
        floating_components=int(np.count_nonzero(~net.component_HV&~net.component_RETURN)),
        one_terminal_components=int(np.count_nonzero(net.component_HV^net.component_RETURN)),
        floating_gauge='arbitrary zero representative; zero energy/current, no epsilon connection',
        G_frozen_S=dissipated,HV_unit_voltage_current_A=ih,RETURN_unit_voltage_current_A=ir,
        terminal_conservation_abs_A=abs(ih-ir),energy_minus_HV_current_abs_W=abs(dissipated-ih),
        max_node_current_residual_A=float(np.max(np.abs(balance))) if n else 0.,
        terminal_conservation_relative=abs(ih-ir)/norm if norm else None,
        energy_current_relative=abs(dissipated-ih)/norm if norm else None,
        max_node_residual_over_G=float(np.max(np.abs(balance)))/norm if norm and n else None,
        potential_min=float(u.min()) if n else None,potential_max=float(u.max()) if n else None,
        conditional_reciprocal_ohm=1/dissipated if dissipated>0 else None,
        reciprocal_identity='frozen discrete reciprocal only, NOT physical spark R')
    if trial is not None:
        trial=np.asarray(trial).ravel()[net.nodes]
        test=energy(net,trial)
        result.update(matched_edge_trial_energy_S=test,variational_difference_S=test-dissipated,
             frozen_over_matched_trial=dissipated/test if test else None,
             variational_identity='same network links and imposed terminal1/0; native psi samples used as an admissible trial, not re-solved harmonic field')
    return result,u
