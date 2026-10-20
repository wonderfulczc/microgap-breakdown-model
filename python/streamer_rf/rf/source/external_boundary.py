"""Source-bound linear electrostatic response, not a new field solver.

The constrained free block A is SPD, B=K_free,terminal, D=K_terminal,terminal.
Charges q are independent free-node sources. Terminal charge is B.T u + D V.
The reference/gauge and physical energy network must be supplied, not inferred
from a mixed-unit Poisson matrix or an assumed extra capacitance.
"""
from dataclasses import dataclass
import hashlib
from pathlib import Path
import numpy as np
from .global_packet import require, finite

VERSION = 'I1_EXTERNAL_BOUNDARY_RESPONSE_V1'


@dataclass(frozen=True)
class I1ExternalBoundaryResponse:
    contract: dict
    A: np.ndarray
    B: np.ndarray
    D: np.ndarray
    local_indices: tuple

    def validate(self):
        c = self.contract
        require(c.get('version') == VERSION, 'BOUNDARY_OPERATOR_VERSION')
        for key in ('operator_id','geometry_id','reference_plane_id','common_boundary_id',
                    'local_control_volume','global_control_volume','gauge', 'terminal_ids'):
            require(c.get(key), 'BOUNDARY_OPERATOR_IDENTITY_MISSING', key)
        require(c.get('units') == {'matrix':'F','potential':'V','charge':'C','energy':'J'},
                'BOUNDARY_OPERATOR_UNIT_MISMATCH')
        require(c.get('energy_network') == 'SYMMETRIC_CONSTRAINED_ELECTROSTATIC',
                'NOT_AN_ENERGY_NETWORK')
        require(c.get('capacitance_owner') in ('TOTAL_CHARGE','SPLIT_GEOMETRIC_CAPACITANCE')
                and c.get('capacitance_copies') == 1, 'CAPACITANCE_DOUBLE_COUNT')
        source = Path(c.get('source_file',''))
        require(source.is_file() and hashlib.sha256(source.read_bytes()).hexdigest() == c.get('source_sha256'),
                'BOUNDARY_OPERATOR_SOURCE_MISMATCH')
        require(c.get('evidence_class') in ('SYNTHETIC_TEST_OPERATOR','SOURCE_BOUND_OPERATOR'),
                'BOUNDARY_OPERATOR_SOURCE_CLASS')
        # Source-bound payloads must bind arrays, not only a descriptive paper/config.
        if c['evidence_class'] == 'SOURCE_BOUND_OPERATOR':
            try:
                with np.load(source,allow_pickle=False) as data:
                    matches=all(np.array_equal(data[k],v) for k,v in
                                (('A',self.A),('B',self.B),('D',self.D),('local_indices',self.local_indices)))
            except (ValueError,KeyError,OSError):
                require(False,'BOUNDARY_OPERATOR_PAYLOAD_UNREADABLE')
            require(matches,'BOUNDARY_OPERATOR_PAYLOAD_MISMATCH')
        n = len(self.A); nt = len(c['terminal_ids'])
        A = finite(self.A,(n,n),'A'); B = finite(self.B,(n,nt),'B'); D = finite(self.D,(nt,nt),'D')
        require(np.array_equal(A,A.T) and np.array_equal(D,D.T), 'NONSYMMETRIC_ENERGY_OPERATOR')
        try: np.linalg.cholesky(A)
        except np.linalg.LinAlgError: require(False,'UNCONSTRAINED_OR_NONPOSITIVE_OPERATOR')
        full = np.block([[A,B],[B.T,D]])
        scale = float(np.linalg.norm(full,2))
        require(scale > 0 and np.linalg.eigvalsh(full).min() >= -1e-12*scale,
                'NEGATIVE_ELECTROSTATIC_ENERGY')
        l = np.asarray(self.local_indices)
        require(l.ndim==1 and np.issubdtype(l.dtype,np.integer) and 0<len(l)<n
                and len(set(l.tolist()))==len(l) and np.all((l>=0)&(l<n)), 'INVALID_DOMAIN_PARTITION')
        return A,B,D,l,np.array([k for k in range(n) if k not in l])

    def evaluate(self, free_charge_C, terminal_voltage_V):
        A,B,D,l,o = self.validate()
        q = finite(free_charge_C,(len(A),),'q'); V = finite(terminal_voltage_V,(len(D),),'V')
        All,Alo,Aoo = A[np.ix_(l,l)],A[np.ix_(l,o)],A[np.ix_(o,o)]
        Bl,Bo = B[l],B[o]; ql,qo=q[l],q[o]
        X=np.linalg.solve(Aoo,Alo.T); Y=np.linalg.solve(Aoo,Bo); z=np.linalg.solve(Aoo,qo)
        Ae=All-Alo@X; Be=Bl-Alo@Y; De=D-Bo.T@Y; qe=ql-Alo@z
        ul=np.linalg.solve(Ae,qe-Be@V)
        uo=np.linalg.solve(Aoo,qo-Alo.T@ul-Bo@V)
        u=np.empty(len(A));u[l]=ul;u[o]=uo
        terminal_q=Be.T@ul+De@V+Bo.T@z
        # Completing the square in outer unknowns retains the q_o energy.
        reduced=np.block([[Ae,Be],[Be.T,De]]); v=np.r_[ul,V]
        local_energy=.5*v@reduced@v
        outer_source_energy=.5*qo@z
        full_v=np.r_[u,V];full_K=np.block([[A,B],[B.T,D]])
        return dict(free_potential_V=u,local_potential_V=ul,outer_potential_V=uo,
                    effective_operator_F=Ae,effective_source_C=qe-Be@V,
                    terminal_charge_C=terminal_q,outer_charge_offset_C=Bo.T@z,
                    reduced_energy_J=float(local_energy),outer_source_energy_J=float(outer_source_energy),
                    total_energy_J=float(local_energy+outer_source_energy),
                    reconstructed_full_energy_J=float(.5*full_v@full_K@full_v),
                    capacitance_owner=self.contract['capacitance_owner'])
