import numpy as np

def bernoulli(x):
    x=np.asarray(x,dtype=float); out=np.empty_like(x); small=np.abs(x)<1e-3
    y=x[small]; out[small]=1-y/2+y*y/12-y**4/720+y**6/30240
    pos=x>50; neg=x<-50; mid=~(small|pos|neg)
    out[pos]=x[pos]*np.exp(-x[pos]);out[neg]=-x[neg];out[mid]=x[mid]/np.expm1(x[mid])
    return float(out) if out.ndim==0 else out

def sg_flux(nl,nr,w,d,h):
    p=w*h/d
    return d/h*(bernoulli(-p)*np.asarray(nl)-bernoulli(p)*np.asarray(nr))

def isg0_flux(nl,nr,wl,wr,d,h,epsilon,n_ref,diagnostics=None):
    if n_ref<=0: raise ValueError("n_ref must be explicit and positive")
    dw=wr-wl; w=.5*(wl+wr)
    if abs(dw)<64*np.finfo(float).eps*max(1,abs(wl),abs(wr)):
        if diagnostics is not None: diagnostics.update(branch="small_velocity_gradient")
        return sg_flux(nl,nr,w,d,h)
    hv=np.sqrt(2*epsilon*d*h/abs(dw))
    if hv>=h:
        if diagnostics is not None: diagnostics.update(branch="virtual_width_ge_cell",h_virtual=hv)
        return sg_flux(nl,nr,w,d,h)
    ul,ur=nl/n_ref,nr/n_ref;a=(np.log1p(ur)-np.log1p(ul))/h
    nmid=n_ref*(np.sqrt((1+ul)*(1+ur))-1)
    if hv/h<1e-6:
        if diagnostics is not None: diagnostics.update(branch="zero_width",h_virtual=hv)
        return nmid*(w-d*a)
    nvl=n_ref*((1+ul)*np.exp(a*(h-hv)/2)-1);nvr=n_ref*((1+ul)*np.exp(a*(h+hv)/2)-1)
    if diagnostics is not None: diagnostics.update(branch="isg0",h_virtual=hv)
    return sg_flux(nvl,nvr,w,d,hv)

def conservative_rk2(n,dt,flux):
    def rhs(q):
        f=flux(q);return -(f-np.roll(f,1))
    q1=n+dt*rhs(n);return .5*n+.5*(q1+dt*rhs(q1))

