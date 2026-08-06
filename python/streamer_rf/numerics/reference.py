import numpy as np
EPS0=8.8541878128e-12
def ring_axis_potential(charge,radius,dz): return charge/(4*np.pi*EPS0*np.sqrt(radius**2+dz**2))
def zheleznyak_g(R,p_o2_torr=150):
    x=p_o2_torr*np.asarray(R)*100 # Torr cm
    return p_o2_torr*100*np.where(x>0,(np.exp(-.035*x)-np.exp(-2*x))/(x*np.log(2/.035)),0)
def sp3_kernel_fit(R,p_o2_torr=150):
    A=np.array([.0067,.0346,.3059]);lam=np.array([.0447,.1121,.5994]);x=p_o2_torr*np.asarray(R)*100
    return p_o2_torr*100*np.sum(A[:,None]*np.exp(-lam[:,None]*x.ravel()),axis=0).reshape(np.asarray(R).shape)
def sp3_constants(p_o2_torr=150):
    s=np.sqrt(6/5);return {"kappa1_sq":3/7-2*s/7,"kappa2_sq":3/7+2*s/7,
      "gamma1":5/7*(1-3*s),"gamma2":5/7*(1+3*s),
      "A_original":[.0067,.0346,.3059],"lambda_original":[.0447,.1121,.5994],
      "a_si":(np.array([.0067,.0346,.3059])*p_o2_torr*100).tolist(),
      "k_si":(np.array([.0447,.1121,.5994])*p_o2_torr*100).tolist()}

