"""Independent transcription of Morrow & Lowke (1997), Appendix A1--A11."""
from math import exp
def evaluate(E,N):
 E=abs(E)
 if E==0:return dict(mobility=6.87e24/N,diffusion=0.,ionization_frequency=0.,attachment_two_body_frequency=0.,attachment_three_body_frequency=0.)
 x=E/N*1e4; nc=N/1e6
 an=(2e-16*exp(-7.248e-15/x) if x>1.5e-15 else 6.619e-17*exp(-5.593e-15/x))
 e2=max(0.,8.889e-5*x+2.567e-19 if x>1.05e-15 else 6.089e-4*x-2.893e-19);e3=4.7778e-59*x**-1.2749
 if x>2e-15:w=7.4e21*x+7.1e6
 elif x>=1e-16:w=1.03e22*x+1.3e6
 elif x>=2.6e-17:
  x0,x1=2.6e-17,1e-16;t=(x-x0)/(x1-x0);y0=6.87e22*x0+3.38e4;y1=1.03e22*x1+1.3e6;m0=6.87e22*(x1-x0);m1=1.03e22*(x1-x0);w=(2*t**3-3*t*t+1)*y0+(t**3-2*t*t+t)*m0+(-2*t**3+3*t*t)*y1+(t**3-t*t)*m1
 else:w=6.87e22*x+3.38e4
 speed=w/100
 return dict(mobility=speed/E,diffusion=.3341e9*x**.54069*w/E*1e-2,ionization_frequency=an*nc*100*speed,attachment_two_body_frequency=e2*nc*100*speed,attachment_three_body_frequency=e3*nc*nc*100*speed)
