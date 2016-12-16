import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from p3 import sync #assumes that fortran module sync has been compiled with f2py as p3.so



def oscillator(Nt,T,N,c,s,mu=1):
    """Simulate fully-coupled network of oscillators
    Compute simulation for Nt time steps between 0 and T (inclusive)
    input: N: number of oscillators
           c: coupling coefficient
           mu,sigma: distribution parameters for omega_i
    """
    sync.c = c
    sync.w = np.random.normal(mu,s,N)
    initialphases = np.random.uniform(0,2*np.pi,N)
    
    t0 = float(0)
    
    
    
    y, order = sync.rk4(float(0),initialphases,float(T)/float(Nt),Nt)
    
    
    time = np.linspace(0,T,Nt)
    return time,sync.w,initialphases,y,order
    

if __name__ == '__main__':
    n,c,m,s = 101,10.0,1.0,0.1
    Nt,T = 500,100
    t,omega,theta0,theta,order = oscillator(Nt,T,n,c,s,m) #i swapped m and s
