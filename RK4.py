import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ





def equa_diff(x,v,t):
    return v-np.exp(t)*x

x0=np.sin(1)
v0=np.cos(1)
time=np.linspace(0,10,500)

def RK4(x0,v0,equa_diff,time):
    
    X=[x0]
    V=[v0]
    
    h=time[1]-time[0]
    
    for t in time[:-1]:
        
        k1=equa_diff(X[-1],V[-1],t)
        k2=equa_diff(X[-1]+h/2*V[-1],V[-1]+h/2*k1,t+h/2)
        k3=equa_diff(X[-1]+h/2*V[-1]+h**2/4*k1,V[-1]+h/2*k2,t+h/2)
        k4=equa_diff(X[-1]+h*V[-1]+h**2/2*k2,V[-1]+h*k3,t+h)
        
        X.append(X[-1]+h*V[-1]+h**2/6*(k1+k2+k3))
        V.append(V[-1]+h/6*(k1+2*k2+2*k3+k4))
    
    return np.array(X)

#plt.plot(RK4(x0,v0,equa_diff,time))


def testRK(a,b):
    H=[]
    err=[]
    eq=lambda x,v,t: -x
    for i in range(a,b+1): 
    
        h=10**(-i)
        time=np.linspace(0,10,10/h+1)
        X=RK4(1,0,eq,time)
        Y=[np.cos(t) for t in time]
        erreur= max(abs(X[t]-Y[t]) for t in range(len(X)))
    
    
        H.append(h)
        err.append(erreur)
    
    plt.loglog(H,err)





    
