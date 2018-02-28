#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created during a boring boring session of projnum

@author: henry
"""

import numpy as np
import matplotlib.pyplot as plt


def equa_diff(x,v,t):
    return v-np.exp(t)*x

x0=np.sin(1)
v0=np.cos(1)
time=np.linspace(0,10,2000)

def RK4(x0,v0,equa_diff,time,n=0):
    assert n>2,"Trop profond"
    X=[x0]
    V=[v0]
    h=time[1]-time[0]
    for t in time[:-1]:
        
        k1=equa_diff(X[-1],V[-1],t)
        k2=equa_diff(X[-1]+h/2*V[-1],V[-1]+h/2*k1,t+h/2)
        k3=equa_diff(X[-1]+h/2*V[-1]+h**2/4*k1,V[-1]+h/2*k2,t+h/2)
        k4=equa_diff(X[-1]+h*V[-1]+h**2/2*k2,V[-1]+h*k3,t+h)
        dx=h*V[-1]+h**2/6*(k1+k2+k3)
        dv=h/6*(k1+2*k2+2*k3+k4)
        if abs(dv)>1:
            X.append(X[-1]+dx)
            V.append(V[-1]+dv)
        else:
            time2=np.linspace(t,t+h,10)
            result = RK4(X[-1],V[-1],equa_diff,time2,n+1)
            X.append(result[0][-1])
            V.append(result[0][-1])
    return np.array(X),np.array(V)

plt.plot(RK4(x0,v0,equa_diff,time)[0])


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





    
