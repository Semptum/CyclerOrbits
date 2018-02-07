#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:03:56 2018

@author: andre
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

#Constants definition
G=4*np.pi**2
SIMULATION_LENGTH=2
STEP=1/365
TIME=np.arange(0,SIMULATION_LENGTH,STEP)
n=TIME.shape[0]
R=np.zeros((n,2))
PR=np.zeros((n,2))
eM=3.003e-6
#Masses are in solar masses 1 sM= 2e30 kg, so 1 kg = 0.5e-30 sM
#Lengths are in AU
#Durations are in AU


def position_func(apogee,perigee,argument,mass,phase=0,):
    if apogee*perigee==0:
        return (lambda t:[np.zeros_like(t),np.zeros_like(t)]),mass
    c=(apogee-perigee)/2
    a=(perigee+apogee)/2
    e=c/a
    b=a*np.sqrt(1-e**2)
    T=a**(2/3)
    angle=np.pi/180
    rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return lambda t:rot.dot(np.array([a*np.cos(2*np.pi*t/T+phase),b*np.sin(2*np.pi*t/T+phase)])),mass

earth=position_func(1.0167,0.98329,288.1,eM)
mars=position_func(1.666,1.3814,286.5,0.107*eM)
sun=position_func(0,0,0,1)

planets=[sun,earth,mars]

init=earth[0](0)+np.array([-0.1,0])
pinit=(earth[0](1)-earth[0](0))/STEP

def acc(planets,t,pos):
    FORCE=np.zeros(2)
    for p in planets:
        r=(p[0](t)-pos)
        D=(r.dot(r))**0.5
        assert D>0.01
        FORCE+=G*p[1]*r/D**3
    return FORCE*STEP

for i in range(1,len(TIME)):
    t=i*STEP
    R[i]=R[i-1]+PR[i-1]
    PR[i]=PR[i-1]+acc(planets,t,R[i])*STEP
R=R.T

for k in planets:
    x,y=k[0](TIME)
    plt.plot(x,y)
    plt.plot(x[-1],y[-1],"o")
x,y=R
plt.plot(x,y)
plt.plot(0,0,"x")
plt.plot(x[-1],y[-1],"o")