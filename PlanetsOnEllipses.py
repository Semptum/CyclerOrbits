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
SIMULATION_LENGTH=10
STEP=1/365/2
TIME=np.arange(0,SIMULATION_LENGTH,STEP)
eM=3.003e-6
#Masses are in solar masses 1 sM= 2e30 kg, so 1 kg = 0.5e-30 sM
#Lengths are in AU
#Durations are in AU


def position_func(apogee,perigee,argument,phase=0):
    if apogee*perigee==0:
        return lambda t:np.array([np.zeros_like(t),np.zeros_like(t),np.zeros_like(t),np.zeros_like(t)])
    c=(apogee-perigee)/2
    a=(perigee+apogee)/2
    e=c/a
    b=a*np.sqrt(1-e**2)
    T=a**(2/3)
    angle=argument*np.pi/180
    rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return lambda t:np.concatenate((rot.dot(np.array([a*np.cos(2*np.pi*t/T+phase)-a+perigee,b*np.sin(2*np.pi*t/T+phase)])),np.expand_dims(np.zeros_like(t),0),np.expand_dims(np.zeros_like(t),0)))

earth=position_func(1.0167,0.98329,288.1)
mars=position_func(1.666,1.3814,286.5)
sun=position_func(0,0,0,1)


PLANETS={sun:1,mars:0.107*eM,earth:eM}

def acceleration(phasevec,t):
    dphasevec=np.zeros(4)
    for p in PLANETS.keys():# p est mtn la fonction qui a t associe la position d'une planete
        phaseplan=p(t)
        diff=phasevec-phaseplan
        D=diff[:2]
        R=(D.dot(D))**(3/2)
        dphasevec[2:]-=G*PLANETS[p]*D/R
    dphasevec[:2]=phasevec[2:]
    return dphasevec

def plot(phase):
    x,y=phase.T[:2]
    plt.plot(x,y)
    plt.plot(x[-1],y[-1],"x")
    for p in PLANETS.keys():
        x,y=p(TIME)[:2]
        plt.plot(x,y)
        plt.plot(x[-1],y[-1],"o")
    plt.axis("equal")
init=earth(0)*(1+0.2*np.random.random())
init[2:]=(earth(STEP)-earth(0))[:2]/STEP

phase=integ.odeint(acceleration,init,TIME)