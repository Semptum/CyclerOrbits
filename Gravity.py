#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:16:30 2018

@author: Andre & Henry 
@title: Gravity framework
@readme: This is a set of functions designed to facilitate simulation of the movement of a spaceship through 
the grav field of a set of planets moving on rails. Multiple integrators are provided.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


#Constants definition
G=4*np.pi**2
eM=3.003e-6
km=1/(1.496e+8)
rearth=6400*km
#moon=position_func()
#Masses are in solar masses 1 sM= 2e30 kg, so 1 kg = 0.5e-30 sM
#Lengths are in AU
#Each celestial body has two invariant parameters: mass (in solar masses) and deathradius (in AU). When the spaceship enters the
#deathradius it ceases to exist, and an error is thrown


#Although inconsistent, it was more appropriate to use 2-vectors for planets and 4-vectors for the spaceship in order
#to perform the integration

def position_func(apogee,perigee,argument,phase=0,refbody=lambda t:np.zeros(np.shape(t)+(2,))):
    """
    This function returns a function which return the position of the planet with the specified parameters
    The parameters are:
        -apogee: self evident
        -perigee: alse self evident
        -argument: the angle between the x-axis and the semi-major axis
        -phase: the moment in the orbit at the wchich the planet starts on its orbit
        -refbody: a predefined reference body which dominates the gravitational force acting on the body. For instance, earth to the moon
    """
    if apogee*perigee==0:
        return lambda t:np.zeros(np.shape(t)+(2,))+refbody(t)
    c=(apogee-perigee)/2
    a=(perigee+apogee)/2
    e=c/a
    b=a*np.sqrt(1-e**2)
    T=a**(2/3)
    angle=argument*np.pi/180
    rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    return lambda t:rot.dot(np.array([a*np.cos(2*np.pi*t/T+phase)-a+perigee,#
                                      b*np.sin(2*np.pi*t/T+phase)])).T+refbody(t)

sun=position_func(0,0,0,1)
earth=position_func(1.0167,0.98329,288.1)
mars=position_func(1.666,1.3814,286.5)
PLANETS={sun:[1,0.05],mars:[0.107*eM,50*km],earth:[eM,0]}
STEP=1/365
TIME=np.arange(0,10,STEP)

def closest_dist(phasevect,PLANETS,t):
    R=10 #pk pas
    for p in PLANETS.keys():
        phaseplan=p(t)
        diff=phasevect[:2]-phaseplan
        D=diff[:2]
        R=min(R,(D.dot(D))**(1/2))
    return R


    
def acceleration(phasevec,t,PLANETS):
    """
    This function simply returns the derivative of the 4-vector position-acceleration of the
    spaceship
    """
    dphasevec=np.zeros(4)
    for p in PLANETS.keys():# p est mtn la fonction qui a t associe la position d'une planete
        phaseplan=p(t)
        D=phasevec[:2]-phaseplan
        R=(D.dot(D))**(1/2)
        assert R<10, "trop loin"
        assert R>PLANETS[p][1],PLANETS[p]
        dphasevec[2:]-=G*PLANETS[p][0]*D/R**3
    dphasevec[:2]=phasevec[2:]
    return dphasevec


def RK4(phasevect0,equa_diff,time,PLANETS,n=0):
    assert n<2,"Trop profond"
    N=time.shape[0]
    PHASE=np.zeros((N,4))
    PHASE[0]=phasevect0
    h=time[1]-time[0]
    for i in range(time[:-1].shape[0]):
        t=time[i]
        k1=equa_diff(PHASE[i],t,PLANETS)*h
        k2=equa_diff(PHASE[i]+k1/2,t+h/2,PLANETS)*h
        k3=equa_diff(PHASE[i]+k2/2,t+h/2,PLANETS)*h
        k4=equa_diff(PHASE[i]+k3,t+h,PLANETS)*h
        dph=1/6*(k1+2*k2+2*k3+k4)
        rmin=closest_dist(PHASE[i],PLANETS,t)/km
        if rmin<10000 and n<1:
            print(rmin)
            PHASE[i+1]=RK4(PHASE[i],equa_diff,np.linspace(t,t+h,10),PLANETS,2)[-1]
        else:
            PHASE[i+1]=PHASE[i]+dph
    return PHASE
        

def plot(phasevec,PLANETS,TIME,refPlanet=None):
    if refPlanet==None:
        ref=phasevec.T[:2]
    else:
        ref=refPlanet(TIME).T
    x,y=phasevec.T[:2]-ref
    plt.plot(x,y)
    plt.plot(x[-1],y[-1],"x")
    for p in PLANETS.keys():
        x,y=p(TIME).T-ref
        plt.plot(x,y)
        plt.plot(x[-1],y[-1],"o")
    plt.axis("equal")


def initial(t0,leo=200):
    earth=position_func(1.0167,0.98329,288.1)
    STEP=1/365
    R=rearth+leo*km
    v=(G*eM/R)**0.5
    period=2*np.pi/(G*eM/R**3)**2
    dt=t0%period*STEP
    init=np.concatenate((earth(t0),(earth(STEP)-earth(0))/STEP))
    init[:2]+=R*np.array([np.cos(dt),np.sin(dt)])
    init[2:]+=v*np.array([-np.sin(dt),np.cos(dt)])
    print((R/km,v/km/(365*24*3600)))
    return init



def evolution_given_adn(adn):
    daytime,corrections=adn[0],adn[1:]
    day=int(daytime)
    time=(daytime-day)*24
    try:
        deltai=10
        dv0=corrections[0]
        init=initial(day+time/24,0)
        TMAX=10
        DAY=1/365
        TIME=np.arange(0,TMAX,DAY)
        PHASE=np.zeros((TIME.shape[0],4))
        PHASE[0]=init
        for i in range(corrections.shape[0]//2):
            corr=corrections[i],corrections[i+corrections.shape[0]//2]
            time=TIME[deltai*i:(i+1)*deltai+1]
            init=PHASE[i*deltai]
            init[2:]+=corr
            PHASE[deltai*i:(i+1)*deltai+1]=RK4(init,acceleration,time,PLANETS)
        return PHASE
    except:
        return None