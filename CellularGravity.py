#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 00:44:39 2018

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


G=4*np.pi**2

N=1000
TMAX=125
DENSITY=np.zeros((N,N))
TIME=np.arange(0,TMAX,1)
k1=np.round(np.linspace(-N//2,N//2,N))
K1=np.meshgrid(k1,k1)

K2=K1[0]**2+K1[1]**2
#K2+=1/(1+K2)

def planet(coord,size,mass):
    result=np.zeros((N,N))
    k1=np.round(np.linspace(0,N,N))
    temp=np.meshgrid(k1-coord[0],k1-coord[1])
    a=1/(1+(temp[0]**2+temp[1]**2)/size**2)
    result[np.where(temp[0]**2+temp[1]**2<size**2)]=mass
    return result*a
DENSITY+=planet((400,400),20,1000)
DENSITY+=planet((450,600),20,1000)
DENSITY+=planet((550,400),20,1000)


FOUR=np.fft.ifftshift(np.fft.fftshift(G*np.fft.fft2(DENSITY))/K2)
FIELD=-np.real(np.fft.ifft2(FOUR))
#FIELD[np.where(DENSITY>0.5)]=0
FORCE=np.array(np.gradient(FIELD))
DENSITY=DENSITY.T

coord=np.array([800,800,0,-0.2])

def acc(coord,t):
    x,y=np.round(coord[:2])
    if DENSITY[x,y]>0.5:
        return np.array([0,0,0,0])
    a=-FORCE[:,x,y]
    return np.concatenate([coord[2:],a])

R=integ.odeint(acc,coord,TIME,printmessg=True)
x,y=R.T[:2]