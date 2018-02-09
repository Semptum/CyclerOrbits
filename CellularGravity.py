#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 00:44:39 2018

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm


G=4*np.pi**2

N=1000
TMAX=1000
DENSITY=np.zeros((N,N))
TIME=np.arange(0,TMAX,0.1)
k1=np.round(np.linspace(-N//2,N//2,N))
K1=np.meshgrid(k1,k1)

K2=K1[0]**2+K1[1]**2
K2+=1/(1+K2)

def planet(coord,size,mass):
    result=np.zeros((N,N))
    k1=np.round(np.linspace(0,N,N))
    temp=np.meshgrid(k1-coord[0],k1-coord[1])
    a=1/(1+(temp[0]**2+temp[1]**2)/size**2)
    result[np.where(temp[0]**2+temp[1]**2<size**2)]=mass
    return result*a
DENSITY+=planet((400,400),20,100)
DENSITY+=planet((450,600),20,1000)
DENSITY+=planet((550,400),20,0)


FOUR=np.fft.ifftshift(np.fft.fftshift(G*np.fft.fft2(DENSITY))/K2)
FIELD=-np.real(np.fft.ifft2(FOUR))
#FIELD[np.where(DENSITY>0.5)]=0
FORCE=np.array(np.gradient(FIELD))
DENSITY=DENSITY
FIELD=FIELD

coord=np.array([600,600,-5,0])
def plot(X,Y):
    t = np.linspace(0,1,X.shape[0])
    dx=X[1:]-X[:-1]
    dy=Y[1:]-Y[:-1]
    dr=(dx**2+dy**2)**0.5
    dr/=dr.max()
    points = np.array([X,Y]).transpose().reshape(-1,1,2)
    segs = np.concatenate([points[:-1],points[1:]],axis=1)
    lc = LineCollection(segs, cmap=plt.get_cmap('jet'))
    lc.set_array(dr)
    plt.gca().add_collection(lc)


def acc(coord,t):
    x,y=np.round(coord[:2])
    if not 0<=x<N or not 0<=y<N:
        return np.array([0,0,0,0])
    if DENSITY[x,y]>0.5:
        return np.array([0,0,0,0])
    a=-FORCE[:,x,y]
    return np.concatenate([coord[2:],a])

R,msg=integ.odeint(acc,coord,TIME,printmessg=True,full_output=1)
y,x=R.T[:2]
plt.imshow(DENSITY)
#plt.contourf(FIELD,10)
plot(x,y)
