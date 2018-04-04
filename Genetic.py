#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:24:42 2018

@author: Andre & Henry 
@title: Genetical algorithm framework
@readme: This is a framework implementing methods for Genetic algorithms.
"""

import numpy as np
import matplotlib.pyplot as plt


# ADN : np.array([t0, dV1x, dV2x, ..., dV365x, dV1y, dV2y, ..., dV365y])

dVmax = 5 #km/s
T0max=365 #jours

nb_dV=365 #simulation de 10, un dv tous les 10 jours

nb_simu=5 #nb cycler / génération

nb_croisement=5 #nb de nouvel ADN à partir des sélectionés

nb_selectionne=3 #np d'ADN séléctionné par génération


def nouvelADN():
    A=np.random.random(nb_dV+1)
    A[0]=0
    A=np.sort(A)
    ADN_dv=A[1:]-A[:-1]
    ADN_dv=ADN_dv*dVmax
    
    phi=np.random.random(nb_dV)*2*np.pi
    
    t0=np.random.random(1)*T0max
    
    ADN=np.concatenate([t0,ADN_dv*np.sin(phi),ADN_dv*np.cos(phi)])
    return ADN


def generation1():
    gen=np.concatenate([nouvelADN() for _ in range(nb_simu)])
    gen=np.reshape(gen,(nb_simu,2*nb_dV+1))
    return gen

def croisement(séléctionés):
    nb_parents=séléctionés.shape[0]
    NouvelleGen=np.concatenate([(séléctionés[np.random.randint(0,nb_parents),:]+séléctionés[np.random.randint(0,nb_parents),:])/2 for _ in range(nb_croisement)])
    NouvelleGen=np.reshape(NouvelleGen,(nb_croisement,2*nb_dV+1))
    return NouvelleGen

def mutation(ADN):
    
    a=np.random.randint(0,2)
    if a==0:
        r=np.random.randint(0,nb_dV)
        mut=0*ADN
        phi=np.random.random(nb_dV)*2*np.pi
        mut[[r+1,2*r+1]]=[dVmax*np.cos(phi),dVmax*np.sin(phi)]
    
    else:
        mut=0*ADN
    
    ADN_muté=0.99*ADN+0.01*mut
    
    dt0=np.random.random(1)*T0max/50-T0max/100
    t0=ADN[0]
    t=max(min(t0+dt0,T0max),0)
    
    ADN_muté[0]=t
    
    return ADN_muté
    

def selection(generation):
    
    dict={}
    
    for i in range(len(generation)):
        phase=evolution_given_adn(generation[i])
        if phase!=None:
            x,y=earth(TIME).T-phase.T[:2]
            d2earth_min=np.min(x**2+y**2)
            
            x,y=mars(TIME).T-phase.T[:2]
            d2mars_min=np.min(x**2+y**2)
            
            dict[i]=d2earth_min+d2mars_min
    
    preselection=np.array([generation[i] for i in sorted(dict.keys(),key=dict.__getitem__)])
    
    n=len(preselection)
    
    return preselection[-nb_selectionne:]
    
    
    
    
    
    
    
    
    
    
    
    
    x,y=mars(TIME).T-ph.T[:2]
    d2mars_min=np.min(x**2+y**2)
    
    x,y=earth(TIME).T-ph.T[:2]
    d2earth_min=np.min(x**2+y**2)
    
    
    


        
    
    
    
    