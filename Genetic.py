#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 20:24:42 2018

@author: Andre & Henry 
@title: Genetical algorithm framework
@readme: This is a framework implementing methods for Genetic algorithms.
"""

from random import randint

# ADN : np.array([t0, dV1x, dV2x, ..., dV365x, dV1y, dV2y, ..., dV365y])

dVmax = 5 #km/s
T0max=365 #jours

nb_dV=5 #simulation de 10, un dv tous les 10 jours

nb_simu=5 #nb cycler / génération

nb_croisement=5 #nb de nouvel ADN à partir des sélectionés


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
    NouvelleGen=np.concatenate([(séléctionés[:,randint(nb_parents)]+séléctionés[:,randint(nb_parents)])/2 for _ in range(np_croisement)])
    NouvelleGen=np.reshape(nb_croisement,(nb_croisement,2*nb_dV+1))