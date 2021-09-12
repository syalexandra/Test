#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 21:25:39 2020

@author: yuesun
"""
from Point import *
import numpy as np
from matplotlib import pyplot
from scipy.special import binom
from quadtree import *


def multipole(sources,center,P):
    #sources is the list of points
    qhat=np.zeros(P,dtype=complex)
    for s in sources:
        qhat[0]+=s.q
        for i in range(1,P):
            qhat[i]-=complex(s.x-center.x,s.y-center.y)**i/i*s.q
    return qhat



#lemma 2.3 in original paper
def multipleShift(qhat,z0):
    newqhat=np.empty_like(qhat)
    newqhat[0]=qhat[0]
    for l in range(1,len(qhat)):
        newqhat[l]=sum([qhat[k]*complex(z0.x,z0.y)**(l-k)*binom(l-1,k-1) for k in range(1,l+1)])-qhat[0]*complex(z0.x,z0.y)**l/l
    return newqhat
    
    
    
def outgoingExpansion(tnode,P=10):
    if tnode.isleaf():
        tnode.outmultipole=multipole(tnode.getPoints(),tnode.center,P)
    else:
        tnode.outmultipole=np.zeros(P,dtype=complex)
        for child in tnode.children:
            outgoingExpansion(child,P)
            z0=child.center-tnode.center
            tnode.outmultipole+=multipleShift(child.outmultipole,z0)
            
            

            
            
            