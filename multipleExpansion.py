#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 13:36:06 2020

@author: yuesun
"""

from Point import *
import numpy as np
from matplotlib import pyplot


def sourceToCenter(sources,center,P):
    qhat=np.zeros(P,dtype=complex)
    for s in sources:
        qhat[0]+=s.q
        for i in range(1,P):
            qhat[i]-=complex(s.x-center.x,s.y-center.y)**i/i*s.q
    return qhat
    

def multipole(targets,center,P,qhat):
    m=len(targets)
    v=np.zeros(m,dtype=complex)
    weight=np.zeros(P,dtype=complex)
    for j in range(m):
        t=targets[j]
        weight[0]=np.log(complex(t.x-center.x,t.y-center.y))
        for p in range(1,P):
            weight[p]=1.0/complex(t.x-center.x,t.y-center.y)**p
        v[j]=np.dot(weight,qhat)
    return v


def direct_sources_to_targets(targets,sources):
    m=len(targets)
    v=np.zeros(m,dtype=complex)
    for i,target in enumerate(targets):
        for source in sources:
            kernel=np.log(complex(target.x-source.x,target.y-source.y))
            v[i]+=source.q*kernel
    return v
            
            
def z_to_z0(z,z0):
    return np.log(complex(z.x-z0.x,z.y-z0.y))


def z_to_z0_expansion(z,z0):
    summ=0
    for k in range(1,10):
        summ+=(complex(z0.x,z0.y)**k)/(complex(z.x,z.y)**k)/k
    return np.log(complex(z.x,z.y))-summ
            
if __name__=='__main__':
    """
    p1=Point(1,2)
    p2=Point(3,4)
    print(distance(p1,p2))
    pp1=Particle(1,2,3)
    pp2=Particle(4,4,5)
    print(distance(pp1,pp2))
    """
    n=100
    sources=[Particle(x,y,1) for x,y in np.random.random((n,2))]
    targets=[Particle(x,y,1) for x,y in -np.random.random((n,2))]
    
    
    #targets=[Particle(5,1,1)]
    #sources=[Particle(1,2,1)]
    
    center=Point(0.5,0.5)
    #print(distance(sources[0],center))
    #print(distance(targets[0],center))
    P=10
    qhat=sourceToCenter(sources,center,P)
    u=multipole(targets,center,P,qhat)
    u1=direct_sources_to_targets(targets,sources)

    #print(qhat)
    error=0
    for i in range(n):
        error+=(u[i]-u1[i])**2
    print(error/n)
    
    """
    z=Particle(5,1,1)
    z0=Particle(1,2,1)
    print(z_to_z0(z,z0))
    print(z_to_z0_expansion(z,z0))
    """