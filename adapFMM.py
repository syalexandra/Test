#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 13:49:32 2020

@author: yuesun
"""

from Point import *
import numpy as np
from matplotlib import pyplot
from scipy.special import binom
from adapQuadtree import *
import time
import random
import math
import copy


def directSource(particles):
    #particles=particles.tolist()
    for i, particle in enumerate(particles):
        for source in particles[:i]+particles[i+1:]:
            
            kernel=np.log(complex(particle.x-source.x,particle.y-source.y))
            particle.u += source.q*kernel
            
            
            
def multipole(sources,center,P=10):
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
    
    
def setAllZero(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        tnode=queue.pop(0)
        
        tnode.inmultipole=np.zeros(P, dtype=complex)
        tnode.outmultipole=np.zeros(P, dtype=complex)
        for tnode in tnode.children:
        
            queue.append(tnode)
            
            
"""   
def outgoingExpansion(tnode,P=10):
    if tnode.isleaf():
        tnode.outmultipole=multipole(tnode.getPoints(),tnode.center,P)
    
    else:
        
        for child in tnode.children:
            outgoingExpansion(child,P)
            z0=child.center-tnode.center
            tnode.outmultipole+=multipleShift(child.outmultipole,z0)
""" 

def outgoingExpansion(tnode,P=10):
    if tnode.isleaf():
        tnode.outmultipole=multipole(tnode.getPoints(),tnode.center,P)
    else:
        tnode.outmultipole=np.zeros(P,dtype=complex)
        for child in tnode.children:
            outgoingExpansion(child,P)
            z0=child.center-tnode.center
            tnode.outmultipole+=multipleShift(child.outmultipole,z0)           
            
#lemma 2.4 in original paper
def inFromOut(qhat,z0):
    z0=complex(z0.x,z0.y)
    newqhat=np.empty_like(qhat)
    P=len(qhat)
    newqhat[0]=qhat[0]*np.log(-z0)+sum([qhat[k]*(-1)**k/z0**k for k in range(1,P)])
    for l in range(1,P):
        newqhat[l]=sum([qhat[k]*binom(l+k-1,k-1)*(-1)**k/(z0**k) for k in range(1,P)])/(z0**l)-qhat[0]/(l*z0**l)
    return newqhat



def shiftInnerExpansion(qhat,z0):
    newqhat=np.empty_like(qhat)
    z0=complex(z0.x,z0.y)
    P=len(qhat)
    for l in range(P):
        
        newqhat[l]=sum([qhat[k]*binom(k,l)*(-z0)**(k-l) for k in range(l,P)])
    
    return newqhat


def directSourceToTargets(targets,sources):
    #targets abd sources are array of particles
    for i,target in enumerate(targets):
        for source in sources:
            kernel=np.log(complex(target.x-source.x,target.y-source.y))
            target.u+=source.q*kernel
            
            
            
def directSourceInOneLeaf(particles):
    particles=particles.tolist()
    for i, particle in enumerate(particles):
        for source in particles[:i]+particles[i+1:]:
            
            kernel=np.log(complex(particle.x-source.x,particle.y-source.y))
            particle.u += source.q*kernel
            

                    
def evaluateFarField(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
            
            tnode.inmultipole=np.zeros(P, dtype=complex)
            
            for tin in tnode.getInteractionList():
                
                z0=tnode.center-tin.center
                IFO=inFromOut(tnode.outmultipole,z0)
                #coeffs=tnode.inmultipole
                z0=tin.center
                for p in tin.getPoints():
                    z=p-z0
                    p.u+=np.polyval(IFO[::-1],complex(z.x,z.y))
                
            if tnode.isleaf():
                
                for nn in tnode.getNeighbors():
                    directSourceToTargets(tnode.getPoints(),nn.getPoints())
                
                directSourceInOneLeaf(tnode.getPoints())    
                


def targetFromOut(tin,tnode):
    z0=tnode.center-tin.center
    IFO=inFromOut(tnode.outmultipole,z0)
    #coeffs=tnode.inmultipole
    z0=tin.center
    for p in tin.getPoints():
        z=p-z0
        p.u+=np.polyval(IFO[::-1],complex(z.x,z.y))
    


def inFromSource(tin,tnode,P=10):
    
    qhat=multipole(tnode.getPoints(),tnode.center,P)
    z0=tin.center-tnode.center
    newqhat=inFromOut(qhat,z0)
    
    return newqhat
    
def targetFromIn(tnode):
    coeffs=tnode.inmultipole
    z0=tnode.center
    for p in tnode.getPoints():
        z=p-z0
        p.u=np.polyval(coeffs[::-1],complex(z.x,z.y))
            
"""               
def incomingExpansionAdap(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
            tnode.inmultipole=np.zeros(P, dtype=complex)
            
            for tin in tnode.interactionList:
                z0=tin.center-tnode.center
                #T^{ifo}*q
                tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
            
            
            for tin in tnode.l4List:
                tnode.inmultipole+=inFromSource(tin,tnode)
                    
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        
        for tnode in node.children:
            if tnode.level>1:
                queue.append(tnode)
                
                
                    
                for child in tnode.children:
                    z0=tnode.center-child.center
                    #T^{ifi}*(Uv)
                    child.inmultipole+=shiftInnerExpansion(tnode.inmultipole,z0)
                        
                        
            if tnode.isleaf():
                targetFromIn(tnode)
                
                
            
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        
        for tnode in node.children: 
            queue.append(tnode)
            
            for tin in tnode.l3List:
                
                targetFromOut(tin,tnode)
            
            if tnode.isleaf():
                for nn in tnode.neighborList:
                    directSourceToTargets(tnode.getPoints(),nn.getPoints())
                
                directSourceInOneLeaf(tnode.getPoints()) 
                              
                   
            
def incomingExpansionAdap(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        tnode=queue.pop(0)
        
        tnode.inmultipole=np.zeros(P, dtype=complex)
        
        for tin in tnode.interactionList:
            z0=tin.center-tnode.center
            #T^{ifo}*q
            tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
        
        
        for tin in tnode.l4List:
            tnode.inmultipole+=inFromSource(tin,tnode)
            
        for tnode in tnode.children:
        
            queue.append(tnode)
            
            
                    
    queue=[]
    queue.append(root)
    while len(queue)>0:
        tnode=queue.pop(0)
        
        
        if tnode.level>1:
                
            for child in tnode.children:
                z0=tnode.center-child.center
                #T^{ifi}*(Uv)
                child.inmultipole+=shiftInnerExpansion(tnode.inmultipole,z0)
                
                        
        if tnode.isleaf():
            targetFromIn(tnode)
            
        for tnode in tnode.children:
            queue.append(tnode)
                
            
    queue=[]
    queue.append(root)
    while len(queue)>0:
        tnode=queue.pop(0)
            
        for tin in tnode.l3List:
            
            targetFromOut(tin,tnode)
        
        if tnode.isleaf():
            for nn in tnode.neighborList:
                directSourceToTargets(tnode.getPoints(),nn.getPoints())
            
            directSourceInOneLeaf(tnode.getPoints())         
                    
        for tnode in tnode.children: 
            queue.append(tnode)           
            
        
"""

def incomingExpansionAdap(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
            
            
            z0=tnode.parent.center-tnode.center
            #T^{ifi}*(Uv)
            tnode.inmultipole=shiftInnerExpansion(tnode.parent.inmultipole,z0)
            
            for tin in tnode.interactionList:
                z0=tin.center-tnode.center
                #T^{ifo}*q
                tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
                
              
            if tnode.isleaf():
                
                
                for tin in tnode.l3List:
                    z0=tin.center-tnode.center
                    #T^{ifo}*q
                    #print(tin.x,tin.y,tin.width,tin.height)
                    #print(tin.outmultipole)
                    tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
                 
                for tin in tnode.l4List:
                    z0=tin.center-tnode.center
                    #T^{ifo}*q
                    #print(tin.x,tin.y,tin.width,tin.height)
                    #print(tin.outmultipole)
                    tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
                    
def summation(root,P):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
                  
            if tnode.isleaf():        
                coeffs=tnode.inmultipole
                z0=tnode.center
                for p in tnode.getPoints():
                    z=p-z0
                    p.u+=np.polyval(coeffs[::-1],complex(z.x,z.y))
                
                for nn in tnode.neighborList:
                    directSourceToTargets(tnode.getPoints(),nn.getPoints())
                
                directSourceInOneLeaf(tnode.getPoints())    
        
    
    
if __name__=='__main__':
    #draw a quad tree, n is number of points, l is the range of the random number, P is the cutoff term
    """
    n=1000
    m=80
    P=10
    
    random.seed(1)
    particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
    """
    n=1000
    m=100
    P=10
    
    
    X=[]
    Y=[]
    for i in range(n):
        # radius of the circle
        circle_r = m/2.0
        # center of the circle (x, y)
        circle_x = m/2.0
        circle_y = m/2.0
    
        # random angle
        alpha = 2 * math.pi * random.random()
        # random radius
        r = circle_r * math.sqrt(random.uniform(0.8,1))
        # calculating coordinates
        x = r * math.cos(alpha) + circle_x
        y = r * math.sin(alpha) + circle_y
        X.append(x)
        Y.append(y)
    #particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
    
    particles=[Particle(x,y,1) for x,y in zip(X,Y)]
    
    particlesCopy1=copy.deepcopy(particles)
    directSource(particlesCopy1) 
    fig = plt.figure(figsize=(12, 8))
    x = [point.x for point in particlesCopy1]
    y = [point.y for point in particlesCopy1]
    
    u = [point.u.real for point in particlesCopy1]
    
    qt=AdapQTree(3,particles,m)
    qt.graph()
    setAllZero(qt.root,P)
    outgoingExpansion(qt.root,P)
    incomingExpansionAdap(qt.root,P)
    summation(qt.root,P)
    u1 = [point.u.real for point in particles]
    error=0
    base=0
    for i in range(len(u)):
        error+=(u[i]-u1[i])**2
        base+=u[i]**2
        #print(x[i],y[i],u[i],u1[i])
    print(error/base)
    
    
    