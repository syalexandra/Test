#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 20:23:53 2020

@author: yuesun
"""

from Point import *
import numpy as np
from matplotlib import pyplot
from scipy.special import binom
from quadtree import *
import time
import copy
import random
import math

def directSource(particles):
    #particles=particles.tolist()
    for i, particle in enumerate(particles):
        for source in particles[:i]+particles[i+1:]:
            
            kernel=np.log(complex(particle.x-source.x,particle.y-source.y))
            particle.u += source.q*kernel

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
            
                    
def incomingExpansion(root,P=10):
    # using bfs to loop over levels
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
            
            if tnode.level==1:
                tnode.inmultipole=np.zeros(P, dtype=complex)
                
            else:
                z0=tnode.parent.center-tnode.center
                #T^{ifi}*(Uv)
                tnode.inmultipole=shiftInnerExpansion(tnode.parent.inmultipole,z0)
                
                for tin in tnode.getInteractionList():
                    z0=tin.center-tnode.center
                    #T^{ifo}*q
                    tnode.inmultipole+=inFromOut(tin.outmultipole,z0)
                    
                  
                if tnode.isleaf():
                    
                    coeffs=tnode.inmultipole
                    z0=tnode.center
                    for p in tnode.getPoints():
                        z=p-z0
                        p.u+=np.polyval(coeffs[::-1],complex(z.x,z.y))
                    
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

                    
def evaluateFarField(root,P=10):
    queue=[]
    queue.append(root)
    while len(queue)>0:
        node=queue.pop(0)
        for tnode in node.children:
            
            queue.append(tnode)
            
            tnode.inmultipole=np.zeros(P, dtype=complex)
            
            for tin in tnode.getInteractionList():
                """
                z0=tnode.center-tin.center
                IFO=inFromOut(tnode.outmultipole,z0)
                #coeffs=tnode.inmultipole
                z0=tin.center
                for p in tin.getPoints():
                    z=p-z0
                    p.u+=np.polyval(IFO[::-1],complex(z.x,z.y))
                """
                targetFromOut(tin,tnode)
                
            if tnode.isleaf():
                
                for nn in tnode.getNeighbors():
                    directSourceToTargets(tnode.getPoints(),nn.getPoints())
                
                directSourceInOneLeaf(tnode.getPoints())    
                
    
if __name__=='__main__':
    #draw a quad tree, n is number of points, l is the range of the random number, P is the cutoff term
    
    n=100
    m=80
    
    
    np.random.seed(seed=1)
    particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
    
    
    particlesCopy1=copy.deepcopy(particles)
    directSource(particlesCopy1) 
    fig = plt.figure(figsize=(12, 8))
    x1 = [point.x for point in particlesCopy1]
    y1 = [point.y for point in particlesCopy1]
    
    u1 = [point.u.real for point in particlesCopy1]
    
    BHtime=[]
    FMMtime=[]
    BHerror=[]
    FMMerror=[]
    
    
    
    for P in range(1,10):
        particlesCopy=copy.deepcopy(particles)
        particlesCopy2=copy.deepcopy(particles)
        
        start=time.clock()
        qt=QTree(3,particlesCopy,m)
        outgoingExpansion(qt.root,P)
        
        evaluateFarField(qt.root,P)
        BHtime.append(time.clock()-start)
        
        u = [point.u.real for point in particlesCopy]
        errorSum=0
        baseSum=0
        for i in range(n):
            if (u[i] != float("inf") and u[i] != float("-inf")):

                errorSum+=(u[i]-u1[i])**2
                baseSum+=u1[i]**2
            else:
                continue
        BHerror.append(np.sqrt(errorSum/baseSum))
        
        #print(u)
        #qt.graph_potential()
        
        #plt.scatter(x1, y1,s=u1)
        #plt.title('Brute Force')
        #plt.show()
        #print(u1) 
        start=time.clock()
        qt=QTree(3,particlesCopy2,m)
        outgoingExpansion(qt.root,P)
        incomingExpansion(qt.root,P)
        FMMtime.append(time.clock()-start)
        u2 = [point.u.real for point in particlesCopy2]
        #print(u)
        #print(u2)
        
        errorSum=0
        baseSum=0
        for i in range(n):
            if (u2[i] != float("inf") and u2[i] != float("-inf")):

                errorSum+=(u2[i]-u1[i])**2
                baseSum+=u1[i]**2
            else:
                continue
        FMMerror.append(np.sqrt(errorSum/baseSum))
        
        
        
        
    #print(FMMerror)
    plist=list(range(1,10))
    fig = plt.figure(figsize=(12, 12))
    plt.plot(plist, FMMerror,label='FMM')
    plt.plot(plist, BHerror,label='BH')
    #plt.yscale('log')
    plt.xlabel('cut off term P')
    plt.ylabel('error')
    plt.title('error vs cut off term')
    plt.legend(loc="upper left")
    plt.show()

    print(FMMtime)
    fig = plt.figure(figsize=(12, 12))
    #plt.plot(glist, FmmTimeList,glist,[originalTime]*9)
    plt.plot(plist, FMMtime,label='FMM')
    plt.plot(plist, BHtime,label='BH')
    plt.yscale('log')
    plt.xlabel('cut off term P')
    plt.ylabel('running time')
    plt.title('running time vs cut off term')
    plt.legend(loc="upper left")
    plt.show()
    
    """
    P=10
    FmmTimeList=[]
    BHTimeList=[]
    errorList=[]
    m=800
    
    
    #pig = plt.figure(figsize=(12, 8))
    #plt.title("Quadtree")
    
    #x1 = [point.x for point in particlesCopy]
    #y1 = [point.y for point in particlesCopy]
    
    
    
    #plt.scatter(x1, y1,s=1)
    #plt.show()
    
    
        
    for g in range(2,8):
        n=(2**(2*g))*10
        print(n)
        #m=2**g
        #particles=[Particle(x,y,1) for x,y in np.random.randint(1,m,(n,2))]
        particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
        particlesCopy=copy.deepcopy(particles)
        
        start_time = time.clock()
        qt=QTree(g,particlesCopy,m)
        outgoingExpansion(qt.root,P)
        
        evaluateFarField(qt.root,P)
        
        
        BHTimeList.append(time.clock() - start_time)
        
        u1 = [point.u.real for point in particlesCopy]
        #print(u1)
        
        start_time = time.clock()
        qt=QTree(g,particles,m)
        outgoingExpansion(qt.root,P)
        incomingExpansion(qt.root,P)
        FmmTimeList.append(time.clock() - start_time)
        #print(time.clock() - start_time)
        
        
        #qt.graph_potential()
        
        #x = [point.x for point in qt.points]
        #y = [point.y for point in qt.points]
        
    
    #print(errorList)
    
    #glist=list(range(2,6))
    #fig = plt.figure(figsize=(12, 8))
    #plt.plot(glist, errorList)
    #plt.yscale('log')
    #plt.xlabel('level L')
    #plt.ylabel('error')
    #plt.show()
   
    #print(FmmTimeList)
    glist=list(range(2,8))
    
    fig = plt.figure(figsize=(12, 12))
    plt.plot(glist, FmmTimeList,label='FMM')
    plt.plot(glist,BHTimeList,label='BH')
    sizen=[(2**(2*g))/(2**(2*2)) for g in glist]
    nsquare=[BHTimeList[0]*s*s for s in sizen]
    nlogn=[BHTimeList[0]*s*np.log(s) for s in sizen]
    n=[BHTimeList[0]*s for s in sizen]
    #nlist=[(2**(2*g))*10 for g in glist]
    plt.plot(glist,nlogn,'o',label='nlogn')
    plt.plot(glist,n,'o',label='n')
    
    plt.legend(loc="upper left")
    plt.xlabel('level L')
    plt.ylabel('running time')
    plt.yscale('log')
    #plt.xscale('log')
    plt.show()
    
    """