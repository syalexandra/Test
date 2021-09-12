#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 22:17:17 2020

@author: yuesun
"""

from quadtree import *
#from adapQuadtree import *
from upwardPass import *
from downwardPass import *
import copy
import matplotlib.pyplot as plt
import numpy as np
import time
import random
import math

def directSource(particles):
    #particles=particles.tolist()
    for i, particle in enumerate(particles):
        for source in particles[:i]+particles[i+1:]:
            
            kernel=np.log(complex(particle.x-source.x,particle.y-source.y))
            particle.u += source.q*kernel
            
            
if __name__=='__main__':
    #draw a quad tree, n is number of points, l is the range of the random number, P is the cutoff term
    n=8
    m=80
    P=10
    np.random.seed(seed=1)
    
    particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
    particlesCopy=copy.deepcopy(particles)
    """
    qt=AdapQTree(1,particles,m)
    outgoingExpansion(qt.root,P)
    qt.root.inmultipole=np.zeros(P,dtype=complex)
    for i in range(4):
        incomingExpansion(qt.root.children[i],P)
    """
    qt=QTree(3,particles,m)
    outgoingExpansion(qt.root,P)
    incomingExpansion(qt.root,P)
    u = [point.u.real for point in particles]
    print(u)
    qt.graph_potential()
    directSource(particlesCopy) 
    fig = plt.figure(figsize=(12, 8))
    x1 = [point.x for point in particlesCopy]
    y1 = [point.y for point in particlesCopy]
    
    u1 = [point.u.real for point in particlesCopy]
    plt.scatter(x1, y1,s=u1)
    plt.title('Brute Force')
    plt.show()
    print(x1)
    print(y1)
    print(u1)
    
    """TEST 3 
    P=10
    start_time = time.clock()
    qt=QTree(3,particles,l)
    outgoingExpansion(qt.root,P)
    incomingExpansion(qt.root,P)
    
    qt.graph_potential()
    
    
    directSource(particlesCopy)
    
    
    fig = plt.figure(figsize=(12, 8))
    #plt.title("Quadtree")
    
    x1 = [point.x for point in particlesCopy]
    y1 = [point.y for point in particlesCopy]
    
    u1 = [point.u.real for point in particlesCopy]
    
    #ax = fig.add_subplot(111)
    #ax.add_patch(plt.Rectangle((0, 0), 80, 80, fill=False))
    plt.scatter(x1, y1,s=u1)
    plt.show()
    """
    
    
    """
    
    #TEST 1
    start_time = time.clock()
    directSource(particlesCopy)
    originalTime=time.clock() - start_time
    print(originalTime)
    u1 = [point.u.real for point in particlesCopy]
    
    
    
    
    for P in range(1,10):
        
        particlesCopy2=copy.deepcopy(particles)
        
        start_time = time.clock()
        qt=QTree(3,particlesCopy2,l)
        outgoingExpansion(qt.root,P)
        incomingExpansion(qt.root,P)
        FmmTimeList.append(time.clock() - start_time)
        
        #qt.graph_potential()
        
        x = [point.x for point in qt.points]
        y = [point.y for point in qt.points]
        
        u = [point.u.real for point in qt.points]
        
        
        
        #plt.scatter(x1, y1,s=1)
        #plt.show()
        errorSum=0
        data=0
        for i in range(n):
            if (u[i] != float("inf") and u[i] != float("-inf")):
                
                errorSum+=(u[i]-u1[i])**2
                data+=1
            else:
                continue
        errorList.append(np.sqrt(errorSum/data))
    
    print(errorList)
    plist=list(range(1,10))
    fig = plt.figure(figsize=(12, 8))
    plt.plot(plist, errorList)
    plt.yscale('log')
    plt.xlabel('')
    plt.show()
    
    print(FmmTimeList)
    glist=list(range(1,10))
    fig = plt.figure(figsize=(12, 8))
    #plt.plot(glist, FmmTimeList,glist,[originalTime]*9)
    plt.plot(glist, FmmTimeList)
    plt.show()
    
    
    """
    
    
    
    """#TEST 2
    P=10
    FmmTimeList=[]
    #OrginalTimeList=[]
    
    particles=[Particle(x,y,1) for x,y in np.random.randint(1,l,(n,2))]
    particlesCopy=copy.deepcopy(particles)
    
    start_time = time.clock()
    directSource(particlesCopy)
    originalTime=time.clock() - start_time
    print(originalTime)
    
    #pig = plt.figure(figsize=(12, 8))
    #plt.title("Quadtree")
    
    #x1 = [point.x for point in particlesCopy]
    #y1 = [point.y for point in particlesCopy]
    
    u1 = [point.u.real for point in particlesCopy]
    
    #plt.scatter(x1, y1,s=1)
    #plt.show()
    
    
        
    for g in range(2,8):    
        particlesCopy2=copy.deepcopy(particles)
        start_time = time.clock()
        

        qt=QTree(g,particlesCopy2,l)
        outgoingExpansion(qt.root,P)
        incomingExpansion(qt.root,P)
        FmmTimeList.append(time.clock() - start_time)
        print(time.clock() - start_time)
        
        #qt.graph_potential()
        
        #x = [point.x for point in qt.points]
        #y = [point.y for point in qt.points]
        
        u = [point.u.real for point in qt.points]
        
        
        errorSum=0
        data=0
        for i in range(n):
            if (u[i] != float("inf") and u[i] != float("-inf")):
                
                errorSum+=(u[i]-u1[i])**2
                data+=1
            else:
                continue
        
        
        errorList.append(np.sqrt(errorSum/data))
    
    print(errorList)
    glist=list(range(2,8))
    fig = plt.figure(figsize=(12, 8))
    plt.plot(glist, errorList)
    plt.show()
    
    print(FmmTimeList)
    glist=list(range(2,8))
    fig = plt.figure(figsize=(12, 8))
    plt.plot(glist, FmmTimeList,glist,[originalTime]*6)
    plt.show()
    
    
    
    X=[]
    Y=[]
    for i in range(10):
        # radius of the circle
        circle_r = 10
        # center of the circle (x, y)
        circle_x = 10
        circle_y = 10
    
        # random angle
        alpha = 2 * math.pi * random.random()
        # random radius
        r = circle_r * math.sqrt(random.uniform(0.8,1))
        # calculating coordinates
        x = r * math.cos(alpha) + circle_x
        y = r * math.sin(alpha) + circle_y
        X.append(x)
        Y.append(y)
    print(X)
    print(Y)
    """