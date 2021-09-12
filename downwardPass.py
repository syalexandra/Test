#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 21:28:22 2020

@author: yuesun
"""
import numpy as np
from scipy.special import binom

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

"""
def direct_sources_to_targets(targets,sources):
    m=len(targets)
    v=np.zeros(m,dtype=complex)
    for i,target in enumerate(targets):
        for source in sources:
            kernel=np.log(complex(target.x-source.x,target.y-source.y))
            v[i]+=source.q*kernel
    return v
"""

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
                 
"""              
def incomingExpansion(tnode,P=10):
    # using bfs to loop over levels
    
    z0=tnode.parent.center-tnode.center    
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
    else:
        for child in tnode.children:
            incomingExpansion(child,P)
"""
                    
            
        


