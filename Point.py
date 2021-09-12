# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

class Point():
    
    def __init__(self,x,y):
        self.x=x
        self.y=y
        
        
    def __sub__(self,other):
        return Point(self.x-other.x,self.y-other.y)

class Particle(Point):
    
    def __init__(self,x,y,charge):
        
        super().__init__(x,y)
        self.q=charge
        self.u=0
        
        
        
def distance(p1,p2):
    #distance between two points
    
    return np.sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2)

def isequal(p1,p2):
    return p1.x==p2.x and p1.y==p2.y
    
    
def direct_sum(particles):
    for i,target in enumerate(particles):
        for source in (particles[:i]+particles[i+1:]):
            if isequal(target,source):
                kernel=0
            else:
                kernel=np.log(distance(target,source))
                
            target.u+=source.q*kernel
            
    

            
if __name__=='__main__':
    
    p1=Particle(1,2,1)
    p2=Particle(3,5,1)
    print((p1-p2).y)
    """
    print(distance(p1,p2))
    pp1=Particle(1,2,3)
    pp2=Particle(4,4,5)
    print(distance(pp1,pp2))
    
    n=50
    particles=[Particle(x,y,1) for x,y in np.random.randint(1,100,(n,2))]
    direct_sum(particles)
    x=[]
    y=[]
    u=[]
    for i in particles:
        x.append(i.x)
        y.append(i.y)
        u.append(i.u)
        
    
    plt.scatter(x,y,s=u)
    plt.show()
    """