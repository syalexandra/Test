#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 20:04:30 2020

@author: yuesun
"""

import numpy as np
from Point import *
import matplotlib.pyplot as plt

class Node:
    #define tree node
    
    def __init__(self,x,y,w,h,points,level):
        self.x=x
        self.y=y
        self.width=w
        self.height=h
        self.level=level
        self.points=points#array of points
        self.children=[]
        #children is a list of nodes
        self.outmultipole=None
        self.inmultipole=None
        self.center=Point(x+w/2.0,y+h/2.0)
        
    def isleaf(self):
        return len(self.children)==0
    
    
    def getPoints(self):
        return self.points
    
    def hasChildren(self):
        return len(self.children)!=0
        
        
        
    def addChildren(self):
        return
    
    
    def contains(self,x,y,w,h,points):
        ptList=[]
        for p in points:
            if p.x>=x and p.x<x+w and p.y>=y and p.y<y+h:
                ptList.append(p)
        return np.array(ptList)
    
    
    
    def split_adaptive(self,threshold):
        if self.hasChildren():
            return
        
        if len(self.points)>threshold:
            w=self.width/2.0
            h=self.height/2.0
            #3,5
            #2,4
            x0=self.x
            y0=self.y
            
            p=self.contains(x0,y0,w,h,self.points)
            n1=Node(x0,y0,w,h,p,self.level+1)
            n1.split(threshold)
            
            p=self.contains(x0,y0+h,w,h,self.points)
            n2=Node(x0,y0+h,w,h,p,self.level+1)
            n2.split(threshold)
            
            p=self.contains(x0+w,y0,w,h,self.points)
            n3=Node(x0+w,y0,w,h,p,self.level+1)
            n3.split(threshold)
            
            p=self.contains(x0+w,y0+h,w,h,self.points)
            n4=Node(x0+w,y0+h,w,h,p,self.level+1)
            n4.split(threshold)
            
            self.children=[n1,n2,n3,n4]
            
        
class QTree:
    def __init__(self,k,points,n):
        
        self.threshold=k
        self.points=points
        self.root=Node(0,0,n,n,self.points,0)
        
        self.buildTree(self.threshold)
        
    def buildTree(self,threshold):
        self.root.split(threshold)
        
    def getRoot(self):
        return self.root
    
    def graph(self):
        fig = plt.figure(figsize=(12, 8))
        plt.title("Quadtree")
        ax = fig.add_subplot(111)
        c = self.find_children(self.root)
        areas = set()
        for el in c:
            areas.add(el.width*el.height)
            
        for n in c:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=False))
            
        x = [point.x for point in self.points]
        y = [point.y for point in self.points]
        plt.plot(x, y, 'ro')
        plt.show()
        return
            
            
    def find_children(self,node):
        if not node.children:
           return [node]
        else:
           children = []
           for child in node.children:
               children += (self.find_children(child))
        return children
        
    
    
    
    
    
if __name__=='__main__':
    #draw a quad tree, n is number of points, l is the range of the random number
    n=500
    l=100
    particles=[Particle(x,y,1) for x,y in np.random.randint(1,l,(n,2))]
    qt=QTree(5,particles,l)
    qt.graph()