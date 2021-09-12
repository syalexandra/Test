#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 16:03:52 2020

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
        self.parent=None
        self.neighbors=None
        #children is a list of nodes
        self.outmultipole=None
        self.inmultipole=None
        self.center=Point(x+w/2.0,y+h/2.0)
        
    def isleaf(self):
        return len(self.children)==0
    
    def setParent(self,parentNode):
        self.parent=parentNode
        
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
    
    
    
        
    def split(self,depth):
        if self.hasChildren():
            return
        
        if self.level<depth:
            w=self.width/2.0
            h=self.height/2.0
            #3,5
            #2,4
            x0=self.x
            y0=self.y
            
            p=self.contains(x0,y0,w,h,self.points)
            n1=Node(x0,y0,w,h,p,self.level+1)
            n1.setParent(self)
            n1.split(depth)
            
            p=self.contains(x0,y0+h,w,h,self.points)
            n2=Node(x0,y0+h,w,h,p,self.level+1)
            n2.setParent(self)
            n2.split(depth)
            
            p=self.contains(x0+w,y0,w,h,self.points)
            n3=Node(x0+w,y0,w,h,p,self.level+1)
            n3.setParent(self)
            n3.split(depth)
            
            p=self.contains(x0+w,y0+h,w,h,self.points)
            n4=Node(x0+w,y0+h,w,h,p,self.level+1)
            n4.setParent(self)
            n4.split(depth)
            
            self.children=[n1,n2,n3,n4]
            
    def getUpNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[0]==self:
            return self.parent.children[1]
        if self.parent.children[2]==self:
            return self.parent.children[3]
        node=self.parent.getUpNeighbor()
        if node is None or node.isleaf():
            return node
        if self.parent.children[1]==self:
            return node.children[0]
        else:
            return node.children[2]
        
    def getDownNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[1]==self:
            return self.parent.children[0]
        if self.parent.children[3]==self:
            return self.parent.children[2]
        node=self.parent.getDownNeighbor()
        if node is None or node.isleaf():
            return node
        if self.parent.children[0]==self:
            return node.children[1]
        else:
            return node.children[3]
        
    def getLeftNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[3]==self:
            return self.parent.children[1]
        if self.parent.children[2]==self:
            return self.parent.children[0]
        node=self.parent.getLeftNeighbor()
        if node is None or node.isleaf():
            return node
        if self.parent.children[1]==self:
            return node.children[3]
        else:
            return node.children[2]
        
        
    def getRightNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[1]==self:
            return self.parent.children[3]
        if self.parent.children[0]==self:
            return self.parent.children[2]
        node=self.parent.getRightNeighbor()
        if node is None or node.isleaf():
            return node
        if self.parent.children[3]==self:
            return node.children[1]
        else:
            return node.children[0]
    
    
    def getNeighbors(self):
        if self.neighbors:
            return self.neighbors
        neighborList=[]
        upNei=self.getUpNeighbor()
        downNei=self.getDownNeighbor()
        rightNei=self.getRightNeighbor()
        leftNei=self.getLeftNeighbor()
        
        if upNei:
            neighborList.append(upNei)
        if downNei:
            neighborList.append(downNei)
        if rightNei:
            neighborList.append(rightNei)
        if leftNei:
            neighborList.append(leftNei)
            
        if upNei and rightNei:
            neighborList.append(upNei.getRightNeighbor())
        if upNei and leftNei:
            neighborList.append(upNei.getLeftNeighbor())
        if downNei and rightNei:
            neighborList.append(downNei.getRightNeighbor())
        if downNei and leftNei:
            neighborList.append(downNei.getLeftNeighbor())
            
        self.neighbors =neighborList
        return neighborList
    

    def getInteractionList(self):
        if not self.neighbors:
            nn=self.getNeighbors()
        else:
            nn=self.neighbors
        if not self.parent.neighbors:
            pn=self.parent.getNeighbors()
        else:
            pn=self.parent.neighbors
            
        interactionList=[]
        for n in pn:
            if n.hasChildren():
                interactionList+=[c for c in n.children if c not in nn]
            elif n not in nn:
                interactionList.append(n)
        return interactionList
        
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
        #plt.title("Quadtree")
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
    
    def graph_neighbors(self,nodes):
        fig = plt.figure(figsize=(12, 8))
        #plt.title("Quadtree")
        ax = fig.add_subplot(111)
        c = self.find_children(self.root)
        areas = set()
        for el in c:
            areas.add(el.width*el.height)
            
        for n in c:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=False))
           
        for n in nodes:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=True))
        
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
        
    
    
    def graph_potential(self):
        fig = plt.figure(figsize=(12, 8))
        #plt.title("Quadtree")
        ax = fig.add_subplot(111)
        c = self.find_children(self.root)
        areas = set()
        for el in c:
            areas.add(el.width*el.height)
            
        for n in c:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=False))
            
        
        x = [point.x for point in self.points]
        y = [point.y for point in self.points]
        u = [point.u.real for point in self.points]
        plt.scatter(x, y,s=u)
        plt.title('FMM')
        plt.show()
        return
    
    
    
if __name__=='__main__':
    #draw a quad tree, n is number of points, l is the range of the random number
    n=0
    l=80
    particles=[Particle(x,y,1) for x,y in np.random.randint(1,l,(n,2))]
    qt=QTree(3,particles,l)
    #qt.graph()
    """
    for i in range(4):
        for j in range(4):
            for k in range(4):
                node=qt.root.children[i].children[j].children[k]
                nnList=node.getNeighbors()
                #print(node.x,node.y)
                #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
                qt.graph_neighbors(nnList)           
    """
    node=qt.root.children[2].children[0].children[0]
    nnList=node.getNeighbors()
    node.parent.getNeighbors()
    print(node.x,node.y)
    #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
    qt.graph_neighbors(nnList)   
    nnList=node.getInteractionList()
    
    #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
    
    qt.graph_neighbors(nnList)