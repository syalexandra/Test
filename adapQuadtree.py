#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 19:31:53 2020

@author: yuesun
"""

import numpy as np
from Point import *
import matplotlib.pyplot as plt
import random
import math
class AdapNode:
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
        self.neighborList=[]
        self.greaterNeighborList=[]
        #children is a list of nodes
        self.interactionList=[]
        self.l3List=[]
        self.l4List=[]
        self.outmultipole=None
        self.inmultipole=None
        self.center=Point(x+w/2.0,y+h/2.0)
        
    def isleaf(self):
        return len(self.children)==0
    
    def isNW(self):
        return self.parent.children[1]==self
    
    def isNE(self):
        return self.parent.children[3]==self

    def isSW(self):
        return self.parent.children[0]==self
    
    def isSE(self):
        return self.parent.children[2]==self
    
    def setParent(self,parentNode):
        self.parent=parentNode
        
    def getPoints(self):
        return self.points
    
    def hasChildren(self):
        return len(self.children)!=0
        
    
    def belongsTo(self,ancestor):
        node=self
        while node.parent:
            node=node.parent
            if node==ancestor:
                return True
        return False
        
    def addChildren(self):
        return
    
    
    def contains(self,x,y,w,h,points):
        ptList=[]
        for p in points:
            if p.x>=x and p.x<x+w and p.y>=y and p.y<y+h:
                ptList.append(p)
        return np.array(ptList)
    
    
    
        
    def split(self,threshold):
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
            n1=AdapNode(x0,y0,w,h,p,self.level+1)
            n1.setParent(self)
            n1.split(threshold)
            
            p=self.contains(x0,y0+h,w,h,self.points)
            n2=AdapNode(x0,y0+h,w,h,p,self.level+1)
            n2.setParent(self)
            n2.split(threshold)
            
            p=self.contains(x0+w,y0,w,h,self.points)
            n3=AdapNode(x0+w,y0,w,h,p,self.level+1)
            n3.setParent(self)
            n3.split(threshold)
            
            p=self.contains(x0+w,y0+h,w,h,self.points)
            n4=AdapNode(x0+w,y0+h,w,h,p,self.level+1)
            n4.setParent(self)
            n4.split(threshold)
            
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
    
    def getUpRightNeighbor(self):
        
        if self.parent is None:
            return None
        if self.parent.children[0]==self:
            return self.parent.children[3]
        
        
        if self.parent.children[3]==self:
            node=self.parent.getUpRightNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[0]
        
        
        
        if self.parent.children[1]==self:
            node=self.parent.getUpNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[2]
            
        if self.parent.children[2]==self:
            node=self.parent.getRightNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[1]
        
        
        
    def getDownRightNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[1]==self:
            return self.parent.children[2]
        
        
        if self.parent.children[0]==self:
            node=self.parent.getDownNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[3]
        
        
        
        if self.parent.children[2]==self:
            node=self.parent.getDownRightNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[1]
            
        if self.parent.children[3]==self:
            node=self.parent.getRightNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[0]
        
        
        
        
    def getDownLeftNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[3]==self:
            return self.parent.children[0]
        
        
        if self.parent.children[1]==self:
            node=self.parent.getLeftNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[2]
        
        
        
        if self.parent.children[0]==self:
            node=self.parent.getDownLeftNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[3]
            
        if self.parent.children[2]==self:
            node=self.parent.getDownNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[1]
        
    
    
    def getUpLeftNeighbor(self):
        if self.parent is None:
            return None
        if self.parent.children[2]==self:
            return self.parent.children[1]
        
        
        if self.parent.children[3]==self:
            node=self.parent.getUpNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[0]
        
        
        
        if self.parent.children[1]==self:
            node=self.parent.getUpLeftNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[2]
            
        if self.parent.children[0]==self:
            node=self.parent.getLeftNeighbor()
            if node is None or node.isleaf():
                return node
            else:
                return node.children[3]
    
    


       
                
                
        
        
    def getSmallerUpNeighbor(self,neighbor):
        #neighbor=self.getUpNeighbor()
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[0])
                candidates.append(current.children[2])
        return neighbors
    
    def getSmallerDownNeighbor(self,neighbor):
        #neighbor=self.getDownNeighbor()
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[1])
                candidates.append(current.children[3])
        return neighbors
    
    
    def getSmallerLeftNeighbor(self,neighbor):
        #neighbor=self.getLeftNeighbor()
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[2])
                candidates.append(current.children[3])
        return neighbors
    
    
    def getSmallerRightNeighbor(self,neighbor):
        #neighbor=self.getRightNeighbor()
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[0])
                candidates.append(current.children[1])
        return neighbors
    
    
    def getSmallerUpRightNeighbor(self,neighbor):
        
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[0])
        return neighbors
    
    
    def getSmallerUpLeftNeighbor(self,neighbor):
        
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[2])
        return neighbors
    
    
    def getSmallerDownLeftNeighbor(self,neighbor):
        
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[3])
        return neighbors
    
    
    def getSmallerDownRightNeighbor(self,neighbor):
        
        candidates=[] if neighbor is None else [neighbor]
        neighbors=[]
        while len(candidates)>0:
            current=candidates.pop(0)
            if current.isleaf():
                neighbors.append(current)
            else:
                candidates.append(current.children[1])
        return neighbors
    
    
    def getGreaterNeighbors(self):
        
        neighborList=[]
        upNei=self.getUpNeighbor()
        downNei=self.getDownNeighbor()
        rightNei=self.getRightNeighbor()
        leftNei=self.getLeftNeighbor()
        upLeftNei=self.getUpLeftNeighbor()
        upRightNei=self.getUpRightNeighbor()
        downLeftNei=self.getDownLeftNeighbor()
        downRightNei=self.getDownRightNeighbor()
        
        if upNei:
            neighborList.append(upNei)
        if downNei:
            neighborList.append(downNei)
        if rightNei:
            neighborList.append(rightNei)
        if leftNei:
            neighborList.append(leftNei)
            
        if upLeftNei:
            neighborList.append(upLeftNei)
        if upRightNei:
            neighborList.append(upRightNei)
        if downLeftNei:
            neighborList.append(downLeftNei)
        if downRightNei:
            neighborList.append(downRightNei)
            
        neighborList = list(dict.fromkeys(neighborList))
        self.greaterNeighborList =neighborList
        return neighborList
    
    def getSmallerNeighbors(self):
        
        neighborList=[]
        upNei=self.getUpNeighbor()
        downNei=self.getDownNeighbor()
        rightNei=self.getRightNeighbor()
        leftNei=self.getLeftNeighbor()
        
        if upNei:
            neighborList+=self.getSmallerUpNeighbor(upNei)
        if downNei:
            neighborList+=self.getSmallerDownNeighbor(downNei)
        if rightNei:
            neighborList+=self.getSmallerRightNeighbor(rightNei)
        if leftNei:
            neighborList+=self.getSmallerLeftNeighbor(leftNei)
          
        
        upright=self.getUpRightNeighbor()
        upleft=self.getUpLeftNeighbor()
        downright=self.getDownRightNeighbor()
        downleft=self.getDownLeftNeighbor()    
        
        if upright:
            neighborList+=self.getSmallerUpRightNeighbor(upright)
            
        if upleft:
            neighborList+=self.getSmallerUpLeftNeighbor(upleft)
            
        if downright:# 
            neighborList+=self.getSmallerDownRightNeighbor(downright)
                
        if downleft:
            neighborList+=self.getSmallerDownLeftNeighbor(downleft)
        neighborList = list(dict.fromkeys(neighborList))
        
        self.neighborList =neighborList
        return neighborList
    
    
    
    
    def getInteractionList(self):
        
        nn=self.greaterNeighborList
        """    
        if not self.parent.neighbors:
            pn=self.parent.getParentNeighbors()
        else:
            pn=self.parent.neighbors
        """
        pn=self.parent.greaterNeighborList
        
        interactionList=[]
        
        
        for n in pn:
            queue=[n]
            while(len(queue)>0):
                node=queue.pop(0)
                if node.level==self.level and node not in nn:
                    interactionList.append(node)    
                
                #if node.isleaf() and node.level<self.level and node not in nn and 
                for tnode in node.children:
            
                    queue.append(tnode)
        
        self.interactionList=interactionList
        return interactionList 
    
    
    def getL3List(self):
        if not self.isleaf():
            return []
        
        pn=self.greaterNeighborList
        
        L3=[]
        
        
        for n in pn:
            queue=[n]
            while(len(queue)>0):
                node=queue.pop(0)
                """
                if (node.level>self.level) and (self in node.parent.neighborList) and (self not in node.neighborList):
                    L3.append(node)   
                    node.l4List.append(self)
                """
                if node.isleaf() and (node.level>self.level) and (node not in self.neighborList):# and (node.parent in self.parent.neighborList):
                    L3.append(node)   
                    node.l4List.append(self)
                
                #if node.isleaf() and node.level<self.level and node not in nn and 
                for tnode in node.children:
            
                    queue.append(tnode)
        self.l3List=L3
        return L3 
    
    
    
    
class AdapQTree:
    def __init__(self,k,points,n):
        
        self.threshold=k
        self.points=points
        self.root=AdapNode(0,0,n,n,self.points,0)
        
        self.buildTree(self.threshold)
        queue=[]
        queue.append(self.root)
        
        while len(queue)>0:
            node=queue.pop(0)
            if node.level!=0:
                
                node.getSmallerNeighbors()
                node.getGreaterNeighbors()
                
                node.getInteractionList()
                node.getL3List()
                
            for tnode in node.children:
                
                queue.append(tnode)
        
    def buildTree(self,threshold):
        self.root.split(threshold)
        
    #def findAllList(self):
        
        
        
        
        
    def getRoot(self):
        return self.root
    
    def graph(self):
        fig = plt.figure(figsize=(12, 12))
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
        fig = plt.figure(figsize=(12, 12))
        #plt.title("Quadtree")
        ax = fig.add_subplot(111)
        c = self.find_children(self.root)
        areas = set()
        for el in c:
            areas.add(el.width*el.height)
            
        for n in c:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=False))
           
        for n in nodes:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(1,0,0,0.5), ec=(0,0,0,1)))
        
        x = [point.x for point in self.points]
        y = [point.y for point in self.points]
        plt.plot(x, y, 'ro')
        plt.show()
        return
    def graph_neighbors_interactions(self,l1,l2,l3,l4,l5):
        fig = plt.figure(figsize=(12, 12))
        #plt.title("Quadtree")
        ax = fig.add_subplot(111)
        c = self.find_children(self.root)
        areas = set()
        for el in c:
            areas.add(el.width*el.height)
            
        for n in c:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height, fill=False))
           
        for n in l1:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(0.5,0.5,0,0.5), ec=(0,0,0,1)))
        
        for n in l2:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(0.5,0,0.5,0.5), ec=(0,0,0,1)))
            
        for n in l3:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(0,0.5,0.5,0.5), ec=(0,0,0,1)))
            
        for n in l4:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(0.5,1,0.5,0.5),ec=(0,0,0,1)))
            
        for n in l5:
            ax.add_patch(plt.Rectangle((n.x, n.y), n.width, n.height,  fc=(0.5,0.5,1,0.5), ec=(0,0,0,1)))
        
        
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
    n=50
    m=80
    P=10
    
    random.seed(1)
    particles=[Particle(random.uniform(0,m),random.uniform(0,m),1) for i in range(n)]
    
    qt=AdapQTree(2,particles,m)
    qt.graph()
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
    """
    node=qt.root.children[2].children[1].children[0]
    nnList=node.getNeighbors()
    node.parent.getNeighbors()
    print(node.x,node.y)
    #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
    qt.graph_neighbors(nnList)   
    nnList=node.getInteractionList()
    
    #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
    
    qt.graph_neighbors(nnList)
    """
    
    #node=qt.root.children[2].children[1]
    
    queue=[]
    queue.append(qt.root)
    while len(queue)>0:
        node=queue.pop(0)
        if node.level!=0:
        #if node.x==7.5 and node.y==15 and node.width==7.5 and node.height==7.5:
            if node.isleaf():
                nn=node.neighborList
                l3=node.l3List
                l4=node.l4List
            else:
                nn=[]
                l3=[]
                l4=[]
            inter=node.interactionList
            print(node.x,node.y)
            qt.graph_neighbors_interactions([node],nn,inter,l3,l4)
            #qt.graph_neighbors(nnList2)
            #qt.graph_neighbors(nnList)
            #qt.graph_neighbors(nnList)
            """
            if nnList1 and nnList:
                qt.graph_neighbors_interactions(nnList,nnList1) 
            else:
                qt.graph_neighbors([node])
            """
            #nnList=node.getInteractionList()
        for tnode in node.children:
            
            queue.append(tnode)
        
    
      
         
    
    #print([str(nn.x)+' '+str(nn.y) for nn in nnList])
    
    #qt.graph_neighbors(nnList)
    