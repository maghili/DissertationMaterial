import numpy as np
import random as r
from math import *
def randomfinder(d,N,nr,link,weight):
	''' A function that generated a random path from the
	minimal to the maximal element '''
	##################
	# converting to list
	fnn=[[i] for i in range(N+2)] # Finding future nearest neighbors
	for e1 in range(N+2):
		for e2 in range(N+2):
			if link[e1][e2]==1:
				fnn[e1].append(e2)
	totpath=weight[0]
	###################
	ls=[] # A list of lengths of chosen paths
	for k in range(nr):
		x=fnn[0] # List of nearest neighbors of current point
		lcounter=0 # A number to count number of links
		while len(x)>1:
			plist=[0]*len(x)
			z=0
			for m in range(1,len(x)):
				if weight[x[m]]!=0:
					plist[m]=z+weight[x[m]] # An increasing list, starting at 1 and end at (sum weights)
				else:
					plist[m]=z+1
				z=plist[m]
			if len(plist)==0: # If there is no points to the future
				plist=[1]
			a=r.randint(1,max(plist)) # Randomly choosing a number in range(1,sum(weights))
			for l in range(len(plist)-1):
				if a>plist[l] and a<=plist[l+1]: # Seeing which point was chosen
					x=fnn[x[l+1]] # Replacing old point with new point
				else: #if there is only one
					x=fnn[x[0]]
			lcounter=lcounter+1
		ls.append(lcounter)
	avglen=sum(ls)/float(len(ls))
	alpha=avglen*1.0/N**(1.0/d)
	return alpha,avglen
