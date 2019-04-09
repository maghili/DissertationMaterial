import numpy as np
import itertools
def pathfinder(link,Relations,PastRelations,N,d):
	''' Finding the path length distribution in a
	causal set, between minimal and maximal points
	and some extra calculations for various methods
	to calculate the dimension '''
	######################
	#MidPoint
	if N>=100:
		minvols = [0]*(N+2);#futures=[0]*(N+2);pasts=[0]*(N+2)
		maxvols = [0]*(N+2)
		for i in range(N+2):
			minvols[i] = min([sum(Relations[i]),sum(PastRelations[i])])
			maxvols[i] = max([sum(Relations[i]), sum(PastRelations[i])])
		HalfminVol = max(minvols)
		HalfmaxVol = maxvols[minvols.index(HalfminVol)]
		MidPoint1 = np.log(N/(HalfminVol-1.))/np.log(2)
		MidPoint2 = np.log(N/(HalfmaxVol-1.))/np.log(2)
	else:
		MidPoint1 = 0
		MidPoint2 = 0
	######################
	#Myrheim-Meyer
	mR=np.matrix(Relations)
	R3=mR**3
	c2=R3[0,-1]
	MyrheimMeyer=c2*1.0/(N*N-1)
	######################
	#Brightwell-Gregory and ABP
	mlinks=np.matrix(link)
	zero=[[0]*len(link) for j in range(len(link))]
	A=np.identity(len(mlinks))
	pathlength=[]
	m=1
	tot=[]
	weights=[0]*(N+2)
	mnk=[]
	while m<=len(link):
		if (A*mlinks).tolist()!=zero:
			A=A*mlinks
			pathlength.append(m)
			tot.append(A[0,-1])
			mnk.append(m*A[0,-1])
			weights=[weights[i]+A[i,-1] for i in range(N+2)]
			m=m+1
		else:
			m=len(link)+1
	totpath=sum(tot)
	Longest=pathlength[-1]
	beta=Longest/(N)**(1.0/d)
	avglen=sum(mnk)*1.0/sum(tot)
	alpha=avglen/(N)**(1.0/d)
	##################################
	return weights, MidPoint1, MidPoint2, MyrheimMeyer, beta, Longest, avglen, alpha, tot
