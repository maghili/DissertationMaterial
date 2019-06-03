import numpy as np
from scipy import sparse
import itertools
def pathfinder(link,Relations,PastRelations,N,d):
	''' Finding the path length distribution in a
	causal set, between minimal and maximal points
	and some extra calculations for various methods
	to calculate the dimension '''
	######################
	#MidPoint
	if N>=100:
		minvols = [ min(sum(Relations[i]),sum(PastRelations[i])) for i in range(N+2) ]
		maxvols = [ max(sum(Relations[i]), sum(PastRelations[i])) for i in range(N+2) ]
		HalfminVol = max(minvols)
		HalfmaxVol = maxvols[np.argmax(minvols)]
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
	mlinks=sparse.csr_matrix( np.matrix(link) )
	zero=sparse.csr_matrix( np.zeros((N+2, N+2)) )
	A=sparse.csr_matrix( np.identity(N+2) )
	pathlength=[]
	m=1
	tot=[]
	weights=[0]*(N+2)
	mnk=[]
	while m <= N+2:
		if ((A*mlinks)!=zero).nnz:
			A=A*mlinks
			pathlength.append(m)
			tot.append(A[0,-1])
			mnk.append(m*A[0,-1])
			weights=[weights[i]+A[i,-1] for i in range(N+2)]
			m = m+1
		else:
			m = N+3
	totpath=sum(tot)
	Longest=pathlength[-1]
	beta=Longest/(N)**(1.0/d)
	avglen=sum(mnk)*1.0/sum(tot)
	alpha=avglen/(N)**(1.0/d)
	##################################
	return weights, MidPoint1, MidPoint2, MyrheimMeyer, beta, Longest, avglen, alpha, tot
