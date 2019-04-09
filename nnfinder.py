import numpy as np
def nnfinder(s,N,d):
	###################
	#function
	def future(xf,xp,d):
		''' Checking if point xf is in the future of xp '''
		R=np.sqrt(sum([(f-p)**2 for f,p in zip(xf[1:], xp[1:])]))
		if xp[0]<xf[0] and R < (xf[0]-xp[0]):
			return 1
		else:
			return 0
	##
	def condition(xf,xm,xp,d):
		''' Checking whether xm is in the midle of xp
		and xf '''
		if future(xf,xm,d)==1 and future(xm,xp,d)==1:
			return 0
		else:
			return 1
	##
	def empty (xf,xp,L,d):
		''' Checking if there is any element in list L
		that is between xf and xp '''
		count=0
		for i in range(len(L)):
			if condition(xf,L[i],xp,d)==0:
				return 0
		if count==0:
			return 1
	###################
	link = np.zeros((N+2, N+2))# Link adjacent matrix
	Relations = np.zeros((N+2, N+2))# Relations Matrix
	PastRelations = np.zeros((N+2, N+2))
	for e1 in range(len(s)):
		for e2 in range(e1, len(s)):
			if future(s[e2],s[e1],d)==1:# if s[e2] in the future of s[e1]
				Relations[e1][e2]=1
				if empty(s[e2],s[e1],s,d)==1: # if interval s[e1] to s[e2] is empty
					link[e1][e2]=1
			if future(s[e1],s[e2],d)==1:
				PastRelations[e1][e2]=1
				if empty(s[e1], s[e2], s, d) == 1:
					link[e2][e1] = 1
	return link,Relations,PastRelations
