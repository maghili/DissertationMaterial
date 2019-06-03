import numpy as np
from pool import empty, future
def nnfinder(s,N,d):

	link = np.zeros((N+2, N+2))# Link adjacent matrix
	Relations = np.zeros((N+2, N+2))# Relations Matrix
	PastRelations = np.zeros((N+2, N+2))

	for e1 in range(len(s)):
		for e2 in range(e1, len(s)):

			if future(s[e2],s[e1],d):# if s[e2] in the future of s[e1]
				Relations[e1][e2]=1
				if empty(s[e2],s[e1],s,d): # if interval s[e1] to s[e2] is empty
					link[e1][e2]=1
			if future(s[e1],s[e2],d):
				PastRelations[e1][e2]=1
				if empty(s[e1], s[e2], s, d):
					link[e2][e1] = 1

	return link,Relations, PastRelations
