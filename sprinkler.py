import numpy as np
from pool import future
def sprinkler(N,d):
	''' Sprinkling the points in the manifold '''
	boxdim=[1]*d
	xmin=[0]*d
	xmax=[1]+[0]*(d-1)
	###################
	points=[xmin]
	while len(points)<N+1:
		point=[np.random.uniform(0,1)]+[np.random.uniform(-1.0/2,1.0/2) for dim in range(1,d)]#a list of length d of random numbers
		if future(xmax,point,d)==1 and future(point,xmin,d)==1:# if it is in the required interval
			points.append(point)
	points.append(xmax)
	return points
###################
