#def funtion_pool():

	
def dot(xf,xp,d):
	'''The Lorentzian dot product in a flat spacetime'''
	eta = [-1]+[1 for dim in range(d-1)]
	return sum(m*(x-y)**2 for m, x, y in zip(eta, xf, xp))

def future(xf, xp, d):
	'''we are trying to check whether the first 
	argument is in the future of the second argument
	for a d-dimensional spacetime'''
	if dot(xf, xp, d) < 0 and xf[0] > xp[0]:
		return 1
	return 0

def past(xf, xp, d):
	'''we are trying to check whether the second 
	argument is in the past of the first argument
	for a d-dimensional spacetime'''
	if dot(xf, xp, d) < 0 and xf[0] > xp[0]:
		return 1
	return 0

def condition(xf,xm,xp,d):
	'''This condition checks whether xm is in the
	intersection of pas lightcone of xf and future
	lightcone of xp'''
	if future(xf,xm,d) and future(xm,xp,d):
		return 0
	return 1

def empty (xf,xp,L,d):
	'''This function checks whether the timelike
	interval between xf and xp is empty'''
	count=0
	for point in L:
		if not condition(xf,point,xp,d):
			return 0
	return 1



class MinkowskiSprinkling():

	def __init__(self, xmin = [0,0], xmax = [1,0], d=2):
		self.dimension = d
		self.xmin = xmin
		self.xmax = xmax
