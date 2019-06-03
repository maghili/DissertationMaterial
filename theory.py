import numpy as np
import math as m
import mpmath as mpm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(i,length):
    l=length;a=m.gamma(i+1)
    while l>1:
        x=[m.gamma(1+i-j)*f(j,length-1) for j in range(i+1)]
        a=sum(x)
        l=l-1
    return a

def alt_f(n):
    func=[[0]*(n+1) for i in range(n+1)]
    #func[0]=[mpm.gamma(i+1) for i in range(n+1)]
    func[0]=[1 for i in range(n+1)]
    for k in range(1,n+1):
        for i in range(n+1):
            func[k][i]=sum([mpm.gamma(j+1)*mpm.gamma(1-j+i)*func[k-1][j]/mpm.gamma(i+1) for j in range(i+1)])
    print('table created')
    return func


def fac(a):
    z=1;c=a
    while c!=0:
        z=a*fac(a-1)
        c=c-1
    return z


def nk(n,A):
    n_k=[]
    for length in range(1,n+1):
        #x=[n**(length-1)*(-1)**i*mpm.gamma(n+1)*A[length-1][i]*mpm.gamma(i+1)/(mpm.gamma(i+length)**2*mpm.gamma(1-i+n)) for i in range(n+1)]
        x=[mpm.gamma(n+1.)*(-1)**i*A[length-1][i]*mpm.gamma(i+1)/(mpm.gamma(i+length)**2*mpm.gamma(2-i-length+n)) 
        for i in range(n-length+2)]
        n_k.append(sum(x))
    knk=sum(k*n_k[k-1] for k in range(1,n+1))
    snk=sum(n_k)
    kavg=knk*1./snk
    a=kavg/n**(1.0/2)
    return n_k, kavg, a

N=[10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

A=alt_f(max(N))

y=[];bar=[];x=[];alpha=[]
for num in N:
    dist,kbar, a=nk(num,A)
    y.append(dist)
    bar.append(kbar)
    x.append(range(1,num+1))
    alpha.append(float(a))
#expol = lambda x,a,b:a/x+b
#N=N+[100];alpha=alpha+[1.225]
def expol(x,a,b,c):
    return [a/j**c+b for j in x]

best_val, cov=curve_fit(expol, N, alpha, [1, 1, 1])

print(len(alpha))
#print best_val

plt.plot(N,alpha)
plt.plot(N,expol(N,best_val[0],best_val[1], best_val[2]))
plt.show()

np.save('coef',best_val)
np.save('talpha',alpha)
