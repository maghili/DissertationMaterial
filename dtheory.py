import numpy as np
import math as m
import mpmath as mpm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sprinkler as sp
import nnfinder as nf
import Calculations as cal
from matplotlib import style

def f(i,length):
    ''' f_{i,k} coefficient in Theoretical Results'''
    l=length;a=m.gamma(i+1)
    while l>1:
        x=[m.gamma(1+i-j)*f(j,length-1) for j in range(i+1)]
        a=sum(x)
        l=l-1
    return a

def alt_f(n,d):
    ''' f_{i,k} coefficient in Theoretical Results'''
    func=[[0]*(n+1) for i in range(n+1)]
    #func[0]=[mpm.gamma(i+1) for i in range(n+1)]
    func[0]=[mpm.gamma((i+1.)*d/2)*mpm.gamma(1.+i*d/2.)/mpm.gamma(i+1.) for i in range(n+1)]
    for k in range(1,n+1):
        for i in range(n-k+1):
            func[k][i]=sum([mpm.gamma((i-j+1.)*d/2.0)*mpm.gamma(1.0+(i-j)*d/2.0)/mpm.gamma(1-j+i)*func[k-1][j] for j in range(i+1)])
    print 'table created'
    return func


def nk(n,A,d):
    ''' path length distribution '''
    n_k=[]
    for length in range(1,n+1):
        x=[(mpm.gamma(d+1)/(2.*mpm.gamma(d/2.)))**(length-1)*(-1)**i*mpm.gamma(n+1)*A[d-2][length-1][i]/
        (mpm.gamma((i+length)*d/2.0)*mpm.gamma(1+(length-1.0+i)*d/2.0)*mpm.gamma(2-i+n-length)) for i in range(n-length+2)]
        n_k.append(sum(x))
    knk=sum([k*n_k[k-1] for k in range(1,n+1)])
    snk=sum(n_k)
    kavg=knk*1./snk
    a=kavg/n**(1.0/d)
    return n_k, kavg, a

N=[10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
Nmax=50
D=5; nsp=1
Nm=100
A=[alt_f(N[-1],d) for d in range(2,D+1)]
dist=[]
for spr in range(nsp):
    s = sp.sprinkler(Nm,2)
    link, R, PR = nf.nnfinder(s,Nm,2)
    weights,MidPoint, MidPoint2,MyrheimMeyer,beta,Longest,avglen,alpha,tot = cal.pathfinder(link, R, PR, Nm, 2)
    dist.append(tot)
print sum(tot)

lens = [len(a) for a in dist]
for a in dist:
    if len(a)<max(lens):
        while len(a)<max(lens):
            a.append(0)

avg=[]
for i in range(max(lens)):
    avg.append(np.average([a[i] for a in dist]))

y=[];bar=[];x=[];Alphat=[]
for d in range(2,D+1):
    alpha=[]
    for num in range(10,Nmax+1):
        dist,kbar, a=nk(num,A,d)
        y.append(dist)
        bar.append(kbar)
        x.append(range(1,num+1))
        alpha.append(float(a))
    Alphat.append(alpha)

print len(y)
def expol(x,a,b,c):
    return [a*1./j**c+b for j in x]

markers = ['o', 's', 'd', '^', '*']
style.use('bmh')
plt.figure()
for d in range(2,D+1):
    plt.plot(range(10,Nmax+1),Alphat[d-2], '--', marker = markers[d-2], mec = 'k', lw = 1, label='%dD'%d)
plt.legend(loc=1,ncol=2)
plt.grid(linestyle='dashed')
plt.xlabel(r'$N$')
plt.ylabel(r'$\alpha_{d,N}$')
plt.ylim(1.05,1.42)
plt.savefig('dthoeory.pdf')
plt.show()
np.save('dtalpha',alpha)
