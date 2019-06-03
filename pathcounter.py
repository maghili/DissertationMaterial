#!/usr/bin/env python
import sprinkler
import Calculations
import nnfinder
import randomfinder as rf
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes, inset_axes
import numpy as np
import matplotlib
import math as m
from scipy.optimize import fsolve
import csv
#####################
nsp=20
N=[10+5*i for i in range(10)]+[100*i for i in range(1, 16)]
D=5; tAlpha = np.load('talpha.npy')
def job(N,nsp,D):
	''' Implementing all the modules (calculation, nnfinder, etc)
	to find all the dimensions and distributions '''
	num=[0]*(nsp*len(N))
	counter=0
	nr=100;
	MMd=[[0]*len(N) for i in range(2,D+1)];MMdStd=[[0]*len(N) for i in range(2,D+1)]
	BGbeta=[[0]*len(N) for i in range(2,D+1)];BGbetaStd=[[0]*len(N) for i in range(2,D+1)]
	BGd=[[0]*len(N) for i in range(2,D+1)];BGdStd=[[0]*len(N) for i in range(2,D+1)]
	ABPd=[[0]*len(N) for i in range(2,D+1)];ABPdStd=[[0]*len(N) for i in range(2,D+1)]
	MidPd=[[0]*len(N) for i in range(2,D+1)];MidPdStd=[[0]*len(N) for i in range(2,D+1)]
	MidPd2 =[[0]*len(N) for i in range(2, D+1)];MidPdStd2 = [[0]*len(N) for i in range(2, D+1)]
	vals=[[0]*len(N) for i in range(2,D+1)]
	stdvals=[[0]*len(N) for i in range(2,D+1)]
	for d in range(2,D+1):
		c=0
		peaksStd=[];widthsStd=[];peaksAvg=[];widthsAvg=[];
		for n in N:
			alpha=[];MP=[];beta=[];M=[];ABP=[];BG=[];AcAVG=[];EsAVG=[];MP2=[]; nct = 0
			for sp in range(nsp):
				print('\r%dD %d points, sprinkling %d'%(d,n,sp+1), end = '')
				s=sprinkler.sprinkler(n,d)
				link, Relations, PastRelations = nnfinder.nnfinder(s,n,d)
				[weights,MidPoint1,MidPoint2,MyrheimMeyer,b,lmax,avglen,a,tot]=Calculations.pathfinder(link,Relations,PastRelations,n,d)
				beta.append(b)
				alpha.append(a)
				####################
				#Midpoint dimension
				MP.append(MidPoint1)
				MP2.append(MidPoint2)
				####################
				#BG dimension
				BG.append(np.log(n)/np.log(lmax/2.0))
				####################
				#ABP dimension
				if n < 56:
					ABP.append(np.log(n)/np.log(avglen/tAlpha[nct]))
				else:
					ABP.append(np.log(n)/np.log(avglen/1.15))
				#####################
				#Myrheim-Meyer dimension calculation
				guess_dim=1.5;
				func =lambda dims: MyrheimMeyer-m.gamma(dims+1)*m.gamma(dims/2.)/(4*m.gamma(3*dims/2.0))
				dim=fsolve(func,guess_dim)
				M.append(dim)
			############################
			MidPd[d-2][c]=np.average(MP); MidPdStd[d-2][c]=np.std(MP)#Midpoint estimated avg dimension and std
			MidPd2[d-2][c] = np.average(MP2); MidPdStd2[d-2][c] = np.std(MP2)
			vals[d-2][c]=np.average(alpha); stdvals[d-2][c]=np.std(alpha)
			BGbeta[d-2][c]=np.average(beta); BGbetaStd[d-2][c]=np.std(beta)
			MMd[d-2][c]=np.average(M); MMdStd[d-2][c]=np.std(M)
			ABPd[d-2][c]=np.average(ABP); ABPdStd[d-2][c]=np.std(ABP)
			BGd[d-2][c]=np.average(BG); BGdStd[d-2][c]=np.std(BG)
			nct += 1
			#######################
			c += 1
	#################################
	#saving
	np.save('Alpha',vals); np.save('StdAlpha',stdvals)
	np.save('MidPoint',MidPd); np.save('MidPointStd',MidPdStd)
	np.save('MidPoint2', MidPd2); np.save('MidPointStd2', MidPdStd2)
	np.save('BGbeta',BGbeta); np.save('BGbetaStd',BGbetaStd)
	np.save('MMd',MMd); np.save('MMdStd',MMdStd)
	np.save('ABPd',ABPd); np.save('ABPdStd',ABPdStd)
	np.save('BGd',BGd); np.save('BGdStd',BGdStd)

def pl(N,D):
	''' Plotting the results'''
	Alpha=np.load('Alpha.npy').tolist(); StdAlpha=np.load('StdAlpha.npy').tolist()
	MidPoint=np.load('Midpoint.npy'); MidPointStd=np.load('MidPointStd.npy')
	MidPoint2=np.load('Midpoint2.npy'); MidPointStd2=np.load('MidPointStd2.npy')
	BGbeta=np.load('BGbeta.npy'); BGbetaStd=np.load('BGbetaStd.npy')
	MMd=np.load('MMd.npy'); MMdStd=np.load('MMdStd.npy')
	ABPd=np.load('ABPd.npy'); ABPdStd=np.load('ABPdStd.npy')
	BGd=np.load('BGd.npy'); BGdStd=np.load('BGdStd.npy')
	tAlpha=np.load('talpha.npy'); C=np.load('coef.npy')
	alpha50=np.load('Alpha50.npy').tolist(); alphastd50=np.load('StdAlpha50.npy').tolist()
	##################################
	cs=['tab:blue','tab:red','tab:green','tab:cyan','tab:pink','b','r','g']
	mk=['o','v','s','d','<','>','8','h']
	def expol(x,a,b,c):
	    return [a/j**c+b for j in x]
	fig, ax = subplots()
	newN=[50]+N[10:];
	for d in range(2,D+1):
		ax.errorbar(np.array(newN)-32+d*8,np.array(alpha50[d-2]+Alpha[d-2][10:]),yerr=alphastd50[d-2]+StdAlpha[d-2][10:],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)

	writer = csv.writer(open('AlphaTheory.csv', 'w'))
	writer.writerow(tAlpha)

	Nt=[10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
	ax.plot(Nt, tAlpha, '--', c='k',label='theory 2D')
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$\alpha_{d,N}$')
	ax.grid(linestyle='dashed',c='grey')
	ax.legend(loc=1)
	savefig('alpha_%dsp_1000.pdf'%nsp)
	##
	fig, ax = subplots()
	for d in range(2,D+1):
		ax.errorbar(np.array(N[10:])-32+d*8,BGbeta[d-2][10:],yerr=BGbetaStd[d-2][10:],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$\beta_{d,N}$')
	ax.grid(linestyle='dashed',c='grey')
	ax.legend(loc=1,ncol=D)
	savefig('beta_%dsp_1500.pdf'%nsp)
	##
	fig, ax = subplots()
	for d in range(2,D+1):
		ax.errorbar(N,MidPoint[d-2],yerr=MidPointStd[d-2],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$d$')
	ax.grid(linestyle='dashed',c='grey')
	ax.legend(loc=1,ncol=D)
	savefig('MidPoint_%dsp_1500_low.pdf'%nsp)
	##
	figure()
	for d in range(2,D+1):
		errorbar(N,MidPoint2[d-2],yerr=MidPointStd2[d-2],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	xlabel(r'$N$')
	ylabel(r'$d$')
	grid(linestyle='dashed',c='grey')
	legend(loc=1,ncol=D)
	savefig('MidPoint_%dsp_1500_up.pdf'%nsp)
	##
	figure()
	for d in range(2,D+1):
		avg = (np.array(MidPoint[d-2][10:]) + np.array(MidPoint2[d-2][10:]))/2.0
		std = (np.array(MidPointStd[d-2][10:]) + np.array(MidPointStd2[d-2][10:]))/2.0
		errorbar(N[10:],avg,yerr=std,c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	xlabel(r'$N$')
	ylabel(r'$d$')
	grid(linestyle='dashed',c='grey')
	legend(loc=1,ncol=D)
	savefig('MidPoint_%dsp_1500_avg.pdf'%nsp)
	##
	figure()
	for d in range(2,2+1):
		avg = (np.array(MidPoint[d-2]) + np.array(MidPoint2[d-2]))/2.0
		std = (np.array(MidPointStd[d-2]) + np.array(MidPointStd2[d-2]))/2.0
		errorbar(N,MidPoint[d-2],yerr=MidPointStd[d-2],c=cs[0],marker=mk[0],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D-MinofMax'%d)
		errorbar(N,MidPoint2[d-2],yerr=MidPointStd2[d-2],c=cs[1],marker=mk[1],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D-MaxofMin'%d)
		errorbar(N, avg, yerr = std, c = cs[2], marker = mk[2],elinewidth=1,capsize=2, ls = 'none', mec = 'k', label = '%d D-Average'%d)
	xlabel(r'$N$')
	ylabel(r'$d$')
	grid(linestyle='dashed',c='grey')
	legend(loc=1,ncol=D)
	savefig('MidPoint_%dsp_1500_compare.pdf'%nsp)
	##
	fig, ax = subplots()
	for d in range(2,D+1):
		ax.errorbar(N,ABPd[d-2],yerr=ABPdStd[d-2],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$d$')
	ax.grid(linestyle='dashed',c='grey')
	ax.set_ylim(ymax = 8)
	legend(loc=1,ncol=D)
	savefig('ABPd_%dsp_1500.pdf'%nsp)
	##
	fig, ax = subplots()
	for d in range(2,D+1):
		ax.errorbar(N,BGd[d-2],yerr=BGdStd[d-2],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$d$')
	ax.grid(linestyle='dashed',c='grey')
	ax.legend(loc=2,ncol=2)
	savefig('BGd_%dsp_1500.pdf'%nsp)
	##
	fig, ax = subplots()
	for d in range(2,D+1):
		ax.errorbar(N[10:],MMd[d-2][10:],yerr=MMdStd[d-2][10:],c=cs[d-2],marker=mk[d-2],mec='k',elinewidth=1,capsize=2,ls='none',label='%d D'%d)
	ax.set_xlabel(r'$N$')
	ax.set_ylabel(r'$d$')
	ax.grid(linestyle='dashed',c='grey')
	ax.set_ylim(ymax=D+0.8)
	ax.legend(loc=1,ncol=D)
	ax.set_ylim(1.5, 6)
	savefig('MMd_%dsp_1500.pdf'%nsp)
	##

	for d in range(2,D+1):
		fig, ax = subplots()
		errorMMd = abs(d - np.array(MMd[d-2]))
		errorMid = abs(d - (np.array(MidPoint[d-2])+np.array(MidPoint2[d-2]))/2.0)
		errorABPd = abs(d - np.array(ABPd[d-2]))
		errorBGd = abs(d - np.array(BGd[d-2]))
		ax.errorbar(N,errorMMd,yerr = MMdStd[d-2], c=cs[0],marker=mk[0],mec='k',ls='none',elinewidth=1,capsize=2,label='%d D-MMd error'%d)
		ax.errorbar(N,errorMid, yerr = MidPointStd[d-2], c=cs[1],marker=mk[1],mec='k',ls='none',elinewidth=1,capsize=2,label='%d D-MidPoint error'%d)
		ax.errorbar(N,errorABPd,yerr = ABPdStd[d-2], c=cs[2], marker=mk[2],mec='k',ls='none',elinewidth=1,capsize=2,label='%d D-ABPd error'%d)
		ax.errorbar(N,errorBGd,yerr = BGdStd[d-2], c=cs[3],marker=mk[3],mec='k',ls='none',elinewidth=1,capsize=2,label='%d D-BGd error'%d)

		ax.set_xlabel(r'$N$')
		ax.set_ylabel(r'$d$')
		ax.set_ylim(-1, 1)
		#title('Myrheim-Meyer dimension')
		ax.grid(linestyle='dashed',c='grey')
		ax.set_ylim(ymax=1)
		ax.legend(loc=4,ncol=1)
		savefig('%dsp_1500_errors_%dD.pdf'%(nsp,d))
##################################
#Running
job(N,nsp,D)
pl(N,D)
