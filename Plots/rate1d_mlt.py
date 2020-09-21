import numpy as np 
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

TYPE='vespaBC_1_gals_SFH6HSsig2_dr7_GR2binR'
Rbnsm_gr=np.loadtxt('Rden_obs1_{}.dat'.format(TYPE),unpack=False)
DAT1=np.loadtxt('{}.dat'.format(TYPE),unpack=False)
ng_gr=DAT1[:,0]
ng_gr=ng_gr/ng_gr.sum()

gr=np.arange(0.0,1.25,0.05)
grcen=gr+0.05/2.
grcen=grcen[:-1]
dgr=0.05

ng_gr=ng_gr/dgr
Rgr=Rbnsm_gr/dgr
Ngr=len(Rgr[:,0])


TYPE='vespaBC_1_gals_SFH6HSsig2_dr7_MR2binR'
Rbnsm_mr=np.loadtxt('Rden_obs1_{}.dat'.format(TYPE),unpack=False)
DAT1=np.loadtxt('{}.dat'.format(TYPE),unpack=False)
ng_mr=DAT1[:,0]
ng_mr=ng_mr/ng_mr.sum()

mr=np.arange(-22.0,-16.375,0.125)
mrcen=mr+0.125/2.
mrcen=mrcen[:-1][::-1]
dmr=0.125

ng_mr=ng_mr/dmr
Rmr=Rbnsm_mr/dmr
Nmr=len(Rmr[:,0])


TYPE='vespaBC_1_gals_SFH6HSsig2_dr7_SM2binR'
Rbnsm_sm=np.loadtxt('Rden_obs1_{}.dat'.format(TYPE),unpack=False)
DAT1=np.loadtxt('{}.dat'.format(TYPE),unpack=False)
ng_sm=DAT1[:,0]
ng_sm=ng_sm/ng_sm.sum()

dsm=0.1
sm=np.arange(7.5,11.5+dsm,dsm)
smcen=sm+dsm/2.
smcen=smcen[:-1]

ng_sm=ng_sm/dsm
Rsm=Rbnsm_sm/dsm
Nsm=len(Rsm[:,0])


TYPE='vespaBC_1_gals_SFH6HSsig2_dr7_sSFR2binR'
Rbnsm_ssfr=np.loadtxt('Rden_obs1_{}.dat'.format(TYPE),unpack=False)
DAT1=np.loadtxt('{}.dat'.format(TYPE),unpack=False)
ng_ssfr=DAT1[:,0]
ng_ssfr=ng_ssfr/ng_ssfr.sum()

dssfr=0.125
ssfr=np.arange(-12.625,-8.5,dssfr)
ssfrcen=ssfr+dssfr/2.
ssfrcen=ssfrcen[:-1]

ng_ssfr=ng_ssfr/dssfr
Rssfr=Rbnsm_ssfr/dssfr
Nssfr=len(Rssfr[:,0])

nobs=[-1.5,-1.1,-0.5]
Nnobs=len(nobs)
tminobs=[0.01,0.035,1.0]
Ntminobs=len(tminobs)
dnmod=0.05
nmod=np.around(np.arange(-1.5,-0.45,dnmod),decimals=2)
tminmod=np.around((np.power(10,np.linspace(1,3,21))*10**6)/10**9,decimals=3)
Nnmod=len(nmod)
Ntminmod=len(tminmod)
NNmod=Nnmod*Ntminmod
NNobs=3


color_idx = np.linspace(0, 1,Ngr)
v=0
LS=['--','-',':']


fig, axs = plt.subplots(2, 2, sharex=False,sharey=False)
fig.subplots_adjust(wspace=0,hspace=0)
(ax1,ax2),(ax3,ax4)=axs


ax1.plot(grcen,ng_gr,lw=1,color='grey',ls='-')
ax2.plot(mrcen,ng_mr,lw=1,color='grey',ls='-')
ax3.plot(smcen,ng_sm,lw=1,color='grey',ls='-')
ax4.plot(ssfrcen,ng_ssfr,lw=1,color='grey',ls='-')

ax1.plot(grcen,Rgr[0,:],lw=5,color='r',ls=LS[0],label=r'n=-1.5, t$_{\rm m}$=0.01 Gyr')
ax2.plot(mrcen,Rmr[0,:],lw=5,color='b',ls=LS[0])
ax3.plot(smcen,Rsm[0,:],lw=5,color='cyan',ls=LS[0])
ax4.plot(ssfrcen,Rssfr[0,:],lw=5,color='orange',ls=LS[0])

ax1.plot(grcen,Rgr[1,:],lw=5,color='r',ls=LS[1],label=r'n=-1.1, t$_{\rm m}$=0.035 Gyr')
ax2.plot(mrcen,Rmr[1,:],lw=5,color='b',ls=LS[1])
ax3.plot(smcen,Rsm[1,:],lw=5,color='cyan',ls=LS[1])
ax4.plot(ssfrcen,Rssfr[1,:],lw=5,color='orange',ls=LS[1])

ax1.plot(grcen,Rgr[2,:],lw=5,color='r',ls=LS[2],label=r'n=-0.5, t$_{\rm m}$=1.0 Gyr')
ax2.plot(mrcen,Rmr[2,:],lw=5,color='b',ls=LS[2])
ax3.plot(smcen,Rsm[2,:],lw=5,color='cyan',ls=LS[2])
ax4.plot(ssfrcen,Rssfr[2,:],lw=5,color='orange',ls=LS[2])

#for tick in ax1.get_xticklabels():
#	tick.set_rotation(90)
#for tick in ax2.get_xticklabels():
#        tick.set_rotation(90)
ax1.tick_params(axis="y", labelsize=20)
ax1.tick_params(axis="x", labelsize=20)
ax1.set_xlabel(r'g-r',fontsize=20)
ax1.set_ylabel(r'PDF',fontsize=20)
ax2.tick_params(axis="y", labelsize=20)
ax2.tick_params(axis="x", labelsize=20)
ax2.invert_xaxis()
ax2.set_xlabel(r'M$_{\rm r}$-$\rm 5$log$\rm h$',fontsize=20)
ax2.set_ylabel(r'PDF',fontsize=20)
#ax2.set_ylabel(r'log(PDF)',fontsize=20,weight='bold')
#for tick in ax3.get_xticklabels():
#        tick.set_rotation(90)
ax3.tick_params(axis="y", labelsize=20)
ax3.tick_params(axis="x", labelsize=20)
ax3.set_xlabel(r'log(M$_*$/M$_\odot$)',fontsize=20)
ax3.set_ylabel(r'PDF',fontsize=20)
#for tick in ax4.get_xticklabels():
#        tick.set_rotation(90)
ax4.tick_params(axis="y", labelsize=20)
ax4.tick_params(axis="x", labelsize=20)
ax4.set_xlabel(r'log(sSFR/yr$^{-1}$)',fontsize=20)
ax4.set_ylabel(r'PDF',fontsize=20)
#ax4.set_ylabel(r'log(PDF)',fontsize=20,weight='bold')
ax1.set_ylim(0.0,3.6)
ax2.set_ylim(0.0,0.5)
ax3.set_ylim(0.0,1.0)
ax4.set_ylim(0.0,1.0)

ax1.set_xlim(grcen[0],grcen[-1])
ax2.set_xlim(mrcen[0],mrcen[-1])
ax3.set_xlim(8,smcen[-1])
ax4.set_xlim(ssfrcen[0],ssfrcen[-1])

ax1.axvspan(0.7,0.74,alpha=0.25,color='k')
ax2.axvspan(-19.93,-20.25,alpha=0.25,color='k')
ax3.axvspan(10.62,10.68,alpha=0.25,color='k')
ax4.axvspan(-13.34,-12.18,alpha=0.25,color='k')

#ax1.set_yticks(np.arange(-2.0,0.5,0.5))
#ax3.set_yticks(np.arange(-2.0,0.5,0.5))
ax1.set_xticks(np.arange(0.1,1.0+0.2,0.2))
ax2.set_xticks(np.arange(-21,-16.,1.0))
ax3.set_xticks(np.arange(8.0,12.0,1.0))
ax4.set_xticks(np.arange(-12.0,-8.0,1.0))

ax1.tick_params('both', length=15, width=2, which='major',labelsize=20)
ax1.tick_params('both', length=10, width=2, which='minor')
ax2.tick_params('both', length=15, width=2, which='major',labelsize=20)
ax2.tick_params('both', length=10, width=2, which='minor')
ax3.tick_params('both', length=15, width=2, which='major',labelsize=20)
ax3.tick_params('both', length=10, width=2, which='minor')
ax4.tick_params('both', length=15, width=2, which='major',labelsize=20)
ax4.tick_params('both', length=10, width=2, which='minor')

#ax1.legend(loc='upper left',fontsize=15)
custom_lines = [Line2D([0], [0], color='k', ls='--', lw=5),Line2D([0], [0], color='k', ls='-', lw=5),Line2D([0], [0], color='k', ls=':', lw=5)]
ax1.legend(custom_lines,['n=-1.5, t$_m$=0.01 Gyr', 'n=-1.1, t$_m$=0.035 Gyr', 'n=-0.5, t$_m$=1.0 Gyr'],loc='upper left')

plt.show()

