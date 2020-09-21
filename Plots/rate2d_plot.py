import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy


def pois(lam,k):
	p=(lam**k*np.exp(-lam))/math.factorial(k)
	return p;

TYPE1='CMD2'
TYPE='vespaBC_1_gals_SFH6HSsig2_dr7_{}binR'.format(TYPE1)
Rbnsm_modden=np.loadtxt('Rden_obs1_{}.dat'.format(TYPE),unpack=False)
magcol=np.loadtxt('magcolVLIM_matrix_KSM3.dat',unpack=False)
magcol[:,1]=np.around(magcol[:,1],decimals=2)
magcol[:,2]=np.around(magcol[:,2],decimals=2)
magcol[:,3]=np.around(magcol[:,3],decimals=2)
magcol[:,4]=np.around(magcol[:,4],decimals=2)
dmr=0.25
dgr=0.1
centerVLIM=np.loadtxt('magcolVLIM_center_KSM3.dat',unpack=False)
MrcenVLIM=centerVLIM[centerVLIM<0.0]
grcenVLIM=centerVLIM[centerVLIM>=0.0]
MrBINvlim=22
grBINvlim=24

DAT1=np.loadtxt('{}.dat'.format(TYPE),unpack=False)
if TYPE1=='CMD2':
        DAT=DAT1[:,5:]
else:
        DAT=DAT1[:,3:]
ng=DAT1[:,0]
ng=ng/ng.sum()


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

contlevel=np.zeros((NNobs,2))
L=len(Rbnsm_modden[0,:])
q=0
sig1=0.682
sig2=0.954
for k in range(Nnobs):
	#for i in range(0,Ntminobs):
	F=np.zeros(L)
	Rtmp=np.sort(Rbnsm_modden[q,:])[::-1]
	for ii in range(0,L):
		F[ii]=Rtmp[:ii].sum()
	contlevel[q,0]=np.log10(Rtmp[F>=sig2][0]/(dgr*dmr))
	contlevel[q,1]=np.log10(Rtmp[F>=sig1][0]/(dgr*dmr))
	q=q+1


contlevelng=np.zeros(2)
L=len(ng)
q=0
sig1=0.682
sig2=0.954
sig3=0.997
        #for i in range(0,Ntminobs):
F=np.zeros(L)
Ngtmp=np.sort(ng)[::-1]
for ii in range(0,L):
	F[ii]=Ngtmp[:ii].sum()
#contlevelng[0]=np.log10(Ngtmp[F>=sig3][0]/(dgr*dmr))
contlevelng[0]=np.log10(Ngtmp[F>=sig2][0]/(dgr*dmr))
contlevelng[1]=np.log10(Ngtmp[F>=sig1][0]/(dgr*dmr))

NG=np.log10(ng/(dgr*dmr))
NG=NG[~np.isinf(NG)]
contlevelNG=np.arange(NG.min(),NG.max(),0.35)
plot6=1
if plot6==1:
	gr=np.arange(0,1.21,0.01)
	Mr=np.around(-1*np.arange(16.5,22.25,0.25)[:-1],decimals=1)
	Y=len(gr)
	X=len(Mr)

	mat_bnsmdencont=np.zeros((NNobs,MrBINvlim,grBINvlim))
        for jj in range(0,NNobs):
                kk=0
                for i in range(0,MrBINvlim):
                        for j in range(0,grBINvlim): 
                                mat_bnsmdencont[jj,i,j]=Rbnsm_modden[jj,kk]/(dgr*dmr)
                                kk=kk+1

        mat_bnsmdencontng=np.zeros((MrBINvlim,grBINvlim)) 
	kk=0
        for i in range(0,MrBINvlim):
        	for j in range(0,grBINvlim):
                	mat_bnsmdencontng[i,j]=ng[kk]/(dgr*dmr)
                        kk=kk+1


	mat_bnsmden=np.zeros((NNobs,X-1,Y-1))
	for jj in range(0,NNobs):
		for i in range(0,Y-1):
			for j in range(0,X-1):
				magtmp=magcol[Mr[j]<=magcol[:,1]]
				magtmp=magtmp[Mr[j]>=magtmp[:,2]]
				magtmp=magtmp[gr[i]>=magtmp[:,3]]
				ind=int(magtmp[gr[i]<=magtmp[:,4]][0][0]) 
				mat_bnsmden[jj,j,i]=Rbnsm_modden[jj,ind]/(dgr*dmr)
			
				
        	
	mat_bnsmLOG=np.log10(mat_bnsmden)
	mat_bnsmcontLOG=np.log10(mat_bnsmdencont)
        mat_bnsmcontngLOG=np.log10(mat_bnsmdencontng)
	levelsavg=np.arange(-5.0,mat_bnsmLOG[~np.isinf(mat_bnsmLOG)].max()+0.1,0.01)
	vmin=levelsavg.min()
	vmax=levelsavg.max()
	MIN=vmin
	#mat_bnsmLOG[~np.isinf(mat_bnsmLOG)].min()
	mat_bnsmLOG[mat_bnsmLOG<MIN]=MIN-0.1
	fig, axs = plt.subplots(1, 3, sharex=True,sharey=True)
	fig.subplots_adjust(wspace=0,hspace=0)
	(ax1,ax2,ax3)=axs
	COL='cool'

	#im1=ax1.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[0,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im2=ax2.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[1,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im3=ax3.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[2,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
        im1=ax1.contourf(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[0,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
        im2=ax2.contourf(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[1,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
        im3=ax3.contourf(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[2,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im4=ax4.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[3,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im5=ax5.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[4,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im6=ax6.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[5,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im7=ax7.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[6,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im8=ax8.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[7,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
	#im9=ax9.contourf(Mr[:-1],gr[:-1],mat_bnsmLOG[8,:,:].T,levelsavg,vmin=vmin,vmax=vmax,cmap=COL)
        ax1.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,contlevelNG,linestyles='--',colors='white')
        ax2.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,contlevelNG,linestyles='--',colors='white')
        ax3.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,contlevelNG,linestyles='--',colors='white')
	ax1.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[0,:,:].T,contlevel[0,:],linestyles=['-','-'],colors='k')
        ax2.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[1,:,:].T,contlevel[1,:],linestyles=['-','-'],colors='k')
        ax3.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[2,:,:].T,contlevel[2,:],linestyles=['-','-'],colors='k') 
        #ax1.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,linestyles='--',colors='white')
        #ax2.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,linestyles='--',colors='white')
        #ax3.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontngLOG[:,:].T,linestyles='--',colors='white')
	#ax4.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[3,:,:].T,contlevel[3,:],linestyles=['--','-'],colors='k')
        #ax5.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[4,:,:].T,contlevel[4,:],linestyles=['--','-'],colors='k')
        #ax6.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[5,:,:].T,contlevel[5,:],linestyles=['--','-'],colors='k')
        #ax7.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[6,:,:].T,contlevel[6,:],linestyles=['--','-'],colors='k')
        #ax8.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[7,:,:].T,contlevel[7,:],linestyles=['--','-'],colors='k')
        #ax9.contour(MrcenVLIM,grcenVLIM,mat_bnsmcontLOG[8,:,:].T,contlevel[8,:],linestyles=['--','-'],colors='k')
	ax1.set_title(r'n=-1.5, t$_{\rm m}$=0.01 Gyr',fontsize=25)
	ax2.set_title(r'n=-1.1, t$_{\rm m}$=0.035 Gyr',fontsize=25)
	ax3.set_title(r'n=-0.5, t$_{\rm m}$=1.0 Gyr',fontsize=25)
	#ax4.set_title(r'n={},t$_m$={}'.format(nobs[1],tminobs[0]),weight="bold",fontsize=20)
	#ax5.set_title(r'n={},t$_m$={}'.format(nobs[1],tminobs[1]),weight="bold",fontsize=20)
	#ax6.set_title(r'n={},t$_m$={}'.format(nobs[1],tminobs[2]),weight="bold",fontsize=20)
	#ax7.set_title(r'n={},t$_m$={}'.format(nobs[2],tminobs[0]),weight="bold",fontsize=20)
	#ax8.set_title(r'n={},t$_m$={}'.format(nobs[2],tminobs[1]),weight="bold",fontsize=20)
	#ax9.set_title(r'n={},t$_m$={}'.format(nobs[2],tminobs[2]),weight="bold",fontsize=20)
	fig.subplots_adjust(right=0.8)
	divider = make_axes_locatable(ax3)
	cax = divider.append_axes('right', size='10%', pad=0)
	colorrange=np.round(np.arange(vmin,vmax+abs(vmax-vmin)/4.,abs(vmax-vmin)/4.),decimals=1)
	cbar=fig.colorbar(im1,cax,ticks=colorrange)
	cbar.set_label(r'log(PDF)',labelpad=30,rotation=270,fontsize=25)
	cbar.ax.tick_params(labelsize=25)
	#divider = make_axes_locatable(ax6)
	#cax = divider.append_axes('right', size='10%', pad=0)
	#cbar=fig.colorbar(im1,cax,ticks=colorrange)
	#cbar.set_label(r'log(PDF)', weight='bold',labelpad=30,rotation=270,fontsize=20)
	#cbar.ax.tick_params(labelsize=20)
	#divider = make_axes_locatable(ax9)
	#cax = divider.append_axes('right', size='10%', pad=0)
	#cbar=fig.colorbar(im1,cax,ticks=colorrange)
	#cbar.set_label(r'log(PDF)', weight='bold',labelpad=30,rotation=270,fontsize=20)
	#cbar.ax.tick_params(labelsize=20)
	ax1.set_ylabel('g-r',fontsize=25)
	ax1.invert_xaxis()
	ax1.errorbar(-20.09,0.72,yerr=0.02,xerr=0.16,lw=3,color='w')
        ax2.errorbar(-20.09,0.72,yerr=0.02,xerr=0.16,lw=3,color='w')
        ax3.errorbar(-20.09,0.72,yerr=0.02,xerr=0.16,lw=3,color='w')
        #ax4.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
        #ax5.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
        #ax6.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
        #ax7.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
        #ax8.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
        #ax9.errorbar(-20.895,0.64,yerr=0.23,xerr=0.165,lw=3)
	ax1.margins(x=0)
	ax1.margins(y=0)
#	for tick in ax1.get_xticklabels():
#		tick.set_rotation(90)
#	for tick in ax2.get_xticklabels():
#		tick.set_rotation(90)
#	for tick in ax3.get_xticklabels():
#		tick.set_rotation(90)
	#ax4.tick_params(axis="y", labelsize=20)
	#ax7.tick_params(axis="y", labelsize=20)



	ax1.set_xlabel(r'M$_{\rm r}$-$\rm 5$log$\rm h$',fontsize=25)
	ax2.set_xlabel(r'M$_{\rm r}$-$\rm 5$log$\rm h$',fontsize=25)
	ax3.set_xlabel(r'M$_{\rm r}$-$\rm 5$log$\rm h$',fontsize=25) 
	ax1.set_ylabel('g-r',fontsize=25)
	ax1.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax1.tick_params('both', length=10, width=2, which='minor')
	ax2.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax2.tick_params('both', length=10, width=2, which='minor')
	ax3.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax3.tick_params('both', length=10, width=2, which='minor')
	#ax4.set_ylabel('g-r',fontsize=20,weight='bold')
	#ax7.set_ylabel('g-r',fontsize=20,weight='bold')
	plt.show()



