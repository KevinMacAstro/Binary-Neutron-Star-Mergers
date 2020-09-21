import numpy as np
import pyfits
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy
import scipy.stats
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF

####### From https://gist.github.com/eteq/4599814
try:
        from scipy.spatial import cKDTree as KDT
except ImportError:
        from scipy.spatial import KDTree as KDT

def spherematch(ra1, dec1, ra2, dec2, tol=None, nnearest=1):
    """
    Finds matches in one catalog to another.
    Parameters
    ra1 : array-like
        Right Ascension in degrees of the first catalog
    dec1 : array-like
        Declination in degrees of the first catalog (shape of array must match `ra1`)
    ra2 : array-like
        Right Ascension in degrees of the second catalog
    dec2 : array-like
        Declination in degrees of the second catalog (shape of array must match `ra2`)
    tol : float or None, optional
        How close (in degrees) a match has to be to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned.
    nnearest : int, optional
        The nth neighbor to find.  E.g., 1 for the nearest nearby, 2 for the
        second nearest neighbor, etc.  Particularly useful if you want to get
        the nearest *non-self* neighbor of a catalog.  To do this, use:
        ``spherematch(ra, dec, ra, dec, nnearest=2)``
    Returns
    -------
    idx1 : int array
        Indecies into the first catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    idx2 : int array
        Indecies into the second catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    ds : float array
        Distance (in degrees) between the matches
    """
    ra1 = np.array(ra1, copy=False)
    dec1 = np.array(dec1, copy=False)
    ra2 = np.array(ra2, copy=False)
    dec2 = np.array(dec2, copy=False)

    if ra1.shape != dec1.shape:
        raise ValueError('ra1 and dec1 do not match!')
    if ra2.shape != dec2.shape:
        raise ValueError('ra2 and dec2 do not match!')

    x1, y1, z1 = _spherical_to_cartesian(ra1.ravel(), dec1.ravel())

    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1

    x2, y2, z2 = _spherical_to_cartesian(ra2.ravel(), dec2.ravel())
    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0] = x2
    coords2[:, 1] = y2
    coords2[:, 2] = z2

    kdt = KDT(coords2)
    if nnearest == 1:
        idxs2 = kdt.query(coords1)[1]
    elif nnearest > 1:
        idxs2 = kdt.query(coords1, nnearest)[1][:, -1]
    else:
        raise ValueError('invalid nnearest ' + str(nnearest))

    ds = _great_circle_distance(ra1, dec1, ra2[idxs2], dec2[idxs2])

    idxs1 = np.arange(ra1.size)

    if tol is not None:
        msk = ds < tol
        idxs1 = idxs1[msk]
        idxs2 = idxs2[msk]
        ds = ds[msk]

    return idxs1, idxs2, ds


def _spherical_to_cartesian(ra, dec):
    """
    (Private internal function)
    Inputs in degrees.  Outputs x,y,z
    """
    rar = np.radians(ra)
    decr = np.radians(dec)

    x = np.cos(rar) * np.cos(decr)
    y = np.sin(rar) * np.cos(decr)
    z = np.sin(decr)

    return x, y, z

def _great_circle_distance(ra1, dec1, ra2, dec2):
    """
    (Private internal function)
    Returns great circle distance.  Inputs in degrees.
    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.
    """
    from numpy import radians, degrees, sin, cos, arctan2, hypot

    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"
    lambs = radians(ra1)
    phis = radians(dec1)
    lambf = radians(ra2)
    phif = radians(dec2)

    dlamb = lambf - lambs

    numera = cos(phif) * sin(dlamb)
    numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
    numer = hypot(numera, numerb)
    denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)
    return degrees(arctan2(numer, denom))
########################################
################# From Z.Z.
def binarySearch(alist, item):
        first = 0
        last = len(alist)-1
        found = False
        while first<=last and not found:
                midpoint = (first + last)//2
                if alist[midpoint] == item:
                    found = True
                else:
                    if item < alist[midpoint]:
                        last = midpoint-1
                    else:
                        first = midpoint+1
        if found:
                return found,midpoint
        else:
                return found, 00
##########################################
#(0)ra,(1)dec,(2)z,(3-19 vespa bins)
Data1=np.loadtxt('vespaBC_1_dr7_gal.dat',unpack=False)
#Data2=np.loadtxt('vespaM05_2_dr7_gal.dat',unpack=False)

N=16
lookstart=np.array([0.002,0.020,0.030,0.048,0.074,0.115,0.177,0.275,0.425,0.658,1.02,1.57,2.44,3.78,5.84,9.04])
lookstop=np.array([0.020,0.030,0.048,0.074,0.115,0.177,0.275,0.425,0.658,1.02,1.57,2.44,3.78,5.84,9.04,14.00])
xbox=np.array([0.002,0.020,0.030,0.048,0.074,0.115,0.177,0.275,0.425,0.658,1.02,1.57,2.44,3.78,5.84,9.04,14.00])

lookwidth=lookstop-lookstart
lookcenter=lookstart+lookwidth/2.
deltlook=lookwidth*10**9
tlook=lookcenter*10**9
xboxcen=np.zeros(N)
for i in range(0,N):
	xboxcen[i]=xbox[i]+(xbox[i+1]-xbox[i])/2.


CMDcut=1
SDSS1=1
if SDSS1==1:
        y07_sdss_ID=np.loadtxt('SDSS7_ID.dat',unpack=False)

        y07_sdss=np.loadtxt('SDSS7.dat',unpack=False)
        y07_sdss_NYU=y07_sdss[:,1].astype(int)
        Y07_sdss_ra=y07_sdss[:,2]
        Y07_sdss_dec=y07_sdss[:,3]
        Y07_sdss_z=y07_sdss[:,4]
        Y07_sdss_Mr=y07_sdss[:,8]
        Y07_sdss_gr=y07_sdss[:,9]
        Y07_sdss_rlim=y07_sdss[:,6]


	if CMDcut==1:
		cut1=Y07_sdss_Mr>=-22.0
		cut2=Y07_sdss_Mr<=-16.5
		y07_sdss_ra=Y07_sdss_ra[cut1*cut2]
        	y07_sdss_dec=Y07_sdss_dec[cut1*cut2]
        	y07_sdss_z=Y07_sdss_z[cut1*cut2]
        	y07_sdss_gr=Y07_sdss_gr[cut1*cut2]
        	y07_sdss_rlim=Y07_sdss_rlim[cut1*cut2]
		y07_sdss_Mr=Y07_sdss_Mr[cut1*cut2]


		cut1=y07_sdss_gr>=0.0
		cut2=y07_sdss_gr<=1.2
        	y07_sdss_ra=y07_sdss_ra[cut1*cut2]
        	y07_sdss_dec=y07_sdss_dec[cut1*cut2]
        	y07_sdss_z=y07_sdss_z[cut1*cut2]
        	y07_sdss_gr=y07_sdss_gr[cut1*cut2]
        	y07_sdss_rlim=y07_sdss_rlim[cut1*cut2]
		y07_sdss_Mr=y07_sdss_Mr[cut1*cut2]

	else:
                y07_sdss_ra=Y07_sdss_ra
                y07_sdss_dec=Y07_sdss_dec
                y07_sdss_z=Y07_sdss_z
                y07_sdss_gr=Y07_sdss_gr
                y07_sdss_rlim=Y07_sdss_rlim
                y07_sdss_Mr=Y07_sdss_Mr



        Olam=0.726
        Om=1-Olam
        h=0.705

	y07_sdss_d=2000*(np.sqrt(Olam+Om*(1+3*y07_sdss_z))-1)/Om

        TOTsa=7984.*(np.pi/180.)**2
        Dmax=np.power(10,(y07_sdss_rlim-y07_sdss_Mr)/5+1)/(10**6) # max luminosity distance with flux limit in Mpc
        Vmax=TOTsa*Dmax**3/3.

        n1v1,n2v1,ds=spherematch(y07_sdss_ra,y07_sdss_dec,Data1[:,0],Data1[:,1],tol=0.00027778)
        data1=np.vstack((y07_sdss_ra[n1v1],y07_sdss_dec[n1v1],y07_sdss_z[n1v1],y07_sdss_d[n1v1],Data1[n2v1,3:19].T,y07_sdss_Mr[n1v1],y07_sdss_gr[n1v1],Dmax[n1v1],Vmax[n1v1])).T


mpa=1
if mpa==1:

# From http://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/sfrs.html
#mpa_SFR=pyfits.open('gal_totsfr_dr7_v5_2.fits')
        mpa_sSFR=pyfits.open('gal_totspecsfr_dr7_v5_2.fits')
# From http://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/Data/stellarmass.html
        mpa_MASS=pyfits.open('totlgm_dr7_v5_2.fit')
# From http://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/raw_data.html
        mpa_galINFO=pyfits.open('gal_info_dr7_v5_2.fit')
	#mpa_oxy=pyfits.open('gal_fiboh_dr7_v5_2.fits')


# Some of the MPA/JHU galaxies have duplicates. This file identifies which galaxies in the cat have no matches.
        dup=open('all_matches_dr7.dat','r').readlines()
        Ndup=len(dup)
        dups=np.zeros(Ndup)
        for i in range(0,Ndup):
                if i < Ndup:
                        dups[i]=dup[i].split()[1]

        mpa_ssfrAVG=mpa_sSFR[1].data['AVG'][dups<0]
        mpa_ssfrMED=mpa_sSFR[1].data['MEDIAN'][dups<0]
        mpa_ssfrMODE=mpa_sSFR[1].data['MODE'][dups<0]
        mpa_massAVG=mpa_MASS[1].data['AVG'][dups<0]
        mpa_massMED=mpa_MASS[1].data['MEDIAN'][dups<0]
        mpa_massMODE=mpa_MASS[1].data['MODE'][dups<0]
        mpa_ra=mpa_galINFO[1].data['RA'][dups<0]
        mpa_dec=mpa_galINFO[1].data['DEC'][dups<0]
        mpa_z=mpa_galINFO[1].data['Z'][dups<0]
        mpa_plateID=mpa_galINFO[1].data['PLATEID'][dups<0]
        mpa_MJD=mpa_galINFO[1].data['MJD'][dups<0]
        mpa_fiberID=mpa_galINFO[1].data['FIBERID'][dups<0]
	#mpa_OXY=mpa_oxy[1].data['AVG'][dups<0]

	cut=mpa_massAVG!=-1.0   
	mpa_mass=mpa_massAVG[cut]
	mpa_ra=mpa_ra[cut]
	mpa_dec=mpa_dec[cut]
	mpa_ssfr=mpa_ssfrAVG[cut]
	#mpa_OXY=mpa_OXY[cut]
	cut=mpa_ssfr!=-99.0 
        mpa_mass=mpa_mass[cut]
        mpa_ra=mpa_ra[cut]
        mpa_dec=mpa_dec[cut]
        mpa_ssfr=mpa_ssfr[cut]
	#mpa_OXY=mpa_OXY[cut]

        n1v1,n2v1,ds=spherematch(data1[:,0],data1[:,1],mpa_ra,mpa_dec,tol=0.00027778)
   #     n1v2,n2v2,ds=spherematch(data2[:,0],data2[:,1],mpa_ra,mpa_dec,tol=0.00027778)

        data1=np.vstack((data1[n1v1,:].T,mpa_mass[n2v1],mpa_ssfr[n2v1])).T
    #    data2=np.vstack((data2[n1v2,:].T,mpa_mass[n2v2],mpa_ssfr[n2v2])).T


	

SFHcalc=1
if SFHcalc==1:
	for k in range(0,N):
        	data1[:,4+k]=data1[:,4+k]/deltlook[k]
     #           data2[:,4+k]=data2[:,4+k]/deltlook[k]

outs1=1
outs2=2
M05=0
BC=1
if outs1==1:
	SIG=6
	dataTMP=data1
        L=len(dataTMP)
	index=np.arange(0,L,1)
        LL=0  
        RA=[]
	DEC=[]
        while L!=LL:
        	L=len(dataTMP)
                for i in range(0,N):
                	out1=[]
                        outs=index[dataTMP[:,4+i]!=0.0][np.log10(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0])>=(np.log10(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0]).mean()+SIG*np.log10(dataTMP[:,+i][dataTMP[:,4+i]!=0.0]).std())] 
                        outs=np.unique(outs)
                        if len(outs)!=0.0:
                        	for j in range(0,len(outs)):
                                	out1.append(outs[j])
                                out1=np.array(out1)
                                dataTMP=np.delete(dataTMP,out1,axis=0)
			LL=len(dataTMP)
                        index=np.arange(0,LL,1)
        N1=len(data1)
        data1=dataTMP

       # dataTMP=data2
       # L=len(dataTMP)
       # index=np.arange(0,L,1)
       # LL=0
       # RA=[]
       # DEC=[]
       # while L!=LL:
       #         L=len(dataTMP)
       #         for i in range(0,N):
       #                 out1=[]
       #                 outs=index[dataTMP[:,4+i]!=0.0][np.log10(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0])>=(np.log10(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0]).mean()+SIG*np.log10(dataTMP[:,+i][dataTMP[:,4+i]!=0.0]).std())]
       #                 outs=np.unique(outs)
       #                 if len(outs)!=0.0:
       #                         for j in range(0,len(outs)):
       #                                 out1.append(outs[j])
       #                         out1=np.array(out1)
       #                         dataTMP=np.delete(dataTMP,out1,axis=0)
       #                 LL=len(dataTMP)
       #                 index=np.arange(0,LL,1)
#	N2=len(data2)
#        data2=dataTMP
      

if outs2==1:
        SIG=20
        dataTMP=data1
        L=len(dataTMP)
        index=np.arange(0,L,1)
        LL=0
        RA=[]
        DEC=[]
        while L!=LL:
                L=len(dataTMP)
                for i in range(0,N):
                        out1=[]
                        outs=index[dataTMP[:,4+i]!=0.0][dataTMP[:,4+i][dataTMP[:,4+i]!=0.0]>=(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0].mean()+SIG*dataTMP[:,+i][dataTMP[:,4+i]!=0.0].std())]
                        outs=np.unique(outs)
                        if len(outs)!=0.0:
                                for j in range(0,len(outs)):
                                        out1.append(outs[j])
                                out1=np.array(out1)
                                dataTMP=np.delete(dataTMP,out1,axis=0)
                        LL=len(dataTMP)
                        index=np.arange(0,LL,1)
        N1=len(data1)
        data1=dataTMP

        dataTMP=data2
        L=len(dataTMP)
        index=np.arange(0,L,1)
        LL=0
        RA=[]
        DEC=[]
        while L!=LL:
                L=len(dataTMP)
                for i in range(0,N):
                        out1=[]
                        outs=index[dataTMP[:,4+i]!=0.0][dataTMP[:,4+i][dataTMP[:,4+i]!=0.0]>=(dataTMP[:,4+i][dataTMP[:,4+i]!=0.0].mean()+SIG*dataTMP[:,+i][dataTMP[:,4+i]!=0.0].std())]
                        outs=np.unique(outs)
                        if len(outs)!=0.0:
                                for j in range(0,len(outs)):
                                        out1.append(outs[j])
                                out1=np.array(out1)
                                dataTMP=np.delete(dataTMP,out1,axis=0)
                        LL=len(dataTMP)
                        index=np.arange(0,LL,1)
        N2=len(data2)
        data2=dataTMP

if outs2==2:
	if M05==1:
		data1=data1[data1[:,5]<=1200]
        	data1=data1[data1[:,8]<=500]
        	data1=data1[data1[:,9]<=800]
        	data1=data1[data1[:,10]<=500]
        	data1=data1[data1[:,11]<=500]
        	data1=data1[data1[:,12]<=600]
        	data1=data1[data1[:,13]<=350]
        	data1=data1[data1[:,14]<=500]
        	data1=data1[data1[:,15]<=400]
        	data1=data1[data1[:,16]<=600]
        	data1=data1[data1[:,17]<=600]
        	data1=data1[data1[:,18]<=600]
        	data1=data1[data1[:,19]<=1000]
	
#	        data2=data2[data2[:,5]<=1200]
#	        data2=data2[data2[:,8]<=250]
#	        data2=data2[data2[:,9]<=100]
#	        data2=data2[data2[:,10]<=70]
#	        data2=data2[data2[:,11]<=640]
#	        data2=data2[data2[:,12]<=30]
#	        data2=data2[data2[:,13]<=300]
#	        data2=data2[data2[:,14]<=200]
#	        data2=data2[data2[:,15]<=300]
#	        data2=data2[data2[:,16]<=500]
#	        data2=data2[data2[:,17]<=400]
#	        data2=data2[data2[:,18]<=600]
#	        data2=data2[data2[:,19]<=600]
	
	if BC==1:
                data1=data1[data1[:,5]<=350]
                data1=data1[data1[:,8]<=250]
                data1=data1[data1[:,9]<=140]
                data1=data1[data1[:,10]<=100]
                data1=data1[data1[:,11]<=150]
                data1=data1[data1[:,12]<=400]
                data1=data1[data1[:,13]<=800]
                data1=data1[data1[:,14]<=700]
                data1=data1[data1[:,15]<=500]
                data1=data1[data1[:,16]<=400]
                data1=data1[data1[:,17]<=500]
                data1=data1[data1[:,18]<=400]
                data1=data1[data1[:,19]<=500]



smooth=0
if smooth==1:
        tot=np.sum(data1[:,4:20],axis=0)

	if M05==1:
        	sfdtot=np.array([tot[7],tot[10]])
        	tlook=np.array([xboxcen[7],xboxcen[10]])
        	f=interp1d(tlook,sfdtot)
        	sfdsk=np.zeros(2)
        	for i in range(0,2):
                	sfdsk[i]=f(xboxcen[i+8])
                	data1[:,4+i+8]=(((data1[:,4+i+8])/tot[i+8])*sfdsk[i])
        if BC==1:
                sfdtot=np.array([tot[9],tot[11]])
                tlook=np.array([xboxcen[9],xboxcen[11]])
                f=interp1d(tlook,sfdtot)
		sfdsk=f(xboxcen[10])
                data1[:,4+10]=(((data1[:,4+10])/tot[10])*sfdsk)
                sfdtot=np.array([tot[11],tot[13]])
                tlook=np.array([xboxcen[11],xboxcen[13]])
                f=interp1d(tlook,sfdtot)
		sfdsk=f(xboxcen[12])
                data1[:,4+12]=(((data1[:,4+12])/tot[12])*sfdsk)


plot=1
if plot==1:	
	#DAT1=np.dot((1/data1[:,23]),data1[:,4:20])
	DAT1=np.sum(data1[:,4:20],axis=0)/len(data1)
#	DAT2=np.sum(data2[:,4:20],axis=0)/len(data2)
        ybox1=np.hstack((DAT1.T,DAT1[-1]))
 #       ybox2=np.hstack((DAT2.T,DAT2[-1]))
	plt.step(xbox,ybox1,where='post',color='k',ls='-')
  #      plt.step(xbox,ybox2,where='post',color='k',ls='--')
	plt.yscale('log')
        plt.xlabel('lookback [Gyr]',fontsize=20) 
        plt.ylabel('<SFR> [M$_\odot$/yr/galaxy]',fontsize=20)
        plt.show()


np.savetxt('vespaBC_1_gals_SFH6HSsig2_dr7.dat',data1)
#np.savetxt('vespaM05_2_gals_SFH6HSsigsmoo_dr7.dat',data2)

