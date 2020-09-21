import numpy as np
import pyfits
import matplotlib.pyplot as plt


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
v1=np.loadtxt('vespaM05_1_dr7_gal.dat',unpack=False)
v2=np.loadtxt('vespaM05_2_dr7_gal.dat',unpack=False)

SDSS1=1
SDSS2=1
xumosaic=0
xubin=0
if SDSS1==1:
        y07_sdss_ID=np.loadtxt('SDSS7_ID.dat',unpack=False)

        y07_sdss=np.loadtxt('SDSS7.dat',unpack=False)
        y07_sdss_NYU=y07_sdss[:,1].astype(int)
        y07_sdss_ra=y07_sdss[:,2]
        y07_sdss_dec=y07_sdss[:,3]
        y07_sdss_z=y07_sdss[:,4]
        y07_sdss_Mr=y07_sdss[:,8]
        y07_sdss_gr=y07_sdss[:,9]
	y07_sdss_rlim=y07_sdss[:,6]


        Ngal=len(y07_sdss[:,0])
        y07_sdss_mjd=np.zeros(Ngal)
        y07_sdss_plate=np.zeros(Ngal)
        y07_sdss_fib=np.zeros(Ngal)
        for i in range(0,Ngal):
                if y07_sdss_ID[i,1]==y07_sdss[i,1]:
                        y07_sdss_mjd[i]=y07_sdss_ID[i,2]
                        y07_sdss_plate[i]=y07_sdss_ID[i,3]
                        y07_sdss_fib[i]=y07_sdss_ID[i,4]


if SDSS2==1:
	# WMAP5 cosmology as used by VESPA
	Olam=0.726
	Om=1-Olam
	h=0.705


	# removed due to bad pixels from deblending, cosmic rays, bleed trails, etc that where interpolated over
	#outsRA=np.array([186.3177,195.1226,181.262,158.877,160.641,166.574,197.4458,120.6555,340.45544,130.49335,170.4634,175.36584,156.0873,205.8393,202.47016,41.902536,29.7265,24.5377,26.6083,126.4852])
	#outsDEC=np.array([-0.4118,0.5988,0.959,-0.4387,-0.037,-0.109,23.3844,40.3529,-9.5470,27.59766,40.8630,11.89559,10.52864,4.63482,41.84246,-0.37038,0.05016,0.12796,0.93717,47.60725])

	outsRA,outsDEC=np.loadtxt('out6sigma_SFH_M05_1_dr7_ksm.dat',unpack='True')
        #in1,in2,ds=spherematch(outsRA1,outsDEC1,outsRA,outsDEC,tol=0.00027778)


	ind1,ind2,ds=spherematch(y07_sdss_ra,y07_sdss_dec,outsRA,outsDEC,tol=0.00027778)
	y07_sdss_NYU[ind1]=-1.
	cut=y07_sdss_NYU>=0.0
	Mr=y07_sdss_Mr[cut]
	rlim=y07_sdss_rlim[cut]
	ra=y07_sdss_ra[cut]
	dec=y07_sdss_dec[cut]
	gr=y07_sdss_gr[cut]
	z=y07_sdss_z[cut]

	TOTsa=7984.*(np.pi/180.)**2
	Dmax=np.power(10,(rlim-Mr)/5+1)/(10**6) # max luminosity distance with flux limit in Mpc
	Vmax=TOTsa*Dmax**3/3.

	n1v1,n2v1,ds=spherematch(ra,dec,v1[:,0],v1[:,1],tol=0.00027778)
        n1v2,n2v2,ds=spherematch(ra,dec,v2[:,0],v2[:,1],tol=0.00027778)

	h=1
	mstar=-0.406+1.097*gr-0.4*(Mr-5*np.log10(h)-4.64)
	
#	dataV1=np.vstack((Mr,gr,ra,dec,z,Vmax)).T
        dataV1=np.vstack((Mr[n1v1],gr[n1v1],ra[n1v1],dec[n1v1],z[n1v1],Vmax[n1v1],mstar[n1v1],v1[n2v1,3:].T)).T
	dataV2=np.vstack((Mr[n1v2],gr[n1v2],ra[n1v2],dec[n1v2],z[n1v2],Vmax[n1v2],mstar[n1v2],v2[n2v2,3:].T)).T

        dataV=np.vstack((Mr[n1v1],gr[n1v1],ra[n1v1],dec[n1v1],z[n1v1],Vmax[n1v1],mstar[n1v1],v1[n2v1,3:].T,v2[n2v2,3:].T)).T
        np.savetxt('vespaSFH6sig_M05_dr7_KSM.dat',dataV)

	#np.savetxt('dr7_KSM.dat',dataV1)
        #np.savetxt('vespaSFH6sig_M05_1_dr7_KSM.dat',dataV1)
	#np.savetxt('vespaSFH6sig_M05_2_dr7_KSM.dat',dataV2)


if xumosaic==1:
	L=79
	perc=np.zeros(L)
	for i in range(0,L):
		TYPE='voldr72brightMosaic{}.dat'.format(i)
		#(0)ra,(1)dec,(2)redshift,(3)comving distance,(4)sector,(5)vmax,(6)fgot, (7)jackknifeID(Haojie's way)
		data=np.loadtxt(TYPE,unpack=False)
		n1v1,n2v1,x=spherematch(v1[:,0],v1[:,1],data[:,0],data[:,1],tol=0.00027778)
		n1v2,n2v2,x=spherematch(v2[:,0],v2[:,1],data[:,0],data[:,1],tol=0.00027778)
		N1=len(data[n2v1,0])
		N2=len(data[n2v2,0])
		if ~np.any(data[n2v1,0]!=data[n2v2,0]):
			perc[i]=float(N1)/len(data[:,0])
			data1=np.vstack((data[n2v1,:].T,v1[n1v1,3:].T,v2[n1v2,3:].T)).T

		else:
			print('Error')
			break
		TYPE='voldr72brightMosaic{}_vespaM05.dat'.format(i)
		np.savetxt(TYPE,data1)

if xubin==1:
	data=np.loadtxt('voldr72brightMosaic_18.00_18.25.dat',unpack=False)




