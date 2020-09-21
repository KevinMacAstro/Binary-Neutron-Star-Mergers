import pyfits
import numpy as np

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

# Saves number density in stellar mass bin calculated from number of objects in MPA catalog in particular bin via 1/Vmax method, limits of SM bin, and average SFR in each bin of available volume limited VESPA objects.

# (0)  ra, (1) dec, (2) d, (3:19) VESPA SFR, (19) Mr, (20) g-r, (21) Dmax, (22) Vamx, (23) Mstar, (24) sSFR
TYPE='SFH6HSsig2'
data1=np.loadtxt('vespaBC_1_gals_{}_dr7R.dat'.format(TYPE),unpack=False)
#data2=np.loadtxt('vespaM05_2_gals_{}_dr7.dat'.format(TYPE),unpack=False)



sm=np.arange(7.5,11.6,0.1)
digit1=np.digitize(data1[:,24],sm)
#digit2=np.digitize(data2[:,24],sm)


SFH1=np.zeros((40,16))
#SFH2=np.zeros((13,16))
ng1=np.zeros(40)
#ng2=np.zeros(13)
SMrange=np.zeros((40,2))
for i in range(1,41):
	dataTMP1=data1[digit1==i]
	ng1[i-1]=np.sum(1/dataTMP1[:,23])	
	dmax1=dataTMP1[:,22][dataTMP1[:,20]==dataTMP1[:,20].max()]
	dataTMP1=dataTMP1[dataTMP1[:,3]<dmax1]
#	dataTMP2=data2[digit2==i]
#        ng2[i-1]=np.sum(1/dataTMP2[:,23])
#        dmax2=dataTMP2[:,22][dataTMP2[:,20]==dataTMP2[:,20].max()]
#        dataTMP2=dataTMP2[dataTMP2[:,3]<dmax2]
	SFH1[i-1,:]=np.sum(dataTMP1[:,4:20],axis=0)/len(dataTMP1)
#	SFH2[i-1,:]=np.sum(dataTMP2[:,4:20],axis=0)/len(dataTMP2)
	SMrange[i-1,0]=sm[i-1]
	SMrange[i-1,1]=sm[i]

DATA1=np.vstack((ng1,SMrange.T,SFH1.T)).T
#DATA2=np.vstack((ng2,SMrange.T,SFH2.T)).T

np.savetxt('vespaBC_1_gals_{}_dr7_SM2binR.dat'.format(TYPE),DATA1)
#np.savetxt('vespaM05_2_gals_{}_dr7_SMbin.dat'.format(TYPE),DATA2)

