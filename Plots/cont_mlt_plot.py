import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable



from scipy.stats import chisqprob
from scipy.stats import chi2
import math
from scipy.stats.mstats import gmean
import copy
from matplotlib.lines import Line2D

def pois(lam,k):
	p=(lam**k*np.exp(-lam))/math.factorial(k)
	return p;
dets1=10

TYPE='vespaBC_1_gals_SFH6HSsig2'
TY='vespaBC_1_gals_SFH6HSsig'

TYPE1='{}_dr7_SM2binR'.format(TYPE)
loglike11=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)
likecont11=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)
loglike21=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)
likecont21=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)
loglike31=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)
likecont31=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE1),unpack=False)


TYPE2='{}_dr7_MR2binR'.format(TYPE)
loglike12=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)
likecont12=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)
loglike22=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)
likecont22=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)
loglike32=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)
likecont32=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE2),unpack=False)


TYPE3='{}_dr7_GR2binR'.format(TYPE)
loglike13=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)
likecont13=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)
loglike23=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)
likecont23=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)
loglike33=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)
likecont33=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE3),unpack=False)


TYPE4='{}_dr7_CMD2binR'.format(TYPE)
loglike14=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)
likecont14=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)
loglike24=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)
likecont24=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)
loglike34=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)
likecont34=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE4),unpack=False)


TYPE5='{}_dr7_sSFR2binR'.format(TYPE)
loglike15=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)
likecont15=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)
loglike25=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)
likecont25=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)
loglike35=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)
likecont35=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE5),unpack=False)


TYPE6='{}_dr7R'.format(TYPE)
loglike16=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)
likecont16=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)
loglike26=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)
likecont26=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)
loglike36=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)
likecont36=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets1,TYPE6),unpack=False)


dets2=100

TYPE7='{}_dr7_SM2binR'.format(TYPE)
loglike17=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)
likecont17=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)
loglike27=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)
likecont27=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)
loglike37=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)
likecont37=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE7),unpack=False)


TYPE8='{}_dr7_MR2binR'.format(TYPE)
loglike18=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)
likecont18=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)
loglike28=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)
likecont28=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)
loglike38=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)
likecont38=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE8),unpack=False)


TYPE9='{}_dr7_GR2binR'.format(TYPE)
loglike19=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)
likecont19=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)
loglike29=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)
likecont29=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)
loglike39=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)
likecont39=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE9),unpack=False)


TYPE10='{}_dr7_CMD2binR'.format(TYPE)
loglike110=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)
likecont110=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)
loglike210=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)
likecont210=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)
loglike310=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)
likecont310=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE10),unpack=False)

TYPE11='{}_dr7_sSFR2binR'.format(TYPE)
loglike111=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)
likecont111=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)
loglike211=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)
likecont211=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)
loglike311=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)
likecont311=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE11),unpack=False)


TYPE12='{}_dr7R'.format(TYPE)
loglike112=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)
likecont112=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)
loglike212=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)
likecont212=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)
loglike312=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)
likecont312=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets2,TYPE12),unpack=False)

dets3=1000

TYPE13='{}_dr7_SM2binR'.format(TYPE)
loglike113=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)
likecont113=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)
loglike213=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)
likecont213=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)
loglike313=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)
likecont313=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE13),unpack=False)


TYPE14='{}_dr7_MR2binR'.format(TYPE)
loglike114=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)
likecont114=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)
loglike214=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)
likecont214=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)
loglike314=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)
likecont314=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE14),unpack=False)

TYPE15='{}_dr7_GR2binR'.format(TYPE)
loglike115=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)
likecont115=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)
loglike215=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)
likecont215=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)
loglike315=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)
likecont315=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE15),unpack=False)

TYPE16='{}_dr7_CMD2binR'.format(TYPE)
loglike116=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)
likecont116=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)
loglike216=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)
likecont216=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)
loglike316=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)
likecont316=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE16),unpack=False)

TYPE17='{}_dr7_sSFR2binR'.format(TYPE)
loglike117=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)
likecont117=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)
loglike217=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)
likecont217=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)
loglike317=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)
likecont317=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE17),unpack=False)

TYPE18='{}_dr7R'.format(TYPE)
loglike118=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike1_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)
likecont118=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf1_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)
loglike218=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike2_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)
likecont218=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf2_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)
loglike318=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/loglike3_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)
likecont318=np.loadtxt('/uufs/astro.utah.edu/common/uuastro/astro_data/zhengzheng/mccarthy/GravBinary_Data/likecontf3_{0}_{1}.dat'.format(dets3,TYPE18),unpack=False)

#loglike1[np.isinf(loglike1)]=loglike1[~np.isinf(loglike1)].min()
#loglike2[np.isinf(loglike2)]=loglike2[~np.isinf(loglike2)].min()
#loglike3[np.isinf(loglike3)]=loglike3[~np.isinf(loglike3)].min()
#loglike4[np.isinf(loglike4)]=loglike4[~np.isinf(loglike4)].min()
#loglike5[np.isinf(loglike5)]=loglike5[~np.isinf(loglike5)].min()
#loglike6[np.isinf(loglike6)]=loglike6[~np.isinf(loglike6)].min()
#loglike7[np.isinf(loglike7)]=loglike7[~np.isinf(loglike7)].min()
#loglike8[np.isinf(loglike8)]=loglike8[~np.isinf(loglike8)].min()
#loglike9[np.isinf(loglike9)]=loglike9[~np.isinf(loglike9)].min()
#loglike10[np.isinf(loglike10)]=loglike10[~np.isinf(loglike10)].min()
#loglike11[np.isinf(loglike11)]=loglike11[~np.isinf(loglike11)].min()
#loglike12[np.isinf(loglike12)]=loglike12[~np.isinf(loglike12)].min()
#loglike13[np.isinf(loglike13)]=loglike13[~np.isinf(loglike13)].min()
#loglike14[np.isinf(loglike14)]=loglike14[~np.isinf(loglike14)].min()
#loglike15[np.isinf(loglike15)]=loglike15[~np.isinf(loglike15)].min()
#loglike16[np.isinf(loglike16)]=loglike16[~np.isinf(loglike16)].min()
#loglike17[np.isinf(loglike17)]=loglike17[~np.isinf(loglike17)].min()
#loglike18[np.isinf(loglike18)]=loglike18[~np.isinf(loglike18)].min()

nobs=[-1.5,-1.1,-0.5]
Nnobs=len(nobs)
tminobs=[0.01,0.035,1.0]
Ntminobs=len(tminobs)
dnmod1=0.05
dnmod2=0.01

nmod11=np.around(np.arange(-2.0,-1.7,dnmod1),decimals=2)
nmod12=np.around(np.arange(-1.7+dnmod2,-1.3+dnmod2,dnmod2),decimals=2)
nmod13=np.around(np.arange(-1.3+dnmod1,0.0+dnmod1,dnmod1),decimals=2)
nmod1=np.hstack((nmod11,nmod12,nmod13))

nmod21=np.around(np.arange(-2.0,-1.3,dnmod1),decimals=2)
nmod22=np.around(np.arange(-1.3,-0.9+dnmod2,dnmod2),decimals=2)
nmod23=np.around(np.arange(-0.9+dnmod1,0.0+dnmod1,dnmod1),decimals=2)
nmod2=np.hstack((nmod21,nmod22,nmod23))

nmod31=np.around(np.arange(-2.0,-0.7,dnmod1),decimals=2)
nmod32=np.around(np.arange(-0.7,-0.3+dnmod2,dnmod2),decimals=2)
nmod33=np.around(np.arange(-0.3+dnmod1,0.0+dnmod1,dnmod1),decimals=2)
nmod3=np.hstack((nmod31,nmod32,nmod33))

dtmin1=0.1
dtmin2=0.01
tminmod11=np.arange(-2.7,-2.2,dtmin1)
tminmod12=np.arange(-2.2,-1.8,dtmin2)
tminmod13=np.arange(-1.8+dtmin1, 0.7+dtmin1,dtmin1)
tminmod1=np.hstack((tminmod11,tminmod12,tminmod13))
tminmod1=np.power(10,tminmod1)


tminmod21=np.arange(-2.7,-1.66,dtmin1)
tminmod22=np.arange(-1.66,-1.26,dtmin2)
tminmod23=np.arange(-1.26+dtmin1, 0.7+dtmin1,dtmin1)
tminmod2=np.hstack((tminmod21,tminmod22,tminmod23))
tminmod2[-1]=0.7
tminmod2=np.power(10,tminmod2)


tminmod31=np.arange(-2.7,-0.2,dtmin1)
tminmod32=np.arange(-0.2,0.2,dtmin2)
tminmod33=np.arange(0.2+dtmin1, 0.7+dtmin1,dtmin1)
tminmod3=np.hstack((tminmod31,tminmod32,tminmod33))
tminmod3=np.power(10,tminmod3)


Nnmod1=len(nmod1)
Ntminmod1=len(tminmod1)
Nnmod2=len(nmod2)
Ntminmod2=len(tminmod2)
Nnmod3=len(nmod3)
Ntminmod3=len(tminmod3)

NNmod1=Nnmod1*Ntminmod1
NNmod2=Nnmod2*Ntminmod2
NNmod3=Nnmod3*Ntminmod3

nmodf=np.around(np.arange(-2.0,0.0+dnmod2,dnmod2),decimals=2)
tminmodf=np.arange(-2.7,0.7+dtmin2,dtmin2)
Nnmodf=len(nmodf)-1
Ntminmodf=len(tminmodf)-1
NNmodf=Nnmodf*Ntminmodf


nmodf_wid=np.zeros(Nnmodf)
nmodf_cen=np.zeros(Nnmodf)
for i in range(0,Nnmodf):
        nmodf_wid[i]=nmodf[i+1]-nmodf[i]
        nmodf_cen[i]=nmodf[i]+nmodf_wid[i]/2.

tminmodf_wid=np.zeros(Ntminmodf)
tminmodf_cen=np.zeros(Ntminmodf)
for i in range(0,Ntminmodf):
        tminmodf_wid[i]=tminmodf[i+1]-tminmodf[i]
        tminmodf_cen[i]=tminmodf[i]+tminmodf_wid[i]/2.

q=0
modelf=np.zeros((NNmodf,3))
for i in range(Nnmodf):
        for j in range(Ntminmodf):
                modelf[q,0]=nmodf_cen[i]
                modelf[q,1]=tminmodf_cen[j]
                q=q+1



nmod1_wid=np.zeros(Nnmod1-1)
nmod1_cen=np.zeros(Nnmod1-1)
for i in range(0,Nnmod1-1):
        nmod1_wid[i]=nmod1[i+1]-nmod1[i]
        nmod1_cen[i]=nmod1[i]+nmod1_wid[i]/2.

tminmod1_wid=np.zeros(Ntminmod1-1)
tminmod1_cen=np.zeros(Ntminmod1-1)
for i in range(0,Ntminmod1-1):
        tminmod1_wid[i]=tminmod1[i+1]-tminmod1[i]
        tminmod1_cen[i]=tminmod1[i]+tminmod1_wid[i]/2.


nmod2_wid=np.zeros(Nnmod2-1)
nmod2_cen=np.zeros(Nnmod2-1)
for i in range(0,Nnmod2-1):
        nmod2_wid[i]=nmod2[i+1]-nmod2[i]
        nmod2_cen[i]=nmod2[i]+nmod2_wid[i]/2.

tminmod2_wid=np.zeros(Ntminmod2-1)
tminmod2_cen=np.zeros(Ntminmod2-1)
for i in range(0,Ntminmod2-1):
        tminmod2_wid[i]=tminmod2[i+1]-tminmod2[i]
        tminmod2_cen[i]=tminmod2[i]+tminmod2_wid[i]/2.

nmod3_wid=np.zeros(Nnmod3-1)
nmod3_cen=np.zeros(Nnmod3-1)
for i in range(0,Nnmod3-1):
        nmod3_wid[i]=nmod3[i+1]-nmod3[i]
        nmod3_cen[i]=nmod3[i]+nmod3_wid[i]/2.

tminmod3_wid=np.zeros(Ntminmod3-1)
tminmod3_cen=np.zeros(Ntminmod3-1)
for i in range(0,Ntminmod3-1):
        tminmod3_wid[i]=tminmod3[i+1]-tminmod3[i]
        tminmod3_cen[i]=tminmod3[i]+tminmod3_wid[i]/2.


NNobs=3


plotlik=1
if plotlik==1:
	tminmod1=np.log10(tminmod1)
        tminmod2=np.log10(tminmod2)
        tminmod3=np.log10(tminmod3)
	tminobs=np.log10(tminobs)
	tminmod1_cen=np.log10(tminmod1_cen)
        tminmod2_cen=np.log10(tminmod2_cen)
        tminmod3_cen=np.log10(tminmod3_cen)


        fig, axs = plt.subplots(3, 3, sharex=True,sharey=True)
        fig.subplots_adjust(wspace=0,hspace=0)
        (ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)=axs

	one=0
	two=1
	three=2
	rr=3


        ax1.contour(tminmod1_cen,nmod1_cen,loglike12.reshape(Nnmod1-1,Ntminmod1-1),[likecont12[1]],linewidths=rr,colors='b')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike22.reshape(Nnmod2-1,Ntminmod2-1),[likecont22[1]],linewidths=rr,colors='b')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike32.reshape(Nnmod3-1,Ntminmod3-1),[likecont32[1]],linewidths=rr,colors='b')
        ax1.contour(tminmod1_cen,nmod1_cen,loglike13.reshape(Nnmod1-1,Ntminmod1-1),[likecont13[1]],linewidths=rr,colors='r')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike23.reshape(Nnmod2-1,Ntminmod2-1),[likecont23[1]],linewidths=rr,colors='r')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike33.reshape(Nnmod3-1,Ntminmod3-1),[likecont33[1]],linewidths=rr,colors='r')
        ax1.contour(tminmod1_cen,nmod1_cen,loglike14.reshape(Nnmod1-1,Ntminmod1-1),[likecont14[1]],linewidths=rr,colors='magenta')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike24.reshape(Nnmod2-1,Ntminmod2-1),[likecont24[1]],linewidths=rr,colors='magenta')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike34.reshape(Nnmod3-1,Ntminmod3-1),[likecont34[1]],linewidths=rr,colors='magenta')
        ax1.contour(tminmod1_cen,nmod1_cen,loglike11.reshape(Nnmod1-1,Ntminmod1-1),[likecont11[1]],linewidths=rr,colors='cyan')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike21.reshape(Nnmod2-1,Ntminmod2-1),[likecont21[1]],linewidths=rr,colors='cyan')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike31.reshape(Nnmod3-1,Ntminmod3-1),[likecont31[1]],linewidths=rr,colors='cyan')
	ax1.contour(tminmod1_cen,nmod1_cen,loglike15.reshape(Nnmod1-1,Ntminmod1-1),[likecont15[1]],linewidths=rr,colors='orange')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike25.reshape(Nnmod2-1,Ntminmod2-1),[likecont25[1]],linewidths=rr,colors='orange')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike35.reshape(Nnmod3-1,Ntminmod3-1),[likecont35[1]],linewidths=rr,colors='orange')
        ax1.contour(tminmod1_cen,nmod1_cen,loglike16.reshape(Nnmod1-1,Ntminmod1-1),[likecont16[1]],linewidths=rr,colors='k')
        ax4.contour(tminmod2_cen,nmod2_cen,loglike26.reshape(Nnmod2-1,Ntminmod2-1),[likecont26[1]],linewidths=rr,colors='k')
        ax7.contour(tminmod3_cen,nmod3_cen,loglike36.reshape(Nnmod3-1,Ntminmod3-1),[likecont36[1]],linewidths=rr,colors='k')

        
        ax2.contour(tminmod1_cen,nmod1_cen,loglike18.reshape(Nnmod1-1,Ntminmod1-1),[likecont18[1]],linewidths=rr,colors='b')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike28.reshape(Nnmod2-1,Ntminmod2-1),[likecont28[1]],linewidths=rr,colors='b')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike38.reshape(Nnmod3-1,Ntminmod3-1),[likecont38[1]],linewidths=rr,colors='b')
        ax2.contour(tminmod1_cen,nmod1_cen,loglike19.reshape(Nnmod1-1,Ntminmod1-1),[likecont19[1]],linewidths=rr,colors='r')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike29.reshape(Nnmod2-1,Ntminmod2-1),[likecont29[1]],linewidths=rr,colors='r')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike39.reshape(Nnmod3-1,Ntminmod3-1),[likecont39[1]],linewidths=rr,colors='r')
        ax2.contour(tminmod1_cen,nmod1_cen,loglike110.reshape(Nnmod1-1,Ntminmod1-1),[likecont110[1]],linewidths=rr,colors='magenta')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike210.reshape(Nnmod2-1,Ntminmod2-1),[likecont210[1]],linewidths=rr,colors='magenta')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike310.reshape(Nnmod3-1,Ntminmod3-1),[likecont310[1]],linewidths=rr,colors='magenta')
        ax2.contour(tminmod1_cen,nmod1_cen,loglike17.reshape(Nnmod1-1,Ntminmod1-1),[likecont17[1]],linewidths=rr,colors='cyan')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike27.reshape(Nnmod2-1,Ntminmod2-1),[likecont27[1]],linewidths=rr,colors='cyan')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike37.reshape(Nnmod3-1,Ntminmod3-1),[likecont37[1]],linewidths=rr,colors='cyan')
	ax2.contour(tminmod1_cen,nmod1_cen,loglike111.reshape(Nnmod1-1,Ntminmod1-1),[likecont111[1]],linewidths=rr,colors='orange')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike211.reshape(Nnmod2-1,Ntminmod2-1),[likecont211[1]],linewidths=rr,colors='orange')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike311.reshape(Nnmod3-1,Ntminmod3-1),[likecont311[1]],linewidths=rr,colors='orange')
        ax2.contour(tminmod1_cen,nmod1_cen,loglike112.reshape(Nnmod1-1,Ntminmod1-1),[likecont112[1]],linewidths=rr,colors='k')
        ax5.contour(tminmod2_cen,nmod2_cen,loglike212.reshape(Nnmod2-1,Ntminmod2-1),[likecont212[1]],linewidths=rr,colors='k')
        ax8.contour(tminmod3_cen,nmod3_cen,loglike312.reshape(Nnmod3-1,Ntminmod3-1),[likecont312[1]],linewidths=rr,colors='k')

        
        ax3.contour(tminmod1_cen,nmod1_cen,loglike114.reshape(Nnmod1-1,Ntminmod1-1),[likecont114[1]],linewidths=rr,colors='b')
        ax6.contour(tminmod2_cen,nmod2_cen,loglike214.reshape(Nnmod2-1,Ntminmod2-1),[likecont214[1]],linewidths=rr,colors='b')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike314.reshape(Nnmod3-1,Ntminmod3-1),[likecont314[1]],linewidths=rr,colors='b')
        ax3.contour(tminmod1_cen,nmod1_cen,loglike115.reshape(Nnmod1-1,Ntminmod1-1),[likecont115[1]],linewidths=rr,colors='r')
        ax6.contour(tminmod2_cen,nmod2_cen,loglike215.reshape(Nnmod2-1,Ntminmod2-1),[likecont215[1]],linewidths=rr,colors='r')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike315.reshape(Nnmod3-1,Ntminmod3-1),[likecont315[1]],linewidths=rr,colors='r')
        ax3.contour(tminmod1_cen,nmod1_cen,loglike116.reshape(Nnmod1-1,Ntminmod1-1),[likecont116[1]],linewidths=rr,colors='magenta')
        ax6.contour(tminmod2_cen,nmod2_cen,loglike216.reshape(Nnmod2-1,Ntminmod2-1),[likecont216[1]],linewidths=rr,colors='magenta')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike316.reshape(Nnmod3-1,Ntminmod3-1),[likecont316[1]],linewidths=rr,colors='magenta')
        ax3.contour(tminmod1_cen,nmod1_cen,loglike113.reshape(Nnmod1-1,Ntminmod1-1),[likecont113[1]],linewidths=rr,colors='cyan')
        ax6.contour(tminmod2_cen,nmod2_cen,loglike213.reshape(Nnmod2-1,Ntminmod2-1),[likecont213[1]],linewidths=rr,colors='cyan')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike313.reshape(Nnmod3-1,Ntminmod3-1),[likecont313[1]],linewidths=rr,colors='cyan')
	ax3.contour(tminmod1_cen,nmod1_cen,loglike117.reshape(Nnmod1-1,Ntminmod1-1),[likecont117[1]],linewidths=rr,colors='orange')
	ax6.contour(tminmod2_cen,nmod2_cen,loglike217.reshape(Nnmod2-1,Ntminmod2-1),[likecont217[1]],linewidths=rr,colors='orange')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike317.reshape(Nnmod3-1,Ntminmod3-1),[likecont317[1]],linewidths=rr,colors='orange')
        ax3.contour(tminmod1_cen,nmod1_cen,loglike118.reshape(Nnmod1-1,Ntminmod1-1),[likecont118[1]],linewidths=rr,colors='k')
        ax6.contour(tminmod2_cen,nmod2_cen,loglike218.reshape(Nnmod2-1,Ntminmod2-1),[likecont218[1]],linewidths=rr,colors='k')
        ax9.contour(tminmod3_cen,nmod3_cen,loglike318.reshape(Nnmod3-1,Ntminmod3-1),[likecont318[1]],linewidths=rr,colors='k')

 
	ax1.set_ylabel(r'n',fontsize=20)
       # ax1.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
       # ax2.set_ylabel(r'n',fontsize=20,weight='bold')
       # ax2.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
       # ax3.set_ylabel(r'n',fontsize=20,weight='bold')
       # ax3.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
        ax4.set_ylabel(r'n',fontsize=20)
       # ax4.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
       # ax5.set_ylabel(r'n',fontsize=20,weight='bold')
       # ax5.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
       # ax6.set_ylabel(r'n',fontsize=20,weight='bold')
       # ax6.set_xlabel(r'log[t$_{min}$/Myr]',fontsize=20,weight='bold')
        ax7.set_ylabel(r'n',fontsize=20)
        ax7.set_xlabel(r'log(t$_{\rm m}$/Gyr)',fontsize=20)
       # ax8.set_ylabel(r'n',fontsize=20,weight='bold')
        ax8.set_xlabel(r'log(t$_{\rm m}$/Gyr)',fontsize=20)
       # ax9.set_ylabel(r'n',fontsize=20,weight='bold')
        ax9.set_xlabel(r'log(t$_{\rm m}$/Gyr)',fontsize=20) 
        ax1.plot(tminobs[0],nobs[0],'o',markersize=5,color='y')
        ax4.plot(tminobs[1],nobs[1],'o',markersize=5,color='y') 
        ax7.plot(tminobs[2],nobs[2],'o',markersize=5,color='y')
        ax2.plot(tminobs[0],nobs[0],'o',markersize=5,color='y')
        ax5.plot(tminobs[1],nobs[1],'o',markersize=5,color='y')
        ax8.plot(tminobs[2],nobs[2],'o',markersize=5,color='y')
        ax3.plot(tminobs[0],nobs[0],'o',markersize=5,color='y')
        ax6.plot(tminobs[1],nobs[1],'o',markersize=5,color='y')
        ax9.plot(tminobs[2],nobs[2],'o',markersize=5,color='y')

        ax1.margins(x=0)
        ax1.margins(y=0)
        ax2.margins(x=0)
        ax2.margins(y=0)
        ax3.margins(x=0)
        ax3.margins(y=0)
        ax4.margins(x=0)
        ax4.margins(y=0)
        ax5.margins(x=0)
        ax5.margins(y=0)
        ax6.margins(x=0)
        ax6.margins(y=0)
        ax7.margins(x=0)
        ax7.margins(y=0)
        ax8.margins(x=0)
        ax8.margins(y=0)
        ax9.margins(x=0)
        ax9.margins(y=0) 

	ax1.tick_params(axis="y", labelsize=20)
	ax1.tick_params(axis="x", labelsize=20)
	ax2.tick_params(axis="y", labelsize=20)
	ax2.tick_params(axis="x", labelsize=20)
	ax3.tick_params(axis="y", labelsize=20)
	ax3.tick_params(axis="x", labelsize=20)
	ax4.tick_params(axis="y", labelsize=20)
	ax4.tick_params(axis="x", labelsize=20)
	ax5.tick_params(axis="y", labelsize=20)
	ax5.tick_params(axis="x", labelsize=20)
	ax6.tick_params(axis="y", labelsize=20)
	ax6.tick_params(axis="x", labelsize=20)
        ax7.tick_params(axis="y", labelsize=20)
        ax7.tick_params(axis="x", labelsize=20)
        ax8.tick_params(axis="y", labelsize=20)
        ax8.tick_params(axis="x", labelsize=20)
        ax9.tick_params(axis="y", labelsize=20)
        ax9.tick_params(axis="x", labelsize=20)

	ax7.set_xticks(np.arange(-2.5,1.5,1.0))
	ax8.set_xticks(np.arange(-2.5,1.5,1.0))
	ax9.set_xticks(np.arange(-2.5,1.5,1.0))

        ax1.set_yticks(np.arange(-2.0,0.5,0.5))
        ax4.set_yticks(np.arange(-2.0,0.5,0.5))
        ax7.set_yticks(np.arange(-2.0,0.5,0.5))

#	for tick in ax7.get_xticklabels():
#	        tick.set_rotation(90)
#	for tick in ax8.get_xticklabels():
#        	tick.set_rotation(90)
#	for tick in ax9.get_xticklabels():
#        	tick.set_rotation(90)

	ax1.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax1.tick_params('both', length=10, width=2, which='minor')
	ax2.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax2.tick_params('both', length=10, width=2, which='minor')
	ax3.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax3.tick_params('both', length=10, width=2, which='minor')
	ax4.tick_params('both', length=15, width=2, which='major',labelsize=20)
	ax4.tick_params('both', length=10, width=2, which='minor')
        ax5.tick_params('both', length=15, width=2, which='major',labelsize=20)
        ax5.tick_params('both', length=10, width=2, which='minor')
        ax6.tick_params('both', length=15, width=2, which='major',labelsize=20)
        ax6.tick_params('both', length=10, width=2, which='minor')
        ax7.tick_params('both', length=15, width=2, which='major',labelsize=20)
        ax7.tick_params('both', length=10, width=2, which='minor')
        ax8.tick_params('both', length=15, width=2, which='major',labelsize=20)
        ax8.tick_params('both', length=10, width=2, which='minor')
        ax9.tick_params('both', length=15, width=2, which='major',labelsize=20)
        ax9.tick_params('both', length=10, width=2, which='minor')
	custom_lines = [Line2D([0], [0], color='b', lw=2),Line2D([0], [0], color='r', lw=2),Line2D([0], [0], color='magenta', lw=2),Line2D([0], [0], color='cyan', lw=2),Line2D([0], [0], color='orange', lw=2),Line2D([0], [0], color='k', lw=2)]
	ax3.legend(custom_lines, ['M$_r$','g-r','(M$_r$,g-r)','M$_*$','sSFR','perGAL'],loc='upper right',fontsize=20)
	ax1.text(-2.45,-0.30,r'N$_{\rm obs}$=10',fontsize=20)
        ax2.text(-2.45,-0.30,r'N$_{\rm obs}$=100',fontsize=20)
        ax3.text(-2.45,-0.30,r'N$_{\rm obs}$=1000',fontsize=20)
        plt.show()
        

