# Paper2_DTDofBNSM https://arxiv.org/abs/2007.15024
Constraining Delay Time Distribution of Binary Neutron Star Mergers from Host Galaxy Properties 
Kevin Spencer McCarthy, Zheng Zheng, and Enrico Ramirez-Ruiz
Submitted MNRAS, MN-20-3071-MJ
----------------------------------------------------------------------------

This project involved the likelihood calculations of particular delay-time distributions (DTD) of binary neutron star mergers (BNSM) models of 3 simulated universes; a slow DTD, a canonical DTD, and a fast DTD. An 'observation' of the truth models (slow, canonical, and fast) are performed and the possible DTDs are constrainted by producing a prediction for the BNSM event rate PDF through a convolution of the star formation history of a galaxy (SFH) or a galaxy property bin (<SFH>) with a power-law form of the DTD with two free variables, the power law index (n) and the minimum delay time (t_min) for an event to occur following the formation of stars from gas, following the procedure employed by Zheng & Ramirez_Ruiz (2007). We also evaluate how these constraints improve with an increase in # of BNSM multi-messenger event detections. We employ the SFHs (VESPA; Tojeiro R. et al. 2009) of ~500k galaxies of the SDSS DR7 main galaxy catalog. 
 
I provide likelihood constraints on these three truth models considering every galaxy in the catalog (convolving each individual galaxy SFH with the potential DTD) and when considering the expected SFH from a particular host galaxy property bin (<SFH(gal prop.)>), evaluating the constraining behavior of 5 different host galaxy properties: r-band absolute magnitude (M_r), color (g-r), color-magnitude (M_r,g-r), stellar mass (M_*), and specific star formation rate (sSFR).


CODE FILES
-----------------------------------------------------------
Fileprep_Galaxies



Plots
- fig1.py (for 2PCF in linear rp,r_pi 2D space)
- fig2.py (for 2PCF (wp,xi_0,xi_2,xi_4) measurement as function of rp or s)
- fig3.py (for fsig_8 estimator)
- fig5.py (for chi^2 and n_g values from HOD/SCAM fits)
 


