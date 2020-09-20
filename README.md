# Paper2_DTDofBNSM https://arxiv.org/abs/2007.15024
Constraining Delay Time Distribution of Binary Neutron Star Mergers from Host Galaxy Properties 
Kevin Spencer McCarthy, Zheng Zheng, and Enrico Ramirez-Ruiz
Submitted MNRAS, MN-20-3071-MJ
----------------------------------------------------------------------------

This project involved the likelihood calculations of particular delay-time distributions (DTD) of binary neutron star mergers (BNSM) models of 3 simulated universes; a slow DTD, a canonical DTD, and a fast DTD. An 'observation' of the truth models (slow, canonical, and fast) are performed and the possible DTDs are constrained by producing a prediction for the BNSM event rate PDF through a convolution of the star formation history of a galaxy (SFH) or a galaxy property bin (<SFH>) with a power-law form of the DTD with two free variables, the power law index (n) and the minimum delay time (t_min) for an event to occur following the formation of stars from gas, following the procedure employed by Zheng & Ramirez_Ruiz (2007). We also evaluate how these constraints improve with an increase in # of BNSM multi-messenger event detections (N_obs). We employ the SFHs (VESPA; Tojeiro R. et al. 2009) of ~500k galaxies of the SDSS DR7 main galaxy catalog. 
 
I provide likelihood constraints on these three truth models considering every galaxy in the catalog (convolving each individual galaxy SFH with the potential DTD) and when considering the expected SFH from a particular host galaxy property bin (<SFH(gal prop.)>), evaluating the constraining behavior of 5 different host galaxy properties: r-band absolute magnitude (M_r), color (g-r), color-magnitude (M_r,g-r), stellar mass (M_*), and specific star formation rate (sSFR). I also provide the figure-of-merit scores for the different constraining methods, and estimates of the errors, both as a function of N_obs.

I show that there is a simple way to calculate the likelihood across DTD parameters that scales with N_obs. This relationship shown in Eq. 9 and 10 of the publication amounts to the expected likelihood being equal to the multiple of N_obs with the relative entropy of the model and observed BNSM event rate PDFs. This arises from considering a Poisson distribution of events in a particular galaxy/galaxy property bin and observing the ensemble average of the N_obs events expected from a particular truth model. This method can also be used in real-world observations of BNSM event rates, once enough detections have been made. Currently, there is only one BNSM GW event with an electromagnetic counterpart detection; GW170817 occurring in NGC4993.

I find that the strongest constraints on DTD occur when we consider the SFH of individual galaxies (perGAL), seeing as we do not marginalize over any detail of the SFH. When marginalizing over the galaxies properties to produce a rate prediction as a function of galaxy property (M_r,g-r,(M_r,g-r),M_*,sSFR), I find that all properties besides M_r perform similarly, i.e., the weakest constraints come from considering the luminosity of the galaxies. We also see that the constraint contours depend on the truth model. In other words, some DTDs may be more difficult to constrain than others given the details of the SFH measurement and inherent degeneracies of the DTD functional form. I find that it is possible to differentiate between physically meaningful DTD truth models with perGAL after O(10) detections. To constrain a DTD with 10% accuracy with perGAL, order O(100) detections are necessary. With the current BNSM GW event rate density estimations from the LIGO-VIRGO collaboration (250–2810 Gpc−3yr−1; Abbott et al. 2020), this could occur within 2-40 years.



CODE FILES
-----------------------------------------------------------
Fileprep_Galaxies



Plots
- fig1.py (for 2PCF in linear rp,r_pi 2D space)
- fig2.py (for 2PCF (wp,xi_0,xi_2,xi_4) measurement as function of rp or s)
- fig3.py (for fsig_8 estimator)
- fig5.py (for chi^2 and n_g values from HOD/SCAM fits)
 


