import numpy as np
import astropy.constants as const
from astropy.table import Table
import astropy.time as time
import astropy.coordinates as coords
import astropy.units as u
import h5py

whichdata="kawa"
h5f_techni = h5py.File("techni-"+str(whichdata)+"-finultihalf.info", 'w')
pi=3.141592653589793
mjd,bjd,rvcor,hjd,airmass=np.loadtxt("wasp33-15.info",unpack=True,dtype="float")


###########################################################################
####################  Exoplanet's parameters ##############################
t0=2456934.77146 #Johnson et al. 2015 BJD Ephemeris
t0_kovacs= 2452950.6724 #HJD Kovacs et al. 2013 HJD
P=1.2198709 #Kovacs et al. 2013 #Period

aRs=3.69 #semi major axis in stellar radius unit
Rs=1.509 #stellar radius in km
a_p=aRs*Rs*const.R_sun.value #semi major axis in km

t_trans=0.1143 #transit duration in days
t_ingress=0.0124 #Ingress duration in days

## Expected V orbital
vorb=2*np.pi*a_p/(P*24.*3600.)/1000. #orbital velocity of wasp 33 b
vstar=-2. #km/s #system radial velocity from cross correlation Collier Cameron et al. 2010
###########################################################################
###########################################################################



phase_john=((bjd+2400000.5)-t0)/P % 1 #phase of the frames
phase_kov=(hjd-t0_kovacs)/P % 1 #phase of the frames

last_ingress_ph= 0.5-(t_trans/2.-t_ingress)/P #phase of the last ingress
start_ingress_ph=0.5-(t_trans/P/2.) #phase of the beginning of ingress


out_ecl_ph_fr=[]
ingress_ph_fr=[]
ecl_ph_fr=[]

for ph in phase_john:
    if ph <start_ingress_ph:
        out_ecl_ph_fr.append(ph)
    elif ph <last_ingress_ph:
        ingress_ph_fr.append(ph)
    else:
        ecl_ph_fr.append(ph)

## Calculating cross correlation parameters


vwasp33b_obs=vorb*np.sin(2.*pi*phase_john)+vstar-rvcor #Positive means red shifted, predicted observed RV

v_width=150.

v_step=0.5

#Kp variation
kpmax=310.
kpmin=150.
kpstep= v_step

#V system variation
vsysmax=80.
vsysmin=-80.
vsysstep= v_step

kp= np.arange(kpmin,kpmax+kpstep,kpstep)
vsys= np.arange(vsysmin,vsysmax+vsysstep,vsysstep)

v1=kpmax*np.sin(2.*pi*phase_john)+(vsysmax)-rvcor #Positive means red shifted, predicted observed RV 
v2=kpmin*np.sin(2.*pi*phase_john)+(vsysmax)-rvcor #Positive means red shifted, predicted observed RV                                                            
v3=kpmax*np.sin(2.*pi*phase_john)+(vsysmin)-rvcor #Positive means red shifted, predicted observed RV                                                            
v4=kpmin*np.sin(2.*pi*phase_john)+(vsysmin)-rvcor #Positive means red shifted, predicted observed RV                                                          


drvs= np.arange(np.min([v1,v2,v3,v4])-30., np.max([v1,v2,v3,v4])+30.+v_step, v_step)


h5f_techni.create_dataset("phase_john", data=phase_john)
h5f_techni.create_dataset("out of eclipse phase", data=out_ecl_ph_fr)
h5f_techni.create_dataset("in ingress phase", data=ingress_ph_fr)
h5f_techni.create_dataset("in eclipse phase", data=ecl_ph_fr)

h5f_techni.create_dataset("kp", data=kp)
h5f_techni.create_dataset("vsys", data=vsys)
h5f_techni.create_dataset("Doppler shifts variation", data=drvs)
h5f_techni.create_dataset("Doppler shifts step", data=v_step)
h5f_techni.close()
