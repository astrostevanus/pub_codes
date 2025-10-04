import numpy as np
import astropy.units as u
import astropy.constants as const
from PyAstronomy import pyasl
import h5py
from astropy.convolution import convolve, Box1DKernel
from  tqdm import tqdm
import scipy.interpolate as sci
import pyfits
import scipy
from ultramodule import intensity

c = 299792.458
pi=3.141592653589793
specdir="/home/stev/fin-pycovfefe/rainbow/"
exotempdir="/home/stev/fin-pycovfefe/spec/wasp33b/"


#### Dynamical Param ###########
#############far################                                                                                              
vorb_obs=237.
vstar_obs=-1.5 #km/s
################################
#vorb_obs=239.5
#vstar_obs=-3.5 #km/s                                                                                                                                                           
################################

whichdata="kawa"
#rot_coef=0.4
inj_mode="addmul" #add or rem
#vmrrem="7"
vmrinj="8"
mol="tio"
#rotcoef=1.
injected_filename= "rainbow-wasp33"+str(whichdata)+"_"+str(mol)+"all-addmul-fineF.h5"

#### Exo Info ###
t0_john=2456934.77146 #Johnson et al. 2015 BJD                                                                                                         
p=1.2198709 # days Kovacs et al. 2013                                                                                                                                           
a_Rst=3.69 #semi major axis in R star                                                                                                                                           
Rp= 1.679 #Radius of the planet in jupiter radius                                                                                                                               
Rp_km=Rp*const.R_jup.value/1000. #Radius of the planet in km                                                                                                                    
Rs=1.509 #Radius of star in solar radius                                                                                                                                        
ap=a_Rst*Rs*const.R_sun.value/1000.
t_trans=0.1143 #days
t_ingress=0.0124 #days
################

info_inj = h5py.File(specdir+"injected/injection-"+str(whichdata)+".info", 'r')
bjd_start=info_inj["bjd_start"][:]
bjd_mid=info_inj["bjd_mid"][:]
bjd_last=info_inj["bjd_last"][:]
rvcor_start=info_inj["rvcor_start"][:]
rvcor_mid=info_inj["rvcor_mid"][:]
rvcor_last=info_inj["rvcor_last"][:]


phase_start=((bjd_start+2400000.5)-t0_john)/p % 1
phase_mid=((bjd_mid+2400000.5)-t0_john)/p % 1
phase_last=((bjd_last+2400000.5)-t0_john)/p % 1
    
v_start=vorb_obs*np.sin(2.*pi*phase_start)+vstar_obs-rvcor_start #Positive means red shifted, predicted observed RV
v_mid=vorb_obs*np.sin(2.*pi*phase_mid)+vstar_obs-rvcor_mid #Positive means red shifted, predicted observed RV
v_last=vorb_obs*np.sin(2.*pi*phase_last)+vstar_obs-rvcor_last #Positive means red shifted, predicted observed RV

drv_exps=v_last-v_start #Orbital-exposure time linked broadening
vp_rotation=2.*np.pi*Rp_km/(p*24.*3600.) #Rotational broadening (km/s)

h5f_techni = h5py.File(specdir+"injected/techni-"+str(whichdata)+".info", 'r')
out_ecl_ph_fr=h5f_techni["out of eclipse phase"][:]
ingress_ph_fr=h5f_techni["in ingress phase"][:]

len_contributed=len(out_ecl_ph_fr)+len(ingress_ph_fr)

rainbow_h5f_inj = h5py.File(specdir+"/injected/"+injected_filename, 'w')
ccd_list=["blue","red"]

if inj_mode=="add":
    varimin=4
    varimax=11
    vari_name="vmr"
    mod=1.
elif inj_mode=="rem":
    varimin=1
    varimax=6
    vari_name="mul"
    mod=-1.
elif inj_mode=="rotfit":
    varimin=1
    varimax=6
    vari_name="rot"
    mul=0.5
elif inj_mode=="addmul":
    varimin=1
    varimax=6
    vari_name="addmul"
    mod=1.
for vari in tqdm (range(varimin,varimax), ascii="True",desc="Doppler Shifting Template"):
    if inj_mode=="add":
        vmr=vari
        mul=1.
    elif inj_mode=="rem":
        vmr=vmrrem
        mul=vari
    elif inj_mode=="addmul":
        vmr=vmrinj
        rot_coef=np.float(0.4)
        mul= vari/5.
    elif inj_mode=="rotfit":
        rot_coef=np.float(vari)/5.
        vmr=vmrinj
        mod=1.
#        print rot_coef
    temp_filename= "spec_wasp33_"+str(mol)+"_inv"+str(vmr)+".dat"
    print temp_filename
    xtemp,ytemp=np.loadtxt(exotempdir+temp_filename,unpack=True,dtype="float") #Loading the template spectrum
    Is=np.array(intensity(xtemp*1.e-10,7400.)) #Calculating blackbody to normalize the template
    ytempnorm=(mul)*ytemp*((Rp*const.R_jup.value)/(Rs*const.R_sun.value))**2/Is #Scaling to the real contrast
    ################
    function_broad_dopp=[]
    broadened_spec=pyasl.fastRotBroad(xtemp,ytempnorm, 0.0, vp_rotation*rot_coef)
    for filenum in range (len(phase_mid)):
        #Doppler Shifting
        fi = sci.interp1d(xtemp*(1.0 + v_mid[filenum]/c), broadened_spec)
        function_broad_dopp.append(fi)    
    ################
#    print vp_rotation*rot_coef
    for ccd in ccd_list:
        if ccd=="blue":
            len_order=19
        else:
            len_order=13 
        h5f_eqwv = h5py.File(specdir+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
        for order in tqdm(range (1,len_order),ascii="True",desc=str(ccd)+" "+str(vari_name)+" = "+str(vari)):
            order_spec = h5f_eqwv[str(ccd)+"-flux-order-"+str(order)][:]
            order_spec=np.array(order_spec,dtype="float")
            wv_order= h5f_eqwv[str(ccd)+"-wv-order-"+str(order)][:]
            if ccd=="blue":
                if order==1 or order==2 or order==3:
                    wv_order=wv_order*(1.0 + -21./c)
            injected_model_order = np.zeros(shape=(len(order_spec),len(wv_order)),dtype="float")
            pix_size=wv_order[1]-wv_order[0]
            for filenum in range (len(order_spec)):
                if filenum < len_contributed:
                    boxsize=(np.abs(drv_exps[filenum])/(const.c.value/1000.))*np.median(wv_order)/(pix_size)
                    if boxsize >1.:
                        boxkernel = Box1DKernel(boxsize)
                        simulated_spec = convolve(function_broad_dopp[filenum](wv_order), boxkernel)
                        injected_spec= (1+mod*simulated_spec)*order_spec[filenum] 
                    else:
                        simulated_spec=function_broad_dopp[filenum](wv_order) 
                        injected_spec= (1+mod*simulated_spec)*order_spec[filenum]
                    injected_model_order[filenum]=injected_spec
                else:
                    injected_model_order[filenum]=order_spec[filenum]
            rainbow_h5f_inj.create_dataset(str(ccd)+"-flux-"+str(vari_name)+"-"+str(vari)+"-order-"+str(order), data=injected_model_order)
            rainbow_h5f_inj.create_dataset(str(ccd)+"-wv-"+str(vari_name)+"-"+str(vari)+"-order-"+str(order), data=wv_order)
print "Wrote it in "+str(injected_filename)
rainbow_h5f_inj.close()
h5f_eqwv.close()
