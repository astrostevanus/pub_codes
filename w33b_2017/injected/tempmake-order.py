import numpy as np
import h5py
from  tqdm import tqdm
from tqdm import tnrange
import os
import scipy.interpolate as sci
import scipy
import h5py
import astropy
import astropy.constants as const
from ultramodule import *
from exoparam import *
c = 299792.458
pi=3.141592653589793
######################################################################                                                                                                                                                  
specdir="./rainbow/"
exotempdir="/home/stev/stdaftersvd/spec/wasp33b/"

mol="tio"
vmr_min=4
vmr_max=11

whichdata= "este"
opt="-fin"
######################################################################                                                                                                                                                  
### Template Making ###                                                                                                                                                                                                 
######################################################################                                                                                                                                                  
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt)+".info", 'r')

drvs=h5f_techni["Doppler shifts variation"][:]


#ccd_list_2breduced=["blue","red"]
ccd_list_2breduced=["red"]
output_tmp= "rainbow-wasp33"+str(whichdata)+"all_"+str(mol)+"_template"+str(opt)+".h5"
h5f_temp = h5py.File(specdir+"template/"+output_tmp, 'w')
vmr=vmr_min
for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    for vmr in tqdm (range(vmr_min,vmr_max), ascii="True",desc="VMR 10^"+str(vmr)):
        xtemp,ytemp=np.loadtxt(exotempdir+"spec_wasp33_"+str(mol)+"_inv"+str(vmr)+".dat",unpack=True,dtype="float")
        ytemp, fwhm = pyasl.instrBroadGaussFast(xtemp,ytemp, 100000, edgeHandling="firstlast", fullout=True, maxsig=5.0)
        Ip=np.array(intensity(xtemp*1.e-10,3610.))
        ytempnorm=ytemp/Ip
        h5f_eqwv = h5py.File(specdir+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
        for order in tqdm(range (1,len_order),ascii="True",desc=str(ccd)+" ccd"):
            wv_order= h5f_eqwv[str(ccd)+"-wv-order-"+str(order)][:]
            template_exo = np.zeros(shape=(len(drvs),len(wv_order)),dtype="float") 
            for i in range(len(drvs)):
                fi = sci.interp1d(xtemp*(1.0 + drvs[i]/c), ytempnorm) #Doppler shifting the template at drvs[i] km/s
                template_exo[i,:]=fi(wv_order) #Extracting the Doppler shifted template with the same wavelength of the observed spectrum in that order and put it in template_exo[i,:]
            h5f_temp.create_dataset(str(ccd)+"-vmr-"+str(vmr)+"-order-"+str(order), data=template_exo)#Writing the template of the specific vmr
    print "Wrote it in "+str(output_tmp)
h5f_temp.close()
h5f_techni.close()
