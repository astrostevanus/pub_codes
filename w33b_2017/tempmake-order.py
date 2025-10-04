import numpy as np
import h5py
from random import shuffle
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
exotempdir="/home/stev/fin-pycovfefe/spec/"

mol="tio"
vmr_min=65
vmr_max=85

whichdata= "kawa"
opt="-fin"

print "Which data: "+str(whichdata)
print "Param option:" +str(opt)
print "Mol:" +str(mol)
print "log VMR= -"+str(vmr_min)+"'till -"+str(vmr_max)
######################################################################                                                                                                                                                  
### Template Making ###                                                                                                                                                                                                 
######################################################################                                                                                                                                                  
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt)+".info", 'r')
kp=h5f_techni["kp"][:]
vsys=h5f_techni["vsys"][:]
print "Kp: "+str(min(kp))+" to "+str(max(kp))+" km.s-1"
print "Vsys: "+str(min(vsys))+" to "+str(max(vsys))+" km.s-1" 
drvs=h5f_techni["Doppler shifts variation"][:]


ccd_list_2breduced=["blue","red"]
#ccd_list_2breduced=["red"]
output_tmp= "rainbow-wasp33"+str(whichdata)+"allhalf_"+str(mol)+"_template"+str(opt)+".h5"
print output_tmp
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
    for vmr in tqdm (range(vmr_min,vmr_max,10), ascii="True",desc="VMR 10^"+str(vmr)):
#        shuf_order=np.arange(0,len_order-1,1)
#        shuffle(shuf_order)  spec_wasp33_tio_all_inv7.dat spec-wasp_tio_allinv7.dat 
        xtemp,ytemp=np.loadtxt(exotempdir+"spec_wasp33_"+str(mol)+"_inv"+str(vmr)+".dat",unpack=True,dtype="float")
#        ytemp, fwhm = pyasl.instrBroadGaussFast(xtemp,ytemp, 100000, edgeHandling="firstlast", fullout=True, maxsig=5.0)
        Ip=np.array(intensity(xtemp*1.e-10,3610.))
        ytempnorm=ytemp/Ip
        h5f_eqwv = h5py.File(specdir+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
        for order in tqdm(range (1,len_order),ascii="True",desc=str(ccd)+" ccd"):
            wv_order= h5f_eqwv[str(ccd)+"-wv-order-"+str(order)][:]
#            wv_shuf_sample=h5f_eqwv[str(ccd)+"-wv-order-"+str(shuf_order[order-1]+1)][:]
#            wv_spc=wv_shuf_sample[1]-wv_shuf_sample[0]
#            new_wv_shuf=np.arange(wv_shuf_sample[0],wv_shuf_sample[0]+wv_spc*len(wv_order),wv_spc)
#            new_wv_shuf=np.linspace(wv_shuf_sample[0],wv_shuf_sample[0]+wv_spc*(len(wv_order)-1),len(wv_order),endpoint=True)
            template_exo = np.zeros(shape=(len(drvs),len(wv_order)),dtype="float") 
#            moddrvs=np.arange(drvs[0]-21.,drvs[-1]-21.+drvs[1]-drvs[0],drvs[1]-drvs[0])
#            print len(moddrvs)
#            print len(drvs)
            for i in range(len(drvs)):
#                if ccd=="blue":
#                    if order==1 or order==2 or order==3:
#                        print moddrvs
#                        fi = sci.interp1d(xtemp*(1.0 + (moddrvs[i])/c), ytempnorm) #Doppler shifting the template at drvs[i] km/s
#                    else:
#                        fi = sci.interp1d(xtemp*(1.0 + (drvs[i])/c), ytempnorm)
#                else:
                fi = sci.interp1d(xtemp*(1.0 + (drvs[i])/c), ytempnorm)
                #sm1 = pyasl.smooth(fi(wv_order), 501, 'flat')
                #sm2 = pyasl.smooth(sm1, 201, 'flat')
                
#                template_exo[i,:]=fi(new_wv_shuf) #Extracting the Doppler shifted template with the same wavelength of the observed spectrum in that order and put it in template_exo[i,:]
                template_exo[i,:]=fi(wv_order)
#            h5f_temp.create_dataset(str(ccd)+"-vmr-"+str(vmr)+"-order-"+str(order), data=template_exo)#Writing the template of the specific vmr
            h5f_temp.create_dataset(str(ccd)+"-vmr-"+str(vmr)+"-order-"+str(order), data=template_exo)
    print "Wrote it in "+str(output_tmp)
h5f_temp.close()
h5f_techni.close()
