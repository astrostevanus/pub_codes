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
exotempdir="./spec/wasp33b/"

#len_sv=10
mol="vo"
segment=10.
vmr_min=5
vmr_max=12
whichdata="kawa"
opt_tech="-fin"
######################################################################
### Template Making ###
######################################################################
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt_tech)+".info", 'r')

drvs=h5f_techni["Doppler shifts variation"][:]


ccd_list_2breduced=["blue","red"]

output_tmp= "rainbow-wasp33"+str(whichdata)+"all_"+str(mol)+"_segment-"+str(segment)+"pix_template"+str(opt_tech)+".h5"
h5f_temp = h5py.File(specdir+"template/"+output_tmp, 'w')
vmr=vmr_min
for vmr in tqdm (range(vmr_min,vmr_max), ascii="True",desc="VMR 10^"+str(vmr)):
    xtemp,ytemp=np.loadtxt(exotempdir+"spec_wasp33_"+str(mol)+"_inv"+str(vmr)+".dat",unpack=True,dtype="float")
    Ip=np.array(intensity(xtemp*1.e-10,3610.))
    ytempnorm=ytemp/Ip
    for ccd in ccd_list_2breduced:
        if ccd=="blue":
            len_order=19
        else:
            len_order=13
        h5f_eqwv = h5py.File(specdir+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
        for order in tqdm(range (1,len_order),ascii="True",desc=str(ccd)+" ccd"):
            wv_order= h5f_eqwv[str(ccd)+"-wv-order-"+str(order)][:]
            for seg_num in range (int(segment)):
                low_lim=int(((seg_num)/segment)*len(wv_order))
                up_lim=int(((seg_num+1.)/segment)*len(wv_order))
                wvsegment=wv_order[low_lim:up_lim]
                masked_temp=(xtemp>wv_order[low_lim]-50.)*(xtemp<wv_order[up_lim-1]+50.)
                masked_xtemp=xtemp[masked_temp]
                masked_ytemp=ytemp[masked_temp]
                template_exo = np.zeros(shape=(len(drvs),len(wvsegment)),dtype="float")
                for i in range(len(drvs)):
                    fi = sci.interp1d(masked_xtemp*(1.0 + drvs[i]/c), masked_ytemp)
                    template_exo[i,:]=fi(wvsegment)
                h5f_temp.create_dataset(str(ccd)+"vmr-"+str(vmr)+"-order-"+str(order)+"-segment-"+str(seg_num), data=template_exo)
    print "Wrote it in "+str(output_tmp)
h5f_temp.close()
h5f_techni.close()
