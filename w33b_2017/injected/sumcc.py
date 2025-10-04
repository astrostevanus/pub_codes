import h5py
import numpy as np
from tqdm import tqdm
#vmrinj=9
vmr_min=4
vmr_max=11
red_mode="svd"

whichdata="kawa"
opt="-fin"
specdir="./rainbow/"
#input_cc="wasp33_tioinj"+str(vmrinj)+"_refarRV_segment_cc.h5"
#input_cc="wasp33_vo_test_segment_cc.h5"

input_cc="wasp33"+str(whichdata)+"_tiocc"+str(opt)+".h5"
output_SA="rainbow-wasp33"+str(whichdata)+"_tioccSA"+str(opt)+".h5"

total_cc_SA_rainbow = h5py.File(specdir+"post_"+red_mode+"/sum_cc/"+output_SA, 'w')
h5f_cc_blue=h5py.File(specdir+"post_"+red_mode+"/cc_result/blue-"+input_cc,"r")
h5f_cc_red=h5py.File(specdir+"post_"+red_mode+"/cc_result/red-"+input_cc,"r")

if red_mode=="svd":
    varname="sv"
    len_var=15
else:
    varname="systematics"
    len_var=15

for vmr in tqdm (range(vmr_min,vmr_max),ascii="True",desc="VMR"):
    cc_map_full_blue_SA= h5f_cc_blue["blue_cc_vmr_"+str(vmr)][:]
    cc_map_full_red_SA= h5f_cc_red["red_cc_vmr_"+str(vmr)][:]
    for var in tqdm(range(len_var),ascii="True",desc="R"+str(varname)):
        SA_blue=np.mean(cc_map_full_blue_SA[var],axis=0)
        SA_red=np.mean(cc_map_full_red_SA[var],axis=0)
#        SA_mean=(np.sum(cc_map_full_blue_SA[var],axis=0)+np.sum(cc_map_full_red_SA[var],axis=0))/30.
#        SA_total_mean=(SA_blue+SA_red)/2.
        SA_blue_norm=(SA_blue.transpose()-np.mean(SA_blue,axis=1)).transpose()
        SA_red_norm=(SA_red.transpose()-np.mean(SA_red,axis=1)).transpose()
#        SA_total_mean_norm=(SA_total_mean.transpose()-np.mean(SA_total_mean,axis=1)).transpose()
        total_cc_SA_rainbow.create_dataset("blue_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_blue_norm)#Writing the template of the specific vmr
        total_cc_SA_rainbow.create_dataset("red_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_red_norm)#Writing the template of the specific vmr
#        total_cc_SA_rainbow.create_dataset("rainbow_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_mean)
#        total_cc_SA_rainbow.create_dataset("mean_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_total_mean_norm)        #Writing the template of the specific vmr
total_cc_SA_rainbow.close()
h5f_cc_blue.close()
h5f_cc_red.close()
print output_SA
