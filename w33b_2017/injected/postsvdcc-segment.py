import numpy as np
from ultramodule import divstd_prime, crossco, cc_multiorde
import h5py
from  tqdm import tqdm
from tqdm import tnrange
###################################################################################################
len_sv=10
vmr_min=4
vmr_max=11
segment=10.
mol="tio"
specdir="./rainbow/"
exotempdir="./spec/"
###################################################################################################

ccd_list_2breduced=["blue","red"]
h5f_techni = h5py.File("techni.info", 'r')
drvs=h5f_techni["Doppler shifts variation"][:]

output_tmp="rainbow-wasp33all_"+str(mol)+"_segment-"+str(segment)+"pix_template.h5"
vmrinj=9
h5f_temp = h5py.File(specdir+"template/"+output_tmp, 'r')
h5f_std = h5py.File(specdir+"Std_fullmatrix.h5", 'r')
h5f_reduced = h5py.File(specdir+"injected/post_svd/rainbow_postsvd_tio_full-vmr_hotspot-re_farRV.h5", 'r')
for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        len_order=13
    output_cc=str(ccd)+"-wasp33_tioinj"+str(vmrinj)+"_refarRV_segment_cc.h5"
    h5f_cc = h5py.File(specdir+"injected/post_svd/cc_result/"+str(output_cc), 'w')
#    h5f_reduced = h5py.File(specdir+"post_svd/"+str(ccd)+"_postsvd_fullmatrix.h5", 'r')
    for vmr in tqdm(range(vmr_min,vmr_max),ascii="True",desc=str(ccd)+" CCD"):
        cc_map_full_SA=[]
        for sv in tqdm(range(len_sv),ascii="True",desc="VMR= 10^-"+str(vmr)+" RSV"):
            cc_svd_collect_SA=[] #Different SV are saved in this matrix
            for order in tqdm (range(1,len_order),ascii="True",desc="Order"):
#                yobs= h5f_reduced[str(ccd)+"-flux-order-"+str(order)+"sv-"+str(sv)][:]
                yobs= h5f_reduced[str(ccd)+"_vmr"+str(vmrinj)+"-flux-order-"+str(order)+"sv-"+str(sv)][:]
                std_frames=[]
                std_pix=[]
                for wvbin in range (len(yobs[1])):
                    std_pix.append(np.std(yobs[:,wvbin]))
                for filenum in range (len(yobs)):
                    std_frames.append(np.std(yobs[filenum,:]))
                std_frames=np.array(std_frames,dtype="float")

                yobs=divstd_prime(yobs,std_pix,std_frames)

                cc_order_SA=np.zeros((len(yobs),len(drvs)),dtype="float")
                for numspec in range(len(yobs)):
                    cc_numspec=[]
                    for seg_num in range (int(segment)):
                        y_temp_order=h5f_temp[str(ccd)+"vmr-"+str(vmr)+"-order-"+str(order)+"-segment-"+str(seg_num)][:] #Loading template of vmr, v= -55., 355
                        low_lim=int(((seg_num)/segment)*len(yobs[0]))
                        up_lim=int(((seg_num+1.)/segment)*len(yobs[0]))

                        fluxsegment=yobs[numspec][low_lim:up_lim]

                        cc_each_seg=[]
                        for rv in range(len(drvs)):
                            cc_rv=crossco(fluxsegment,y_temp_order[rv])
                            cc_each_seg.append(cc_rv)
                        cc_numspec.append(np.array(cc_each_seg))
                    cc_numspec=np.array(cc_numspec)[~np.isnan(np.mean(cc_numspec,axis=1))]
                    cc_numspec_SA= np.mean(cc_numspec,axis=0)
                    cc_order_SA[numspec]=cc_numspec_SA
                cc_svd_collect_SA.append(cc_order_SA)
            cc_map_full_SA.append(cc_svd_collect_SA)
        h5f_cc.create_dataset("SA"+str(ccd)+"_cc_vmr_"+str(vmr), data= cc_map_full_SA)#Writing the template of the specific vmr
                
    h5f_cc.close()
h5f_reduced.close()
h5f_temp.close()
h5f_std.close()
h5f_techni.close()
