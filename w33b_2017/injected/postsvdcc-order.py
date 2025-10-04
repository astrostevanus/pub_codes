import numpy as np
from ultramodule import divstd_prime, crossco, cc_multiorde
import h5py
from  tqdm import tqdm
from tqdm import tnrange
###################################################################################################                                                                            
len_sv=15
vmr_min=4
vmr_max=11
whichdata="kawa"
opt="-fin"
mol="tio"
specdir="./rainbow/"

###################################################################################################                                                                            


ccd_list_2breduced=["blue","red"]
#ccd_list_2breduced=["red"]
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt)+".info", 'r')
drvs=h5f_techni["Doppler shifts variation"][:]

output_tmp="rainbow-wasp33"+str(whichdata)+"all_"+str(mol)+"_template"+str(opt)+".h5"
h5f_temp = h5py.File(specdir+"template/"+output_tmp, 'r')
for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    output_cc=str(ccd)+"-wasp33"+str(whichdata)+"_"+str(mol)+"cc"+str(opt)+".h5" 
    h5f_cc = h5py.File(specdir+"post_svd/cc_result/"+str(output_cc), 'w')
    h5f_reduced = h5py.File(specdir+"post_svd/"+str(ccd)+str(whichdata)+"_postsvd_fullmatrix.h5", 'r')         
    for vmr in tqdm(range(vmr_min,vmr_max),ascii="True",desc=str(ccd)+" CCD"):
        cc_map_full=[]
        for sv in tqdm(range(len_sv),ascii="True",desc="VMR= 10^-"+str(vmr)+" RSV"):
            cc_svd_collect=[] #Different SV are saved in this matrix                                                                                                       
            for order in tqdm (range(1,len_order),ascii="True",desc="Order"):
                y_temp_order=h5f_temp[str(ccd)+"-vmr-"+str(vmr)+"-order-"+str(order)][:]
                yobs= h5f_reduced[str(ccd)+"-flux-order-"+str(order)+"sv-"+str(sv)][:]                                                                                        
#                yobs= h5f_reduced[str(ccd)+"_vmr"+str(vmrinj)+"-flux-order-"+str(order)+"sv-"+str(sv)][:]
                std_frames=[]
                std_pix=[]
                for wvbin in range (len(yobs[1])):
                    std_pix.append(np.std(yobs[:,wvbin]))
                for filenum in range (len(yobs)):
                    std_frames.append(np.std(yobs[filenum,:]))
                std_frames=np.array(std_frames,dtype="float")
                yobs=divstd_prime(yobs,std_pix,std_frames)
                cc_order=np.zeros((len(yobs),len(drvs)),dtype="float")
                for numspec in range(len(yobs)):
                    for rv in range(len(drvs)):
                        cc_order[numspec][rv]=crossco(yobs[numspec],y_temp_order[rv])
                cc_svd_collect.append(cc_order)
            cc_map_full.append(cc_svd_collect)
        h5f_cc.create_dataset(str(ccd)+"_cc_vmr_"+str(vmr), data=cc_map_full)#Writing the template of the specific vmr
    h5f_cc.close()
    h5f_reduced.close()
