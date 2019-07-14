import h5py
import numpy as np
from tqdm import tqdm
#vmrinj=9
vmr_min=7
vmr_max=8
#segment=10.
red_mode="sysrem"
mol="tio"
whichdata="kawa"
opt_tech="-fin"
specdir="./rainbow/"
#red-wasp33kawa_vo_cc_pearson-uncor-fin.h5
#input_cc="wasp33_tioinj"+str(vmrinj)+"_refarRV_segment_cc.h5"
#red-wasp33kawa_tiocc_pearson-usual-fin.h5
#/blue-wasp33kawa_tio_allinv_cc_pearson-uncor-fin.h5
#-wasp33kawa_tiocc_pearson--finultihalf-superorder.h
#red-wasp33kawa_tioallinv_cc_pearson--finultihalf-superorder.h5#
#input_cc="wasp33kawa_tiocc_pearson--fin-superorder-hp.h5"
input_cc="wasp33kawa_tionew_allnoinv_cc_pearson--fin-superorder-hp.h5"

#input_cc="wasp33kawa_tiohay_cc_pearson--fin-superorder-hp.h5"
#input_cc="wasp33kawa_tionew_allinv_cc_pearson--fin-superorder-hp.h5"
#input_cc="wasp33kawa_tiocc_pearson--finultihalf-usualorder.h5"
#input_cc="wasp33kawa_tio_cc_pearson--finultihalf-fulbre.h5"
#input_cc="wasp33"+str(whichdata)+"_"+str(mol)+"_allinv_cc_pearson-uncor"+str(opt_tech)+".h5"
#input_cc="wasp33"+str(whichdata)+"_"+str(mol)+"cc"+str(opt)+".h5"
#input_cc="rainbow-wasp33kawaall_vo_segment-10.0pix_cc-fin.h5"
#input_cc="rainbow-wasp33kawaall_tio_segment-10.0pix_cc-fin.h5"
#input_cc="wasp33"+str(whichdata)+"_tiocc"+str(opt)+".h5"
#input_cc="wasp33"+str(whichdata)+"_"+str(mol)+"cc_spearman"+str(opt)+"test.h5"
#input_cc="rainbow-wasp33"+str(whichdata)+"all_"+str(mol)+"_segment-"+str(segment)+"pix_cc"+str(opt_tech)+".h5"

output_SA="rainbow-"+input_cc[:-3]+"_SA-masked.h5"
total_cc_SA_rainbow = h5py.File(specdir+"post_"+red_mode+"/sum_cc/"+output_SA, 'w')

h5f_cc_blue=h5py.File(specdir+"post_"+red_mode+"/cc_result/blue-"+input_cc,"r")
h5f_cc_red=h5py.File(specdir+"post_"+red_mode+"/cc_result/red-"+input_cc,"r")
#h5f_cc=h5py.File(specdir+"post_"+str(red_mode)+"/cc_result/"+input_cc,"r")

if red_mode=="svd":
    varname="sv"
    len_var=10
else:
    varname="systematics"
    len_var=10

for vmr in tqdm (range(vmr_min,vmr_max),ascii="True",desc="VMR"):
    blue_usual= h5f_cc_blue["blue_cc_vmr_"+str(vmr)][:]
    red_usual= h5f_cc_red["red_cc_vmr_"+str(vmr)][:]
#    bluebl_sv_col=[]
#    for sv in range (len(blue_usual)):
#        bluebl_order_col=[]
#        for order in range (0,18):
#            bluebl_order=np.zeros((len(blue_usual[0][0]),len(blue_usual[1][1][1][42:-42])))
#            for filenum in range (len(blue_usual[0][0])):
#                if order==0 or order==1 or order==2:
#                    bluebl_order[filenum]=blue_usual[sv][order][filenum][:-84]
#               else:
#                    bluebl_order[filenum]=blue_usual[sv][order][filenum][42:-42]
#            bluebl_order_col.append(bluebl_order)
#        bluebl_sv_col.append(bluebl_order_col)
        
#    redbl_sv_col=[]
#    for sv in range (len(red_usual)):
#        redbl_order_col=[]
#        for order in range (0,12):
#            redbl_order=np.zeros((len(red_usual[0][0]),len(red_usual[1][1][1][42:-42])))
#            for filenum in range (len(red_usual[0][0])):
#                    redbl_order[filenum]=red_usual[sv][order][filenum][42:-42]
#            redbl_order_col.append(redbl_order)
#        redbl_sv_col.append(redbl_order_col)

    for var in tqdm(range(len_var),ascii="True",desc="R"+str(varname)):
#        SA_blue= np.sum(bluebl_sv_col[var],axis=0)
#        SA_red= np.sum(redbl_sv_col[var],axis=0)
        blue_ccf=[]
        for i in range (len(blue_usual[var])):
            if i==0 or i== 1 or i==2 or i==17:
                0
            else:
                blue_ccf.append(blue_usual[var][i])
        red_ccf=[]
        for i in range (len(red_usual[var])):
            if i==0 or i==1  or i== 6 or i ==7 or i==11:
                0
            else:
                red_ccf.append(red_usual[var][i])
        SA_blue=np.sum(blue_ccf,axis=0)
        SA_red=np.sum(red_ccf,axis=0)
#        SA_blue= np.sum(blue_usual[var],axis=0)
#        SA_red= np.sum(red_usual[var],axis=0)
#        SA_blue_norm=(SA_blue.transpose()-np.mean(SA_blue,axis=1)).transpose()
#        SA_red_norm=(SA_red.transpose()-np.mean(SA_red,axis=1)).transpose()
        SA_blue_norm=(SA_blue.transpose()-np.mean(SA_blue,axis=1)).transpose()                                                         
        SA_red_norm=(SA_red.transpose()-np.mean(SA_red,axis=1)).transpose() 
#        SA_mean=(np.mean(bluebl_sv_col[var],axis=0)+np.mean(redbl_sv_col[var],axis=0))/2.
        SA_mean=(np.mean(blue_usual[var],axis=0)+np.mean(red_usual[var],axis=0))/2.
        SA_total=SA_blue_norm+SA_red_norm
        SA_total_mean_norm=(SA_mean.transpose()-np.mean(SA_mean,axis=1)).transpose()
#        total_cc_SA_rainbow.create_dataset("blue_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_blue_norm)#Writing the template of the specific vmr
#        total_cc_SA_rainbow.create_dataset("red_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_red_norm)#Writing the template of the specific vmr
        total_cc_SA_rainbow.create_dataset("rainbow_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_total)
#        print np.shape(SA_total_mean)
        total_cc_SA_rainbow.create_dataset("mean_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var), data=SA_total_mean_norm)        #Writing the template of the specific vmr
total_cc_SA_rainbow.close()
h5f_cc_blue.close()
h5f_cc_red.close()
#h5f_cc.close()
print output_SA
