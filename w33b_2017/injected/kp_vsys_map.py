import numpy as np
import h5py
from tqdm import tqdm
import scipy.interpolate as sci
specdir="./rainbow/"
red_mode="svd"
vmr_min= 4
vmr_max= 11
len_sv=15
whichdata="kawa"
opt="-fin"
##################################################################
# Importing technical parameters from techni.info and wasp33.info
##################################################################
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt)+".info", 'r')
drvs=h5f_techni["Doppler shifts variation"][:]
kp=h5f_techni["kp"][:]
vsys=h5f_techni["vsys"][:]
phase_john=h5f_techni["phase_john"][:]
out_ecl_ph_fr=h5f_techni["out of eclipse phase"][:]
ingress_ph_fr=h5f_techni["in ingress phase"][:]
mjd,bjd,rvcor,hjd,airmass=np.loadtxt("wasp33-"+str(whichdata)+".info",unpack=True,dtype="float")

##################################################################
# only out of transit frames = len(out_ecl_ph_fr)
# out of transit and in ingress frames= len(out_ecl_ph_fr)+ len(ingress_ph_fr)
# in transit = len(ecl_ph_fr)
len_frames=len(out_ecl_ph_fr)+len(ingress_ph_fr) 


##################################################################
input_cc="wasp33"+str(whichdata)+"_tioccSA"+str(opt)+".h5"
output_SA="rainbow-wasp33"+(whichdata)+"_tioccSA"+str(opt)+".h5"
output_KpVsys= "rainbow-wasp33"+(whichdata)+"_tioccSA_john_ingress"+str(opt)+".h5"
##################################################################

if red_mode=="svd":
    varname="sv"
    len_var=15
else:
    varname="systematics"
    len_var=15
if whichdata=="este":
    ccd_list=["red"]
else:
    ccd_list=["blue","red"]

total_cc_SA = h5py.File(specdir+"post_"+red_mode+"/sum_cc/"+output_SA, "r")
Kp_Vsys_SA = h5py.File(specdir+"post_"+red_mode+"/kp_vsys_result/"+output_KpVsys, 'w')
for ccd in ccd_list:
    print str(ccd)
    for vmr in tqdm (range(vmr_min,vmr_max),ascii="True",desc="For Various VMR n Sys"):
        for var in tqdm (range(len_sv),ascii="True",desc="VMR= 10^-"+str(vmr)+" Sys"):
            norm_mat_cc_SA=total_cc_SA[str(ccd)+"_SA_cc_vmr_"+str(vmr)+str(varname)+"_"+str(var)][:]
            cc_Kp_Vsys=np.zeros((len(kp),len(vsys)))
            rv_int_function=[]
            for filenum in range(len_frames):
                rv_int_function.append(sci.interp1d(drvs,norm_mat_cc_SA[filenum],kind="cubic"))
            for i_kp in range(len(kp)):
                for j_vsys in range(len(vsys)):
                    cc_drv=np.zeros(len_frames)
                    for filenum in range(len_frames):
                        v_predicted=kp[i_kp]*np.sin(2.*np.pi*phase_john[filenum])+vsys[j_vsys]-rvcor[filenum] #Positive means red shifted, predicted observed RV
                        cc_drv[filenum]=rv_int_function[filenum](v_predicted)
                    cc_Kp_Vsys[i_kp,j_vsys]=np.mean(cc_drv)
            Kp_Vsys_SA.create_dataset(str(ccd)+"vmr-"+str(vmr)+"_"+str(varname)+"-"+str(var), data=cc_Kp_Vsys)
total_cc_SA.close()
Kp_Vsys_SA.close()
h5f_techni.close()
print output_KpVsys
