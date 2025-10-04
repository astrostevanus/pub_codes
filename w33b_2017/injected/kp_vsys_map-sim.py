import numpy as np
import h5py
from tqdm import tqdm
import scipy.interpolate as sci
specdir="./post_sysrem/"
red_mode="sysrem"
len_sv=10
whichdata="kawa"
vari="addmul" #add or rem
opt_tech="-fin"
#opt="tio_finalmodel-7addlow5far2ex"
opt="tio-modx9-atre-cor2tempuncor"

##################################################################
# Importing technical parameters from techni.info and wasp33.info
##################################################################
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt_tech)+".info", 'r')
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

#"rainbow-wasp33"+str(whichdata)+"_cc_pearson-SA"+str(opt)+".h5"
##################################################################
output_SA="rainbow-wasp33"+str(whichdata)+"_cc_pearson-SA"+str(opt)+".h5"
output_KpVsys= "rainbow-wasp33"+str(whichdata)+"_cc_pearson-SA_john_ingress"+str(opt)+".h5"

##################################################################
if vari=="addmul":
    varimin=1
    varimax=6
    variname="addmul"
elif vari=="rem":
    varimin=1
    varimax=6
    variname="mul"


if red_mode=="svd":
    varname="sv"
    len_var=15
else:
    varname="systematics"
    len_var=10
if whichdata=="este":
    ccd_list=["red"]
else:
    ccd_list=["rainbow","mean"]

total_cc_SA = h5py.File(specdir+"sum_cc/"+output_SA, "r")
Kp_Vsys_SA = h5py.File(specdir+"kp_vsys_result/"+output_KpVsys, 'w')
for ccd in ccd_list:
    print str(ccd)
    for va in tqdm (range(varimin,varimax),ascii="True",desc=str(variname)+" n "+str(varname)):
        for var in tqdm (range(len_sv),ascii="True",desc=str(variname)+"= "+str(va)+" "+str(varname)):
            norm_mat_cc_SA=total_cc_SA[str(ccd)+"_SA_cc_"+str(variname)+"_"+str(va)+str(varname)+"_"+str(var)][:]
            cc_Kp_Vsys=np.zeros((len(kp),len(vsys)))
            rv_int_function=[]
            for filenum in range(len_frames):
                rv_int_function.append(sci.interp1d(drvs,norm_mat_cc_SA[filenum]))
            for i_kp in range(len(kp)):
                for j_vsys in range(len(vsys)):
                    cc_drv=np.zeros(len_frames)
                    for filenum in range(len_frames):
                        v_predicted=kp[i_kp]*np.sin(2.*np.pi*phase_john[filenum])+vsys[j_vsys]-rvcor[filenum] #Positive means red shifted, predicted observed RV
                        cc_drv[filenum]=rv_int_function[filenum](v_predicted)
                    cc_Kp_Vsys[i_kp,j_vsys]=np.mean(cc_drv)
            Kp_Vsys_SA.create_dataset(str(ccd)+str(variname)+"-"+str(va)+"_"+str(varname)+"-"+str(var), data=cc_Kp_Vsys)
total_cc_SA.close()
Kp_Vsys_SA.close()
h5f_techni.close()
print output_KpVsys
