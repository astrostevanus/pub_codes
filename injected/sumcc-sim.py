import h5py
import numpy as np
from tqdm import tqdm

red_mode="sysrem"
vari="rot" #add or rem
whichdata="kawa"
opt_tech="-fin"
specdir="./post_"+str(red_mode)+"/"
varimin=1
varimax=10

#opt="tio_finalmodel-7addlow5far2ex"
#opt="tio-modx9-atre-cor2tempuncor"
opt="tio-9-farre-rotficor2tempuncor"
#rainbow-wasp33kawa_cc_pearsontio-modx9-atre-cor2tempuncor.h5
output_SA="rainbow-wasp33"+str(whichdata)+"_cc_pearson-SA"+str(opt)+".h5"

total_cc_SA_rainbow = h5py.File(specdir+"sum_cc/"+output_SA, 'w')
h5f_cc=h5py.File(specdir+"cc_result/rainbow-wasp33"+str(whichdata)+"_cc_pearson"+str(opt)+".h5","r")
#h5f_cc=h5py.File(specdir+"cc_result/rainbow-wasp33kawa_cc_pearsontio-modx9-atre-cor2tempuncor.h5","r")
#print total_cc_SA_rainbow
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
    len_var=10

else:
    varname="systematics"
    len_var=10

for va in tqdm (range(varimin,varimax),ascii="True",desc=str(variname)):
    cc_map_full_blue_SA= h5f_cc["blue_cc_"+str(variname)+"_"+str(va)][:]
    cc_map_full_red_SA= h5f_cc["red_cc_"+str(variname)+"_"+str(va)][:]
    for var in tqdm(range(len_var),ascii="True",desc="R"+str(varname)):
        SA_blue= np.sum(cc_map_full_blue_SA[var],axis=0)
        SA_red= np.sum(cc_map_full_red_SA[var],axis=0)
        SA_blue_norm=(SA_blue.transpose()-np.mean(SA_blue,axis=1)).transpose()
        SA_red_norm=(SA_red.transpose()-np.mean(SA_red,axis=1)).transpose()
        SA_mean=(np.mean(cc_map_full_blue_SA[var],axis=0)+np.mean(cc_map_full_red_SA[var],axis=0))/2.
        SA_total=SA_blue_norm+SA_red_norm
        SA_total_mean_norm=(SA_mean.transpose()-np.mean(SA_mean,axis=1)).transpose()
#        total_cc_SA_rainbow.create_dataset("blue_SA_cc_"+str(variname)+"_"+str(va)+str(varname)+"_"+str(var), data=SA_blue_norm)#Writing the template of the specific vmr
#        total_cc_SA_rainbow.create_dataset("red_SA_cc_"+str(variname)+"_"+str(va)+str(varname)+"_"+str(var), data=SA_red_norm)#Writing the template of the specific vmr
        total_cc_SA_rainbow.create_dataset("rainbow_SA_cc_"+str(variname)+"_"+str(va)+str(varname)+"_"+str(var), data=SA_total)
        total_cc_SA_rainbow.create_dataset("mean_SA_cc_"+str(variname)+"_"+str(va)+str(varname)+"_"+str(var), data=SA_total_mean_norm)
total_cc_SA_rainbow.close()
h5f_cc.close()
print output_SA
