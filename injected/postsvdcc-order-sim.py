import numpy as np
from ultramodule import divstd_prime, crossco, cc_multiorde, print_item, print_structure
import h5py
from  tqdm import tqdm
from tqdm import tnrange
from PyAstronomy import pyasl
###################################################################################################                                        #rainbow-wasp33"+str(whichdata)+"_postsvd_"+str(opt)+".h5"                                    
len_sv=10
#vari="mul"
#varimin=1
#varimax=6
vari="addmul"
varimin=1
varimax=6
whichdata="kawa"
mol="tio"
#opt="vo_8x10-farre"
#opt="tio-7x7-farre-corr2-fin"
#opt="tio_finalmodel-7addlow5far2ex"
#opt="tio-9-atre" # rainbowkawa_postsysrem-superorder_tio-9-atre.h5
#rainbowkawa_postsysrem-usualorder_tioall-rotfit-atre.h5 
opt="tioall-addmul-fineF"
#opt="tioall-rotfit-atre"
#opt="tioall-rotfit-boostkp4"
#opt="tioall-rotfit-boostatre"
redmode="sysrem"
opt_tech="-fin"
specdir="./"
print "vari= "+str(vari)
print "varimin= "+str(varimin)
print "varimax= "+str(varimax)
print "whichdata= "+str(whichdata)
print "opt= "+str(opt)
print "redmode= "+str(redmode)
print "opt_tech="+str(opt_tech)

def rankmat(mat):
    temp = mat.argsort()
    rank = temp.argsort()
    return rank
def cc_sp(mat1,mat2):
    r1=rankmat(mat1)
    r2=rankmat(mat2)
    sp1=r1-np.mean(r1)
    sp2=r2-np.mean(r2)
    up=np.sum(sp1*sp2)
    below=np.sqrt(np.sum(sp1*sp1)*np.sum(sp2*sp2))
    cc=up/below
    return cc

###################################################################################################                                                                            
ccd_list_2breduced=["blue","red"]
#ccd_list_2breduced=["red"]
h5f_techni = h5py.File("techni-"+str(whichdata)+str(opt_tech)+".info", 'r')
drvs=h5f_techni["Doppler shifts variation"][:]
output_tmp="rainbow-wasp33"+str(whichdata)+"all_"+str(mol)+"_template"+str(opt_tech)+".h5"
#rainbow-wasp33kawaall_tio_template-corr2-finmod.h5
h5f_temp = h5py.File("/home/stev/fin-pycovfefe/rainbow/template/"+output_tmp, 'r')

output_cc="rainbow-wasp33"+str(whichdata)+"_cc_pearson"+str(opt)+"cor2tempuncor"+str(opt_tech)+".h5"#rainbowkawa_postsysrem-superorder_tio-9-atre-finultihalf.h5
h5f_cc = h5py.File(specdir+"post_"+str(redmode)+"/cc_result/"+str(output_cc), 'w') 
#print_structure(specdir+"post_"+str(redmode)+"/rainbow"+str(whichdata)+"_post"+str(redmode)+"-superorder_"+str(opt)+".h5")
h5f_reduced = h5py.File(specdir+"post_"+str(redmode)+"/rainbow"+str(whichdata)+"_post"+str(redmode)+"-usualorder_"+str(opt)+".h5", 'r')

print "Template= "+str(output_tmp)
print "Output=" +str(output_cc)
print "#############################################"
print ""
if redmode=="sysrem":
    varname="systematics"
elif redmode=="svd":
    varname="sv"

for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
#    output_cc=str(ccd)+"-wasp33"+str(whichdata)+"_cc_"+str(opt)+".h5"
#    h5f_cc = h5py.File(specdir+"cc_result/"+str(output_cc), 'w')
    for va in tqdm(range(varimin,varimax),ascii="True",desc=str(ccd)+" CCD"):
        cc_map_full=[]
        for var in tqdm(range(len_sv),ascii="True",desc=str(vari)+"="+str(va)+" R"+str(varname)):
            cc_svd_collect=[] #Different SV are saved in this matrix                       
            for order in tqdm (range(1,len_order),ascii="True",desc="Order"):
#                y_temp_order=h5f_temp[str(ccd)+"-vmr-"+str(va)+"-order-"+str(order)][:]
                yobs=h5f_reduced[str(ccd)+"-flux"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sysrem"][:][var]
#                yobs=h5f_reduced[str(ccd)+"-flux-"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sysrem"][:][var]
                y_temp_order=h5f_temp[str(ccd)+"-vmr-9-order-"+str(order)][:]
#                yobs= h5f_reduced[str(ccd)+"-flux-"+str(vari)+"-7-order-"+str(order)+str(varname)+"-"+str(var+1)][:]
                for filenum in range (len(yobs)):
                    sm = pyasl.smooth(yobs[filenum], 25, 'flat')
                    sm1 = pyasl.smooth(sm, 51, 'flat')
                    yobs[filenum]=(yobs[filenum]+1.)/(sm1+1.)-1.
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
#                        cc_order[numspec][rv]=cc_sp(yobs[numspec],y_temp_order[rv])
                cc_svd_collect.append(cc_order)
            cc_map_full.append(np.array(cc_svd_collect))
        h5f_cc.create_dataset(str(ccd)+"_cc_"+str(vari)+"_"+str(va), data= np.array(cc_map_full)) #Writing the template of the specific vmr
h5f_cc.close()
h5f_reduced.close()
