import numpy as np
import h5py
from  tqdm import tqdm
from tqdm import tnrange, tqdm_notebook
specdir="./"

def sysrem_iter(order_spec,error,airmass):
    # SYSREM Iteration sequence
    a_fake=[] # 'airmass'
    c_fake=[] # 'color'
#    fitt=[]
    a_fake.append(airmass)
    for iter in range (0,40):
        c_i=[]
        a_j=[]
        for wvbin in range (len(std_pix)):
            c_wvbin=np.sum(order_spec[:,wvbin]*a_fake[iter][:]/error[:,wvbin]**2)/np.sum((a_fake[iter][:]/error[:,wvbin])**2)
            c_i.append(c_wvbin)
        c_fake.append(c_i)
        for filenum in range (len(std_frames)):
            a_filenum=np.sum(order_spec[filenum,:]*c_fake[iter][:]/error[filenum,:]**2)/np.sum((c_fake[iter][:]/error[filenum,:])**2)
            a_j.append(a_filenum)
        a_fake.append(a_j)
    sysrem_iter=np.zeros((len(std_frames),len(std_pix)))
    for wvbin in range (len(std_pix)):
        for filenum in range (len(std_frames)):
            sysrem_iter[filenum,wvbin]=c_fake[-1][wvbin]*a_fake[-1][filenum]
    return sysrem_iter

def sysrem_full(spec,error,n_sys,airmass):
    reduced_collect=[]
    sysrem_collect=[]
    data=spec
    for systematics in tqdm (range (n_sys),ascii="True",desc="Systematics:"):
        sysrem_systematics=sysrem_iter(data,error,airmass)
        data=data-sysrem_systematics
        reduced_collect.append(data)
        sysrem_collect.append(sysrem_systematics)
    reduced_collect=np.array(reduced_collect,dtype="float")
    sysrem_collect=np.array(sysrem_collect,dtype="float")  
    return reduced_collect,sysrem_collect

n_sys=10
whichdata="kawa"
opt="tio-9-farre-rotfit-fin"
vari="rot"
varmin=1
varmax=16

mjd,bjd,rvcor,hjd,airmass=np.loadtxt("wasp33-15.info",unpack=True,dtype="float")
h5f_eqwv = h5py.File(specdir+"rainbow-wasp33"+str(whichdata)+"_"+str(opt)+".h5", 'r')
h5f_sysrem = h5py.File(specdir+"post_sysrem/rainbow-wasp33"+str(whichdata)+"_postsysrem_"+str(opt)+".h5", 'w')

ccd_list_2breduced=["blue","red"]
for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    for va in tqdm(range(varmin,varmax),ascii="True",desc=str(vari)):
        for order in tqdm (range(1,len_order),ascii="True",desc=str(ccd)):
            order_spec = h5f_eqwv[str(ccd)+"-flux-"+str(vari)+"-"+str(varmin)+"-order-"+str(order)][:]
            for wvbin in range (len(order_spec[1])):
                order_spec[:,wvbin]=order_spec[:,wvbin]-np.mean(order_spec[:,wvbin])
            std_frames=[]
            std_pix=[]
            for wvbin in range (len(order_spec[1])):
                std_pix.append(np.std(order_spec[:,wvbin]))
            for filenum in range (len(order_spec)):
                std_frames.append(np.std(order_spec[filenum,:]))
            std_frames=np.array(std_frames,dtype="float")
            std_pix=np.array(std_pix,dtype="float") 
            error=np.zeros((len(std_frames),len(std_pix)))
            for wvbin in range (len(std_pix)):
                for filenum in range (len(std_frames)):
                    error[filenum,wvbin]=np.sqrt(std_frames[filenum]**2+std_pix[wvbin]**2)
#            for wvbin in range (len(order_spec[1])):
#                order_spec[:,wvbin]=order_spec[:,wvbin]-np.mean(order_spec[:,wvbin])
            reduced_collect,sysrem_collect=sysrem_full(order_spec,error,n_sys,airmass)
            h5f_sysrem.create_dataset(str(ccd)+"-flux-"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sysrem", data=reduced_collect)
            h5f_sysrem.create_dataset(str(ccd)+"-sys_fr-"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sysrem", data=sysrem_collect)
h5f_sysrem.close()
