import numpy as np
import h5py
from  tqdm import tqdm
len_sv=20  #How many first singular values to be removed?
whichdata="kawa"
specdir="./"
opt="tio-4x10-atre"
#opt="tio_finalmodel-farre"
#opt="tio_finalmodel-7addlow5far2ex"

#vari="mul" #mul or vmr
#varmin=1  #1-6 or 4-11
#varmax=6
vari="vmr"
varmin=10
varmax=11
ccd_list_2breduced=["blue","red"]


h5f_eqwv = h5py.File(specdir+"rainbow-wasp33"+str(whichdata)+"_"+str(opt)+".h5", 'r')
h5f_svd = h5py.File(specdir+"post_svd/rainbow-wasp33"+str(whichdata)+"_postsvd_"+str(opt)+".h5", 'w')

for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    for order in tqdm (range(1,len_order),ascii="True",desc="SVD for Order"):
        wv_order= h5f_eqwv[str(ccd)+"-flux-"+str(vari)+"-"+str(varmin)+"-order-"+str(order)][:]
        h5f_svd.create_dataset(str(ccd)+"-wv-order-"+str(order), data=wv_order)
        for va in tqdm(range(varmin,varmax),ascii="True",desc=str(vari)):
            order_spec = h5f_eqwv[str(ccd)+"-flux-"+str(vari)+"-"+str(va)+"-order-"+str(order)][:]
#            order_spec=np.array(order_spec,dtype=float)[~np.isnan(np.mean(order_spec,axis=1))]
            h5f_svd.create_dataset(str(ccd)+"-flux-"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sv-0", data= order_spec)
            u,w,v=np.linalg.svd(order_spec,full_matrices=False)
            for sv in range (len_sv):
                w[sv]=0
                rm=np.dot(u, np.dot(np.diag(w),v))
                h5f_svd.create_dataset(str(ccd)+"-flux-"+str(vari)+"-"+str(va)+"-order-"+str(order)+"sv-"+str(sv+1), data=rm)
h5f_svd.close()
h5f_eqwv.close()
