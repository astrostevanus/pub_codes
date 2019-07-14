import numpy as np
import h5py
from  tqdm import tqdm
len_sv=20  #How many first singular values to be removed?
whichdata="kawa"
specdir="./"
ccd_list_2breduced=["blue","red"]

for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    h5f_eqwv = h5py.File(specdir+"rainbow/"+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
    h5f_svd = h5py.File(specdir+"rainbow/post_svd/"+str(ccd)+str(whichdata)+"_postsvd_fullmatrix.h5", 'w')

    for order in tqdm (range(1,len_order),ascii="True",desc="SVD for Order"):
        order_spec = h5f_eqwv[str(ccd)+"-flux-order-"+str(order)][:]
        wv_order= h5f_eqwv[str(ccd)+"-wv-order-"+str(order)][:]
        h5f_svd.create_dataset(str(ccd)+"-flux-order-"+str(order)+"sv-0", data= order_spec)
        u,w,v=np.linalg.svd(order_spec,full_matrices=False)
        for sv in range (len_sv):
            w[sv]=0
            rm=np.dot(u, np.dot(np.diag(w),v))
            h5f_svd.create_dataset(str(ccd)+"-flux-order-"+str(order)+"sv-"+str(sv+1), data=rm)
        h5f_svd.create_dataset(str(ccd)+"-wv-order-"+str(order), data=wv_order)
    h5f_svd.close()
    h5f_eqwv.close()
