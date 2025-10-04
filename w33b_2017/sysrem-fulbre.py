import numpy as np
import h5py
from  tqdm import tqdm
specdir="./rainbow/"

def sysrem_iter(order_spec,error,airmass):
    # SYSREM Iteration sequence
    a_fake=[] # 'airmass'
    c_fake=[] # 'color'
    a_fake.append(airmass)
    a_fake_std=[]
    c_fake_std=[]
    for iter in tqdm(range (0,400),ascii="True",desc="Iteran"):
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
        a_fake_std.append(np.std(a_j))
        c_fake_std.append(np.std(c_i))
    sysrem_iter=np.zeros((len(std_frames),len(std_pix)))
    for wvbin in range (len(std_pix)):
        for filenum in range (len(std_frames)):
            sysrem_iter[filenum,wvbin]=c_fake[-1][wvbin]*a_fake[-1][filenum]

    return sysrem_iter, a_fake[-1], c_fake[-1],a_fake_std,c_fake_std

def sysrem_full(spec,error,n_sys,airmass):
    reduced_collect=[]
    sysrem_collect=[]
    data=spec
    a_sys=[]
    c_sys=[]
    a_std_sys=[]
    c_std_sys=[]
    for systematics in tqdm (range (n_sys),ascii="True",desc="Systematics:"):
        sysrem_systematics,a_fake,c_fake,a_fake_std,c_fake_std=sysrem_iter(data,error,airmass)
        data=data-sysrem_systematics
        reduced_collect.append(data)
        sysrem_collect.append(sysrem_systematics)
        a_sys.append(a_fake)
        c_sys.append(c_fake)
        a_std_sys.append(a_fake_std)
        c_std_sys.append(c_fake_std)
    reduced_collect=np.array(reduced_collect,dtype="float")
    sysrem_collect=np.array(sysrem_collect,dtype="float")  
    return reduced_collect,sysrem_collect, a_sys, c_sys,a_std_sys,c_std_sys


n_sys=20
whichdata="kawa"
mode="order"
mjd,bjd,rvcor,hjd,airmass=np.loadtxt("wasp33-15.info",unpack=True,dtype="float")
h5f_sysrem = h5py.File(specdir+"/post_sysrem/rainbow"+str(whichdata)+"_postsysrem-usual"+str(mode)+".h5", 'w')

ccd_list_2breduced=["blue","red"]



print "##########################################################################"
print "SYSREM-Order/Fulbre v1.0 (2017 August 3)"
print "\n"
print "SYSREM Algorithm-> Tamuz et al. 2014"
print "Order           -> Doing SYSREM order by order "
print "Fulbre          -> SYSREM Modification to use full cross dispersed echelle"
print "                   spectroscopy order then break it again to each order"
print "Provided by Stevanus K. Nugroho"
print "Tohoku University"
print "sknugroho@astr,tohoku,ac.jp"
print "##########################################################################"
print "\n"
print "Input file: blue and red "+str(whichdata)+"_eqwv_fullmatrix.h5"
print "# Systematics to be removed: "+str(n_sys)
print "Mode: "+str(mode)
print "Iteran: 250x"
print "\n"
for ccd in ccd_list_2breduced:
    if ccd=="blue":
        len_order=19
    else:
        if whichdata=="este":
            len_order=16
        else:
            len_order=13
    print "\n"
    print str(ccd)+" CCD Reduction"
    print "Making Combined Order Matrix"
    h5f_eqwv = h5py.File(specdir+str(ccd)+str(whichdata)+"_eqwv_fullmatrix.h5", 'r')
    if mode=="fulbre":
        len_spec=[]
        len_spec.append(0)
        for order in range(1,len_order):
            order_spec = h5f_eqwv[str(ccd)+"-flux-order-"+str(order)][:]
            len_spec.append(len(order_spec[0]))
            lenspec=np.cumsum(np.array(len_spec))
        spec=np.zeros((len(order_spec),np.sum(len_spec)))
        print "Checking Full Array Size:"
        print "Length per order: "+str(np.array(len_spec))
        print "Cummulative length:"+str(lenspec)
        print "Shape of Full Array:"+str(np.shape(spec))
        print "\n"
        for order in range(1,len_order):
            order_spec = h5f_eqwv[str(ccd)+"-flux-order-"+str(order)][:]
            for filenum in range (len(order_spec)):
                spec[filenum][lenspec[order-1]:lenspec[order]]=order_spec[filenum]
#        for wvbin in range (len(spec[1])):
#            spec[:,wvbin]=spec[:,wvbin]-np.mean(spec[:,wvbin])
        print "Calculating error for each pixel"
        std_frames=[]
        std_pix=[]
        for wvbin in range (len(spec[1])):
            std_pix.append(np.std(spec[:,wvbin]))
        for filenum in range (len(spec)):
            std_frames.append(np.std(spec[filenum,:]))
        std_frames=np.array(std_frames,dtype="float")
        std_pix=np.array(std_pix,dtype="float") 
        error=np.zeros((len(std_frames),len(std_pix)))
        for wvbin in range (len(std_pix)):
            for filenum in range (len(std_frames)):
                error[filenum,wvbin]=np.sqrt(std_frames[filenum]**2+std_pix[wvbin]**2)
        for wvbin in range (len(spec[1])):
            spec[:,wvbin]=spec[:,wvbin]-np.mean(spec[:,wvbin])
        print "Finished calculating error"
        print "\n"
        print "Doing SYSREM-Fulbre"
        reduced_collect,sysrem_collect,a_sys,c_sys,a_std_sys,c_std_sys=sysrem_full(spec,error,n_sys,airmass)
        print "SYSREM-Fulbre Finished"

    for order in tqdm (range(1,len_order),ascii="True",desc=str(ccd)):
        if mode=="fulbre":
            sys_flux_order=[]
            sys_sys_order=[]
            for systematics in tqdm (range (n_sys),ascii="True",desc="Systematics:"):
                sysrem_flux_order=np.zeros((len(order_spec),len_spec[order]))
                sysrem_sys_order=np.zeros((len(order_spec),len_spec[order]))
                for filenum in range (len(order_spec)):
                    sysrem_flux_order[filenum]=reduced_collect[systematics][filenum][lenspec[order-1]:lenspec[order]]
                    sysrem_sys_order[filenum]=sysrem_collect[systematics][filenum][lenspec[order-1]:lenspec[order]]
                sys_flux_order.append(sysrem_flux_order)
                sys_sys_order.append(sysrem_sys_order)
        elif mode=="order":
            order_spec = h5f_eqwv[str(ccd)+"-flux-order-"+str(order)][:]
#            for wvbin in range (len(order_spec[1])):
#                order_spec[:,wvbin]=order_spec[:,wvbin]-np.mean(order_spec[:,wvbin])
            print "Calculating error for each pixel"
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
            for wvbin in range (len(order_spec[1])):
                order_spec[:,wvbin]=order_spec[:,wvbin]-np.mean(order_spec[:,wvbin])
            print "Finished calculating error"
            print "Doing SYSREM-Order"
            sys_flux_order,sys_sys_order,a_sys,c_sys,a_std_sys,c_std_sys=sysrem_full(order_spec,error,n_sys,airmass)
            h5f_sysrem.create_dataset(str(ccd)+"-a_std_sys-order-"+str(order),data=a_std_sys)
            h5f_sysrem.create_dataset(str(ccd)+"-c_std_sys-order-"+str(order),data=c_std_sys)
        h5f_sysrem.create_dataset(str(ccd)+"-flux-order-"+str(order)+"sysrem", data=sys_flux_order)
        h5f_sysrem.create_dataset(str(ccd)+"-sys_fr-order-"+str(order)+"sysrem", data= sys_sys_order)
        h5f_sysrem.create_dataset(str(ccd)+"-a_sys-order-"+str(order)+"sysrem", data=a_sys)   
        h5f_sysrem.create_dataset(str(ccd)+"-c_sys-order-"+str(order)+"sysrem", data=c_sys)
    if mode=="fulbre":
        h5f_sysrem.create_dataset(str(ccd)+"-a_std_sys",data=a_std_sys)
        h5f_sysrem.create_dataset(str(ccd)+"-c_std_sys",data=c_std_sys)
    elif mode=="order":
        print "SYSREM-Order Finished"
h5f_sysrem.close()
print "\n"
print "Output file: rainbow"+str(whichdata)+"_postsysrem-"+str(mode)+".h5"
print "Thank you. Now you can proceed to the cross correlation step or inspecting the SYSREM result"
