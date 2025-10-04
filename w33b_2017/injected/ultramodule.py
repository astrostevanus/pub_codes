import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
from astropy.table import Table
import astropy.time as time
import astropy.coordinates as coordsimport
import astropy.units as u
import astropy.constants as const
import h5py
import scipy.interpolate as sci
import pyfits
import scipy
from matplotlib.ticker import NullFormatter
import jplephem
import de423

nullfmt = NullFormatter()
plt.rcParams['font.family'] = 'Heltevica'
fontsize = 18
plt.rcParams['font.size'] = fontsize
c = 299792.458
pi=3.141592653589793

hplanck= 6.62607004e-34 # m2 kg / s
c_light = 299792458 #m/s
kb= 1.38064852e-23 #m2 kg s-2 K-1

specdir="/Volumes/SaySomething/Observation-Data/PrimeWASP33/rainbow/"
exotempdir="/Volumes/SaySomething/py4cats/data/spec/wasp33b/"

longitude = -155.47333
latitude = 19.82444
altitude = 4205.

def crossco(f,fi):
    up=np.sum((f-np.mean(f))*(fi-np.mean(fi)))
    bottom=np.sum((f-np.mean(f))**2)*np.sum((fi-np.mean(fi))**2)
    cc = up/(bottom**0.5)
    return cc


def intensity(wavelength,temperature):
    #Wavelength input in Angstrom, Temperature in K
    upper=2*hplanck*c_light**2*wavelength**-5
    lower=wavelength*kb*temperature
    I_lamda= upper/(np.exp(hplanck*c_light/lower)-1)
    return I_lamda

def divstd_prime(rm,std_pix,std_frames):
#Extract standard deviation (std) per wavelength bin and frame, then weight wavelength bin and frame by their std
    for wvbin in range (len(rm[1])):
        rm[:,wvbin]=rm[:,wvbin]/std_pix[wvbin]
    for filenum in range (len(rm)):
        rm[filenum,:]=rm[filenum,:]/std_frames[filenum]
    return rm

def timestamp(objectname):
    mjd,airmass,dom_tmp,dom_hum,dom_prs,adc=np.loadtxt("/Volumes/SaySomething/Observation-Data/WASP-33/obs/red/sr-red/combined/header.info",unpack=True)
    hjd_list=[]
    rvcor=[]
    bjd=jd_corr(mjd,objectname).value #the BJD of the frames
    for i in range (len(mjd)):
        correct=auto_barrcor(objectname,mjd[i])
        hjd_list.append(correct[1])
        rvcor.append(correct[0])
    hjd_list=np.array(hjd_list,dtype="float")

    file=open(str(objectname)+".info","w")
    for i in range (len(bjd)):
        file.write(str.format("{0:.16f}",mjd[i])+" "+str.format("{0:.16f}",bjd[i])+" "+str.format("{0:.16f}",rvcor[i])+" "+str.format("{0:.16f}",hjd_list[i])+"\n")
    file.close()
    print str(objectname)+".info created"


def cc_multiorde(cc):
    #Combining cross correlation from different orders or frames
    #Input => cc : row= specnum; column= rv
    #(Zucker et al. 2003)
    cc_mo=[]
    for rv in range (len(cc[1])):
        cc_mo.append(np.sqrt(1.-(np.prod(1.-cc[:,rv]**2))**(1./len(cc))))
    return np.array(cc_mo,dtype="float")

def sysrem_iter(order_spec,error):
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

def sysrem_full(spec,error,n_sys):
    reduced_collect=[]
    sysrem_collect=[]
    data=spec
    for systematics in tnrange (n_sys,ascii="True",desc="Systematics:"):
        sysrem_systematics=sysrem_iter(data,error)
        data=data-sysrem_systematics
        reduced_collect.append(data)
        sysrem_collect.append(sysrem_systematics)
    reduced_collect=np.array(reduced_collect,dtype="float")
    sysrem_collect=np.array(sysrem_collect,dtype="float")
    return reduced_collect,sysrem_collect
def print_structure(file_name) :
    """Prints the HDF5 file structure"""
    file = h5py.File(file_name, 'r') # open read-only
    item = file #["/Configure:0000/Run:0000"]
    print_item(item)
    file.close()
 
def print_item(g, offset='    ') :
    """Prints the input file/group/dataset (g) name and begin iterations on its content"""
    if   isinstance(g,h5py.File) :
        print g.file, '(File)', g.name
 
    elif isinstance(g,h5py.Dataset) :
        print '(Dataset)', g.name, '    len =', g.shape #, g.dtype
 
    elif isinstance(g,h5py.Group) :
        print '(Group)', g.name
 
    else :
        print 'WORNING: UNKNOWN ITEM IN HDF5 FILE', g.name
        sys.exit ( "EXECUTION IS TERMINATED" )
 
    if isinstance(g, h5py.File) or isinstance(g, h5py.Group) :
        for key,val in dict(g).iteritems() :
            subg = val
            print offset, key, #,"   ", subg.name #, val, subg.len(), type(subg),
            print_item(subg, offset + '    ')


import jplephem
import de423
hplanck= 6.62607004e-34 # m2 kg / s
c_light = 299792458 #m/s
kb= 1.38064852e-23 #m2 kg s-2 K-1

def intensity(wavelength,temperature):
    #Wavelength input in Angstrom, Temperature in K
    upper=2*hplanck*c_light**2*wavelength**-5
    lower=wavelength*kb*temperature
    I_lamda= upper/(np.exp(hplanck*c_light/lower)-1)
    return I_lamda
