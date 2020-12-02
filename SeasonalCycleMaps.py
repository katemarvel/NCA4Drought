#!/usr/bin/env python
# coding: utf-8

# In[39]:


import cdms2 as cdms
import MV2 as MV
import cdtime,cdutil,genutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import string
import glob
import scipy.stats as stats
# Local solution
# If running remotely, uncomment the following code:
# %%bash
# git clone https://github.com/katemarvel/CMIP5_tools
# import CMIP5_tools as cmip5
import sys,os
sys.path.append("../python-utils")

import CMIP5_tools as cmip5
import DA_tools
import Plotting


# cartopy stuff
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry.polygon import LinearRing
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from eofs.cdms import Eof
from eofs.multivariate.cdms import MultivariateEof


rootdirec="/home/kdm2144/"
vcert=stats.norm.interval(.99)[1]

#Region locations
NCA4regions={}
#Northwest (NW): (125°W–111°W, 42°N–49°N)
NCA4regions["NW"]=cdutil.region.domain(longitude=(-125,-111),latitude=(42,49))
#Southwest (SW): (124°W–102°W, 31°N–42°N)
NCA4regions["SW"]=cdutil.region.domain(longitude=(-124,-102),latitude=(31,42))
#Upper Great Plains (GPu): (116°W–95°W, 40°N–49°N)
NCA4regions["GPu"]=cdutil.region.domain(longitude=(-116,-95),latitude=(40,49))
#Lower Great Plains (GPl): (107°W–93°W, 26°N–40°N)
NCA4regions["GPl"]=cdutil.region.domain(longitude=(-107,-93),latitude=(26,40))
#Midwest (MW): (97°W–80°W, 36°N–50°N)
NCA4regions["MW"]=cdutil.region.domain(longitude=(-97,-80),latitude=(36,50))
#Northeast (NE): (82°W–67°W, 37°N–48°N)
NCA4regions["NE"]=cdutil.region.domain(longitude=(-82,-67),latitude=(37,48))
#Southeast (SE): (95°W–76°W, 25°N–39°N)
NCA4regions["SE"]=cdutil.region.domain(longitude=(-95,-76),latitude=(25,39))


external_drive="/home/kdm2144/"

import DroughtHelper2 as dh


# In[40]:


from scipy import signal

def FourierPlot(tas):
    detrend = signal.detrend(tas)
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    tas_fft = np.fft.fft(detrend)
    R = tas_fft.real
    Im = tas_fft.imag
    mag = np.sqrt(R**2+Im**2)
    plt.plot(1/freqs,mag)
def annual_cycle_dominant(tas):
    """Check to see whether the annual cycle is dominant"""
    detrend = signal.detrend(tas)
    L = len(tas)
    freqs = np.fft.fftfreq(L)
    tas_fft = np.fft.fft(detrend)
    R = tas_fft.real
    Im = tas_fft.imag
    mag = np.sqrt(R**2+Im**2)
    the_period = 1./np.abs(freqs[np.argmax(mag)])
    return the_period

def get_dominant_cycle(tas):
    nt,nlat,nlon = tas.shape
    to_mask = MV.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            to_mask[i,j]=annual_cycle_dominant(tas[:,i,j])
    to_mask.setAxisList(tas.getAxisList()[1:])
    return to_mask
def get_cycle(tas,period=12,return_complex=False):

    L = len(tas)
    freqs = np.fft.fftfreq(L)
    closest = np.abs(freqs-1./period)
    
#    i = np.where(freqs == 1./period)[0]
    i = np.argmin(closest)
   # print(1/freqs[i])
    tas_fft = np.fft.fft(tas)/L
    R = tas_fft.real
    Im = tas_fft.imag
    if return_complex:
        return R[i],Im[i]
    else:
        mag = np.sqrt(R**2+Im**2)
        phase = np.arctan2(Im,R)
        return mag[i],phase[i]
def day_of_peak(phase):
    #make sure phase lies between [0,2pi]
    phase=stats.circmean(phase)
    days=np.arange(365)
    func=np.cos(2*np.pi/12*np.linspace(0,12,365)+phase)
    #return the day of max variable if it's between spring and fall
    #if (phase>np.pi/4) and (phase<7*np.pi/4):
    return float(days[np.argmax(func)])
   # #otherwise winter
   # else:
       # return float(days[np.argmin(func)])
vec_day_of_peak=np.vectorize(day_of_peak)

def day_of_trough(phase):
    #make sure phase lies between [0,2pi]
    phase=stats.circmean(phase)
    days=np.arange(365)
    func=np.cos(2*np.pi/12*np.linspace(0,12,365)+phase)
    #return the day of max variable if it's between spring and fall
    #if (phase>np.pi/4) and (phase<7*np.pi/4):
    return float(days[np.argmin(func)])
   # #otherwise winter
   # else:
       # return float(days[np.argmin(func)])
vec_day_of_trough=np.vectorize(day_of_trough)


# In[38]:


def regrid_models(variable,experiment):

    #get the shape from CESM2


    model = "CESM2"

    rip="r1i1p1f1"
    allfiles=sorted(glob.glob("/home/kdm2144/DROUGHT/DOWNLOADED_RAW/"+"/"+variable+"/"+model+"/*."+experiment+".*."+rip+".*"))

    f=cdms.open(allfiles[0])
    data=f(variable)
    #for model in get_ok_models("SW"):

    grid=data.getGrid()
    nyears=86 #historical runs begin in 2015-2100
    model_shape=(12*nyears,)+grid.shape

    f.close()

    models=dh.get_ok_models("SW")
    nmodels=len(models)
    bigshape=(nmodels,)+model_shape


    allmodels=MV.zeros(bigshape)
    for modeli in range(nmodels):
        model=models[modeli]
        #print(model)
        allfiles_rips=sorted(glob.glob("/home/kdm2144/DROUGHT/DOWNLOADED_RAW/"+"/"+variable+"/"+model+"/*."+experiment+".*"))
        rips=sorted(np.unique([x.split(".")[3] for x in allfiles_rips]))
        if len(rips)>0:
            rip=rips[0]
        else:
            continue
        allfiles=sorted(glob.glob("/home/kdm2144/DROUGHT/DOWNLOADED_RAW/"+"/"+variable+"/"+model+"/*."+experiment+".*."+rip+".*"))
        if len(allfiles ) >= nyears:
            for i in range(nyears):
                f=cdms.open(allfiles[i])
                data=f(variable)
                data_regrid=data.regrid(grid,regridTool='regrid2')
                if i==0:
                    bigdata=data_regrid
                else:
                    bigdata=MV.concatenate((bigdata,data_regrid))
                f.close()
            allmodels[modeli]=bigdata
        else:
            allmodels[modeli]=1.e20
    allmodels=MV.masked_where(np.abs(allmodels)>1.e10,allmodels)
    modax=cmip5.make_model_axis(models)
    tax=cdms.createAxis(np.arange(12*86))
    tax.units="months since 2015-1-1"
    tax.designateTime()
    tax.id='time'
    allmodels.setAxis(0,modax)
    allmodels.setAxis(1,tax)
    cdutil.setTimeBoundsMonthly(allmodels)
    allmodels.setAxis(2,data_regrid.getLatitude())
    allmodels.setAxis(3,data_regrid.getLongitude())    
    allmodels.name=variable
    allmodels.id=variable
    return allmodels


# In[27]:


def mma_amp_phase(mma):
    nmonths,nlat,nlon=mma.shape
    nyears=int(nmonths/12)
    A=MV.zeros((nyears,nlat,nlon))
    P=MV.zeros((nyears,nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            if not mma.mask[0,i,j]:
                data=mma[:,i,j]
                if annual_cycle_dominant(data)==12:
                    for year in range(nyears):
                        amp,phase=get_cycle(data[year*12:(year+1)*12])
                        A[year,i,j]=amp
                        P[year,i,j]=phase
    P=MV.masked_where(P==0,P)
    yav=cdutil.YEAR(mma)
    P.setAxisList(yav.getAxisList())
    P.id="phase"

    A=MV.masked_where(A==0,A)
    P=MV.masked_where(A==0,P)
    A.setAxisList(yav.getAxisList())
    A.id="amp"
    return A,P


#Write the data


variables=["prsn","tas","evspsbl","mrro","mrros","mrso","mrsos","pr"]
experiments=["ssp585","historical"]
for experiment in experiments:
    for variable in variables:

        alldata=regrid_models(variable,experiment)
        fw=cdms.open("/home/kdm2144/DROUGHT/Maps/"+variable+"."+experiment+".MMA.2015_2100.nc","w")
        fw.write(alldata)
        fw.close()

        mma=MV.average(alldata,axis=0)

        A,P=mma_amp_phase(mma)

        fw=cdms.open("/home/kdm2144/DROUGHT/Maps/"+variable+"."+experiment+".AmpPhase.2015_2100.nc","w")
        fw.write(A)
        fw.write(P)
        fw.close()







