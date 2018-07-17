# INPUT DATA INTERPOLATION FOR ESIM3
# RIESNA R AUDH
##################################################################################
# RUN FIRST 

# PATH TO BE SET BY USER
import numpy as np
import pandas as pd
import scipy
from scipy import interpolate

datadir='/home/riesna/MSc2018/ESIM_AntTest/NCEP/2/'

fname=[(datadir + 'Cl.txt'),(datadir + 'Fsd_cloud.txt'),(datadir + 'P_rate.txt'),
	(datadir + 'qa.txt'),(datadir + 'qs.txt'),(datadir + 'Ta.txt'),
	(datadir + 'U.txt'),(datadir + 'V.txt')]

x=np.arange(1,5840,4)
xi=np.arange(1,5840,1)

a=fname[0]
Cl=np.loadtxt(a)
Cl=np.reshape(Cl,(1,1460))
Cl1=scipy.interpolate.interp1d(Cl,Cl,kind='linear')

   




