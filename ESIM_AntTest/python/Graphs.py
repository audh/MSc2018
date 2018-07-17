#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:43:15 2018

@author: riesna
loop to plot CTD and nutrient values on one figure 

"""

from seabird.cnv import fCNV
from matplotlib import pyplot as plt
import ctd 
from ctd import DataFrame
from seawater import satO2

cast = [DataFrame.from_cnv('')] #import .cnv file (create list for multiple CTD files)

a = [[var-n-],[var-n2-]]#read in variables (nutrients)
b = [[var-n-],[var-n2-]]
c = [[var-n-],[var-n2-]]
d = [[var-n-],[var-n2-]]

for j in range (0,6): #ref to the cast range
    a1 = cast[j]
    downcast, upcast = a1.split()    #splits up and down cast of CTD file
    plt.figure()
   
    # Temperature
    ax1 = plt.subplot(231)
    plt.plot(downcast['t090C'],downcast['depSM'],'k') #CTD file variables
    plt.setp(ax1.get_yticklabels(), fontsize=6)
    plt.xlabel('Temperature [degC]')
    plt.ylabel('Depth')
    plt.gca().invert_yaxis()
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    # Salinity
    ax2 = plt.subplot(232,sharey=ax1)
    plt.plot(downcast['sal00'],downcast['depSM'],'r')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.xlabel('Salinity')
    plt.gca().invert_yaxis()
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    # oxygen
    ax3 = plt.subplot(233, sharey=ax1)
    plt.plot(downcast['sbeox0ML/L'],downcast['depSM'],'g')
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.xlabel('Oxygen')
    plt.gca().invert_yaxis()
    ax3.xaxis.tick_top()
    ax3.xaxis.set_label_position('top')
    
    #Nutrients
    #Nitrate
    ax4 = plt.subplot(234)
    plt.plot(nitrate[j],depth[j],'.k-')
    plt.setp(ax4.get_yticklabels(), fontsize=6) #set y axis
    plt.xlabel('Nitrate Concentration')
    plt.ylabel('Depth')
    plt.gca().invert_yaxis()
    ax4.xaxis.tick_top()
    ax4.xaxis.set_label_position('top')
    # Silicate
    ax5 = plt.subplot(235, sharey=ax4)
    plt.plot(silicate[j],depth[j],'.r-')
    plt.setp(ax5.get_yticklabels(), visible=False) #hide y axis since shared
    plt.xlabel('Silicate Concentration')
    plt.gca().invert_yaxis()
    ax5.xaxis.tick_top()
    ax5.xaxis.set_label_position('top')
    # Phosphate
    ax6 = plt.subplot(236,sharey=ax4)
    plt.plot(phosphate[j],depth[j],'.g-')
    plt.setp(ax6.get_yticklabels(), visible=False)
    plt.xlabel('Phosphate Concentration')
    plt.gca().invert_yaxis()
    ax6.xaxis.tick_top()
    ax6.xaxis.set_label_position('top')
    if j==0 or j==1 or j==2 or j==3:
        plt.suptitle('Transect A CTD '+str(j+7), y=1.09,x=-0.8,fontsize=15)
    else:
        plt.suptitle('Transect B CTD '+str(j-2), y=1.09,x=-0.8,fontsize=15)

