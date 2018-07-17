#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 01:08:04 2017

@author: riesna
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import xarray as xr
import matplotlib
import matplotlib.patches as mpatches
import pandas as pd

x = (1,2,3,4,5,6)
x1= (7,8,9,10)
y = (2,90,15,10,80,50)
y1= (20,75,5,95)
y2= ()
y3= ()

raw_data = {'CTD #': ['1', '2', '3', '4', '5','6','7','8','9','10'],
        'Observed SIC': [0,80,10,0,70,20,10,50,10,90],
        'Satellite SIC': [67.099,56.7099,73.2299,61.9099,76.6599,76.9899,
                      9.5899,75.9499,93.4599,2.2899]}
df = pd.DataFrame(raw_data, columns = ['CTD #', 'Observed SIC', 'Satellite SIC'])
df

pos = list(range(len(df['Observed SIC']))) 
width = 0.25 

# Plotting the bars
fig, ax = plt.subplots(figsize=(10,5))
plt.fill([-0.5,6,6,-0.5], [0,0,100,100], 'b', alpha=0.2)
plt.fill([10.2,6,6,10.2],[0,0,100,100],'g',alpha=0.2)

# Create a bar with pre_score data,
# in position pos,
plt.bar([p + width for p in pos],df['Observed SIC'],width,alpha=0.8,color='blue',
        label=df['CTD #'][0]) 

# Create a bar with mid_score data,
# in position pos + some width buffer,
plt.bar([p + 2*width for p in pos],df['Satellite SIC'],width,alpha=0.8,
        color='#F78F1E') 

# Set the y axis label
ax.set_ylabel('Sea Ice Concentration [%]',fontsize=15)
ax.set_xlabel('CTD number',fontsize=15)

# Set the chart's title
ax.set_title('Observed and satellite derived sea ice concentration',fontsize=20)

# Set the position of the x ticks
ax.set_xticks([p + 1.5 * width for p in pos])

# Set the labels for the x ticks
ax.set_xticklabels(df['CTD #'],fontsize=15)
ax.set_yticklabels([0,20,40,60,80,100],fontsize=15)
# Setting the x-axis and y-axis limits
plt.xlim(min(pos)-width, max(pos)+width*4)
plt.ylim([0,100])

# Adding the legend and showing the plot
red_patch = mpatches.Patch(color='blue',alpha=0.8, label='Observed SIC')
blue_patch = mpatches.Patch(color='#F78F1E',alpha=0.8, label='Satellite SIC')
box1 = mpatches.Patch(color = 'blue', alpha=0.2,label='Transect B')
box2 = mpatches.Patch(color = 'green', alpha=0.2,label='Transect A')
plt.legend(handles=[red_patch,blue_patch,box1,box2],loc='upper left',fontsize=15)

plt.grid()
plt.show()




## add some text for labels, title and axes ticks
#ax.set_ylabel('Sea Ice Concentration [%]')
#ax.set_xlabel('CTD Number')
#ax.set_title('Sea ice concentration during CTD cast')
#ax.set_yticks([0,10,20,30,40,50,60,70,80,90,100])
#ax.set_xticks([1,2,3,4,5,6,7,8,9,10])
#ax.set_xticklabels(['1','2','3','4','5','6','7','8','9','10'])

#red_patch = mpatches.Patch(color='red', label='Transect A')
#blue_patch = mpatches.Patch(color='blue', label='Transect B')
#ax.legend(handles=[red_patch,blue_patch])
#ax.legend(handles=[blue_patch])
