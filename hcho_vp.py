#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 15:16:25 2025

@author: cnaught
"""

import datetime
import sys
import numpy as np
import matplotlib.pyplot as plt
import pytz
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

pan_dir = r'/projectnb/atmchem/naughton'
sys.path.insert(0, pan_dir)
from pan_hcho_profiles import plot_hcho_profiles_proper_vv, convert_profile_df_hcho_standard, plot_interp_profiles, plot_raw_hcho_profiles, return_profile_quantiles

pan_dir = r'/projectnb/atmchem/shared/for_claire' #move file into /project/atmchem/shared/ 
sys.path.insert(0, pan_dir)
from load_l2h_1_8 import pan_header, header_dict, l2h_to_df

sys.path.append("/usr3/graduate/cnaught/Desktop/pandora files")
from pan_functions import filter_pan_hcho_skyscan, subset_col, mol_per_m3_to_ppb, flatten_extend, filter_pan_hcho_ds

def convert_tz(df):
    boston_tz = pytz.timezone('America/New_York')
    df.index = df.index.tz_convert(boston_tz)
    df.index = df.index.tz_localize(None)
    df = df[(df.index.month >= 5) & (df.index.month <=9)]
    return(df)
#%%
output_files = 'rfuh5p1-8.txt' #HCHO SS
Pan26_fp = '/projectnb/atmchem/shared/pandora/Pandora26s1_CambridgeMA_L2_' + output_files
header_content, col_n = pan_header(Pan26_fp)
Pan26_Data = l2h_to_df(Pan26_fp, header_content)
col_dict = header_dict(header_content, col_n) 
Pan26_Filtered_ss = filter_pan_hcho_skyscan(Pan26_Data)
cam_ss = convert_tz(Pan26_Filtered_ss)
pan_sum = cam_ss.loc[(cam_ss.index.month >=5) & (cam_ss.index.month <= 9)]
#%%
dict2021 = pd.read_pickle(r'/project/atmchem/naughton/met2021.pkl')
dict2022 = pd.read_pickle(r'/project/atmchem/naughton/met2022.pkl')
ds = [dict2021, dict2022]
date_dict = {}
for k in dict2021.keys():
    date_dict[k] = tuple(date_dict[k] for date_dict in ds)
westerly = flatten_extend(date_dict.get("Westerly"))
west = pd.DataFrame(pan_sum[pan_sum.index.floor("D").isin(westerly)])
#%%
alt_int = np.arange(0.10, 2.5, 0.1)

profiles_dict = {}
mean_dict = {}
for i in range(7,20):
    dat_subset = west.loc[west.index.hour == i]
    dat_subset_interp = convert_profile_df_hcho_standard(dat_subset, alt_int)
    dat_quantiles = return_profile_quantiles(dat_subset_interp)
    dat_median = dat_quantiles['q_50']
    profiles_dict[i] = dat_median
    mean_dict[i] = dat_subset_interp.T.mean()

profiles_df = pd.DataFrame(profiles_dict)
alt = profiles_df.index.astype(float)
hour = profiles_df.columns
time_grid, alt_grid = np.meshgrid(hour, alt)
hcho = profiles_df.values
#%%
colors = ['white', 'skyblue', 'mediumseagreen', 'yellow', 'darkorange', 'red', 'darkred']
white_rainbow = LinearSegmentedColormap.from_list("white_rainbow", colors)
    
vmin = 0
vmax = 5

plt.figure(figsize=(10, 3))
dat = plt.pcolormesh(time_grid, alt_grid, hcho,vmin = vmin, vmax = vmax, cmap = white_rainbow)
cb = plt.colorbar(dat, label = 'HCHO (ppb)', pad = 0.02)
plt.xlabel("Hour of Day (LT)")
plt.ylabel("Altitude (km)")

# Show plot
plt.xticks(np.arange(7,19,1))
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("westerly_vp_median", dpi = 600)
plt.show()
#%%
surf = profiles_df.loc[0.1]
surf2 = profiles_df.loc[0.2]
surf3 = profiles_df.loc[0.30000000000000004]
surf4 = profiles_df.loc[0.4]
surf5 = profiles_df.loc[0.5]

plt.rcParams['font.family'] = 'Nimbus Sans'
plt.figure(figsize = (10,6))
plt.plot(surf, color = 'red', marker = 'o', label = '0.1km')
plt.plot(surf2, color = 'blue', marker = 'o', label = '0.2km')
plt.plot(surf3, color = 'green', marker = 'o', label = '0.3km')
plt.plot(surf4, color = 'orange', marker = 'o', label = '0.4km')
plt.plot(surf5, color = 'purple', marker = 'o', label = '0.5km')
plt.xlabel("Hour of Day (LT)")
plt.ylabel("HCHO (ppb)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("vp_surface_layer", dpi = 600)
plt.show()

#%%
mean_prof_df = pd.DataFrame(mean_dict)
alt_mean = mean_prof_df.index.astype(float)
hour_mean = mean_prof_df.columns
time_grid_mean, alt_grid_mean = np.meshgrid(hour_mean, alt_mean)
hcho_mean = mean_prof_df.values

vmin = 0
vmax = 5

plt.figure(figsize=(12, 3))
dat = plt.pcolormesh(time_grid_mean, alt_grid_mean, hcho_mean ,vmin = vmin, vmax = vmax, cmap = white_rainbow)
cb = plt.colorbar(dat, label = 'HCHO (ppb)', pad = 0.02)
plt.xlabel("Hour of Day (LT)")
plt.ylabel("Altitude (km)")

# Show plot
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("westerly_vp_mean", dpi = 600)
plt.show()