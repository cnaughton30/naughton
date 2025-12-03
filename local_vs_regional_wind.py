#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 13:27:42 2025

@author: cnaught
"""

import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
import glob
import datetime
import pytz
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.gridspec import GridSpec

pan_dir = r'/projectnb/atmchem/shared/for_claire' #move file into /project/atmchem/shared/ 
sys.path.insert(0, pan_dir)
from avg_wind_fn import wind_avg
from load_l2h_1_8 import pan_header, header_dict, l2h_to_df

sys.path.append("/usr3/graduate/cnaught/Desktop/pandora files")
from pan_functions import filter_pan_hcho_skyscan, subset_col, mol_per_m3_to_ppb, flatten_extend, filter_pan_hcho_ds

def load_isd_simple(directory, station_short_id):
    
    # Provide column names
    col_list = ['Year', 'Month', 'Day', 'Hour', 'airtemp', 'dewtemp', 
                'slp', 'winddir', 'windspd', 'sky',  'ppt1hr', 'ppt6hr']
    file_path = directory+station_short_id
    isd_df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=col_list)

    # Apply scale factors
    isd_df['airtemp'] = isd_df['airtemp']/10.
    isd_df['windspd'] = isd_df['windspd']/10.
    isd_df['dewtemp'] = isd_df['dewtemp']/10.
    isd_df['slp'] = isd_df['slp']/10.
    isd_df['ppt1hr'] = isd_df['ppt1hr']/10.

    # Clean up any spurious wind direction observations
    isd_df.loc[isd_df['winddir'] < 0, 'winddir'] = np.nan
    isd_df.loc[isd_df['winddir'] > 360, 'winddir'] = np.nan
    isd_df.loc[isd_df['windspd'] <= 0, 'winddir'] = np.nan

    isd_df.loc[isd_df['windspd'] > 20, 'windspd'] = np.nan
               
    # Create datetime column and make it the index
    cols = ['Year', 'Month', 'Day', 'Hour']
    isd_df['date'] = pd.to_datetime(isd_df[cols])
    isd_df.set_index('date', inplace=True)

    # Drop the redundant columns
    isd_df.drop(['Year', 'Month', 'Day', 'Hour'], axis=1, inplace=True)
    
    return isd_df

def convert_tz(df):
    boston_tz = pytz.timezone('America/New_York')
    df.index = df.index.tz_convert(boston_tz)
    df.index = df.index.tz_localize(None)
    df = df[(df.index.month >= 5) & (df.index.month <=9)]
    return(df)

from matplotlib.gridspec import GridSpec

def plot_windrose(df, site, label_pos, title, ax = None):
    label_pos = label_pos
    label_ang = 345
    
    # Set up a grid of plots
    gs = GridSpec(1, 4, width_ratios=[1,1,1,0.1])  # Add an extra column for the color bar
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(gs[0,0], projection='polar') 

    df['winddir_rad'] = np.deg2rad(df['winddir'])
    df['windspd'] = df['windspd']
    
    #Polar Plot:
    ax1.set_theta_offset(np.pi/2.0)
    ax1.set_theta_direction(-1)
    ax1.grid(True)
    rose1 = ax1.scatter(df['winddir_rad'], df['windspd'], c=df[site],
                        cmap='jet', alpha=0.8 , s=20, vmin = 0, vmax = 6)
    #r = np.arange(0,11,2)
    #ax1.set_rticks(r)
    ax1.text(345*np.pi/180, label_pos, title, ha='right')
    ax1.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])

    cbar = plt.colorbar(rose1, shrink = 0.3) #shrink function so the cbar fits
    cbar.set_label('HCHO (ppb)')
    #cbar.set_ticks(np.arange(0, 0.00065, 0.0001))
#%%
directory = '/usr3/graduate/cnaught/Desktop/wind data/'
station_short_id = '725090-14739-2021'
isd_df = load_isd_simple(directory, station_short_id)
isd_df.index = isd_df.index + pd.DateOffset(hours=-4)

directory = '/usr3/graduate/cnaught/Desktop/wind data/'
station_short_id = '725090-14739-2022'
isd_df22 = load_isd_simple(directory, station_short_id)
isd_df22.index = isd_df22.index + pd.DateOffset(hours=-4)

met = pd.concat([isd_df, isd_df22])
met = met.drop(columns = {'airtemp', 'dewtemp', 'slp', 'sky', 'ppt1hr','ppt6hr'})
met = met.dropna()
met = met.sort_index()
met['winddir'] = met['winddir'].astype(int)
met['windspd'] = met['windspd'].astype(int)
#%%
output_files = 'rfuh5p1-8.txt' #HCHO SS
Pan26_fp = '/projectnb/atmchem/shared/pandora/Pandora26s1_CambridgeMA_L2_' + output_files
header_content, col_n = pan_header(Pan26_fp)
Pan26_Data = l2h_to_df(Pan26_fp, header_content)
col_dict = header_dict(header_content, col_n) 
Pan26_Filtered_ss = filter_pan_hcho_skyscan(Pan26_Data)
cam_ss = convert_tz(Pan26_Filtered_ss)
cam_surf = subset_col(cam_ss, 45, 'surf')
t_index = pd.date_range(start=datetime.datetime(2021, 5, 1), end=datetime.datetime(2022, 9, 30), freq='60min')
df = pd.DataFrame(index = t_index )
df['surf'] = cam_surf.resample("60min").mean().reindex(t_index)
df['surf'] = mol_per_m3_to_ppb(df['surf']) #convert to ppb
df = df[(df['surf'] > 0) & (df['surf'] < 25)]
df = df.sort_index()
#%% 
met = met[(met.index.month > 4) & (met.index.month < 10)]
winddir_avgs = []
windspd_avgs = []

for i in met.index:
    #https://pandas.pydata.org/docs/reference/api/pandas.Timedelta.html
    start_time = i - pd.Timedelta(hours=3)
    end_time = i - pd.Timedelta(hours=1)
    print(i, start_time in met.index, end_time in met.index)
    if start_time in  met.index:
        if end_time in met.index:
            window = met.loc[start_time:end_time]
            print(window)
            spd = window.windspd
            direction = window.winddir
            dir_avg, spd_avg = wind_avg(direction, spd)
            winddir_avgs.append((i, dir_avg))
            windspd_avgs.append((i, spd_avg))
    else:
        continue
                          
#%%  
dir_df = pd.DataFrame(winddir_avgs, columns = ['time', 'winddir'])
dir_df = dir_df.set_index(dir_df['time']).drop(columns = 'time')
spd_df = pd.DataFrame(windspd_avgs, columns = ['time', 'windspd'])
spd_df = spd_df.set_index(spd_df['time']).drop(columns = 'time')
avg_wind = pd.merge(dir_df, spd_df, left_index=True, right_index=True, how='inner')
#%%
#winddir and windspd are an average from 3 hours prior
#surface conc is  the measurement from the hour
fdf = pd.merge(df, avg_wind, left_index=True, right_index=True)
elevated_hcho = fdf[fdf['surf'] >= 4]
local = elevated_hcho[(elevated_hcho['windspd'] <= 1.5)|((elevated_hcho['winddir'] > 90) & (elevated_hcho['winddir'] <135))]
long_range = elevated_hcho[
    (elevated_hcho['windspd'] > 1.5) &
    ((elevated_hcho['winddir'] <= 90) | (elevated_hcho['winddir'] >= 135))]
#%%
dir_vec = np.arange(0, 361, 20)
spd_vec = np.arange(0, 12, 2)


pv = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0]])*0.0
kv = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0]])*0.0


# Want to keep track of the count in each coordinate for averaging
count_pv = np.zeros([np.shape(spd_vec)[0], np.shape(dir_vec)[0]]).astype(int)
count_kv = np.zeros([np.shape(spd_vec)[0], np.shape(dir_vec)[0]]).astype(int)


# Keep a 3-d array of all the matching points to do population statistics later
pv_all = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0], 100])*0.0
kv_all = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0], 200])*0.0

ss = local['windspd']
dd = local['winddir']
pan = local['surf']

#Go through each data point one at a time, and assign it a gridded coordinate in "direction" x "speed" space
for i in range(np.shape(dd)[0]):
    # What is the wind direction and speed of the data point?
    wd = dd[i]
    ws = ss[i]
    pl = pan[i]
    # If the Pandora and wind data are real, locate appropriate direction x speed grid coordinate of the data point.
    if (~np.isnan(pl)) & (~np.isnan(wd)):
        ii = np.argmin(np.abs(dir_vec - wd))
        jj = np.argmin(np.abs(spd_vec - ws))
        # At those coordinates, store the pollutant data and the count
        pv[jj, ii] = pv[jj, ii] + pl
        pv_all[jj, ii, count_pv[jj, ii]] = pl
        count_pv[jj, ii] = count_pv[jj, ii] + 1
 
pv_all_2019 = pv_all
pv_all_2019[pv_all_2019 == 0] = np.nan
pv_avg_2019 = pv/count_pv

#%%
# Do a test plot of Pandora data?
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_offset(np.pi/2.0)
ax.set_theta_direction(-1)
rose = ax.pcolormesh(dir_vec*np.pi/180.0, spd_vec, pv_avg_2019, cmap='jet', vmin=0, vmax=6, shading = 'auto')
fig.colorbar(rose, label = 'HCHO (ppb)')
plt.savefig("local_hcho_wr", dpi = 600)
plt.show()
#%%
dir_vec = np.arange(0, 361, 20)
spd_vec = np.arange(0, 12, 2)

pv = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0]])*0.0
kv = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0]])*0.0

# Want to keep track of the count in each coordinate for averaging
count_pv = np.zeros([np.shape(spd_vec)[0], np.shape(dir_vec)[0]]).astype(int)
count_kv = np.zeros([np.shape(spd_vec)[0], np.shape(dir_vec)[0]]).astype(int)

# Keep a 3-d array of all the matching points to do population statistics later
pv_all = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0], 100])*0.0
kv_all = np.empty([np.shape(spd_vec)[0], np.shape(dir_vec)[0], 200])*0.0

ss = long_range['windspd']
dd = long_range['winddir']
pan = long_range['surf']

#Go through each data point one at a time, and assign it a gridded coordinate in "direction" x "speed" space
for i in range(np.shape(dd)[0]):
    # What is the wind direction and speed of the data point?
    wd = dd[i]
    ws = ss[i]
    pl = pan[i]
    # If the Pandora and wind data are real, locate appropriate direction x speed grid coordinate of the data point.
    if (~np.isnan(pl)) & (~np.isnan(wd)):
        ii = np.argmin(np.abs(dir_vec - wd))
        jj = np.argmin(np.abs(spd_vec - ws))
        # At those coordinates, store the pollutant data and the count
        pv[jj, ii] = pv[jj, ii] + pl
        pv_all[jj, ii, count_pv[jj, ii]] = pl
        count_pv[jj, ii] = count_pv[jj, ii] + 1
 
pv_all_2019 = pv_all
pv_all_2019[pv_all_2019 == 0] = np.nan

pv_avg_2019 = pv/count_pv
#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_offset(np.pi/2.0)
ax.set_theta_direction(-1)
rose = ax.pcolormesh(dir_vec*np.pi/180.0, spd_vec, pv_avg_2019, cmap='jet', vmin=0, vmax=6, shading = 'auto')
fig.colorbar(rose, label = 'HCHO (ppb)')
#plt.savefig("local_hcho_wr", dpi = 600)
plt.show()
#%%
# Means, medians, std
print("Local stats:\n", local['surf'].describe())
print("Long-range stats:\n", long_range['surf'].describe())


print("Local stats:\n", local['winddir'].describe())
print("Long-range stats:\n", long_range['winddir'].describe())


print("Local stats:\n", local['windspd'].describe())
print("Long-range stats:\n", long_range['windspd'].describe())
#%%
plt.hist(long_range['surf'], bins=30, alpha=0.5, label='Long-range')
plt.hist(local['surf'], bins=30, alpha=0.5, label='Local')
plt.xlabel('HCHO (ppb)')
plt.ylabel('Frequency')
plt.legend()
plt.show()

plt.boxplot([local['surf'], long_range['surf']], labels=['Local', 'Long-range'])
plt.ylabel('HCHO')
plt.show()

from scipy.stats import ttest_ind

t_stat, p_val = ttest_ind(local['surf'], long_range['surf'], equal_var=False)
print("t-stat:", t_stat, "p-value:", p_val)

from scipy.stats import ks_2samp

stat, p_val = ks_2samp(local['surf'], long_range['surf'])
print("KS stat:", stat, "p-value:", p_val)

