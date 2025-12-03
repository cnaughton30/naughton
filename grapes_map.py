#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 20:45:57 2025

@author: cnaught
"""

import xarray as xr
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
import geopandas as gpd
import numpy as np
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import LinearSegmentedColormap

def plot_individual_map(df, sector_name, fig_title):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    df = df.sum(dim = 'Time')

    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    lat = df['XLAT']
    lon = df['XLONG']
    
    #Create a custom colormap
    #https://matplotlib.org/stable/gallery/color/named_colors.html
    colors = ['white', 'skyblue', 'mediumseagreen', 'yellow', 'darkorange', 'red', 'darkred']
    #https://stackoverflow.com/questions/67605719/displaying-lowest-values-as-white
    white_rainbow = LinearSegmentedColormap.from_list("white_rainbow", colors)

    mesh = ax.pcolormesh(lon, lat, df, transform=ccrs.PlateCarree(),
    cmap= white_rainbow, vmin = 0, vmax = 45)

    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, orientation="vertical", shrink=0.7, pad=0.05)
    cbar.set_label(f"HCHO {sector_name} Emissions (mole km$^{-2}$ hr$^{-1}$)")
    ax.set_extent([-73.5, -69.5, 41.0, 43.5], crs=ccrs.PlateCarree())

    # Define the Boston region (same as you used for subsetting)
    lat_min, lat_max = 42.15, 42.55
    lon_min, lon_max = -71.3, -70.9

    # Add rectangle to map
    rect = Rectangle(
        (lon_min, lat_min),            # lower-left corner
        lon_max - lon_min,             # width
        lat_max - lat_min,             # height
        linewidth=2, edgecolor='red', facecolor='none', zorder=3, linestyle="--")

    ax.add_patch(rect)
    plt.savefig(fig_title, dpi = 1000)
    plt.show()

#%% Load Major GRAPES data Pt. 1
#onroad diesel
path = "/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_DSL_202107/202107/weekdy/"
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_weekdy_12to23Z.nc")]
or_diesel = xr.open_mfdataset(files, combine='by_coords')

#onroad gas
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_GAS_202107/202107/weekdy'
# List of your two files
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_weekdy_12to23Z.nc")]

or_gas = xr.open_mfdataset(files, combine='by_coords')

#Residental
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RES_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_RES_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RES_202107_weekdy_12to23Z.nc")]
res_e = xr.open_mfdataset(files, combine='by_coords')

#Offroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_OFFROAD_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_weekdy_12to23Z.nc")]
offroad_e = xr.open_mfdataset(files, combine='by_coords')

#Commercial
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_COMM_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_COMM_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_COMM_202107_weekdy_12to23Z.nc")]
comm = xr.open_mfdataset(files, combine='by_coords')

#Industrial Fuel
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_INDF_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_INDF_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_INDF_202107_weekdy_12to23Z.nc")]
indf = xr.open_mfdataset(files, combine='by_coords')

#Railroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RAIL_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_weekdy_12to23Z.nc")]
rail = xr.open_mfdataset(files, combine='by_coords')

#%% Load Minor GRAPES data Pt2
#Aviation
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_AVIATION_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_AVIATION_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_AVIATION_202107_weekdy_12to23Z.nc")]
aviation = xr.open_mfdataset(files, combine='by_coords')

#Shipping
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_SHIPPING_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_SHIPPING_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_SHIPPING_202107_weekdy_12to23Z.nc")]
shipping = xr.open_mfdataset(files, combine='by_coords')

#EGU
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_EGU_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_EGU_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_EGU_202107_weekdy_12to23Z.nc")]
egu = xr.open_mfdataset(files, combine='by_coords')

#Oil and Gas
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_OG_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_OG_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_OG_202107_weekdy_12to23Z.nc")]
og = xr.open_mfdataset(files, combine='by_coords')

#FUG
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_FUG_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_FUG_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_FUG_202107_weekdy_12to23Z.nc")]
fug = xr.open_mfdataset(files, combine='by_coords')

#Waste
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_WASTE_202107/202107/weekdy'
files = [
    os.path.join(path, "GRA2PESv1.0_WASTE_202107_weekdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_WASTE_202107_weekdy_12to23Z.nc")]
waste = xr.open_mfdataset(files, combine='by_coords')
#%% Subset data to the surface layer
ship_hcho = shipping.HC14
ship_hcho_surf = ship_hcho.isel(bottom_top=0)
ord_hcho = or_diesel.HC14
ord_hcho_surf = ord_hcho.isel(bottom_top=0)
org_hcho = or_gas.HC14
org_hcho_surf = org_hcho.isel(bottom_top=0)
res_hcho = res_e.HC14
res_hcho_surf = res_hcho.isel(bottom_top=0)
offroad_hcho = offroad_e.HC14
offroad_hcho_surf = offroad_hcho.isel(bottom_top=0)
aviation_hcho = aviation.HC14
aviation_hcho_surf = aviation_hcho.isel(bottom_top=0)
comm_hcho = comm.HC14
comm_hcho_surf = comm_hcho.isel(bottom_top=0)
egu_hcho = egu.HC14
egu_hcho_surf = egu_hcho.isel(bottom_top=0)
fug_hcho = fug.HC14
fug_hcho_surf = fug_hcho.isel(bottom_top=0)
og_hcho = og.HC14
og_hcho_surf = og_hcho.isel(bottom_top=0)
indf_hcho = indf.HC14
indf_hcho_surf = indf_hcho.isel(bottom_top=0)
rail_hcho = rail.HC14
rail_hcho_surf = rail_hcho.isel(bottom_top=0)
waste_hcho = waste.HC14
waste_hcho_surf = waste_hcho.isel(bottom_top=0)

#Add major weekday emissions
total_wkdy = ord_hcho_surf + org_hcho_surf + res_hcho_surf + offroad_hcho_surf + comm_hcho_surf + indf_hcho_surf + rail_hcho_surf
total_hcho_wkdy = total_wkdy.sum(dim = 'Time')
#%% Pie Chart of Emissions in Pandora Grid Box
#lat, lon = 42.38, -71.11
#lat, lon = 42.36, -71.06 # upper right
#lat, lon = 42.33, -71.08 # lower right
lat, lon = 42.33, -71.12 # lower left

#https://stackoverflow.com/questions/58758480/xarray-select-nearest-lat-lon-with-multi-dimension-coordinates?utm_source=chatgpt.com
abslat = np.abs(ord_hcho_surf.XLAT - lat)
abslon = np.abs(ord_hcho_surf.XLONG - lon)
c = np.maximum(abslon, abslat)
([iy], [ix]) = np.where(c == np.min(c))

# Select the nearest grid cell
ord_point = ord_hcho_surf.isel(south_north=iy, west_east=ix)
org_point = org_hcho_surf.isel(south_north=iy, west_east=ix)
res_point = res_hcho_surf.isel(south_north=iy, west_east=ix)
offroad_point = offroad_hcho_surf.isel(south_north=iy, west_east=ix)
comm_point = comm_hcho_surf.isel(south_north=iy, west_east=ix)
indf_point = indf_hcho_surf.isel(south_north=iy, west_east=ix)
rail_point = rail_hcho_surf.isel(south_north=iy, west_east=ix)

emissions_totals = {
    "Onroad Gas": float(org_point.sum(dim="Time").values),
    "Onroad Diesel": float(ord_point.sum(dim="Time").values),
    "Residential": float(res_point.sum(dim="Time").values),
    "Offroad": float(offroad_point.sum(dim="Time").values),
    "Commercial": float(comm_point.sum(dim="Time").values),
    "Industrial": float(indf_point.sum(dim="Time").values),
    "Rail": float(rail_point.sum(dim="Time").values),
}

# Pie chart
plt.figure(figsize=(6,6))
plt.pie(emissions_totals.values(), labels=emissions_totals.keys(),
        autopct='%1.1f%%', startangle=90)
plt.title(f"HCHO Surface Emission Breakdown near ({lat}, {lon})")
plt.savefig("wkdy_emissions_pie_lower_left", dpi = 600)
plt.show()
#%% Saturday Major Emissions
#onroad diesel
path = "/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_DSL_202107/202107/satdy/"
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_satdy_12to23Z.nc")]
or_diesel = xr.open_mfdataset(files, combine='by_coords')

#onroad gas
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_GAS_202107/202107/satdy'
# List of your two files
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_satdy_12to23Z.nc")]
or_gas = xr.open_mfdataset(files, combine='by_coords')

#Residental
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RES_202107/202107/satdy'
files = [
    os.path.join(path, "GRA2PESv1.0_RES_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RES_202107_satdy_12to23Z.nc")]
res_e = xr.open_mfdataset(files, combine='by_coords')

#Offroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_OFFROAD_202107/202107/satdy'
files = [
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_satdy_12to23Z.nc")]
offroad_e = xr.open_mfdataset(files, combine='by_coords')

#Commercial
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_COMM_202107/202107/satdy'
files = [
    os.path.join(path, "GRA2PESv1.0_COMM_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_COMM_202107_satdy_12to23Z.nc")]
comm = xr.open_mfdataset(files, combine='by_coords')

#Industrial Fuel
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_INDF_202107/202107/satdy'
files = [
    os.path.join(path, "GRA2PESv1.0_INDF_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_INDF_202107_satdy_12to23Z.nc")]
indf = xr.open_mfdataset(files, combine='by_coords')

#Railroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RAIL_202107/202107/satdy'
files = [
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_satdy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_satdy_12to23Z.nc")]
rail = xr.open_mfdataset(files, combine='by_coords')

ord_hcho = or_diesel.HC14
ord_hcho_surf = ord_hcho.isel(bottom_top=0)
org_hcho = or_gas.HC14
org_hcho_surf = org_hcho.isel(bottom_top=0)
res_hcho = res_e.HC14
res_hcho_surf = res_hcho.isel(bottom_top=0)
offroad_hcho = offroad_e.HC14
offroad_hcho_surf = offroad_hcho.isel(bottom_top=0)
comm_hcho = comm.HC14
comm_hcho_surf = comm_hcho.isel(bottom_top=0)
indf_hcho = indf.HC14
indf_hcho_surf = indf_hcho.isel(bottom_top=0)
rail_hcho = rail.HC14
rail_hcho_surf = rail_hcho.isel(bottom_top=0)

total_sat = ord_hcho_surf + org_hcho_surf + res_hcho_surf + offroad_hcho_surf + comm_hcho_surf + indf_hcho_surf + rail_hcho_surf
total_sat_hcho = total_sat.sum(dim = 'Time')
#%% Sunday Major Emissions
#onroad diesel
path = "/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_DSL_202107/202107/sundy/"
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_DSL_202107_sundy_12to23Z.nc")]
or_diesel = xr.open_mfdataset(files, combine='by_coords')

#onroad gas
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_ONROAD_GAS_202107/202107/sundy'
# List of your two files
files = [
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_ONROAD_GAS_202107_sundy_12to23Z.nc")]
or_gas = xr.open_mfdataset(files, combine='by_coords')

#Residental
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RES_202107/202107/sundy'
files = [
    os.path.join(path, "GRA2PESv1.0_RES_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RES_202107_sundy_12to23Z.nc")]
res_e = xr.open_mfdataset(files, combine='by_coords')

#Offroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_OFFROAD_202107/202107/sundy'
files = [
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_OFFROAD_202107_sundy_12to23Z.nc")]
offroad_e = xr.open_mfdataset(files, combine='by_coords')

#Commercial
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_COMM_202107/202107/sundy'
files = [
    os.path.join(path, "GRA2PESv1.0_COMM_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_COMM_202107_sundy_12to23Z.nc")]
comm = xr.open_mfdataset(files, combine='by_coords')

#Industrial Fuel
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_INDF_202107/202107/sundy'
files = [
    os.path.join(path, "GRA2PESv1.0_INDF_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_INDF_202107_sundy_12to23Z.nc")]
indf = xr.open_mfdataset(files, combine='by_coords')

#Railroad
path = '/projectnb/atmchem/naughton/GRAPES/GRA2PESv1.0_RAIL_202107/202107/sundy'
files = [
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_sundy_00to11Z.nc"),
    os.path.join(path, "GRA2PESv1.0_RAIL_202107_sundy_12to23Z.nc")]
rail = xr.open_mfdataset(files, combine='by_coords')

ord_hcho = or_diesel.HC14
ord_hcho_surf = ord_hcho.isel(bottom_top=0)
org_hcho = or_gas.HC14
org_hcho_surf = org_hcho.isel(bottom_top=0)
res_hcho = res_e.HC14
res_hcho_surf = res_hcho.isel(bottom_top=0)
offroad_hcho = offroad_e.HC14
offroad_hcho_surf = offroad_hcho.isel(bottom_top=0)
comm_hcho = comm.HC14
comm_hcho_surf = comm_hcho.isel(bottom_top=0)
indf_hcho = indf.HC14
indf_hcho_surf = indf_hcho.isel(bottom_top=0)
rail_hcho = rail.HC14
rail_hcho_surf = rail_hcho.isel(bottom_top=0)

total_sun = ord_hcho_surf + org_hcho_surf + res_hcho_surf + offroad_hcho_surf + comm_hcho_surf + indf_hcho_surf + rail_hcho_surf
total_sun_hcho = total_sun.sum(dim = 'Time')

plot_individual_map(total_sun, 'Sunday', 'sunday_hcho_emissions')
#%% Total Weekday, Sat and Sun Emissions
weekdy_july = total_hcho_wkdy * 22 #22 weekdays in July 2021
sat_july = total_sat_hcho * 5 #5 Saturdays in July 2021
sun_july = total_sun_hcho * 4 #4 Sundays in July 2021
total_july = (weekdy_july + sat_july + sun_july)/24
#%% Plot Regional Emissions Map - July 2021
plt.rcParams['font.family'] = 'Nimbus Sans'
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.STATES, linewidth=0.75, edgecolor='black')
ax.add_feature(cfeature.COASTLINE,linewidth=0.75, edgecolor='black')

roadshp = '/projectnb/atmchem/naughton/Road shapefiles/tl_2021_us_primaryroads.shp'
gdf = gpd.read_file(roadshp)
gdf.plot(ax=ax, transform=ccrs.PlateCarree(), linewidth=0.4, color="black", alpha = 0.75 )

lat = total_july['XLAT']
lon = total_july['XLONG']

#Create a custom colormap
#https://matplotlib.org/stable/gallery/color/named_colors.html
colors = ['white', 'skyblue', 'mediumseagreen', 'yellow', 'darkorange', 'red', 'darkred']
#https://stackoverflow.com/questions/67605719/displaying-lowest-values-as-white
white_rainbow = LinearSegmentedColormap.from_list("white_rainbow", colors)

mesh = ax.pcolormesh(lon, lat, total_july, transform=ccrs.PlateCarree(),
    cmap=white_rainbow , vmin = 0, vmax = 45)

# Colorbar
cbar = plt.colorbar(mesh, ax=ax, orientation="vertical", shrink = 0.825, aspect = 15)
cbar.set_label("Total HCHO Emissions (mole km$^{-2}$ hr$^{-1}$)")
cbar.set_ticks(np.arange(0,46,5))
ax.set_extent([-73.5, -69.5, 41.0, 43.5], crs=ccrs.PlateCarree())

lat_min, lat_max = 42.2, 42.55
lon_min, lon_max = -71.3, -70.9

# Add rectangle around Boston
rect = Rectangle(
    (lon_min, lat_min),            # lower-left corner
    lon_max - lon_min,             # width
    lat_max - lat_min,             # height
    linewidth=2, edgecolor='red', facecolor='none', zorder=3, linestyle="--"
)
lon, lat = -71.11, 42.38  
ax.add_patch(rect)
plt.savefig("total_july_hcho", dpi = 600)
plt.show()

#%% Plot Local Emissions Map - July 2021
plt.rcParams['font.family'] = 'Nimbus Sans'
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

ax.add_feature(cfeature.GSHHSFeature(scale='f'), linewidth = 1)

roadshp = '/projectnb/atmchem/naughton/Road shapefiles/tl_2021_us_primaryroads.shp'
gdf = gpd.read_file(roadshp)
gdf.plot(ax=ax, transform=ccrs.PlateCarree(), linewidth=0.4, color="black", alpha = 0.75 )

lat = total_july['XLAT']
lon = total_july['XLONG']

mesh = ax.pcolormesh(lon, lat, total_july, transform=ccrs.PlateCarree(),
    cmap= white_rainbow, vmin = 0, vmax = 45)

# Colorbar
cbar = plt.colorbar(mesh, ax=ax, orientation="vertical", aspect = 15)
cbar.set_label("Total HCHO Emissions (mole km$^{-2}$ hr$^{-1}$)")
cbar.set_ticks(np.arange(0,46,5))
ax.set_extent([-71.3, -70.9, 42.2, 42.55], crs=ccrs.PlateCarree())

lon, lat = -71.11, 42.38  # Cambridge Pandora
ax.plot(lon, lat,
        marker='o', color='white', 
        markeredgecolor='black', markeredgewidth=2, 
        markersize=10,
        transform=ccrs.PlateCarree())  

#lon, lat = -71.1087, 42.37  #Center of Pandora Gridbox

#https://stackoverflow.com/questions/17086847/box-around-text-in-matplotlib
ax.annotate('Cambridge Pandora', (-71.11, 42.395), ha='center', weight='bold',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.9))

#ax.set_title("GRAPES Total Weekday Emissions ")
plt.savefig("total_hcho_boston", dpi = 1000)
plt.show()
#%%

#Plot Diurnal Averages in the Box area

#%% Subset data to only Boston area
lat_min, lat_max = 42.2, 42.55
lon_min, lon_max = -71.3, -70.9

# Subset grid cells in this box
boston_box = ord_hcho_surf.where(
    (ord_hcho_surf.XLAT >= lat_min) & (ord_hcho_surf.XLAT <= lat_max) &
    (ord_hcho_surf.XLONG >= lon_min) & (ord_hcho_surf.XLONG <= lon_max),
    drop=True
)
boston_box_ship = ship_hcho_surf.where(
    (ord_hcho_surf.XLAT >= lat_min) & (ord_hcho_surf.XLAT <= lat_max) &
    (ord_hcho_surf.XLONG >= lon_min) & (ord_hcho_surf.XLONG <= lon_max),
    drop=True
)
boston_box_org = org_hcho_surf.where(
    (ord_hcho_surf.XLAT >= lat_min) & (ord_hcho_surf.XLAT <= lat_max) &
    (ord_hcho_surf.XLONG >= lon_min) & (ord_hcho_surf.XLONG <= lon_max),
    drop=True
)
boston_box_res = res_hcho_surf.where(
    (ord_hcho_surf.XLAT >= lat_min) & (ord_hcho_surf.XLAT <= lat_max) &
    (ord_hcho_surf.XLONG >= lon_min) & (ord_hcho_surf.XLONG <= lon_max),
    drop=True
)
boston_box_offroad = offroad_hcho_surf.where(
    (ord_hcho_surf.XLAT >= lat_min) & (ord_hcho_surf.XLAT <= lat_max) &
    (ord_hcho_surf.XLONG >= lon_min) & (ord_hcho_surf.XLONG <= lon_max),
    drop=True
)

av_box = aviation_hcho_surf.where(
    (aviation_hcho_surf.XLAT >= lat_min) & (aviation_hcho_surf.XLAT <= lat_max) &
    (aviation_hcho_surf.XLONG >= lon_min) & (aviation_hcho_surf.XLONG <= lon_max),
    drop=True
)
comm_box = comm_hcho_surf.where(
    (comm_hcho_surf.XLAT >= lat_min) & (comm_hcho_surf.XLAT <= lat_max) &
    (comm_hcho_surf.XLONG >= lon_min) & (comm_hcho_surf.XLONG <= lon_max),
    drop=True
)
egu_box = egu_hcho_surf.where(
    (egu_hcho_surf.XLAT >= lat_min) & (egu_hcho_surf.XLAT <= lat_max) &
    (egu_hcho_surf.XLONG >= lon_min) & (egu_hcho_surf.XLONG <= lon_max),
    drop=True
)
og_box = og_hcho_surf.where(
    (og_hcho_surf.XLAT >= lat_min) & (og_hcho_surf.XLAT <= lat_max) &
    (og_hcho_surf.XLONG >= lon_min) & (og_hcho_surf.XLONG <= lon_max),
    drop=True
)
fug_box = fug_hcho_surf.where(
    (fug_hcho_surf.XLAT >= lat_min) & (fug_hcho_surf.XLAT <= lat_max) &
    (fug_hcho_surf.XLONG >= lon_min) & (fug_hcho_surf.XLONG <= lon_max),
    drop=True
)
indf_box = indf_hcho_surf.where(
    (indf_hcho_surf.XLAT >= lat_min) & (indf_hcho_surf.XLAT <= lat_max) &
    (indf_hcho_surf.XLONG >= lon_min) & (indf_hcho_surf.XLONG <= lon_max),
    drop=True)
rail_box = rail_hcho_surf.where(
    (rail_hcho_surf.XLAT >= lat_min) & (rail_hcho_surf.XLAT <= lat_max) &
    (rail_hcho_surf.XLONG >= lon_min) & (rail_hcho_surf.XLONG <= lon_max),
    drop=True)
waste_box = waste_hcho_surf.where(
    (waste_hcho_surf.XLAT >= lat_min) & (waste_hcho_surf.XLAT <= lat_max) &
    (waste_hcho_surf.XLONG >= lon_min) & (waste_hcho_surf.XLONG <= lon_max),
    drop=True)

#%% Calculate diurnal averages
diurnal_cycle = boston_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
ship_cycle = boston_box_ship.groupby("Time.hour").mean(dim=["south_north","west_east"])
org_cycle = boston_box_org.groupby("Time.hour").mean(dim=["south_north","west_east"])
res_cycle = boston_box_res.groupby("Time.hour").mean(dim=["south_north","west_east"])
offroad_cycle = boston_box_offroad.groupby("Time.hour").mean(dim=["south_north","west_east"])
comm_cycle = comm_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
av_cycle = av_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
og_cycle = og_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
egu_cycle = egu_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
fug_cycle = fug_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
indf_cycle = indf_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
rail_cycle = rail_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
waste_cycle = waste_box.groupby("Time.hour").mean(dim=["south_north","west_east"])
#%%Plot diurnal averages
plt.figure(figsize=(8,5))
#Major emissions
plt.plot(diurnal_cycle, marker = 'o', linestyle = '-', label = "Onroad Diesel")
plt.plot(org_cycle, marker = 'o', linestyle = '-', label = "Onroad Gas")
plt.plot(res_cycle, marker = 'o', linestyle = '-', label = "Residential")
plt.plot(offroad_cycle, marker = 'o', linestyle = '-', label = "Offroad")
plt.plot(comm_cycle, marker = 'o', linestyle = '-', label = "Commercial")
plt.plot(indf_cycle, marker = 'o', linestyle = '-', label = "Industrial Fuel")
plt.plot(rail_cycle, marker = 'o', linestyle = '-', label = "Railroad")

#Minor emissions
#plt.plot(ship_cycle, marker = 'o', linestyle = '-', label = "Shipping")
#plt.plot(av_cycle, marker = 'o', linestyle = '-', label = "Aviation")
#plt.plot(og_cycle, marker = 'o', linestyle = '-', label = "Oil and Gas")
#plt.plot(egu_cycle, marker = 'o', linestyle = '-', label = "Electricity Gen")
#plt.plot(fug_cycle, marker = 'o', linestyle = '-', label = "FUG")
#plt.plot(waste_cycle, marker = 'o', linestyle = '-', label = "Waste")

plt.ylabel("HCHO Emissions (mole km$^{-2}$ hr$^{-1}$)")
plt.xlabel("UTC Hour")
plt.title("Diurnal Primary HCHO Emissions")
plt.legend(loc = 'upper left')
plt.grid(True)
plt.savefig('diurnal_grapes_minor', dpi = 600)
plt.show()
#%% Sum of Diurnal average hourly emissions
emissions_totals = {
    "Onroad Gas": float(org_cycle.sum(dim="Time").compute().values),
    "Onroad Diesel": float(diurnal_cycle.sum(dim="Time").compute().values),
    "Residential": float(res_cycle.sum(dim="Time").compute().values),
    "Offroad": float(offroad_cycle.sum(dim="Time").compute().values),
    "Comm": float(comm_cycle.sum(dim="Time").compute().values),
    "Industrial": float(indf_cycle.sum(dim="Time").compute().values),
    "Rail": float(rail_cycle.sum(dim="Time").compute().values),
}
#%% Pie Chart of Major Emissions
plt.figure(figsize=(7,7))
plt.pie(
    emissions_totals.values(),
    labels=emissions_totals.keys(),
    autopct="%1.1f%%",
    startangle=90,
)
plt.title("Total Diurnal Emissions Contribution by Sector")
plt.savefig("pie_major_emissions", dpi = 600)
plt.show()