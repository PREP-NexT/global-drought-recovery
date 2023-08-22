# %%
import netCDF4 as nc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors



# %%
# prob_seasonal_recovery probability
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)),
    )
    return new_cmap


def add_south(data):
    data_nan = np.zeros((120, 1440))
    data_nan[:] = np.nan
    datanew = np.r_[data, data_nan]
    return datanew


def reverse_mapdata(inputdata):
    prob_fli = np.flipud(inputdata[:, :, :, :])
    data_map = np.concatenate((prob_fli[:, 720:1440, :, :], prob_fli[:, 0:720, :, :]), axis=1)
    return data_map


def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i + 1), va="center", ha="center")
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)



# %%
filename_delta = 'J:\\output\\map_delta_sign.npy'
filename_his = 'J:\\output\\map_his_climatology.npy'
filename_pres = 'J:\\output\\map_pres_climatology.npy'
region = nc.Dataset(r'region.id.nc', 'r')
RID = region.variables['rid'][:].data
RID = np.flipud(RID)
cbformat = mticker.ScalarFormatter()
cbformat.set_powerlimits((-2, 2))  # Set size thresholds for scientific notation.（10^(-2)~10^2)
lat = np.arange(89.875, -90, -0.25)
lon = np.arange(-179.875, 180, 0.25)
data_delta_sign = np.load(filename_delta)
data_his = np.load(filename_his)
data_pres = np.load(filename_pres)

RT = ['1-7', '8-14', '15-21', '22-28']
sea = ['MAM', 'JJA', 'SON', 'DJF']
# for i_rt in range(4):
i_rt = 1
title = 'RT=' + RT[i_rt]

# Define the figure and each axis for the 3 rows and 3 columns
fig, axs = plt.subplots(nrows=2, ncols=3,
                        subplot_kw={'projection': ccrs.Robinson()},  # Robinson()},  # PlateCarree()},
                        figsize=(11, 9))

axs = axs.flatten()
num = -1
# Loop over all of the models
for i_sea in range(4):
    for col in range(3):
        num += 1
        cmap = plt.get_cmap('GnBu')
        new_cmap = truncate_colormap(cmap, minval=0.2, maxval=1.0, n=100)  
        bounds = np.linspace(0, 80, 9)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

        if col == 0:
            data = data_his[:, :, i_rt, i_sea]
        elif col == 1:
            data = data_pres[:, :, i_rt, i_sea]
        elif col == 2:

            data = data_delta_sign[:, :, i_rt, i_sea]
            colors_below = plt.cm.RdBu(np.linspace(0, 0.45, 256))
            colors_over = plt.cm.RdBu(np.linspace(0.55, 1, 256))
            all_colors = np.vstack((colors_below, colors_over))
            new_cmap = colors.LinearSegmentedColormap.from_list(
                'RdBu_map', all_colors)
            bounds = np.linspace(-60, 60, 7)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256) #, extend='both')
        cmap_nan = ListedColormap(['#FFFFFF'])
        norm_nan = mpl.colors.Normalize(vmin=0, vmax=2)
        np_condition=np.where((np.isnan(data)) & (RID > 0) & (RID < 11))
        data_nan=np.empty((600,1440))
        data_nan[:]=np.nan
        data_nan[np_condition]=1

        datanew = add_south(data) * 100
        data_nan_new=add_south(data_nan)
        if col == 2:
            cs = axs[num].contourf(lon, lat, datanew, cmap=new_cmap, levels=bounds, extend='both',
                                   transform=ccrs.PlateCarree())
        else:
            cs = axs[num].contourf(lon, lat, datanew, cmap=new_cmap, levels=bounds, extend='max',
                                   transform=ccrs.PlateCarree())

        # Draw the coastines for each subplot
        axs[num].coastlines(lw=0.2, color='#696969')
        axs[num].add_feature(cfeat.LAND.with_scale('110m'), facecolor='#CDCDCD', alpha=0.8, edgecolor='#CDCDCD',
                             linewidth=0.01)  # ,edgecolor='')#color='#D3D3D3'
        gl = axs[num].gridlines(draw_labels=False, lw=0.00007, color='#C0C0C0', linestyle='--', alpha=0.5)

        if i_sea == 3 & col == 1:
            cbar_ax = fig.add_axes([0.154, 0.12, 0.45, 0.02])
            cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal', extend='max')

            for l in cbar.ax.xaxis.get_ticklabels():
                l.set_family('Myriad Pro')
            for l in cbar.ax.yaxis.get_ticklabels():
                l.set_family('Myriad Pro')
            cbar.ax.tick_params(labelsize=10)
            plt.rcParams['font.family'] = 'Myriad Pro'
            plt.rcParams['font.size'] = 11

        elif (i_sea == 3) & (col == 2):
            cbar_ax2 = fig.add_axes([0.6972, 0.12, 0.23, 0.02])
            cbar2 = fig.colorbar(cs, cax=cbar_ax2, orientation='horizontal', extend='both')
            for l in cbar2.ax.xaxis.get_ticklabels():
                l.set_family('Myriad Pro')
            for l in cbar2.ax.yaxis.get_ticklabels():
                l.set_family('Myriad Pro')
            cbar2.ax.tick_params(labelsize=10)
            plt.rcParams['font.family'] = 'Myriad Pro'
            plt.rcParams['font.size'] = 11

        axs[num].set_global()
# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.18, left=0.1, right=0.95,
                    wspace=0.03, hspace=0.05)  # wspace: width

# Add a big title at the top
plt.suptitle(title)
plt.savefig(
    'J:\\output\\results_1\\ ' + title + 'legend.jpg',
    dpi=1000, bbox_inches='tight')

plt.close("all")


#%%
#plot_bar
filename1 = 'J:\\output\\global_clima_his_clima_mean.npy'
filename2 = 'J:\\output\\global\\global_clima_pres_clima_mean.npy'
filename3 = 'J:\\output\\global\\global_clima_pvalue_ks_mask.npy'

data_his0 = np.load(filename1)
data_his0[(data_his0 == -111) | (data_his0 == -999)] = np.nan
data_his0[(data_his0 == 0)] = 0.001
data_his0[(data_his0 == 1)] = 0.999

data_pres0 = np.load(filename2)
data_pres0[(data_pres0 == -111) | (data_pres0 == -999)] = np.nan
data_pres0[(data_pres0 == 0)] = 0.001
data_pres0[(data_pres0 == 1)] = 0.999

data_sign0 = np.load(filename3)

data_his = reverse_mapdata(data_his0)
data_pres = reverse_mapdata(data_pres0)
data_sign = reverse_mapdata(data_sign0)

data_delta = data_pres-data_his

data_delta_sign = np.zeros((600, 1440, 4, 4))

data_delta_sign[:] = np.nan
data_delta_sign[data_sign == 1] = data_delta[data_sign == 1]
region = nc.Dataset(r'E:\phd\data\climate region\region.id.nc', 'r')
RID = region.variables['rid'][:].data
RID = np.flipud(RID)
share_rid = np.zeros((10, 2, 4, 4))  # fraction
share_rid[:] = np.nan
median_rid = np.zeros((10, 2, 4, 4))  # medium
median_rid[:] = np.nan
median_rid_ratio = np.zeros((10, 2, 4, 4))
median_rid_ratio[:] = np.nan

share_global = np.zeros((4,2))
share_global[:] = np.nan
median_global= np.zeros((4,2))
median_global[:]=np.nan
for i_rt in range(4):
    data_rid_all=np.zeros(0)
    for i_sea in range(4):
        data_delta_sign_1 = data_delta_sign[:, :, i_rt, i_sea]
        data_delta_1=data_delta[:,:,i_rt,i_sea]

        for ii in range(1, 11):  # 1-10
            data_rid = data_delta_sign_1[RID == ii]
            data_less = data_rid[data_rid < 0]  # harder to recover
            data_over = data_rid[data_rid > 0]  # easier to recover
            data_has_change = data_delta_1[RID==ii]
            data_has_change_no_nan = np.delete(data_has_change, np.where(np.isnan(data_has_change)))
            data_rid_all=np.concatenate([data_rid_all,data_rid])
            share_rid[ii - 1, 0, i_rt, i_sea] = len(data_less) / len(data_has_change_no_nan)
            share_rid[ii - 1, 1, i_rt, i_sea] = len(data_over) / len(data_has_change_no_nan)
            median_rid[ii - 1, 0, i_rt, i_sea] = np.median(data_less)
            median_rid[ii - 1, 1, i_rt, i_sea] = np.median(data_over)

    data_delta_no_nan=np.delete(data_delta_1.reshape(-1,1), np.where(np.isnan(data_delta_1.reshape(-1,1))))
    data_less = data_rid_all[data_rid_all < 0]
    data_over = data_rid_all[data_rid_all > 0]
    share_global[i_rt, 0] = len(data_less)/len(data_delta_no_nan)/4
    share_global[i_rt, 1] = len(data_over) / len(data_delta_no_nan)/4
    median_global[i_rt,0] = np.median(data_less)
    median_global[i_rt,1] = np.median(data_over)
median_rid_rt_8_14_harder=median_rid[:,0,1,:]
median_rid_rt_8_14_easier=median_rid[:,1,1,:]

import proplot as pplt
RID_neworder=[7,8,3,2,4,10,9,6,1,5]  # 1-10['1SAF', '2NAS', '3EUR', '4SNA', '5OCE', '6SSA', '7NNA', '8CAF', '9AMZ', '10SAS']
imate_region = ['SAF', 'NAS', 'EUR', 'SNA', 'OCE', 'SSA', 'NNA', 'CAF', 'AMZ', 'SAS']
sea = ['MAM', 'JJA', 'SON', 'DJF']
RT = ['1-7', '8-14', '15-21', '22-28']
color_red = '#fc8982'
color_blue = '#6880ca'
font_color = '#525252'
hfont = {'fontname': 'Myriad Pro'}
facecolor = '#eaeaf2'
index=sea
array = [  # the "picture" (1 == subplot A, 2 == subplot B, etc.)
    [1, 2, 3, 4, 5, 6, 7, 8],
    [9, 10, 21, 21, 21, 21, 11, 12],
    [13, 14, 15, 16, 17, 18, 19, 20],
]

for i_rt in range(4):
    fig, axs = pplt.subplots(
        array, refwidth=1.1, refaspect=0.6, span=False,
        bottom='5em', right='5em', top='7em', # margin spacing overrides
        wspace=(0, 2, 0, 2, 0, 2, 0), hspace=(4, 4))
    plt.delaxes(axs[20])
    for ii in range(10):
        column0 = share_rid[RID_neworder[ii]-1, 1, i_rt,:] * 100
        column1 = share_rid[RID_neworder[ii]-1, 0, i_rt,:] * 100
        width0 = median_rid[RID_neworder[ii]-1, 1, i_rt,:] * 3
        width1 = median_rid[RID_neworder[ii]-1, 0, i_rt,:] * 3

        ax0 = axs[ii * 2]
        ax1 = axs[ii * 2 + 1]
        ax0.barh(index, column0, width=width0, align='center', color='#457B9D',  # color=(122/255, 162/255, 170/255),
                 zorder=10)  # barh() would draw horizontal bar plots.

        ax1.barh(index, column1, width=width1, align='center', color='#C47C70', zorder=10)        
        ax0.set_xticks([0, 25, 50])
        ax0.set_xlim(0, 65)
        ax1.set_xticks([0, 25, 50])
        ax1.set_xlim(0, 65)
        ax0.set_ylim(-0.5, 3.5)
        ax0.axvline(share_global[i_rt, 1] * 100, color='grey', linestyle='--')
        ax1.axvline(share_global[i_rt, 0] * 100, color='grey', linestyle='--')
        # If you have positive numbers and want to invert the x-axis of the left plot
        ax0.invert_xaxis()
        ax0.invert_yaxis()
        ax1.axes.yaxis.set_visible(False)

        if ii == 0 or ii == 4 or ii == 6:

            for label in (ax0.get_yticklabels()):
                label.set(fontsize=14, color=font_color, **hfont)
        else:
            ax0.axes.yaxis.set_visible(False)

        for label in (ax0.get_xticklabels()):
            label.set(fontsize=15, color=font_color, **hfont)
        for label in (ax1.get_xticklabels()):
            label.set(fontsize=15, color=font_color, **hfont)

    # title = sea[i_sea]
    title = RT[i_rt]
    plt.suptitle(title)
    plt.savefig('J:\\output\\Results_2\\' + title + '.jpg') #, dpi=600)
    plt.close()

#%%
def classify_region(RID, coef, neworder):
    coef_25_50_75 = np.empty((3, 11))
    coef_25_50_75[:] = np.nan
    coef_all = coef[~np.isnan(coef)]
    coe_sort_all = np.sort(coef_all)
    coef_75 = coe_sort_all[round(len(coef_all) * 0.75)]
    coef_25 = coe_sort_all[round(len(coef_all) * 0.25)]
    coef_25_50_75[:, 10] = [coef_25, np.median(coef_all), coef_75]

    for i_rid in range(10):
        coef_region = coef[RID == neworder[i_rid]]
        coef_region2 = coef_region[~np.isnan(coef_region)]
        len_significant = len(coef_region2)
        coe_sort = np.sort(coef_region2)
        coef_75 = coe_sort[round(len_significant * 0.75)]
        coef_25 = coe_sort[round(len_significant * 0.25)]
        coef_25_50_75[:, i_rid] = [coef_25, np.median(coef_region2), coef_75]

    return coef_25_50_75


cbformat = mticker.ScalarFormatter()
cbformat.set_powerlimits((-2, 2))  # Set size thresholds for scientific notation.
lat = np.arange(89.875, -90, -0.25)
lon = np.arange(-179.875, 180, 0.25)
RT = ['1-7', '8-14', '15-21', '22-28']
sea = ['MAM', 'JJA', 'SON', 'DJF']

# %%
# typical scenario map
ncolor = ['#C47530', '#8c510a', '#bf812d', '#ffbf80', '#776483', '#222B5F', '#776483', '#3e99a7', '#327355',
          '#003c30']
neworder = [7, 3, 2, 4, 8, 10, 9, 6, 1, 5]
newcolor = []
climate_region = ['SAF', 'NAS', 'EUR', 'SNA', 'OCE', 'SSA', 'NNA', 'CAF', 'AMZ',
                  'SAS']
region_new = []
period=['1951-1983','1984-2016']
for ii in range(10):
    newcolor.append(ncolor[neworder[ii] - 1])
    region_new.append(climate_region[neworder[ii] - 1])

region = nc.Dataset(r'E:\phd\data\climate region\region.id.nc', 'r')
RID = region.variables['rid'][:].data
RID = np.flipud(RID)
sea = ['MAM', 'JJA', 'SON', 'DJF']
cmap = plt.get_cmap('YlOrRd')
new_cmap = truncate_colormap(cmap, minval=0.05, maxval=1.0, n=100)
bounds = np.linspace(0, 10, 11)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)


for i_rt in range(4):
    data_his = np.load(
        'J:\\output\\regression\\his_grid_global_coef_rt(' + RT[i_rt] + ').npy')
    p_his = np.load(
        'J:\\output\\regression\\his_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy')
    data_his[p_his > 0.05] = np.nan

    data_pres = np.load(
        'J:\\output\\regression\\pres_grid_global_coef_rt(' + RT[i_rt] + ').npy')
    p_pres = np.load(
        'J:\\output\\regression\\pres_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy')
    data_pres[p_pres > 0.05] = np.nan

    grid_his = np.load('J:\\output\\regression\\his_grid_need_cal.npy')
    grid_pres = np.load('J:\\output\\regression\\pres_grid_need_cal.npy')
    lat_his = grid_his[0, :].astype(int)
    lon_his = grid_his[1, :].astype(int)

    coef_map_his = np.empty([600, 1440, 4, 6])
    coef_map_his[:] = np.nan
    num_grid_his = np.zeros((600, 1440))
    num_grid_his[lat_his, lon_his] = 1
    coef_map_his[lat_his, lon_his, :, :] = data_his * 100

    lat_pres = grid_pres[0, :].astype(int)
    lon_pres = grid_pres[1, :].astype(int)
    coef_map_pres = np.empty([600, 1440, 4, 6])
    coef_map_pres[:] = np.nan
    num_grid_pres = np.zeros((600, 1440))
    num_grid_pres[lat_pres, lon_pres] = 1
    coef_map_pres[lat_pres, lon_pres, :, :] = data_pres * 100

    fig, axs = plt.subplots(nrows=2, ncols=2,
                            subplot_kw={'projection': ccrs.Robinson()},  # Robinson()},  # PlateCarree()},
                            figsize=(20, 13))
    num = -1
    axs=axs.flatten()
    i_sea = 1
    for coef_map in [coef_map_his, coef_map_pres]:
        for i_interval in [1, 4]:
            num+=1
            ax1 = axs[num]
            data1 = coef_map[:, :, i_sea, i_interval]
            np_condition = np.where((np.isnan(data1)) & (RID > 0) & (RID < 11))
            data_nan = np.empty((600, 1440))
            data_nan[:] = np.nan
            data_nan[np_condition] = 1
            data_nan_new = add_south(data_nan)
            coef_25_50_75_oridinary = classify_region(RID, data1, neworder)
            order_descend=np.flipud(np.argsort(coef_25_50_75_oridinary[1, :10]))
            coef_25_50_75=np.zeros((3,11))
            region_new_order=[]
            newcolor_order=[]

            for id in range(10):
                coef_25_50_75[:,id]=coef_25_50_75_oridinary[:,order_descend[id]]
                coef_25_50_75[:,10]=coef_25_50_75_oridinary[:,10]
                region_new_order.append(region_new[order_descend[id]])
                newcolor_order.append(newcolor[order_descend[id]])

            x = np.arange(10)
            earth=np.zeros((720,1440))
            # plot map
            data1_south = add_south(data1)
            cs1 = ax1.pcolormesh(lon, lat, data1_south, shading='auto', norm=norm, cmap=new_cmap,
                                 transform=ccrs.PlateCarree())

            cs2 = ax1.pcolormesh(lon, lat, data_nan_new, shading='auto', norm=norm_nan, cmap=cmap_nan,
                                 transform=ccrs.PlateCarree())
            # Draw the coastines for each subplot
            ax1.coastlines(lw=0.2, color='#696969')
            ax1.add_feature(cfeat.LAND.with_scale('110m'), facecolor='#CDCDCD', alpha=0.8, edgecolor='#CDCDCD',
                            linewidth=0.01)  # ,edgecolor='')#color='#D3D3D3'

            gl = ax1.gridlines(draw_labels=False, lw=0.00007, color='#C0C0C0', linestyle='--', alpha=0.5)
            ax1.patch.set_facecolor('white')
            ax1.patch.set_alpha(1)
            if num==3:
                cbar_ax = fig.add_axes([0.27, 0.12, 0.5, 0.025])  # [0.153, 0.12, 0.45, 0.02])
                cbar = fig.colorbar(cs1, cax=cbar_ax, orientation='horizontal', extend='max')
                for l in cbar.ax.xaxis.get_ticklabels():
                    l.set_family('Myriad Pro')
                cbar.ax.tick_params(labelsize=18)

    fig.subplots_adjust(bottom=0.18, left=0.1,  # right=0.95, top=0.95,
                        wspace=0.06, hspace=0.02)  # wspace: width

    # Add a colorbar axis at the bottom of the graph
    title = 'RT=' + RT[i_rt] + ' ' + sea[i_sea]
    # Add a big title at the top
    plt.suptitle(title)
    plt.savefig('J:\\output\\Results_3\\map_' + title + '.jpg')  # , dpi=600)
    plt.close()


# %%
# box plot
x=np.arange(10)
ncolor = ['#C47530', '#8c510a', '#bf812d', '#ffbf80', '#776483', '#222B5F', '#776483', '#3e99a7', '#327355',
          '#003c30']
neworder = [7, 3, 2, 4, 8, 10, 9, 6, 1, 5]
newcolor = []
climate_region = ['SAF', 'NAS', 'EUR', 'SNA', 'OCE', 'SSA', 'NNA', 'CAF', 'AMZ',
                  'SAS']
region_new = []
period = ['1951-1983', '1984-2016']
for ii in range(10):
    newcolor.append(ncolor[neworder[ii] - 1])
    region_new.append(climate_region[neworder[ii] - 1])

region = nc.Dataset(r'E:\phd\data\climate region\region.id.nc', 'r')
RID = region.variables['rid'][:].data
RID = np.flipud(RID)
sea = ['MAM', 'JJA', 'SON', 'DJF']
cmap = plt.get_cmap('YlOrRd')
new_cmap = truncate_colormap(cmap, minval=0.05, maxval=1.0, n=100)
bounds = np.linspace(0, 10, 11)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap_nan = ListedColormap(['#FFFFFF'])
norm_nan = mpl.colors.Normalize(vmin=0, vmax=2)
norm_sea = mpl.colors.Normalize(vmin=-1, vmax=1)
for i_rt in range(4):
    data_his = np.load(
        'J:\\output\\prob_season_0811\\global\\regression\\his_grid_global_coef_rt(' + RT[i_rt] + ').npy')
    p_his = np.load(
        'J:\\output\\prob_season_0811\\global\\regression\\his_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy')
    data_his[p_his > 0.05] = np.nan

    data_pres = np.load(
        'J:\\output\\prob_season_0811\\global\\regression\\pres_grid_global_coef_rt(' + RT[i_rt] + ').npy')
    p_pres = np.load(
        'J:\\output\\prob_season_0811\\global\\regression\\pres_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy')
    data_pres[p_pres > 0.05] = np.nan

    grid_his = np.load('J:\\output\\prob_season_0811\\global\\regression\\his_grid_need_cal.npy')
    grid_pres = np.load('J:\\output\\prob_season_0811\\global\\regression\\pres_grid_need_cal.npy')
    lat_his = grid_his[0, :].astype(int)
    lon_his = grid_his[1, :].astype(int)

    coef_map_his = np.empty([600, 1440, 4, 6])
    coef_map_his[:] = np.nan
    num_grid_his = np.zeros((600, 1440))
    num_grid_his[lat_his, lon_his] = 1
    coef_map_his[lat_his, lon_his, :, :] = data_his * 100

    lat_pres = grid_pres[0, :].astype(int)
    lon_pres = grid_pres[1, :].astype(int)
    coef_map_pres = np.empty([600, 1440, 4, 6])
    coef_map_pres[:] = np.nan
    num_grid_pres = np.zeros((600, 1440))
    num_grid_pres[lat_pres, lon_pres] = 1
    coef_map_pres[lat_pres, lon_pres, :, :] = data_pres * 100

    fig, axs = plt.subplots(nrows=2, ncols=2,figsize=(29, 9))
    num = -1
    axs=axs.flatten()
    # axs is a 2 dimensional array of `GeoAxes`.  We will flatten it into a 1-D array
    i_sea = 1
    # color = sns.color_palette("BuPu", 11)
    for coef_map in [coef_map_his, coef_map_pres]:
        for i_interval in [1, 4]:
            num+=1
            ax1 = axs[num]

            data1 = coef_map[:, :, i_sea, i_interval]

            np_condition = np.where((np.isnan(data1)) & (RID > 0) & (RID < 11))
            data_nan = np.empty((600, 1440))
            data_nan[:] = np.nan
            data_nan[np_condition] = 1
            data_nan_new = add_south(data_nan)
            # data2 = coef_map_pres[:, :, i_sea, i_interval]
            coef_25_50_75_oridinary = classify_region(RID, data1, neworder)
            order_descend=np.flipud(np.argsort(coef_25_50_75_oridinary[1, :10]))
            coef_25_50_75=np.zeros((3,11))
            region_new_order=[]
            newcolor_order=[]

            for id in range(10):
                coef_25_50_75[:,id]=coef_25_50_75_oridinary[:,order_descend[id]]
                coef_25_50_75[:,10]=coef_25_50_75_oridinary[:,10]
                region_new_order.append(region_new[order_descend[id]])
                newcolor_order.append(newcolor[order_descend[id]])

            yerror = abs(coef_25_50_75[::2, :10] - np.tile(coef_25_50_75[1, :10], [2, 1]))
            for pos, y, err, err1, colors_bar in zip(x, coef_25_50_75[1, :10], yerror[0, :], yerror[1, :], newcolor_order):
                ax1.errorbar(pos, y, np.array([err, err1]).reshape(2, 1), fmt='o',
                                lw=2, capsize=5, capthick=3, markersize=17,elinewidth=2,
                                color=colors_bar)

            for patch in ax1.artists:
                fc = patch.get_facecolor()
                patch.set_facecolor(mpl.colors.to_rgba(fc, 0.85))
            ax1.axhline(coef_25_50_75[1, 10], color='#A52502', linestyle='--', linewidth=1, zorder=0)
            ax1.set_xticks(x)
            ax1.set_xticklabels(region_new_order, rotation=0, fontsize=32, fontfamily='Myriad Pro')

            if i_interval == 1:
                ax1.set(ylim=(0, 7))
                ax1.set_yticks([0, 7])
                ax1.set_yticklabels([0, 7], fontsize=27, fontfamily='Myriad Pro')
            else:
                ax1.set(ylim=(0, 3))
                ax1.set_yticks([0, 3])
                ax1.set_yticklabels([0, 3], fontsize=27, fontfamily='Myriad Pro')
            # axins2.set(ylim=(0, 0.25))

            ax1.set(xlabel=None)
            ax1.set(ylabel=None)
            ax1.patch.set_alpha(0.0)
            [ax1.spines[loc_axis].set_visible(False) for loc_axis in ['top', 'right']]
            ax1.tick_params(bottom=False, top=False, left=True, right=False)
            # axins1.tick_params(top=False, bottom=False, left=False, right=False)
            ax1.grid(False)
            ax1.patch.set_alpha(0)

            plt.rcParams['font.family'] = 'Myriad Pro'
    title = 'RT=' + RT[i_rt] + ' ' + sea[i_sea]
    # Add a big title at the top
    plt.suptitle(title)
    fig.subplots_adjust(bottom=0.018, left=0.01,
                        wspace=0.1, hspace=0.3)
    plt.savefig('J:\\output\\Results_3\\box_' + title + '.jpg')  # , dpi=600)
    plt.close()
