import numpy as np

# find the most severe drought event in the past 66 years.
filename_his = 'global_his_main_charc.npy'
main_charc_his = np.load(filename_his)
main_charc_his = np.c_[main_charc_his, np.ones((main_charc_his.shape[0], 1))]
filename_pres = 'global_pres_main_charc.npy'
main_charc_pres = np.load(filename_pres)
main_charc_pres = np.c_[main_charc_pres, np.ones((main_charc_pres.shape[0], 1)) * 2]


# %%
def find_max_dm(main_charc_his, main_charc_pres, lat_lon_his, indices_his):
    max_serious = np.zeros((0, 7))
    for ii in range(indices_his.shape[0] - 1):
        print(ii)
        lat = lat_lon_his[indices_his[ii], 0]
        lon = lat_lon_his[indices_his[ii], 1]
        main_charc_his_grid = main_charc_his[indices_his[ii]:indices_his[ii + 1], :]
        main_charc_pres_grid = main_charc_pres[(main_charc_pres[:, 0] == lat) & (main_charc_pres[:, 1] == lon), :]
        if main_charc_pres_grid.shape[0] >= 20:
            main_charc_all_grid = np.concatenate((main_charc_his_grid, main_charc_pres_grid))
            dm_max = np.max(main_charc_all_grid[:, 2])
            max_serious = np.concatenate(
                (max_serious, main_charc_all_grid[main_charc_all_grid[:, 2] == dm_max, :].reshape(1, 7)))  #[lat,lon,dm,rt,prec,date_num,his or pres]    return max_serious
    return max_serious

# %%
lat_lon_his = main_charc_his[:, 0:2]
lat_lon_his_uni, indices_his = np.unique(lat_lon_his, return_index=True, axis=0)
indices_his = np.r_[indices_his, lat_lon_his.shape[0]]
indices_his = np.sort(indices_his)
max_serious = find_max_dm(main_charc_his, main_charc_pres, lat_lon_his, indices_his)
np.save('max_serious.npy',arr=max_serious)


# %%
max_serious = np.load('J:\\output\\main charc_add rand\\max_serious.npy')
max_serious_global = np.zeros((600, 1440, 3)) # The latitude and longitude corresponding to the upper-left corner are 59.875°S and 0.125°E respectively.
max_serious_global[:] = None
lat = (59.875 + max_serious[:, 0]) / 0.25
lon = (max_serious[:, 1] - 0.125) / 0.25
lat = lat.astype(int)
lon = lon.astype(int)
max_serious_global[lat, lon, 0] = max_serious[:, 2]
max_serious_global[lat, lon, 1] = max_serious[:, 3]
max_serious_global[lat, lon, 2] = max_serious[:, 4]
np.save('max_serious_global.npy',arr=max_serious_global)
