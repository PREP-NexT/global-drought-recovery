import numpy as np
import netCDF4 as nc
import statsmodels.api as sm


#%%
def reverse_mapdata(inputdata):
    """Convert data from the coordinate reference of "-59.875°N,0.125°E" to "89.875°N,179.875°W" for the first row and first column"""
    prob_fli = np.flipud(inputdata[:, :, :, :])
    data_map = np.concatenate((prob_fli[:, 720:1440, :, :], prob_fli[:, 0:720, :, :]), axis=1)
    return data_map

def regression_step(rp):
    rp[rp<0.001]=0.001
    ln_rp=np.log(rp)
    x=np.linspace(0,27.5,12).reshape(-1,1)
    coef=np.empty([6])
    coef[:]=np.nan
    p_value = np.empty([6])
    p_value[:] = np.nan
    for i in range(6):
        y_interval=ln_rp[i*12:i*12+12].reshape(-1,1)
        X2=sm.add_constant(x)
        est = sm.OLS(y_interval, X2)
        est2=est.fit()
        coef[i]=est2.params[1]
        p_value[i]=est2.pvalues[1]
    return coef, p_value


path1 = 'J:\\output\\global\\'
phase='his'
RT=['1-7','8-14','15-21','22-28']

for i_rt in range(4):
    filename1 = path1 + 'global_sce_' + phase + '_rt(' + RT[i_rt] + ').npy'  # ratio_all=np.linspace(0.025, 2, 80)#1-1:0.025-1.0（0-39）,2_2:1.025~2.0（40~79）
    data_input_global = np.load(filename1)
    data_input_map = reverse_mapdata(data_input_global)
    map_one_sce = data_input_map[:, :, 0, 0]
    map_one_sce[(map_one_sce == -111) | (map_one_sce == -999)] = np.nan
    grid_need_cal=np.array(np.where(~np.isnan(map_one_sce)))
    coef_all=np.empty([grid_need_cal.shape[1],4,6])  # grid_id, season, interval
    coef_all[:]=np.nan
    pvalue_all = np.empty([grid_need_cal.shape[1], 4, 6])  # grid_id, season, interval
    pvalue_all[:] = np.nan
    print(i_rt)
    for grid_i in range(grid_need_cal.shape[1]):
        for i_sea in range(4):
            data_one_scenario = data_input_map[grid_need_cal[0,grid_i], grid_need_cal[1,grid_i], i_sea, 4:76]
            coef_all[grid_i,i_sea,:],pvalue_all[grid_i,i_sea,:]=regression_step(data_one_scenario)
    path = 'J:\\output\\regression\\' + phase
    np.save(path+'_grid_global_coef_rt(' + RT[i_rt] + ').npy',arr=coef_all)
    np.save(path + '_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy', arr=pvalue_all)
np.save(path+'_grid_need_cal.npy', arr=grid_need_cal)

#%%
# Calculate the regression coefficients are for each region of the grid。
region = nc.Dataset(r'E:\phd\data\climate region\region.id.nc', 'r')
RID = region.variables['rid'][:].data
RID = np.flipud(RID)
for phase in ['his','pres']:

    data_coef=np.load('J:\\output\\prob_season_0811\\global\\regression\\'+phase+'_grid_global_coef_rt('+RT[i_rt]+').npy')
    p_his = np.load(
        'J:\\output\\prob_season_0811\\global\\regression\\'+phase+'_grid_global_coef_pvalue_rt(' + RT[i_rt] + ').npy')
    data_coef[p_his > 0.05] = np.nan

    grid=np.load('J:\\output\\regression\\'+phase+'_grid_need_cal.npy')

    lat = grid[0,:].astype(int)
    lon = grid[1,:].astype(int)
    coef_map=np.empty([600,1440,4,6])
    coef_map[:]=np.nan
    coef_map[lat,lon,:,:]=data_coef
    coef_region_all=np.empty((10,50000,4,6))
    coef_region_all[:]=np.nan
    for i_rid in range(1,11):
        coef_region=coef_map[RID==i_rid, :, :]
        coef_region_all[i_rid-1,:coef_region.shape[0],:,:]=coef_region
    np.save('J:\\output\\regression\\region_coef'+phase+RT[i_rt]+'.npy', arr=coef_region_all)
