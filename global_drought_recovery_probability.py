import warnings
import numpy as np
import scipy.stats as st
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2
pandas2ri.activate()
importr('VineCopula')
importr('copula')
importr('CDVineCopulaConditional')


# %%
def best_fit_distribution(data, jj):
    """Model data by finding best fit distribution to data"""
    if jj == 0:
        Distribution = [st.beta]
    elif jj == 1:
        Distribution = [st.genpareto]
    else:
        Distribution = [st.genpareto]
    # Best holders based on KS
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_KS = np.inf
    for distribution in Distribution:
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate the KS statistic
                KS = st.kstest(data, distribution.name, params)
                # identify if this distribution is better
                if best_KS > KS[0] > 0:
                    best_distribution = distribution
                    best_params = params
                    best_KS = KS[0]

        except Exception:
            pass
    return (best_distribution.name, best_params)


# %%
def make_cdf(dist, params, data):
    """Generate distributions's Probability Distribution Function """
    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]
    y = dist.cdf(data, loc=loc, scale=scale, *arg)
    return y


# %%
def get_percentile(ux, dist, params, logtransformed=False):
    """ Get the percentile of variable x based on its best marginal distribution"""
    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get the inversed magnitude
    x = dist.ppf(ux, *arg, loc=loc, scale=scale) if arg else dist.ppf(ux, loc=loc, scale=scale)

    if logtransformed == True:
        x = 1 - 10 ** (-x)
    return x


# %%
def cal_prob(prec_500, prec_seas_grid, rt_s, rt_e):

    prec_pro = np.zeros((4, 80))
    prec_pro[:] = np.nan
    ratio_all = np.linspace(0.025, 2, 80)
    prec_sample=prec_500[:,1]
    for season in range(4):
        date_s=season * 28
        date_s = np.int32(date_s)
        date_e = season * 28+28
        date_e = np.int32(date_e)
        prec_se = prec_seas_grid[date_s:date_e]
        for i_ratio in range(80):
            num_recovery = 0
            for rt in range(rt_s, rt_e):  # rt:1~7；8~14；15~21；22~28
                prec_sce = prec_se * ratio_all[i_ratio]
                np_condition_rt = np.where((prec_500[:, 0] > rt) & (prec_500[:, 0] <= rt + 1))
                prec_sim_rt = prec_sample[np_condition_rt]
                recovery_rt = prec_sim_rt[np.where(prec_sim_rt <= prec_sce[rt - 1])]  # Number of grids recovering from drought at this RT.
                num_recovery = num_recovery + recovery_rt.shape[0]
            prec_pro[season, i_ratio] = num_recovery / 501
    return prec_pro


# %%
def simvinecopula(drought_grid_simvine, max_serious_grid, prec_seas_grid, num_file, num_grid):
    """
    # fit copula
    :param drought_grid_simvine: [lat,lon,dm,rt,prec]
    :param max_serious_grid: [dm,rt,prec]
    :param prec_seas_grid
    :param num_file
    :param num_grid
    :return:
    """
    # output: # -111 indicates that the number of samples meeting the criteria is less than 50, while -999 indicates that the grid cannot undergo vine copula fitting.
    best_fit_name_DM, best_params_DM = best_fit_distribution(drought_grid_simvine[:, 2], 0)
    best_fit_name_RT, best_params_RT = best_fit_distribution(drought_grid_simvine[:, 3], 1)
    best_fit_name_PREC, best_params_PREC = best_fit_distribution(
        drought_grid_simvine[:, 4], 2)

    # if np.min(np.array([best_KS_DM, best_KS_RT, best_KS_PREC])) > 0.05:

    best_dist_DM = getattr(st, best_fit_name_DM)
    best_dist_RT = getattr(st, best_fit_name_RT)
    best_dist_PREC = getattr(st, best_fit_name_PREC)
    cdf_DM = make_cdf(best_dist_DM, best_params_DM, drought_grid_simvine[:, 2])
    cdf_RT = make_cdf(best_dist_RT, best_params_RT, drought_grid_simvine[:, 3])
    cdf_PREC = make_cdf(best_dist_PREC, best_params_PREC, drought_grid_simvine[:, 4])
    DM_max = max_serious_grid[0]
    cdf_DM_serious = make_cdf(best_dist_DM, best_params_DM, DM_max)  # Return the cumulative distribution function (CDF) of the most severe historical event.
    condition_DM = np.ones((400000, 1)) * cdf_DM_serious
    np.random.seed(num_file * 20000 + num_grid)
    RT_ALL = np.r_[np.random.uniform(1, 8, [100000, 1]),
                   np.random.uniform(8, 15, [100000, 1]),
                   np.random.uniform(15, 22, [100000, 1]),
                   np.random.uniform(22, 29, [100000, 1])]
    condition_RT = make_cdf(best_dist_RT, best_params_RT, RT_ALL)
    U2 = np.c_[np.ones((400000, 1)), condition_RT, condition_DM]
    try:
        U = np.array([cdf_PREC, cdf_RT, cdf_DM]).T
        r.assign('U', U)
        r.assign('U2', U2)
        r('RVM = CDVineCondFit(U,Nx=2,type="CVine", c(1:6), treecrit="BIC", selectioncrit="BIC", rotations=TRUE)')  # selectioncrit = 'AIC'(default)
        seed_num = num_file * 10000 + num_grid
        r.assign('seed_num', seed_num)
        r('set.seed(seed_num)')
        r('d=dim(RVM$Matrix)[1]')
        r('cond1 <- U2[,RVM$Matrix[(d+1)-1,(d+1)-1]]')
        r('cond2 <- U2[,RVM$Matrix[(d+1)-2,(d+1)-2]]')
        r('condition_C <- cbind(cond1,cond2)')
        Usim = r('usim=CDVineCondSim(RVM,condition_C)')
        # Sim_time = get_percentile(Usim[:, 1], best_dist_RT, best_params_RT, logtransformed=False)
        Sim_prec = get_percentile(Usim[:, 0], best_dist_PREC, best_params_PREC, logtransformed=False)
        datasim = np.c_[RT_ALL, Sim_prec]

        RT_sim = np.array([1, 8, 15, 22, 29])
        pro_grid_clima_sce = np.zeros((4, 4, 80, 200))  # season, ratio, and 200 repetitions of experiments, respectively
        pro_grid_clima_sce[:] = None

        for ii in range(4):
            rt_s = RT_sim[ii]
            rt_e = RT_sim[ii + 1]
            select_data1 = datasim[ii * 100000:(ii + 1) * 100000,:]

            for sample_num in range(200):
                prec_500 = select_data1[sample_num * 500:(sample_num + 1) * 500] #Sample size is 500, repeated 200 times.
                prec_pro = cal_prob(prec_500, prec_seas_grid, rt_s, rt_e)
                pro_grid_clima_sce[ii, :, :, sample_num] = prec_pro

        pro_grid_climatology = pro_grid_clima_sce[:, :, 39, :]  # the mean for conducting significance testing.
        pro_grid_sce = np.mean(pro_grid_clima_sce, axis=3)

    except rpy2.rinterface_lib.embedded.RRuntimeError:
        pro_grid_climatology = np.ones((4, 4, 200)) * (-999)
        pro_grid_sce = np.ones((4, 4, 80)) * (-999)

    return pro_grid_sce, pro_grid_climatology

# %%
def sim_prob(data_charc, max_serious_global, season_prec, num_file, phase):
    """
    calculate the probability of drought recovery
    :param data_charc:[lat,lon,dm,rt,prec]
    :param max_serious_global:lat*lon*[dm,rt,prec]
    :param season_prec:lat*lon*[3~5 months, 1-28 days moving average precipitation, 6~8 months, 1-28 days moving average precipitation, 9~11 months, ..., 12~1 months, ...]
    :param num_file
    :return:
    """
    
    lat_lon = data_charc[:, 0:2]
    lat_lon_uni, indices = np.unique(lat_lon, return_index=True, axis=0)
    indices = np.r_[indices, lat_lon.shape[0]]
    indices = np.sort(indices)
    pro_file_climatology = np.zeros((0, 4, 4, 200))
    pro_file_sce = np.zeros((0, 4, 4, 80))
    lat_lon_cal = np.zeros((0, 2))
    for ii in range(indices.shape[0] - 1):
        lat_num = (data_charc[indices[ii], 0] + 59.875) / 0.25
        lon_num = (data_charc[indices[ii], 1] - 0.125) / 0.25
        lat_num = lat_num.astype(int)
        lon_num = lon_num.astype(int)
        prec_seas_grid = season_prec[lat_num, lon_num, :]

        if (np.isnan(max_serious_global[lat_num, lon_num, 0])) or np.isnan(prec_seas_grid[0]):  # Find grids that meet the criteria in both historical and current periods.
            print(num_file, ii, indices.shape[0] - 1, 'Nan')
        else:
            print(num_file, ii, indices.shape[0] - 1)
            max_serious_grid = max_serious_global[lat_num, lon_num, :]
            drought_grid_simvine = data_charc[indices[ii]:indices[ii + 1], :]  # all events for a specific grid
            pro_grid_sce, pro_grid_climatology = simvinecopula(drought_grid_simvine, max_serious_grid, prec_seas_grid,
                                                               num_file, ii)
            pro_grid_climatology = pro_grid_climatology.reshape(1, 4, 4, 200)
            pro_grid_sce = pro_grid_sce.reshape(1, 4, 4, 80)
            pro_file_climatology = np.concatenate((pro_file_climatology, pro_grid_climatology), axis=0)
            pro_file_sce = np.concatenate((pro_file_sce, pro_grid_sce), axis=0)

            lat_lon_cal = np.r_[
                lat_lon_cal, np.array([data_charc[indices[ii], 0], data_charc[indices[ii], 1]]).reshape(1, 2)]

    return (pro_file_climatology, pro_file_sce, lat_lon_cal)


# %%
# main code
filename_input1 = 'global_prec_seas.npy'
season_prec = np.load(filename_input1)
filename_input2 = './data/main_charc_add_rand/max_serious_global.npy'
max_serious_global = np.load(filename_input2)
for num in range(420):
    for phase in ['his','pres']:
        # load data
        filename_input3 = './main_charc_add_rand/' + phase + '_reslice/' + phase + '_' + str(num) + '.npy'
        print(phase, num)
        data_charc = np.load(filename_input3)  # [lat,lon,dm,rt,prec,dm_date]
        # calculate
        pro_file_climatology, pro_file_sce, lat_lon_cal = sim_prob(data_charc, max_serious_global, season_prec, num,
                                                                   phase)
        # save data
        filename_output1 = 'J:/output/grid_' + phase + '/lat_lon_' + phase + '_' + str(
            num) + '.npy'
        filename_output2 = 'J:/output/clima_' + phase + '/prob_clima_' + phase + '_' + str(
            num) + '.npy'
        filename_output3 = 'J:/output/sce_' + phase + '/prob_sce_' + phase + '_' + str(
            num) + '.npy'

        np.save(filename_output1, arr=lat_lon_cal)
        np.save(filename_output2, arr=pro_file_climatology)
        np.save(filename_output3, arr=pro_file_sce)
