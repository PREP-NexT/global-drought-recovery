# Load the data and extract characteristics of drought events.
import numpy as np
import time
import os


# %%
def judge_phase(SMPct_grid, Prec_grid, phase):
    """judge 1951-1983 or 1984-2016"""
    # input
    # SMPct_grid and Prec_grid include data series from 1950 to 2016, length:24472
    # phase'0：history' or '1：present'
    # output
    # history or present period SMPct and Prec data
    if phase == 0:
        SMPct_i = SMPct_grid[365:12418]  # 365-12417:1951.1.1~1983.12.31
        Prec_i = Prec_grid[365:12418]
    else:
        SMPct_i = SMPct_grid[12418:24472]  # 12418-24471:1951.1.1~1983.12.31
        Prec_i = Prec_grid[12418:24472]
    return SMPct_i, Prec_i


# %%
def judge_drought_SMPct(SMPct, SM_threshold, min_duration):
    """Record the date when the drought occurred and ended"""
    # input
    # SMPct: t * 1 grid SMPct
    # SM_threshold: If continuously below this value for 14 days or more, it is considered a drought event
    # output
    # d_event: a matrix to record drought events.
    #          [Drought start date: counting from 0, the first day below the threshold, Day 0 is 1951.1.1 or 1984.1.1]
    #          [Drought end date: counting from 0, the last day below the threshold date]
    #          [Duration of the drought, including both start and end dates]

    t = SMPct.shape[0]
    nn = 0
    num_d = 0  # number of drought events
    one_event = np.zeros(10000)
    d_event = np.zeros((1000, 3))
    for ii in range(t):
        if SMPct[ii] < SM_threshold:
            nn += 1
            one_event[nn - 1] = SMPct[ii]
        elif SMPct[ii] >= SM_threshold and nn >= min_duration:
            num_d += 1
            d_event[num_d - 1, :] = np.array([ii - nn, ii - 1, nn])
            # The start date is the first day below 20% SMPct,
            # and the end date is the last day below 20% SMPct.
            nn = 0
            one_event[:] = np.zeros(10000)
        else:
            nn = 0
            one_event[:] = np.zeros(10000)  # If the end date of the last event is the last day of a long sequence,
            # drought recovery cannot be determined, hence no output.
    d_event = d_event[:num_d, :]
    return d_event


# %%
def drought_character(lat, lon, Prep, SMPct, d_event, delay_time, SM_threshold):
    """Calculate drought characteristics"""
    # input：
    # Prep: a column matrix of size t*1 used to store grid daily precipitation sequences.
    # SMPct: a column matrix of size t*1 used to store grid daily smpct sequences.
    # d_event: three columns correspond to: [Drought start date, Drought end date, Duration of the drought].
    # output：
    # grid_charct:0: Grid latitude, 1: Grid longitude, 2: Index of the event in the grid,
    # 3: Drought start date (counted from 0, starting from 1951.1.1 for 'his' and 1984.1.1 for 'pres'),
    # 4: Drought end date (counted in the same way as 3), 5: Duration of the drought (inclusive of start and end dates),
    # 6: Day of drought peak within the drought period (counted from Day 0),
    # 7: Difference between drought peak and the threshold value,
    # 8: Duration of drought recovery (starting from the first day after 'DM' until the first day when SMPct>=20),
    # 9: Precipitation during the drought recovery period (corresponding to the same period as 8),
    # 10: Precipitation on the day following the entire drought period, indicating cumulative precipitation from
    # the day after the drought start date to the day after the drought end date,
    # 11: Soil Moisture on the day after the drought end.
    # grid_charct, columns=['lat', 'lon', ii,'start_date_num', 'end_date_num',
    #                       'duration', 'min_sm_date_num', 'max_deficit',
    #                       'recovery_duration', 'recovery_prcp',
    #                       'drought_prcp', 'SM_after_drought'])

    n = d_event.shape[0]  # n represents the total count of drought events in this grid.
    grid_charct1 = np.zeros((n, 12), dtype="float32")
    for ii in range(n):
        date_start = d_event[ii, 0]
        date_end = d_event[ii, 1]
        SM = SMPct[date_start:date_end + 1]
        min_SM = np.min(SM)
        min_date = np.where(SM == min_SM)[0]
        P = Prep[
            date_start + delay_time:date_end + delay_time + 1]  # Considering the delayed effect between precipitation and soil moisture, a lag of 1 day.
        sum_p1 = np.sum(P[min_date[
                              -1]:])  # cumulative precipitation from the day after the drought peak utill the drought ends.
        sum_p2 = np.sum(
            P)  # Represents cumulative precipitation during the period from the day after the drought start date to the day after the drought end date.
        grid_charct1[ii, :] = np.array([lat, lon, ii, d_event[ii, 0], d_event[ii, 1], d_event[ii, 2],
                                        min_date[-1], SM_threshold - min_SM, d_event[ii, 2] - min_date[-1],
                                        sum_p1, sum_p2, SMPct[date_end + 1]])
    return grid_charct1, n


# %%
def region_drought_charc(data_SMPct, data_Prcp, SM_threshold, min_duration, lat_area, lon_area, event_num_his,
                         event_num_pres):
    """all events in the region"""
    num_empty = 1
    grid_charct_all_his = np.empty((1, 12))
    grid_charct_all_his[:] = None
    grid_charct_all_pres = np.empty((1, 12))
    grid_charct_all_pres[:] = None
    for lat_num in range(60):
        for lon_num in range(72):
            SMPct_grid = data_SMPct[:, lat_num, lon_num]
            if np.isnan(SMPct_grid).all():
                num_empty += 1
            else:
                Prec_grid = data_Prcp[:, lat_num, lon_num]
                phase = 0  # his
                SMPct_i, Prec_i = judge_phase(SMPct_grid, Prec_grid, phase)
                d_event = judge_drought_SMPct(SMPct_i, SM_threshold, min_duration)
                lat = (lat_area * 60 + lat_num) * 0.25 - 59.875
                lon = (lon_area * 72 + lon_num) * 0.25 + 0.125
                d_event = np.array(d_event, dtype=int)
                grid_charct, n = drought_character(lat, lon, Prec_i, SMPct_i, d_event, delay_time, SM_threshold)
                event_num_his[lat_area * 60 + lat_num, lon_area * 72 + lon_num] = n
                grid_charct_all_his = np.concatenate((grid_charct_all_his, grid_charct), axis=0)

                phase = 1  # pres
                SMPct_i, Prec_i = judge_phase(SMPct_grid, Prec_grid, phase)
                d_event = judge_drought_SMPct(SMPct_i, SM_threshold, min_duration)
                lat = (lat_area * 60 + lat_num) * 0.25 - 59.875
                lon = (lon_area * 72 + lon_num) * 0.25 + 0.125
                d_event = np.array(d_event, dtype=int)
                grid_charct, n = drought_character(lat, lon, Prec_i, SMPct_i, d_event, delay_time, SM_threshold)
                event_num_pres[lat_area * 60 + lat_num, lon_area * 72 + lon_num] = n
                grid_charct_all_pres = np.concatenate((grid_charct_all_pres, grid_charct), axis=0)
    return grid_charct_all_his, grid_charct_all_pres, event_num_his, event_num_pres


# %%
SM_threshold = 20
min_duration = 14
delay_time = 1
start_time = time.time()
event_num_his = np.empty((600, 1440))
event_num_his[:] = None
event_num_pres = np.empty((600, 1440))
event_num_pres[:] = None

path1 = 'J:\\output\\his'
path2 = 'J:\\output\\pres'

num_empty = 0
for lat_area in range(10):
    for lon_area in range(20):
        data_SMPct = np.load(r"J:\GDFC\SMPct\numpy2\SMPct" + str(lat_area) + "-" + str(lon_area) + ".npy")
        if np.isnan(data_SMPct).all():
            print(lat_area, lon_area, 'Nan')
        else:
            data_Prcp = np.load(r"J:\\GDFC\\Prcp\\numpy2\\Prcp" + str(lat_area) + "-" + str(lon_area) + ".npy")
            print(lat_area, lon_area)
            print((time.time() - start_time) / 60)
            grid_charct_all_his, grid_charct_all_pres, event_num_his, event_num_pres = region_drought_charc(data_SMPct,
                                                                                                            data_Prcp,
                                                                                                            SM_threshold,
                                                                                                            min_duration,
                                                                                                            lat_area,
                                                                                                            lon_area,
                                                                                                            event_num_his,
                                                                                                            event_num_pres)

            filename1 = os.path.join(path1, 'his_charc' + str(lat_area) + "-" + str(lon_area))
            filename2 = os.path.join(path2, 'pres_charc' + str(lat_area) + "-" + str(lon_area))
            np.save(filename1, arr=grid_charct_all_his[1:, :])
            np.save(filename2, arr=grid_charct_all_pres[1:, :])
            grid_charct_all_his = np.empty((1, 12))
            grid_charct_all_his[:] = None
            grid_charct_all_pres = np.empty((1, 12))
            grid_charct_all_pres[:] = None

np.save('J:\\output\\his\\his_event.npy', arr=event_num_his)
np.save('J:\\output\\pres\\pres_event.npy', arr=event_num_pres)
