import numpy as np 
import matplotlib.pyplot as plt 
import datetime



times = np.array([datetime.datetime(2017, 2, 13, 14, 27, 0), datetime.datetime(2017, 2, 13, 14, 28, 0), datetime.datetime(2017, 2, 13, 14, 29, 0), datetime.datetime(2017, 2, 13, 14, 30, 0), datetime.datetime(2017, 2, 13, 14, 31, 0), datetime.datetime(2017, 2, 13, 14, 32, 0), datetime.datetime(2017, 2, 13, 14, 33, 0), datetime.datetime(2017, 2, 13, 14, 34, 0), datetime.datetime(2017, 2, 13, 14, 35, 0), datetime.datetime(2017, 2, 13, 14, 36, 0), datetime.datetime(2017, 2, 13, 14, 37, 0), datetime.datetime(2017, 2, 13, 14, 38, 0), datetime.datetime(2017, 2, 13, 14, 39, 0), datetime.datetime(2017, 2, 13, 14, 40, 0), datetime.datetime(2017, 2, 13, 14, 41, 0), datetime.datetime(2017, 2, 13, 14, 42, 0), datetime.datetime(2017, 2, 13, 14, 43, 0), datetime.datetime(2017, 2, 13, 14, 44, 0), datetime.datetime(2017, 2, 13, 14, 45, 0), datetime.datetime(2017, 2, 13, 14, 46, 0), datetime.datetime(2017, 2, 13, 14, 47, 0)])

integrated_beam_cur = np.array([0, 37, 75, 112, 150, 188, 226, 264, 301, 339, 377, 414, 451, 489, 527, 566, 604, 642, 680, 718, 756])


def average_beam_cur_func(times, integrated_beam_cur):
    average_beam_cur_list = []
    start_time = times[0]
    times_sec = []
    for time, cur in zip(times[1:-1], integrated_beam_cur[1:-1]):
        delta_time = time-start_time
        delta_time_sec = delta_time.total_seconds()
        average_beam_cur_list.append(cur/delta_time_sec)
        times_sec.append(delta_time_sec)
    return np.array(times_sec), np.array(average_beam_cur_list)*2e-07



times_sec, average_beam_cur = average_beam_cur_func(times, integrated_beam_cur)

plt.plot(times_sec, average_beam_cur*1e9, color = 'hotpink', marker = 'o')
plt.xlabel('Time since start of irradiation (s)')
plt.ylabel('Average beam current (nA)')
# plt.savefig('./Figures/beam_current.pdf')
plt.show()
