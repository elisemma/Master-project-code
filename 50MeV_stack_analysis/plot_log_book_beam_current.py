import numpy as np 
import matplotlib.pyplot as plt 
import datetime



times = np.array([datetime.datetime(2017, 2, 12, 19, 1, 0), datetime.datetime(2017, 2, 12, 19, 2, 0), datetime.datetime(2017, 2, 12, 19, 3, 0), 
                  datetime.datetime(2017, 2, 12, 19, 4, 0), datetime.datetime(2017, 2, 12, 19, 5, 0), datetime.datetime(2017, 2, 12, 19, 6, 0), 
                  datetime.datetime(2017, 2, 12, 19, 7, 0), datetime.datetime(2017, 2, 12, 19, 8, 0), datetime.datetime(2017, 2, 12, 19, 9, 0), 
                  datetime.datetime(2017, 2, 12, 19, 10, 0), datetime.datetime(2017, 2, 12, 19, 11, 0), datetime.datetime(2017, 2, 12, 19, 12, 0), 
                  datetime.datetime(2017, 2, 12, 19, 13, 0), datetime.datetime(2017, 2, 12, 19, 14, 0), datetime.datetime(2017, 2, 12, 19, 15, 0), 
                  datetime.datetime(2017, 2, 12, 19, 16, 0), datetime.datetime(2017, 2, 12, 19, 17, 0), datetime.datetime(2017, 2, 12, 19, 18, 0), 
                  datetime.datetime(2017, 2, 12, 19, 19, 0), datetime.datetime(2017, 2, 12, 19, 20, 0), datetime.datetime(2017, 2, 12, 19, 21, 0)])

integrated_beam_cur = np.array([0, 31, 61, 91, 121, 150, 179, 209, 240, 268, 297, 326, 356, 385, 414, 443, 472, 501, 531, 560, 588])


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
plt.ylabel('Average beam current since start of irradiation (nA)')
plt.ylim(0, 130)
plt.xlim(0)
# plt.savefig('./Figures/beam_current.pdf')
plt.show()
