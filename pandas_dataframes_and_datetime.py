import pandas as pd
import datetime



path = '/Users/elisemma/Library/CloudStorage/OneDrive-Personal/Dokumenter/Master/PlanBCode/MyGeneratedFiles/Ti_foils/CJ010317_Ti01_18cm_30MeV/'
file = 'CJ010317_Ti01_18cm_30MeV_peak_data.csv'

#Importing the data as a pd data frame
df = pd.read_csv(path+file,
        header=0,
        usecols=['isotope', 'decays', 'unc_decays', 'start_time', 'live_time']) 

#Converting the dates to datetime-objects
df['start_time'] = pd.to_datetime(df['start_time'])

#Choosing which isotopes/ rows I want to save as a new data frame
df_48V = df[df['isotope'] == '48V']

# printing dataframes
print('df:__________________________________ \n', df)
print('df_48V:______________________________ \n', df_48V)


#Saving the start time as a variable by choosing row 0 and column 3 in df_48V
start_time = df_48V.iloc[0,3]

#Saving a date as a datetime-object: datetime(year, month, day, hour, minute, second)
time_EOB = datetime.datetime(2017, 2, 13, 12, 50, 3)


#Calculating the difference between the EOB time and start time
test_delta_t = start_time - time_EOB

#Printing the number of seconds between the EOB and start time
print(test_delta_t.total_seconds())

