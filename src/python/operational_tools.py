
import os, datetime as dt , numpy as np

# define epoch time

deltat = 600 

t0 = dt.datetime.utcfromtimestamp(0)

current_time = dt.datetime.utcnow()

path='./'


d1 = deltat*( np.round( ( current_time - t0 ).total_seconds() / deltat ) )

sec_max = ( dt.timedelta( seconds = d1 ) ).total_seconds()
sec_min = sec_max - deltat 

date1= t0 + dt.timedelta(seconds=sec_min)
date2= t0 + dt.timedelta(seconds=sec_max)

file_list=[]

for (dirpath, dirnames, filenames) in os.walk(path):


    for filename in filenames:
        f = '/'.join([dirpath,filename])
        ctime = os.stat(f).st_ctime 
        if ctime>=sec_min and ctime <=sec_max: 
           file_list.append(f)


print(file_list) 
