#!/home/qcradar/.conda/envs/da/bin/python
qc_path = "/home/qcradar/scripts/radar_qc/"

import sys
sys.path.append( qc_path + '/src/python/' )
sys.path.append( qc_path + '/src/fortran/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.

import conf_defaults as conf

datapath = '/ms-36/mrugna/RMA/datos/'  #Main data path.
datapath_out = '/home/qcradar/data'    #Out data path
deltat = 600.0                         #Time window (seconds)
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA3','RMA4','RMA5','RMA6',
                 'RMA7','RMA8','PAR','PER','ANG']         #Instrument list.

file_type_list = ['.h5','.vol','.nc']

current_time = datetime.datetime.now()

ref_date=datetime(current_time.year, current_time.month, 1, 0, 0, 0)
freqtimes = freq*floor((current_date-ref_date).total_seconds()/freq)
c_end_date = current_date + timedelta( seconds=freqtimes )
c_ini_date = current_date + timedelta( seconds=freqtimes-deltat )

print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

print(file_list)
