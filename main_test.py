##!/home/qcradar/.conda/envs/da/bin/python
#qc_path = "/home/qcradar/scripts/radar_qc/"
qc_path = "/home/jruiz/share/radar_qc/"

import sys
sys.path.append( qc_path + '/src/python/' )
sys.path.append( qc_path + '/src/fortran/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.

import conf_defaults as conf
import numpy as np

datapath = '/home/jruiz/test_qc/'  #Main data path.


instrument_list = ['RMA1',# 'RMA2', 'RMA3', 'RMA4', 'RMA5', 'RMA6',
                 # 'RMA7', 'RMA8',
                 # 'PAR',
                 'PER', 'ANG']         #Instrument list.

file_type_list = ['.H5','.VOL','.nc']

init_date='20180824200000'
end_date ='20180824210000'

#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , init_date , end_date , time_search_type='filename' , file_type_list = file_type_list )

print(file_list)

#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )

for radar in radar_list :

## llama a qc

    options = conf.options #This is the default configuration.

    options['name_ref'] ='ZH'         #Reflectivity
    options['name_v']   ='VRAD'       #Dopper velocity
    options['name_rho'] ='RHOHV'      #Rho HV
    options['toporawdatapath']= qc_path + "/data/terrain_data/raw/"
    options['toporadardatapath']= qc_path + "/data/terrain_data/radar/" 

    [ radar , qc_output ] = rqc.main_qc( options , radar )

