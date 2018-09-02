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
deltat = 600.0                         #Time window (seconds)
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1',# 'RMA2', 'RMA3', 'RMA4', 'RMA5', 'RMA6',
                 # 'RMA7', 'RMA8',
                 # 'PAR',
                 'PER', 'ANG']         #Instrument list.

file_type_list = ['.h5','.vol','.nc']

#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , deltat , time_offset , file_type_list )

#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )

for radar in radar_list :

# llama a qc

    options = conf.options #This is the default configuration.

    options['name_ref'] ='ZH'         #Reflectivity
    options['name_v']   ='VRAD'       #Dopper velocity
    options['name_rho'] ='RHOHV'      #Rho HV
    options['toporawdatapath']= qc_path + "/data/terrain_data/raw/"
    options['toporadardatapath']= qc_path + "/data/terrain_data/radar/" 

    [ radar , qc_output ] = rqc.main_qc( options , radar )


# manda qc_ref a superobbing??
