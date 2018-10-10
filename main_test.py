#!/home/qcradar/.conda/envs/da/bin/python
qc_path = "/home/jruiz/Dropbox/DATA/"

import sys
sys.path.append( qc_path + '/radar_qc/' )
sys.path.append( qc_path + '/radar_qc/src/python/' )
sys.path.append( qc_path + '/radar_qc/src/fortran/')
sys.path.append( qc_path + '/radar_so/fortran/')
sys.path.append( qc_path + '/radar_so/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import radar_so as so                #Superobbing module
import matplotlib.pyplot as plt     
import conf_defaults_test as conf

datapath = '/ms-36/mrugna/RMA/datos/'  #Main data path.
datapath_out = '/home/jruiz/Dropbox/DATA/radar_qc/'    #Out data path
deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA3','RMA4','RMA5','RMA6',
                 'RMA7', 'RMA8','PAR','PER', 'ANG']         #Instrument list.

file_type_list = ['.h5','.vol','.nc']

current_date = dt.datetime.utcnow()

ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) ).strftime('%Y%m%d%H%M%S')
c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat ).strftime('%Y%m%d%H%M%S')


print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

#file_list=['./RMA2_0200_02_VRAD_20180921T133347Z.H5']
file_list=['./RMA1_0200_02_TH_20181010T114222Z.H5']

#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )

for radar in radar_list :

# llama a qc

      options = conf.options #This is the default configuration.

      print('=============================================================================')
      print(' QC')
      print('=============================================================================')
      print('')

      options['name_ref'] ='ZH'         #Reflectivity
      options['name_v']   ='VRAD'       #Dopper velocity
      options['name_rho'] ='RHOHV'      #Rho HV
      options['toporawdatapath']= qc_path + "/data/terrain_data/raw/"
      options['toporadardatapath']= qc_path + "/data/terrain_data/radar/"

      #print( np.max( radar.fields['VRAD']['data'] ) , np.min( radar.fields['VRAD']['data'] ) )

      #plt.pcolor( radar.fields['VRAD']['data'] )
      #plt.show()

      [ radar , qc_output ] = rqc.main_qc( options , radar )

      print( np.max( radar.fields['CVRAD']['data'] ) , np.min( radar.fields['CVRAD']['data'] ) )

      print('')
      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')
      print('')

      output_freq = 300
      #        dx    dz   zmax  rmax
      grid = [10000, 1000, 15e3, 240e3]
      opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
      outputpath = datapath_out

      
      #plt.pcolor( radar.fields['CVRAD']['data'] )
      #plt.show()

      letkf_filelist = so.main_radar_so(radar, output_freq, grid, opts, outputpath)
      print(letkf_filelist)

      print('Uploading files to RELAMPAGO FTP')
      for ifile in letkf_filelist:
         # ot.upload_to_ftp(ifile, ftp_host, ftp_user, ftp_passwd)   
         pass

      print('')
      print('=============================================================================')
      print(' END OF SUPEROBBING')
      print('=============================================================================')
      print('')


