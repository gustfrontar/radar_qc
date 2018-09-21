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
                 'RMA7', 'RMA8','PAR','PER', 'ANG']         #Instrument list.

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

      [ radar , qc_output ] = rqc.main_qc( options , radar )

      print('')
      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')
      print('')

      output_freq = 300
      #        dx    dz   zmax  rmax
      grid = [2000, 1000, 15e3, 240e3]
      opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
      outputpath = '/home/qcradar/data/'

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


