##!/home/qcradar/.conda/envs/da/bin/python
#qc_path = "/home/qcradar/scripts/radar_qc/"
qc_path = "/home/jruiz/share/radar_qc/"
so_path = "/home/jruiz/share/radar_so/"

import sys
sys.path.append( qc_path + '/src/python/' )
sys.path.append( qc_path + '/src/fortran/')
sys.path.append( so_path )
sys.path.append( so_path + '/fortran/')


import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
import radar_so as so                #Superobbing python classes


import conf_defaults as conf
import numpy as np
import datetime as dt
import os

os.system('export OMP_NUM_THREADS=12')

#datapath = '/home/jruiz/datosmate/RRA_radar/'  #Min data path.
datapath='/home/jruiz/share/DATA/DATOS_RADAR/RRA_radar/'

instrument_list = ['RMA1','RMA3','RMA4','RMA5','RMA6',
                 'RMA7', 'RMA8',
                 'PAR',
                 'PER', 'ANG']         #Instrument list.

file_type_list = ['cfrad']

time_delta = 300 #Time delta in seconds.

init_date='20170923060000'
end_date ='20170923060000'
#end_date ='20171003000000'

idate=dt.datetime.strptime(init_date,"%Y%m%d%H%M%S")
edate=dt.datetime.strptime(end_date ,"%Y%m%d%H%M%S")
cdate=idate
dated=dt.timedelta( seconds = time_delta )

while cdate <= edate    : 

   c_ini_date = cdate.strftime('%Y%m%d%H%M%S')
   c_end_date = (cdate + dated).strftime('%Y%m%d%H%M%S')

   #Obtenemos la lista de archivos.
   file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

   #Obtenemos la lista de objetos radares.
   radar_list = ot.read_multiple_files(  file_list , instrument_list )

   for radar in radar_list :

   ### llama a qc

      options = conf.options #This is the default configuration.

      options['name_ref'] ='ZH'         #Reflectivity
      options['name_v']   ='VRAD'       #Dopper velocity
      options['name_rho'] ='RHOHV'      #Rho HV
      options['toporawdatapath']= qc_path + "/data/terrain_data/raw/"
      options['toporadardatapath']= qc_path + "/data/terrain_data/radar/" 

      [ radar , qc_output ] = rqc.main_qc( options , radar )

      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')

      output_freq = 300
      #        dx    dz   zmax  rmax
      grid = [2000, 1000, 15e3, 240e3]
      opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
      outputpath = './superobbing'

      letkf_filelist = so.main_radar_so(radar, output_freq, grid, opts, outputpath)
      print(letkf_filelist)

      print('Uploading files to RELAMPAGO FTP')
      for ifile in letkf_filelist:
         # ot.upload_to_ftp(ifile, ftp_host, ftp_user, ftp_passwd)   
         pass

      print('=============================================================================')
      print(' END OF SUPEROBBING')
      print('=============================================================================')



   cdate = cdate + dated 

