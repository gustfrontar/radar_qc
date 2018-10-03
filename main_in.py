#!/home/qcradar/.conda/envs/da/bin/python
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy

#=========================================================================================================
# CONFIGURATION SECTION
#=========================================================================================================

#General section
qc_path = "/home/qcradar/scripts/"

datapath = '/ms-36/mrugna/RMA/datos/'  #Main data path.
datapath_out = '/home/qcradar/data'    #Out data path
deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
deltat_archive = dt.timedelta( seconds=86400 ) #Time window that will be kept in the remote ftp server.
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA4','RMA6','RMA8','RMA11','PAR','PER','ANG']  #Instrument list.

file_type_list = ['.H5','.vol','.nc']

remove_local_pkl = True                #Remove intermediate gridded data in pkl format.
remove_local_dat = False               #Remove gridded data in letkf format.
remove_remote_dat = True               #Remove remote letkf files.
write_cfradial=True                    #Write a cfradial file with the qced data
#Qc section
name_ref='ZH'                                               #Reflectivity name
name_v  ='VRAD'                                             #Doppler velocity name
name_rho='RHOHV'                                            #RhoHV name
toporawdatapath=qc_path + "/data/terrain_data/raw/"         #Raw topography data path (tif format)
toporadardatapath=qc_path + "/data/terrain_data/radar/"     #Interpolated topography data path (unformatted)

#Superobbing section 
output_freq = 600
#        dx    dz   zmax  rmax
grid = [10000, 1000, 15e3, 240e3]
opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
outputpath = datapath_out

#Ftp section
ftp_host='ftp.smn.gob.ar'
ftp_user='rra'
ftp_pass='TwV/27gm7YQsw'
ftp_path='radar'
compress=True

#=========================================================================================================
# END OF CONFIGURATION SECTION
#=========================================================================================================

import sys
sys.path.append( qc_path + '/radar_qc/' )
sys.path.append( qc_path + '/radar_qc/src/python/' )
sys.path.append( qc_path + '/radar_qc/src/fortran/')
sys.path.append( qc_path + '/radar_so/fortran/')
sys.path.append( qc_path + '/radar_so/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
import radar_so as so                #Superobbing module

options = conf.options #This is the default configuration.
options['name_ref'] = name_ref         #Reflectivity
options['name_v']   = name_v           #Dopper velocity
options['name_rho'] = name_rho         #Rho HV
options['toporawdatapath']= toporawdatapath
options['toporadardatapath']= toporadardatapath 

#Set the dates that will be processed. 
current_date = dt.datetime.utcnow()

ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat   ).strftime('%Y%m%d%H%M%S')
c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat*2 ).strftime('%Y%m%d%H%M%S')

#Set the dates that will be archived. 
a_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) + deltat         ).strftime('%Y%m%d%H%M%S')
a_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat_archive ).strftime('%Y%m%d%H%M%S')



print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

#Set the dates that will be archived in the remote server.
a_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) ).strftime('%Y%m%d%H%M%S')
a_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat_archive ).strftime('%Y%m%d%H%M%S')

print('')
print('=============================================================================')
print(' GETTING FILE LIST ')
print('=============================================================================')
print('')


#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

#file_list = ['/ms-36/mrugna/RMA/datos/RMA6/2018/09/25/06/2737/RMA6_0200_01_TH_20180925T062737Z.H5']

print(file_list)

print('')
print('=============================================================================')
print(' READING FILE LIST ')
print('=============================================================================')
print('')

#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )

for radar in radar_list :

      print('=============================================================================')
      print(' QC')
      print('=============================================================================')
      print('')
 
      #Call QC routine

      [ radar , qc_output ] = rqc.main_qc( options , radar )

      #Save data in cfradial format

      if write_cfradial  :

         print('=============================================================================')
         print(' WRITING QC OUTPUT IN CFRADIAL FORMAT')
         print('=============================================================================')
         print('')


         ot.save_cfradial( datapath_out + '/cfradial/' , radar )

      print('')
      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')
      print('')

      #Call SO routine 

      letkf_filelist = so.main_radar_so(radar, output_freq, grid, opts, datapath_out  )

      print('')
      print('=============================================================================')
      print(' UPLOADING FILES TO REMOTE FTP SERVER ' + ftp_host )
      print('=============================================================================')
      print('')

      #Call remote server uploading routine

      ot.upload_to_ftp( letkf_filelist , ftp_host, ftp_user, ftp_pass , ftp_path , compress=compress )


print('')
print('=============================================================================')
print('We will keep all the files within the following dates:' )
print( a_ini_date )
print( a_end_date )
print('=============================================================================')
print('')


print('')
print('=============================================================================')
print(' REMOVE OLD FILES FROM REMOTE SERVER ' + ftp_host )
print('=============================================================================')
print('')

#Call remote server deleting routine

if remove_remote_dat :

   ot.remove_from_ftp_timebased( ftp_host, ftp_user, ftp_pass , ftp_path , a_ini_date , a_end_date , file_format_list = ['letkf'] ) 

print('')
print('=============================================================================')
print(' REMOVE OLD FILES FROM LOCAL SERVER ' + datapath_out )
print('=============================================================================')
print('')

#Call local server deleting routine 
tmp_format_list=[]
if remove_local_pkl : 
   tmp_format_list.append('pickle')
if remove_local_dat :
   tmp_format_list.append('letkf')

ot.remove_from_localpath_timebased( datapath_out , a_ini_date , a_end_date , file_format_list = tmp_format_list )




