

def read_multiple_files(  filelist , instrument_list = None )

   import numpy as np
   #filelist is a list of filenames corresponding to radar files in different formats.
   #instrument_list (if present) a list of instrument names. Only names in instrument_list will be incorporated to the 
   #radar object list.

   #Generate a list of radar objects
   
   radar_list   = [] 
   radar_used   = []
   radar_remove = []
   
   final_radar_list = []

   used_file = np.zeros( len( filelist ) ).astype(bool)

   for ifile , filename in enumerate( filelist )    :

      file_format = get_file_format( filename )
 
      file_time   = get_time_from_filename( filename )

      file_instrument = get_instrument_type_from_filename( filename )

      if ( ( file_instrument in instrument_list ) or ( instrument_list == None ) ) and ( not used_file[ifile] ) :
         #Read the radar 
         my_radar = read_file( filename , file_format )  

         if len( tmp_radar ) != 0  :
            #Add the current radar to the list. 
            #radar_list.append = my_radar
            radar_used[ifile] = True
            my_radar.files = [ filename ]

            #Check if we can associate other radars in the list to the current radar.
  
            for tmp_ifile , tmp_filename  in enumerate( file_list )   :
  
               if ( not radar_used[tmp_ifile] ) and ( tmp_ifile != ifile )       :
                  if ( file_format     ==  get_file_format( tmp_filename ) ) or 
                     ( file_time       ==  get_time_from_filename( tmp_filename ) ) or
                     ( file_instrument ==  get_instrument_type_from_filename( tmp_filename ) )

                     #Read the data
                     tmp_radar = read_file( filename , file_format )
                     if len( tmp_radar ) != 0  :
                        
                        #Check if we can merge current_radar and my_radar objects.
                        if np.shape( my_radar.azimuth['data'] )[0] == np.shape( tmp_radar.azimuth['data'] )[0]   :
                           if np.sum( my_radar.azimuth['data'] - tmp_radar.azimuth['data'] ) == 0                :
                            #my_radar corresponds to the same instrumet and initial time as current_radar
                            #merge tmp_radar into my_radar
                            my_radar.files.append( tmp_filename )
                            used_file[ tmp_ifile ] = True
                            for my_key in my_radar.fields   :
                               if not my_key in my_radar.fields   :
                                  my_radar.fields[my_key] = tmp_radar.fields.pop(my_key)
                        else                                                                                     :
                           print('Warning: Inconsistent shapes found for ' + file_instrument + ' ' + file_time )
                           print('This volume will be processed separately')  

         #So far my_radar contains all the variables corresponding to this instrument and initial time.
         radar_list.append( my_radar )       

         print('RADAR : ' + my_radar.metadata['instrument_name'] + ' ' + my_radar.metadata['start_datetime'] )
         for ifile in my_radar.files   :
            print('   FILE: ' + ifile )
         
   return radar_list

def get_file_format( filename )            :

      if ('.h5' in filename ) or ( '.H5' in filename )   :
         file_format = 'h5'
      if ( '.vol' in filename ) or ( 'VOL' in filename ) :
         file_format = 'vol' 
      if ( '.nc'  in filename ) or ( 'NC' in filename )  :
         file_format = 'nc' 

   return file_format

def read_file( filename , format_file )    :
   import pyart
   from pyart.aux_io.sinarame_h5 import read_sinarame_h5
   from pyart.aux_io.rainbow_wrl import read_rainbow_wrl


   try   :
      if format_file == 'h5'   :
         radar = read_sinarame_h5(filename, file_field_names=True)
      if format_file == 'vol'  :
         radar = read_rainbow_wrl(filename, file_field_names=True)
      if format_file == 'nc'   :
         radar = pyart.io.read(filename)

      print( 'Reading file:' + filename )
      for my_key in radar.fields         
         print( 'Found the following variable: ' + my_key )

   except ValueError   :
         print('Warning: Could not read file ' + filename )
         radar = dict()
   except KeyError     :
         print('Warning: Could not read file ' + filename )
         radar = dict()

   return radar

def rename_fields ( radar )  :

     # Reflectividad RMA
     if 'TH' in radar.fields
        radar.fields['ZH'] = radar.fields.pop('TH')

     if 'V' in radar.fields
        radar.fields['VRAD'] = radar.fields.pop('V')

     if 'ZH' in radar.fields
        radar.fields['ZH'] = radar.fields.pop('dBZ')

     if 'W' in radar.fields
        radar.fields['WRAD'] = radar.fields.pop('W')

   return radar   

def get_strat ( filename , radar )  :

#Include some parameters that are required by the QC module.

    import numpy as np
    import numpy.ma as ma
    import os

    local_fill_value = -9999.0
    levels=np.unique(radar.elevation['data'])

    #Set the instrument name 
    radar.metadata['instrument_name'] = get_instrument_name_from_filename( filename )

    #Add missing structures to the radar object.
    if radar.altitude_agl == None :
       radar.altitude_agl = dict()
    if radar.metadata == None :
       radar.metadata = dict()
    if radar.instrument_parameters == None :
       radar.instrument_parameters = dict()

    if (not  'nyquist_velocity' in radar.instrument_parameters ) or ( radar.instrument_parameters['nyquist_velocity'] == None )  :

       radar.instrument_parameters['nyquist_velocity']=dict()
       radar.instrument_parameters['nyquist_velocity']['long_name']='unambiguous_doppler_velocity'
       radar.instrument_parameters['nyquist_velocity']['units']='meters per second'
       radar.instrument_parameters['nyquist_velocity']['_FillValue']= local_fill_value
       radar.instrument_parameters['nyquist_velocity']['meta_group']='instrument_parameters'

    if (not 'radar_beam_width_v' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_v'] == None ) :
       radar.instrument_parameters['radar_beam_width_v']=dict()
       radar.instrument_parameters['radar_beam_width_v']['long_name']='half_power_radar_beam_width_v_channel'
       radar.instrument_parameters['radar_beam_width_v']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_v']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_v']['meta_group']='instrument_parameters'

    if (not 'radar_beam_width_h' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_h'] == None ) :
       radar.instrument_parameters['radar_beam_width_h']=dict()
       radar.instrument_parameters['radar_beam_width_h']['long_name']='half_power_radar_beam_width_h_channel'
       radar.instrument_parameters['radar_beam_width_h']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_h']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_h']['meta_group']='instrument_parameters'

    if radar.range == None :
       radar.range = dict()
    if radar.ray_angle_res == None :
       radar.ray_angle_res = dict()
       radar.ray_angle_res['long_name']='angular_resolution_between_rays'
       radar.ray_angle_res['units']='degrees'
       radar.ray_angle_res['_FillValue']= local_fill_value

    #Get the corresponding radar strategy depending on the filename.
    radar.altitude_agl['data'] = 0.0


    if 'RMA' in filename  :


       if '0117_01' in filename  :  #9005-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0117_02' in filename  :  #9005-2 STRATEGY
          nyquist_velocity     = 33.04
       if '0117_03' in filename  :  #9005-3 STRATEGY
          nyquist_velocity     = 3.98

       if '0117_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0117_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0121_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0121_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0122_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 8.28
       if '0122_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 39.79
       if '0122_03' in filename  :  #122-3 STRATEGY
          nyquist_velocity     = 13.35

       if '0123_01' in filename  :  #123-1 STRATEGY
          nyquist_velocity     = 8.28
       if '0123_02' in filename  :  #123-2 STRATEGY
          nyquist_velocity     = 39.79
       if '0123_03' in filename  :  #123-3 STRATEGY
          nyquist_velocity     = 13.25
       if '0123_04' in filename  :  #123-4 STRATEGY
          nyquist_velocity     = 8.28

       if '0200_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0200_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0201_01' in filename  :  #201-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0201_02' in filename  :  #201-2 STRATEGY
          nyquist_velocity     = 13.25
       if '0201_03' in filename  :  #201-3 STRATEGY
          nyquist_velocity     = 8.28

       if '0202_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0202_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25


    if ( 'PAR' in filename ) or ( 'ANG' in filename ) or ( 'PER' in filename )  :

       if '/120/' in filename   : 
          nyquist_velocity = 39.8  #120

       if '/240/' in filename   :
          nyquist_velocity = 6.63  #240

    #Correct instrument altitude.
    if 'PAR' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(122.0)

    if 'ANG' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(190.0)

    if 'PER' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(100.0)

    if 'RMA1' in filename  :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(476.0)

    if 'RMA2' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array( 47.0)

    if 'RMA3' in filename  :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(197.0)

    if 'RMA4' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(119.0)

    if 'RMA5' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(841.0)

    if 'RMA6' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array( 80.0)

    if 'RMA7' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(388.0)

    if 'RMA8' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(131.0)

    #Some common parameters
    ray_angle_res        = 1.0
    radar_beam_width_h   = 1.0
    radar_beam_width_v   = 1.0
    meters_between_gates  = radar.range['data'][1]-radar.range['data'][0]

    #Apply the missing parameters to the radar structure.
    radar.instrument_parameters['nyquist_velocity']['data'] = ma.array(np.ones( np.shape(radar.azimuth['data']) )*nyquist_velocity , mask = np.zeros( np.shape(radar.azimuth['data']) , dtype=bool ) , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_h']['data'] = ma.array( radar_beam_width_h , mask =False , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_v']['data'] = ma.array( radar_beam_width_v , mask =False , fill_value = local_fill_value )
    radar.ray_angle_res['data'] = ma.array( np.ones( np.shape( levels ) )*ray_angle_res , mask = np.zeros( np.shape( levels ) , dtype=bool ) , fill_value = local_fill_value )
    radar.range['meters_between_gates']= meters_between_gates
    radar.range['meters_to_center_of_first_gate']= radar.range['data'][0]  #meters_between_gates / 2.0

    return radar


def get_file_list( datapath , init_time , end_time , time_search_type = None , file_types_list = None )

   #datapath : base path of radar data
   #init time: [yyyymmddhhMMss] beginning of the time window
   #end time : [yyyymmddhhMMss] end of the time window
   #time_search_type : [filename] or [timestamp]
   #file_types_list  : a list with file extensions that will be included in the file_list

   import os
   import datetime as dt 
   import numpy as np

   if time_search_type == None :
      time_search_type = 'timestamp'

   date_min = datetime.strptime( init_time , '%Y%m%d%H%M%S')
   date_max = datetime.strptime( end_time  , '%Y%m%d%H%M%S')

   file_list=[]

      for (dirpath, dirnames, filenames) in os.walk( datapath ):

         for filename in filenames            :
            f = '/'.join([dirpath,filename])

            if time_search_type == 'filename'   :
               date_c = get_time_from_filename( filename )
            if time_search_type == 'timestamp'  :
               date_c = dt.fromtimestamp( os.stat(f).st_ctime )

            if date_c >= date_min and date_c <= date_max  :
               file_list.append(f)
   
   #Keep only some file names and some paths.

   final_file_list[]

   if file_type_list != None :

      for my_file in file_list  :

         filename = os.path.basename( my_file )

         if any(ft in filename for ft in file_type_list ):
 
            final_file_list.append( my_file )

   else                      :

      final_flie_list  = file_list 


return file_list



def get_time_from_file_name( filename )    :

   import datetime as dt
   import os 

   filename = os.path.basename( filename )
   file_time = None

   if ( 'PAR' in filename ) or ( 'PER' in filename ) or ( 'ANG' in filename )  :
      file_time  = datetime.strptime(filename[:14], '%Y%m%d%H%M%S')

   if ( 'RMA' in filename )  :
      file_time  = datetime.strptime(filename.split('_')[-1][:15], '%Y%m%dT%H%M%S')  

return file_time


def get_instrument_name_from_file_name( filename ) :
   import os

   filename = os.path.basename( filename )
   instrument_name = None
   if 'RMA' in filename    :
      instrument_name = ( os.path.basename( filename ) ).split('_')[0]
   if 'ANG' in filename    :
      instrument_name = 'ANG'
   if 'PAR' in filename    :
      instrument_name = 'PAR'
   if 'PER' in filename    :
      instrument_name = 'PER'

