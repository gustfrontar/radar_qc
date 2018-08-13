#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

def main_qc( filename , options ) :

   import sys
   import time
   sys.path.append('../fortran')

   import numpy as np
   import numpy.ma as ma
   import pyart
   from common_qc_tools  import qc  #Fortran code routines.
   import netCDF4
   import os
   from netCDF4 import Dataset


   start=time.time()

   output=dict() #Initialize output dictionary.

   computed_etfilter=False   #Flag to indicate if echo top has been computed already.

   #Read the data
   radar = pyart.io.read(filename)

   if options['is_rma']  :

      radar = get_rma_strat( filename , radar ) 

   if radar.altitude_agl['data'] == 0.0    :
      #Radar is assumed to be at ground level. Add 20 meters to account for the height of the tower.
      radar.altitude['data']=radar.altitude['data']+options['radar_altitude_agl']


   #Get the nyquist velocity
   if not options['is_rma']  :
      nyquistv=radar.get_nyquist_vel(0,check_uniform=True)
      

   output['maxw_ref']=0.0
   output['maxw_v']  =0.0


   #===================================================
   # SOME PRELIMINARY PROCEDURES.
   #===================================================

   #From cfradial format to azimuth,range,elevation order

   print('')
   print('-------------------------------------------')
   print('Reshaping the variables')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = reshape_variables( radar , output , options )
   #Compute X,Y,lat and lon for each pixel

   print('')
   print('-------------------------------------------')
   print('Georeferencing the data')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = georeference_data( radar , output , options )

   print('')
   print('-------------------------------------------')
   print('Interpolating the terrain')
   print('-------------------------------------------')
   print('')

   #Compute the terrain height corresponding to each pixel
   [ radar , output ] = get_topography( radar , output , options )

   #===================================================
   # FILTER APPLICATION
   #===================================================

   #This functions runs the selected filters in the selected
   #order.

   print('')
   print('-------------------------------------------')
   print('Quality control routines')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = run_filters( radar , output , options )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("the entire QC",end-start) )

   return radar , output


#===================================================
# FUNCIONES GENERALES ------------------------------
#===================================================


#===================================================
# RESHAPE VARIABLES
#===================================================

def reshape_variables( radar , output , options )    :

   import numpy as np
   import time
   #From time,range -> azimuth,range,elevation

   start=time.time()

   if options['name_ref'] in radar.fields :

        output['undef_ref']=radar.fields[options['name_ref']]['_FillValue']

        [ output['ref'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact'] ]=order_variable( radar , options['name_ref'] , output['undef_ref'] )
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cref'] = np.copy(output['ref'])           #Initialize the corrected reflectivity array.
        output['ref_input'] = np.copy(output['ref'])      #Initialize the input reflectivity array (for ploting only)
 
        output['cref'][ output['cref'] == output['undef_ref'] ]=options['norainrefval']
        output['ref'] [ output['ref']  == output['undef_ref'] ]=options['norainrefval']

        output['qcref'] = np.zeros(output['cref'].shape)  #Set the qc flag array to 0.
        output['wref']  = np.zeros(output['cref'].shape)  #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])

   if options['name_v'] in radar.fields  : 

        output['undef_v']=radar.fields[ options['name_v'] ]['_FillValue']

        [ output['v'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , options['name_v'] , output['undef_v']  )
 
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne
 
        output['cv'] = np.copy(output['v'])              #Initialize the corrected velocity array
        output['input_v'] = np.copy(output['v'])         #Initialize the input velocity array (for ploting only)

        output['qcv'] = np.zeros(output['v'].shape)      #Set the qc flag array to 0.

        output['wv'] = np.zeros(output['v'].shape)       #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])

   end=time.time()
   print("The elapsed time in {:s} is {:2f}".format("Data reshape",end-start) )

   return radar , output

#===================================================
# GEOREFERENCE RADAR DATA
#===================================================

def georeference_data( radar , output , options )    :

   import numpy as np
   import time
   #Use pyart rutines to compute x, y and z at each grid point
   #Reorder data
   start=time.time()

   #dm is a dummy variable

   [ output['altitude'] , dm , dm , dm , dm , dm ] = order_variable( radar , 'altitude' , options['undef'] )
   [ output['x']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'x' , options['undef'] ) 
   [ output['y']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'y' , options['undef'] )


   #Compute distance to radar for the first azimuth (we assume that distance to radar will be the same.
   #for all other azimuths.
   output['distance']=np.power( np.power(output['x'][0,:,:],2)+np.power(output['y'][0,:,:],2) , 0.5 )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("Data georeferentiation",end-start) )

   return radar , output 

#===================================================
# READ AND INTERPOLATE TOPOGRAPHY DATA
#===================================================

def get_topography( radar , output , options )  :

   import time
   import numpy as np
   import os
   from scipy.interpolate import interp2d

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   start=time.time()

   #The idea is to store the topography interpolated to the earth surface in a polar geometry 
   #sourrounding each radar into a binary file. If this binary file is already created we read it
   #and interpolate the topography to the 3D radar geometry. If not we read the raw raster data.

   #Name of the file containing the topography interpolated the earth surface surrounding the radar.
   polar_coord_topo_file=options['toporadardatapath'] + '/' + radar.metadata['instrument_name'] + '_' +  \
             str( radar.range['meters_between_gates'] ) + '_' + str( np.size( radar.range['data'] ) * radar.range['meters_between_gates'] / 1000 ) \
             + '.dat'
   raster_topo_file=options['toporadardatapath'] + '/' + radar.metadata['instrument_name'] + '_' +  \
             str( radar.range['meters_between_gates'] ) + '_' + str( np.size( radar.range['data'] ) * radar.range['meters_between_gates'] / 1000 ) \
             + '.tif'

   if ( os.path.isfile( polar_coord_topo_file ) )  :
      read_raster = False
      print('Using a previously generated topography file')
   else                                :
      read_raster = True
      print('Topography file not found. We will generate a new file from raw raster data')

   if read_raster    :   #We read the original data and interpolate it to a polar grid centered at the radar.

      #Generate topo file ( inputs are , radar lat and lon, range and azimuth )
      my_topo=generate_topo_file( radar.longitude['data'][0] , radar.latitude['data'][0] , radar.range['data'] , output['az'] , 
                                  options['toporawdatapath'] , polar_coord_topo_file )

   else   :

      print('Reading polar coordinate topography from a file')

      my_topo=read_topo( polar_coord_topo_file )

   #Interpolate topography data to the 3D radar structure.

   topo_r = my_topo['range'][0,:]
   topo_a = my_topo['azimuth'][:,0]

   #Define a "interpolator"
   interpolator = interp2d(topo_r,topo_a,my_topo['mean'], kind='linear')

   output['topo'] = np.zeros((output['na'],output['nr'],output['ne']))

   #Loop over the vertical levels.
   for ie in range( 0 , ne )   :

      ppi_r = output['distance'][:,ie]
      ppi_a = output['az']

      output['topo'][:,:,ie] =  interpolator( ppi_r, ppi_a )

   end=time.time()

   #Compute AGL height for each pixel.
   output['altitude_agl']=output['altitude']-output['topo'] 

   print("The elapsed time in {:s} is {:2f}".format("topography interpolation",end-start) )

   return radar , output

#===================================================
# COMPUTE THE FINAL WEIGHT AND APPLY FUZZY LOGIC QC
#===================================================
 
def weight_computation( output , options )   :

   import numpy as np

   if output['maxw_ref'] == 0    :
      output['wref'][:]=0.0
   else                          :
      output['wref']=output['wref'] / output['maxw_ref']

   if output['maxw_v'] == 0      :
      output['wv'][:]=0.0
   else                          :
      output['wv']=output['wv'] / output['maxw_v']

   #All the grid points where the probability of non-meteorological echoes
   #is greather than 0.5 will be flagged as missing.
   output['cref'][ output['wref'] > options['w_tr'] ] = output['undef_ref']
   output['cv'][output['wv'] > options['w_tr'] ] = output['undef_v']

   return output

#===================================================
# ADD CORRECTED DATA TO RADAR OBJECT
#===================================================

def  update_radar_object( radar , output , options )   :

   import numpy as np

   if options['name_ref'] in radar.fields :

      tmp=order_variable_inv( radar , output['cref'] , output['index'] , output['undef_ref'] )

      radar.fields[ options['name_cref'] ] = dict()

      radar.fields[ options['name_cref'] ] = radar.fields[ options['name_ref'] ]

      radar.fields[ options['name_cref'] ]['data']=np.ma.masked_array(tmp , tmp==output['undef_ref'] )

   if options['name_v'] in radar.fields :

      tmp=order_variable_inv( radar , output['cv'] , output['index'] , output['undef_v'] )

      radar.fields[ options['name_cv'] ]= dict()

      radar.fields[ options['name_cv'] ]= radar.fields[ options['name_v'] ]

      radar.fields[ options['name_cv'] ]['data']=np.ma.masked_array(tmp , tmp==output['undef_v'] )

   return radar , output 


#====================================================
# UPDATE OUTPUT BASED ON FILTER DATA
#====================================================

def output_update( output , qc_index , options , filter_name )  :

   from common_qc_tools  import qc  #Fortran code routines.
   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   #Compute the corresponding weigth.
   weigth = qc.multiple_1d_interpolation( field=qc_index , nx=na , ny=nr , nz=ne
                                         , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                         , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

   if 'ref' in options[filter_name]['var_update_list']   : 
      #Update reflectivity
      if options[filter_name]['sequential']   :
         #Preserve the data used as input for the filter for ploting pourposes.
         output['input_ref']=np.copy(output['ref'])

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + weigth * options[filter_name]['w']
         output['qcref'][ weigth > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         if options[filter_name]['fill_value']  == 'undef'       :
            output['cref'][ weigth > options[filter_name]['force_value'] ]=output['undef_ref']
            if options[filter_name]['sequential']   :
               output['ref'][ weigth > options[filter_name]['force_value'] ]=output['undef_ref']
         elif options[filter_name]['fill_value']  == 'min_ref'   :
            output['cref'][ weigth > options[filter_name]['force_value'] ]=options['norainrefval'] 
            if options[filter_name]['sequential']   :
               output['ref'][ weigth > options[filter_name]['force_value'] ]=options['norainrefval']
         else                                                    :
            output['cref'][ weigth > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']
            if options[filter_name]['sequential']   :
               output['ref'][ weigth > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']


         output['qcref'][ weigth > options[filter_name]['force_value'] ] = options[filter_name]['code']
   
   if 'v' in options[filter_name]['var_update_list']   :
      #Update radial velocity
      if options[filter_name]['sequential']   :
         #Preserve the data used as input for the filter for ploting pourposes.
         output['input_v']=np.copy(output['v'])

      if not options[filter_name]['force']   :
         output['wv']=output['wv'] + weigth * options[filter_name]['w']
         output['qcv'][ weigth > 0.5 ] = options[filter_name]['code']
         output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
      else                                   :
         if (options[filter_name]['fill_value']  == 'undef') or  (options[filter_name]['fill_value']  == 'min_ref')  :
            output['cv'][ weigth > options[filter_name]['force_value'] ]=output['undef_v']
            if options[filter_name]['sequential']   :
               output['v'][ weigth > options[filter_name]['force_value'] ]=output['undef_v']
         else                                                    :
            output['cv'][ weigth > options[filter_name]['force_value'] ]= options[filter_name]['fill_value']
            if options[filter_name]['sequential']   :
               output['v'][ weigth > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']


         output['qcv'][ weigth > options[filter_name]['force_value'] ] = options[filter_name]['code']


   #If "save" is true, then store the qc_index and the weigth computed for this particular 
   #filter. This will be stored in a dictionary named after the filter and which is part of the 
   #output dictionary. 
   if  options[filter_name]['save']      :
      output[filter_name]=dict()
      output[filter_name]['qc_index']=qc_index
      output[filter_name]['weight']=weigth

   #Plot filter diagnostics if this option is available.
   if  options['plot']['Enable']    :
      plot_filter( output , qc_index , weigth , options , filter_name )

   return output

#====================================================
# PLOT THE FILTER
#====================================================

def plot_filter( output , qc_index , weigth , options , filter_name )  :

   import numpy as np
   import matplotlib.pyplot as plt

   if 'ref' in options[filter_name]['var_update_list']   :

    tmp_ref=np.ma.masked_array( output['ref_input'] , output['ref_input'] == output['undef_ref'] )
    tmp_cref=np.ma.masked_array( output['cref'] , output['cref'] == output['undef_ref'] )
    tmp_qc_index=np.ma.masked_array( qc_index , qc_index == output['undef_ref'] )

    for ilev in options['plot']['Elevs']  :

       plt.figure(figsize=(8, 8))
       plt.subplot(2,2,1)

       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, tmp_ref[:,:,ilev],vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
       plt.title('Original Reflectivity')
       plt.colorbar()

       plt.subplot(2,2,2)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, tmp_cref[:,:,ilev],vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
       plt.title('Corrected Reflectivity')
       plt.colorbar()

       plt.subplot(2,2,3)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, ( output['qcref'][:,:,ilev]==options[filter_name]['code'] ).astype(float) )
       plt.title('Pixels corrected by ' + filter_name)
       plt.colorbar()

       plt.subplot(2,2,4)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, ( tmp_qc_index[:,:,ilev] ) )
       plt.title('QC Index')
       plt.colorbar()
                                                                                                    
       if options['plot']['Show']  :

           plt.show()

       figname=options['plot']['FigNamePrefix'] + '_elev_' + str(ilev) + '_' + filter_name + '_ref_' + options['plot']['FigNameSufix']
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)
        
       plt.close()

   if 'v' in options[filter_name]['var_update_list']   :

    tmp_v=np.ma.masked_array( output['input_v'] , output['input_v'] == output['undef_v'] )
    tmp_cv=np.ma.masked_array( output['cv'] , output['cv'] == output['undef_v'] )
    tmp_qc_index=np.ma.masked_array( qc_index , qc_index == output['undef_v'] )

    for ilev in options['plot']['Elevs']   :

       plt.figure(figsize=(8, 8))
       plt.subplot(2,2,1)

       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3,tmp_cv[:,:,ilev],vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Corrected Doppler Velocity')
       plt.colorbar()

       plt.subplot(2,2,2)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3,tmp_v[:,:,ilev],vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Original Doppler Velocity')
       plt.colorbar()

       plt.subplot(2,2,3)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, ( output['qcv'][:,:,ilev]==options[filter_name]['code'] ).astype(float) )
       plt.title('Pixels corrected by ' + filter_name)
       plt.colorbar()

       plt.subplot(2,2,4)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, tmp_qc_index[:,:,ilev] )
       plt.title('QC index')
       plt.colorbar()

       if options['plot']['Show']  :

           plt.show()

       figname=options['plot']['FigNamePrefix'] + '_elev_' + str(ilev) + '_' + filter_name + '_v_' + options['plot']['FigNameSufix']
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)

       plt.close()

    return

#====================================================
# RUN THE SELECTED FILTERS IN THE CORRESPONDING ORDER
#====================================================

def run_filters( radar , output , options )   :

   import time
   import numpy as np

   my_filters       = list()

   my_filters_order = list()

   #First get the list of enabled filters from the options structure.
   #Note that 'order' can include several instances so filters can be 
   #Applied several times to the same data set
  
   for key in options  :
       
       if type( options[key] ) == dict  :

          if 'flag' in options[key]  :
           
              if options[key]['flag']  :
                 
                 for ii in range(len(options[key]['order'])) :
 
                     #Filters can be applied several times to the same data.
                     my_filters.append( key )
                     my_filters_order.append( options[key]['order'][ii] )

   #Sort the filters according to the order associated with each filter.   
 
   my_sorted_filters=[x for _,x in sorted(zip(my_filters_order,my_filters))]
   print('-------------------------------------------')
   print('The following filters will be applied')
   print('In the following order')
   print('-------------------------------------------')
   print('')
   for ifilter in my_sorted_filters  :
       print( ifilter )

   print('')
   print('-------------------------------------------')
              
   print('')
   print('-------------------------------------------')
   print('Running the filters')
   print('-------------------------------------------')
   print('')

   #Run the filtesr in the selected order.
   #If ploting is enabled then plot diagnostics.

   for ifilter in my_sorted_filters :

       start=time.time()

       if globals().get(ifilter) is None  :
          #If I can not find the corresponding function raise an error.
          print("Error: No filter called {}".format(ifilter))
       else                                      :
          #Run the filter (call the function named after the filter name.
          print('')
          print('-------------------------------------------')
          print('Running filter ' + ifilter )
          print('-------------------------------------------')
          print('')
          exec("[ radar , output ] = " + ifilter + "( radar , output , options )")

          end=time.time()
          print('')
          print("The elapsed time in {:s} is {:2f}".format(ifilter,end-start) )
          print('')

#       if globals().get('plot_' + ifilter) is None  :
#          #If I can not find the corresponding function raise an error.
#          print("Error: No plotting function for {}".format(ifilter))
#       else                                      :
#          #Run the filter (call the function named after the filter name.
#          print('')
#          print('-------------------------------------------')
#          print('Plotting diagnostics for ' + ifilter )
#          print('-------------------------------------------')
#          print('')
#          exec("plot_" + ifilter + "( output , options )")
#          print('')
#          end=time.time()
#          print("The elapsed time in ploting {:s} is {:2f}".format(ifilter,end-start) )
#          print('')

   print('')
   print('-------------------------------------------')
   print('Finish running the filters')
   print('-------------------------------------------')
   print('')


   return radar , output 

#===================================================
# FUNCIONES PARA LA IMPLEMENTACION DE DIVERSOS 
# FILTROS QUE COMPONEN EL QC.
#===================================================
 
#===================================================
# DEALIASING 
#===================================================

def Dealiasing( radar , output , options )   :

   import numpy as np
   import pyart

   filter_name='Dealiasing'
   if options[filter_name]['flag']  and ( options['name_v'] in radar.fields ) :
      #Define a new instance of the radar strcture containing the wind (potentially affected by previous filters).
      radar.fields[ options['name_cv'] ] = dict()
      radar.fields[ options['name_cv'] ] = radar.fields[ options['name_v'] ]
      tmp = order_variable_inv(  radar , output['v'] , output['index'] , output['undef_v'] )
      radar.fields[ options['name_cv'] ]['data']= np.ma.masked_array( tmp , tmp==output['undef_v'] )

      #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
      winddealias=pyart.correct.region_dealias.dealias_region_based(radar,interval_splits=options[filter_name]['interval_split'],interval_limits=None, 
                 skip_between_rays=options[filter_name]['skip_between_ray'],skip_along_ray=options[filter_name]['skip_along_ray'],centered=True,nyquist_vel=None,
                 check_nyquist_uniform=True,gatefilter=False,rays_wrap_around=True,keep_original=False,set_limits=True,
                 vel_field=options['name_cv'] ,corr_vel_field=None)

      #Replace cv wind by dealiased winds.
      radar.fields[ options['name_cv'] ]['data'] = winddealias['data']

      #Re-order dealiased wind data.
      [ output['cv'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , options['name_cv'] , output['undef_v']  )
    
      mask=np.logical_and( output['cv'] != output['v'] , output['cv'] != output['undef_v'] )
      output['qcv'][ mask ]=options[filter_name]['code']

      #Plot the data (this filter do not call update_output so we need to plot it here).
      tmp_index  = output['cv'] - output['v']
      tmp_weigth = np.zeros( np.shape( tmp_index ) )
      plot_filter( output , tmp_index , tmp_weigth , options , filter_name )   

      if options[filter_name]['sequential']     :
         #Following filters will be computed using the corrected velocity.
         output['v'] = output['cv'] 

   
   return radar , output 

#===================================================
# MODEL FILTER 
#===================================================

def ModelFilter( radar , output , options)   :

   import numpy as np

   #TODO
   #TODO
   print('This filter is not ready yet')
     
   return radar , output 

#===================================================
# DEALIASING BORDER FILTER
#===================================================

def DealiasingBorderFilter( radar , output , options )    :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='DealiasingBorderFilter'  

   if options['Dealiasing']['flag'] :

      #Use a simple edge detection routine to detect the borders between dealiased and non-dealiased regions.
      #Flag this grid points and eliminate all the non-dealiased corrected pixels near by.
      mask =  ( output['qcv'] == options['Dealiasing']['code'] ).astype(float)  
    
      [ edge_intensity , edge_mask ]=qc.simple_edge_filter( field=mask , nx=na,ny=nr,nz=ne,undef=output['undef_v'],
                                                                               nboxx=options[filter_name]['nx'],nboxy=options[filter_name]['ny'],
                                                                               nboxz=options[filter_name]['nz'],edge_tr=0.5 )

      #Find the pixels which are close to a dealiased region but which has not been corrected by the dealiasing.
      tmp_index = np.logical_and( edge_mask , mask == 0 ).astype(float) 

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# ECHO TOP FILTER  
#===================================================

def EchoTopFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='EchoTopFilter'
   if options['name_ref'] in radar.fields  :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']
                                   
      tmp_max_z=np.zeros((na,nr,ne))

      for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
         tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]
 
      [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

      tmp_index=tmp_data_3d[:,:,:,0]
        
      computed_etfilter = True  #In case we need any of the other variables computed in this routine.

      tmp_index[ tmp_max_z < options[filter_name]['heigthtr'] ] = 1.0e6  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)

      output = output_update( output , tmp_index , options , filter_name ) 
    
   return radar , output 

#===================================================
# ECHO DEPTH FILTER 
#===================================================

def EchoDepthFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='EchoDepthFilter'
   if options['name_ref'] in radar.fields  :

     na=output['na']
     nr=output['nr']
     ne=output['ne']
     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     tmp_max_z=np.zeros((na,nr,ne))

     for ii in range(0,ne)     :       #Estimate the maximum radar data height assoicated with each gate.
        tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]

     [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

     tmp_index=tmp_data_3d[:,:,:,2]


     tmp_index[ tmp_max_z < options[filter_name]['heigthtr'] ] = 1.0e6  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)

     output = output_update( output , tmp_index , options , filter_name ) 


   return radar , output 

#===================================================
# LOW ELEVATION ANGLE REFLECTIVITY FILTER
#===================================================

def LowElevFilter( radar , output , options)  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   #Remove echos which are present at low elevation angles but not at higher elevation 
   #angles. This will help to remove clutter produced by anomalous propagation. 
   #This can also help to remove second trip echoes which are usually observed only at low elevation
   #angles.
   #This can also eliminate distant convective cells that are only seen in low 
   #elevation angles near the radar range edge. However these cells are of little interest
   #for data assimilation since that information is usually not enough to adequatelly constrain
   #the evolution of the convective cells.

   filter_name='LowElevFilter'

   if  options['name_ref'] in radar.fields :

     na=output['na']
     nr=output['nr']
     ne=output['ne']
     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     #Get the angles that will be used based on the selected threshold.
     tmp_angles= output['elevations'][ output['elevations'] < options[filter_name]['min_angle']]
     tmp_n_angles = np.size( tmp_angles )

     output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
     
     tmp_index=np.zeros([na,nr,ne]) 
     for ie in range( 0 , tmp_n_angles )  :
        tmp_index[:,:,ie]=np.logical_and( output['ref'][:,:,ie] > options['norainrefval'] , output['smooth_ref'][:,:,tmp_n_angles] <= options['norainrefval'] )
        tmp_index[:,:,ie][ output['altitude'][:,:,tmp_n_angles] > options[filter_name]['height_thr'] ] = 0.0

     output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# RHO HV FILTER
#===================================================

def RhoFilter( radar , output , options )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'RhoFilter'

   na=output['na']
   nr=output['nr']
   ne=output['ne']
   nx=options[filter_name]['nx']
   ny=options[filter_name]['ny']
   nz=options[filter_name]['nz']

   if options['name_rho'] in radar.fields  :

      output['undef_rho']=radar.fields[ options['name_rho'] ]['_FillValue']

      [ rhohv , dm , dm , dm , dm , dm  ]=order_variable( radar , options['name_rho'] , output['undef_rho']  )

      #Compute the filter parameter
      tmp_index=qc.box_functions_2d(datain=rhohv,na=na,nr=nr,ne=ne,undef=output['undef_rho']
                                              ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

      output = output_update( output , tmp_index , options , filter_name ) 
      

   return radar , output
 
#===================================================
# REFLECTIVITY SPECKLE FILTER
#===================================================

def RefSpeckleFilter( radar , output , options )    :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='RefSpeckleFilter'

   if options['name_ref'] in radar.fields :

       #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
       na=output['na']
       nr=output['nr']
       ne=output['ne']

       nx=options[filter_name]['nx']
       ny=options[filter_name]['ny']
       nz=options[filter_name]['nz']
       tr=options[filter_name]['reftr']

       tmp_index=qc.box_functions_2d(datain=output['ref'].data,na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr) 

       output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# DOPPLER SPECKLE FILTER
#===================================================

def DopplerSpeckleFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='DopplerSpeckleFilter'

   if  options['name_v'] in radar.fields :

       na=output['na']
       nr=output['nr']
       ne=output['ne']
       nx=options[filter_name]['nx']
       ny=options[filter_name]['ny']
       nz=options[filter_name]['nz']
       tr=options[filter_name]['dvtr']

       tmp=np.abs(output['v'].data)

       tmp[ output['v'] == output['undef_v'] ] = output['undef_v']

       tmp_index=qc.box_functions_2d(datain=tmp,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr)

       output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# ATTENUATION FILTER
#===================================================

def AttenuationFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='AttenuationFilter'

   if options['name_ref'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']

      beaml=radar.range['data'][1]-radar.range['data'][0] #Get beam length

      tmp_index=qc.get_attenuation( var=output['cref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,beaml=beaml,cal_error=options[filter_name]['attcalerror']
                                                ,is_power=options[filter_name]['is_power']
                                                ,coefs=options[filter_name]['att_coefs'] )
      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# TOPOGRAPHY BLOCKING FILTER
#===================================================

def BlockingFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='BlockingFilter'  #This is not included in the fuzzy logic algorithm

   tmp_index=qc.compute_blocking( radarz=output['altitude'] , topo=output['topo'] , na=na , nr=nr , ne=ne      ,
                                              undef=output['undef_ref'],
                                              radar_beam_width_v=radar.instrument_parameters['radar_beam_width_v']['data'] , 
                                              beam_length=radar.range['meters_between_gates']                              , 
                                              radarrange=radar.range['data'] , radarelev=output['elevations'] )  

   #Compute correction 
   if options[filter_name]['blocking_correction'] & ( options['name_ref'] in radar.fields )  :
      #Correct partially blocked precipitation echoes.
      mask=np.logical_and( tmp_index > 0.1 , tmp_index <= 0.3 )
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 1.0

      mask=np.logical_and( tmp_index > 0.3 , tmp_index <= 0.4 )
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 2.0

      mask= tmp_index > 0.4 
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 3.0

   output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# DOPPLER NOISE FILTER
#===================================================

def DopplerNoiseFilter( radar , output , options )   :

   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='DopplerNoiseFilter'

   if options[filter_name]['flag'] & ( options['name_v'] in radar.fields ) :

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     nx2=options[filter_name]['nx2']
     ny2=options[filter_name]['ny2']
     nz2=options[filter_name]['nz2']

     tr_1=options[filter_name]['threshold_1']
     tr_2=options[filter_name]['threshold_2']

     tmp_dv_1=np.copy(output['v'])
     tmp_dv_2=np.copy(output['v'])

     for ip in range(0,options[filter_name]['n_filter_pass']) :

       output['smooth_v']=qc.box_functions_2d(datain=tmp_dv_2,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

       undef_mask = np.logical_or( output['smooth_v'] == output['undef_v'] , output['v'] == output['undef_v'] )

       output['distance_1'] = np.abs( output['smooth_v'] - output['v'] )

       output['distance_1'][undef_mask] = 0.0


       #output['distance_1']=qc.compute_distance(tmp_dv_1,tmp_dv_2,nx=na,ny=nr,nz=ne,
       #                                         undef=output['undef_v'],nx_box=nx,ny_box=ny,nz_box=nz)

       #tmp_dv_2=np.copy(output['v']) 

       tmp_dv_2[ output['distance_1'] > tr_1 ]=output['undef_v']

       print( np.max( output['distance_1'] ), np.min( output['distance_1'] ) )

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['distance_1'], nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx_1']
                                           , yy=options[filter_name]['ify_1'] , nxx=np.size(options[filter_name]['ifx_1']) )


     if not options[filter_name]['force']   :
       output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
       output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
     else                                   :
       output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     #tmp_dv_1[ np.logical_or( output['distance_1'] > tr_1 , output['distance_1']==output['undef_v'] ) ]=output['undef_v']

     #tmp_dv_2=np.copy(tmp_dv_1)

     for ip in range(0,options[filter_name]['n_filter_pass']) :

        output['smooth_v']=qc.box_functions_2d(datain=tmp_dv_2,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                              ,boxx=nx2,boxy=ny2,boxz=nz2,operation='MEAN',threshold=0.0)

        undef_mask = np.logical_or( output['smooth_v'] == output['undef_v'] , output['v'] == output['undef_v'] )

        output['distance_2'] = np.abs( output['smooth_v'] - output['v'] )

        output['distance_2'][undef_mask] = 0.0

        tmp_dv_2[ output['distance_2'] > tr_2 ]=output['undef_v']


     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['distance_2'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx_2']
                                           , yy=options[filter_name]['ify_2'] , nxx=np.size(options[filter_name]['ifx_2']) )

     if not options[filter_name]['force']   :
       output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
       output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
       output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
       if options[filter_name]['fill_value']  == 'undef'   :
          output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       else                                                :
          output['cv'][ tmp_w > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']

       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

     if  not options[filter_name]['save']  :
          output.pop('distance_1')
          output.pop('distance_2')

   return radar , output 

#===================================================
#  DOPPLER LOCAL STD FILTER
#===================================================

def DopplerLocalStdFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'DopplerLocalStdFilter'

   if options['name_v'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']

      #Compute the filter parameter
      tmp_index=qc.box_functions_2d(datain=output['v'],na=na,nr=nr,ne=ne,undef=output['undef_v']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='SIGM',threshold=0.0)

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# INTERFERENCE FILTER
#===================================================

def InterferenceFilter( radar , output , options )   :

   import numpy as np
   #Clean rays affected by interference.
   filter_name = 'InterferenceFilter'

   if options['name_ref'] in radar.fields :

      tmp_ref=np.copy( output['ref'] )

      #This function returns 1 if the pixel is sopossed to be contaminated by interference and
      # 0 otherwise.
      tmp_index = interference_filter ( tmp_ref , output['undef_ref'] , options['norainrefval'] 
                                            , radar.range['data'] , options[filter_name] ) 

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# DOPPLER SPATIAL COHERENCE FILTER
#===================================================

def DopplerSpatialCoherenceFilter( radar , output , options )   :

   import numpy as np
   #Clean pixels which shows no coherence with their neighbors

   filter_name = 'DopplerSpatialCoherenceFilter'

   if options['name_v'] in radar.fields  :

      tmp_v=np.copy( output['v'] )

      tmp_index = dopplerspatialcoherence_filter( tmp_v , output['undef_v'] , options[filter_name] )

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# POWER FILTER
#===================================================

def PowerFilter( radar , output , options )   :

   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name = 'PowerFilter' 

   if  options['name_v'] in radar.fields  :

      print('This filter has not been coded yet')

      output = output_update( output , tmp_index , options , filter_name )

   #TODO

   return radar , output 

#===================================================
# DETECT MISSING REFLECTIVITY VALUES
#===================================================

def MissingRefFilter( radar , output , options )    :

   import numpy as np

   filter_name='MissingRefFilter'

   if options[filter_name]['flag'] & ( options['name_ref'] in radar.fields ) :

      options['missing_mask'] = qc.detect_missing(  output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                   ,min_ref=options['norainrefval'],threshold=options[filter_name]['threshold'] 
                                                   ,nmissing_max=options[filter_name]['nmissing_max'] )

      tmp_w = options['missing_mask'].astype(int)

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
         output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         if options[filter_name]['fill_value']  == 'undef'       :
             output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
         elif options[filter_name]['fill_value']  == 'min_ref'   :
            output['cref'][ tmp_w > options[filter_name]['force_value'] ]=options['norainrefval']
         else                                                    :
            output['cref'][ tmp_w > options[filter_name]['force_value'] ]=options['filter_name']['fill_value']

         output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      if  not options[filter_name]['save']  :
           output.pop('missing_mask')

   return radar , output 

#===================================================
# REFLECTIIVTY TEXTURE FILTER
#===================================================

def ReflectivityTextureFilter( radar , output , options )   :

   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='ReflectivityTextureFilter'

   if  options[filter_name]['flag'] & ( options['name_ref'] in radar.fields ) : 

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']
     output['texture_ref']=qc.compute_texture(var=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)


     if options[filter_name]['use_smooth_ref'] :
        output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],
                                                 boxx=0,boxy=0,boxz=0,operation='MEAN',threshold=0.0)

        #High reflectivity cores will not be affected by texture filter.
        output['texture_ref'][ output['smooth_ref'] >= options[filter_name]['smooth_ref_tr']  ]= 0.0

     output['texture_ref'][output['ref']==output['undef_ref']] = 0.0

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['texture_ref'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )


     if not options[filter_name]['force']   :
        output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
        output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
        output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
        if options[filter_name]['fill_value']  == 'undef'       :
           output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
        elif options[filter_name]['fill_value']  == 'min_ref'   :
           output['cref'][ tmp_w > options[filter_name]['force_value'] ]=options['norainrefval']
        else                                                    :
           output['cref'][ tmp_w > options[filter_name]['force_value'] ]=options['filter_name']['fill_value']

        output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

     if  not options[filter_name]['save']  :
          output.pop('texture_v')

   return radar , output 

#===================================================
# DOPPER VELOCITY TEXTURE FILTER
#===================================================

def DopplerTextureFilter( radar , output , options )   :

   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='DopplerTextureFilter'

   if  options[filter_name]['flag'] & ( options['name_ref'] in radar.fields ) :

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']
     output['texture_v']=qc.compute_texture(var=output['v'],na=na,nr=nr,ne=ne,undef=output['undef_v'],nx=nx,ny=ny,nz=nz)

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['texture_v'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

     if not options[filter_name]['force']   :
        output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
        output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
        output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
     else                                   :
       if options[filter_name]['fill_value']  == 'undef'   :
          output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       else                                                :
          output['cv'][ tmp_w > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']
                               
       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     if  not options[filter_name]['save']  :
          output.pop('texture_v')

   return radar , output

#===================================================
# LOW DOPPLER VOLOCITY FILTER
#===================================================

def LowDopplerFilter( radar , output , options )   :
 
   #This filter removes reflectivity values which are close to the surface and 
   #have a low associated doppler velocity. 
   import numpy as np

   filter_name='LowDopplerFilter'

   if  options['name_v'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']

      tmp_index=np.abs(output['v'])
      tmp_index[ (output['altitude'] - output['topo']) > options[filter_name]['height_thr'] ]=10.0
    
      output = output_update( output , tmp_index , options , filter_name )
 
   return radar , output 


#===========================================================================================================
# OTRAS FUNCIONES CONTENIDAS EN ESTE MODULO
#===========================================================================================================   

#From Pyart order to ARE (Azmimuth , range , elevation )

def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   import numpy.ma as ma
   import warnings 

   #order_var es la variable ordenada con los azimuths entre 0 y 360 (si hay rayos repetidos se promedian).
   #order_azimuth es el azimuth "aproximado" utilizando 0 como azimuth inicial y avanzando en intervalos regulares e iguales a la resolucion
   #del azimuth en grados.
   #levels son los angulos de elevacion.
   #azimuth_exact es un array que contiene para cada nivel el valor exacto del azimuth que corresponde a cada rayo. 

   ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos ')
   ray_angle_res=np.nanmean( ray_angle_res )
   #print ('The resolution in azimuth is: %5.3f' % ( ray_angle_res ) )


   levels=np.unique(radar.elevation['data'])
   azimuth=radar.azimuth['data']
   time=radar.time['data']

   order_azimuth=np.arange(0.0,360.0-ray_angle_res,ray_angle_res) #Asuming quasi regular azimuth location
   naz=np.size(order_azimuth)
   nel=np.size(levels)

   if ( var_name == 'altitude' ) :
      var=radar.gate_altitude['data']
   elif( var_name == 'longitude' ) :
      var=radar.gate_longitude['data'] 
   elif( var_name == 'latitude'  ) :
      var=radar.gate_latitude['data'] 
   elif( var_name == 'x' )         :
      var=radar.gate_x['data']
   elif( var_name == 'y' )         : 
      var=radar.gate_y['data']
   else  :
      var=radar.fields[var_name]['data'].data


      var[ var == undef ] = np.nan

   nr=var.shape[1]

   #Allocate arrays
   order_var    = undef + np.zeros((naz,nr,nel))
   order_time   =np.zeros((naz,nel)) 
   order_index  =undef + np.zeros((naz,nel))   #This variable can be used to convert back to the azimuth - range array
   azimuth_exact=undef + np.zeros((naz,nel))

   order_var[:]     = undef 
   order_time[:]    = undef 
   azimuth_exact[:] = undef

   

   for ilev in range(0, nel) :

      levmask= radar.elevation['data'] == levels[ilev] 

      min_index = np.min( np.where( levmask ) )


      #Find the azimuths corresponding to the current elevation.
      azlev=azimuth[ levmask ]
      timelev=time[ levmask ]
      #Get variabile values corresponding to the current elevation
      varlev=var[ levmask , : ]

      #For the first azimuth which is a special case because it contains zero.
      az_index=np.logical_or( azlev <= ray_angle_res/2.0 , azlev >= 360 - ray_angle_res/2.0 )
     
      if ( np.sum(az_index) > 0 ) : 
         with warnings.catch_warnings():
              #Run time warnings resulting from all nan array in nanmean are expected
              #and supressed in this block.
              warnings.simplefilter("ignore", category=RuntimeWarning)

              order_var[0,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
              order_time[0,ilev] = np.nanmean( timelev[ az_index ] )
              azimuth_exact[0,ilev] = np.nanmean( azlev[ az_index ] )
         order_index[0,ilev]   = np.where( az_index )[0][0] + min_index

      #Para los que vienen despues.
      for iaz in range(1,naz) :
         #Search for all the rays that are close to the current azimuth and level.
         az_index=np.logical_and( azlev <= order_azimuth[iaz] + ray_angle_res/2.0 , azlev >= order_azimuth[iaz] - ray_angle_res/2.0 )
         if( np.sum( az_index ) > 0 ) :
            with warnings.catch_warnings():
                 #Run time warnings resulting from all nan array in nanmean are expected
                 #and supressed in this block.
                 warnings.simplefilter("ignore", category=RuntimeWarning)
                 order_var[iaz,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
                 order_time[iaz,ilev] = np.nanmean( timelev[ az_index ] )
                 azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )
            order_index[iaz,ilev] = np.where(az_index)[0][0] + min_index #If multiple levels corresponds to a single azimuth / elevation chose the first one.

   order_var[ np.isnan( order_var ) ]= undef
   order_index[ np.isnan( order_index ) ]=undef

   return order_var , order_azimuth , levels , order_time , order_index , azimuth_exact

#From ARE order to Pyart order 
def order_variable_inv (  radar , var , order_index , undef )  :

    import numpy as np
   
    #Esta funcion es la inversa de la funcion order variable. Es decir que toma un array ordenado como azimuth , range y elevation y lo vuelve
    #a ordenar como azimuth-elevation y range. Es decir el orden original que se encuentra en los archivos con formato cfradial y que heredan los objetos radar de pyart.

    #var es la variable ordenada como var(azimuth,range,elevation)
    #order_index (azimuth,elevation) contiene la posicion original de los haces que fueron asignados a cada azimuth y elevacion por la funcion order_variable.
    #nr numero de puntos en la direccion del rango.
    #nb numero de beams. 

    na=var.shape[0]
    nr=var.shape[1]
    ne=var.shape[2]

    nb=radar.azimuth['data'].shape[0]

    output_var=np.ones((nb,nr)) * undef 
    
    for ia in range(0,na)  :

       for ie in range(0,ne)  :

          if ( not order_index[ia,ie] == undef )  :
     
              output_var[int(order_index[ia,ie]),:]=var[ia,:,ie] 

    return output_var


def get_rma_strat ( filename , radar )  :

    import numpy as np
    import numpy.ma as ma

    local_fill_value = -9999.0
    levels=np.unique(radar.elevation['data'])

    #Add missing structures to the radar object.
    if radar.altitude_agl == None :
       radar.altitude_agl = dict()
    if radar.metadata == None :
       radar.metadata = dict()
    if radar.instrument_parameters == None :
       radar.instrument_parameters = dict()
       radar.instrument_parameters['nyquist_velocity']=dict()
       radar.instrument_parameters['nyquist_velocity']['long_name']='unambiguous_doppler_velocity'
       radar.instrument_parameters['nyquist_velocity']['units']='meters per second'
       radar.instrument_parameters['nyquist_velocity']['_FillValue']= local_fill_value
       radar.instrument_parameters['nyquist_velocity']['meta_group']='instrument_parameters'
 
       radar.instrument_parameters['radar_beam_width_v']=dict()
       radar.instrument_parameters['radar_beam_width_v']['long_name']='half_power_radar_beam_width_v_channel'
       radar.instrument_parameters['radar_beam_width_v']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_v']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_v']['meta_group']='instrument_parameters'

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

    radar.metadata['instrument_name'] = 'RMA' + filename[ filename.find('RMA') + 3 ]

    if '0122_03' in filename  :  #122-3 STRATEGY
       nyquist_velocity     = 13.35
       ray_angle_res        = 1.0
       meters_between_gates = 300
       radar_beam_width_h   = 1.0
       radar_beam_width_v   = 1.0


    #Apply the missing parameters to the radar structure.
    radar.instrument_parameters['nyquist_velocity']['data'] = ma.array(np.ones( np.shape(radar.azimuth['data']) )*nyquist_velocity , mask = np.zeros( np.shape(radar.azimuth['data']) , dtype=bool ) , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_h']['data'] = ma.array( radar_beam_width_h , mask =False , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_v']['data'] = ma.array( radar_beam_width_v , mask =False , fill_value = local_fill_value )
    radar.ray_angle_res['data'] = ma.array( np.ones( np.shape( levels ) )*ray_angle_res , mask = np.zeros( np.shape( levels ) , dtype=bool ) , fill_value = local_fill_value )
    radar.range['meters_between_gates']= meters_between_gates
    radar.range['meters_to_center_of_first_gate']= meters_between_gates / 2.0


    return radar


def dopplerspatialcoherence_filter( v , undef , my_conf ) : 

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   from sklearn import linear_model, datasets
   import matplotlib.pyplot as plt

   #Compute a robust correlation between a ray and its neighbors.
   #Detect outliers based on the correlation. 

   na=v.shape[0]
   nr=v.shape[1]
   ne=v.shape[2]

   #Coherence mask starts flagging the undef values
   available_mask = ( v != undef )

   coherence_index = np.zeros( np.shape( v ) )

   ransac = linear_model.RANSACRegressor()

   if my_conf['compute_horizontal_coherence']   :
   #CONSIDER THE COHERENCE WITH THE NEIGHBOR BEAMS IN AZIMUTH DIRECTION

      for k in range(ne)  :

         for i in range(na)  :

              iprev = i - 1
              if i == 0    :
                   iprev = na-1


              #First get the number of undef values and generate a common undef mask.
              local_valid_mask = np.logical_and( available_mask[i,:,k] , available_mask[iprev,:,k] )

              n_valid = np.sum( local_valid_mask.astype(int) ) 

              coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ]=coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ] + 1.0

              if ( n_valid / nr ) > my_conf['threshold_undef']   :

                 tmp_vi = np.copy( v[i,:,k][local_valid_mask] )
                 tmp_viprev = np.copy( v[iprev,:,k][local_valid_mask] )

                 if my_conf['consistency_metric'] == 'ransac'   :
                    #Use a robust correlation method to decide which
                    #elements are spatially coherent.

                    ransac.fit( tmp_vi.reshape(-1, 1) , tmp_viprev.reshape(-1, 1) )
                    inlier_mask = ransac.inlier_mask_
                    outlier_mask = np.logical_not( inlier_mask )

                 if my_conf['consistency_metric'] == 'constant'  :
                     #Use a constant threshold to define spatial coherence.

                     inlier_mask = np.abs( tmp_vi - tmp_viprev ) < my_conf['constant_consistency_threshold'] 
                     outlier_mask = np.logical_not( inlier_mask )


                 corrcoef=np.corrcoef( tmp_vi[inlier_mask] ,  tmp_viprev[inlier_mask] )[0,1]

                 tmp=coherence_index[i,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,k][local_valid_mask]=tmp

                 tmp=coherence_index[iprev,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[iprev,:,k][local_valid_mask]=tmp

                 if corrcoef < my_conf['threshold_corr']   :

                    coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0
                    coherence_index[iprev,:,k] = coherence_index[iprev,:,k] + 3.0
              else  :

                 coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0

   if my_conf['compute_vertical_coherence']   :
   #CONSIDER THE COHERENCE WITH THE NEIGHBOR BEAMS IN VERTICAL DIRECTION

      for k in range(1,ne)  :
 
         kprev = k - 1

         for i in range(na)  :


              #First get the number of undef values and generate a common undef mask.
              local_valid_mask = np.logical_and( available_mask[i,:,k] , available_mask[i,:,kprev] )

              n_valid = np.sum( local_valid_mask.astype(int) )

              coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ]=coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ] + 1.0

              if ( n_valid / nr ) > my_conf['threshold_undef']   :

                 tmp_vi = np.copy( v[i,:,k][local_valid_mask] )
                 tmp_viprev = np.copy( v[i,:,kprev][local_valid_mask] )

                 if my_conf['consistency_metric'] == 'ransac'   :
                    #Use a robust correlation method to decide which
                    #elements are spatially coherent.

                    ransac.fit( tmp_vi.reshape(-1, 1) , tmp_viprev.reshape(-1, 1) )
                    inlier_mask = ransac.inlier_mask_
                    outlier_mask = np.logical_not( inlier_mask )

                 if my_conf['consistency_metric'] == 'constant'  :
                     #Use a constant threshold to define spatial coherence.

                     inlier_mask = np.abs( tmp_vi - tmp_viprev ) < my_conf['constant_consistency_threshold']
                     outlier_mask = np.logical_not( inlier_mask )

                 #Compute the correlation between the adjacent rays for the inlier mask. 
                 corrcoef=np.corrcoef( tmp_vi[inlier_mask] ,  tmp_viprev[inlier_mask] )[0,1]

                 tmp=coherence_index[i,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,k][local_valid_mask]=tmp

                 tmp=coherence_index[i,:,kprev][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,kprev][local_valid_mask]=tmp

                 #coherence_index[i,:,k][local_valid_mask][outlier_mask] = coherence_index[i,:,k][local_valid_mask][outlier_mask] + 3.0 
                 #coherence_index[iprev,:,k][local_valid_mask][outlier_mask] = coherence_index[iprev,:,k][local_valid_mask][outlier_mask] + 3.0 

                 if corrcoef < my_conf['threshold_corr']   :

                    coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0
                    coherence_index[i,:,kprev] = coherence_index[iprev,:,kprev] + 3.0

   for ifilter in range( my_conf['npass_filter'] )  :

      v[ coherence_index > my_conf['threshold_coherence_index'] ] = undef

      for k in range(ne) :

         for i in range(na) :
            if  my_conf['azimuthfilter']         :  #DETECT ISOLATED PIXELS IN AZIMUTH

               if ( i > 1 ) & ( i < na-2 ) :

                  #If we have in only one ray but not in the neighbors this suggest an interference pattern.
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef  ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==na-1   :
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[0,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==na-2   :
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[i,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[0,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==0      :

                  tmp_mask = np.logical_and( v[na-1,:,k] == undef , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[na-2,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==1      :

                  tmp_mask = np.logical_and( v[i-1,:,k] == undef  , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[na-1,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

            if  my_conf['elevationfilter']       :   #DETECT ISOLATED PIXELS IN ELEVATION

               if ( k > 0 ) & ( k < ne-1 ) :

                  tmp_mask = np.logical_and( v[i,:,k-1] == undef , v[i,:,k+1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               if ( k == 0 )                :

                  tmp_mask = np.logical_and( v[i,:,k+2] == undef , v[i,:,k+1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               if ( k == ne-1 )            :

                  tmp_mask = np.logical_and( v[i,:,k-2] == undef  , v[i,:,k-1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

         for i in range(nr) :

            if  my_conf['rangefilter']           :   #DETECT ISOLATED PIXELS IN RANGE

               if ( i  > 0 ) & ( i < nr-1 ) :

                  tmp_mask = np.logical_and( v[:,i-1,k] == undef , v[:,i+1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

               if ( i == 0 )                :

                  tmp_mask = np.logical_and( v[:,i+2,k] == undef , v[:,i+1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

               if ( i == nr-1 )            :

                  tmp_mask = np.logical_and( v[:,i-2,k] == undef  , v[:,i-1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

   if my_conf['enable_speckle']   :

      tr=0.0
      nx=my_conf['nx']
      ny=my_conf['ny']
      nz=my_conf['nz']
      tmp=np.abs( v )
      tmp[ v == undef ] = undef

      speckle_v=qc.box_functions_2d(datain=tmp,na=na,nr=nr,ne=ne,undef=undef
                                 ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr)

      coherence_index[ speckle_v < my_conf['speckle_threshold'] ] = 10.0

   #v[ coherence_index > my_conf['threshold_coherence_index'] ] = undef

   return coherence_index

def interference_filter ( ref , undef , min_ref , r , my_conf )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   from sklearn import linear_model, datasets
   import matplotlib.pyplot as plt

   #ref a 3D array with the reflectivity in the na,nr,ne grid.
   #undef: value corresponding to the missing value in the reflectivity field.
   #min_ref: value corresponding to the no-rain value in the reflectivity field.
   #r : radar range (vector)
   #my_conf : configuration for this filter.

   offset= my_conf['offset']

   na=ref.shape[0]
   nr=ref.shape[1]
   ne=ref.shape[2]

   nx=my_conf['nx']
   ny=my_conf['ny']
   nz=my_conf['nz']

   if my_conf['Smooth_Ref']  :

      tmp_ref = qc.box_functions_2d(datain=ref,na=na,nr=nr,ne=ne,undef=undef
                                 ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
   else                      :
      tmp_ref = np.copy( ref )

   tmp_index=np.zeros( np.shape( ref ) ) 

   tmp_z = np.power( 10.0, tmp_ref / 10.0  )

   tmp_z[tmp_ref == undef ] = undef
  
   
   #offset=my_conf['offset']
   att=my_conf['att']
   AzimuthFilter=my_conf['AzimuthFilter']
   ElevationFilter=my_conf['ElevationFilter']
   npass_filter=my_conf['npass_filter']
   percent_valid_threshold=my_conf['percent_valid_threshold']
   corr_threshold=my_conf['corr_threshold']
   ref_threshold=my_conf['ref_threshold']
   percent_ref_threshold=my_conf['percent_ref_threshold']

   Power_Regression = my_conf['Power_Regression']

   #Main loops

   for k in range(ne)  :

      for i in range(na)  :

           local_ref = np.copy( tmp_ref[i,:,k] )
           local_ref[0:offset]=undef

           undef_mask = local_ref != undef 

           tmp_count = np.sum( (undef_mask).astype(int) )/nr

           power = local_ref  - 20.0 * np.log10( r ) - 2.0 * att * r

           if tmp_count > percent_valid_threshold   :
              ransac = linear_model.RANSACRegressor()     

              if Power_Regression   :
                 ransac.fit( r[undef_mask].reshape(-1, 1) , np.power( 10.0 , power[undef_mask].reshape(-1, 1)/10.0 ) )
                 inlier_mask = ransac.inlier_mask_
                 outlier_mask = np.logical_not( inlier_mask )

                 tmppower=ransac.predict( r.reshape(-1,1) )[:,0]
                 tmppower[ tmppower < 0.0 ]=10e-10

                 powerrayo = 10.0*np.log10( tmppower ) + 20.0 * np.log10( r ) + 2.0 * att * r

              else                  :
                 ransac.fit( r[undef_mask].reshape(-1, 1) , power[undef_mask].reshape(-1, 1) )
                 inlier_mask = ransac.inlier_mask_
                 outlier_mask = np.logical_not(inlier_mask)

                 powerrayo = ransac.predict( r.reshape(-1,1) )[:,0]  + 20.0 * np.log10( r ) + 2.0 * att * r

           else:

              powerrayo = np.zeros( nr ) + 20.0 * np.log10( r ) + 2.0 * att * r

           #print( np.shape( powerrayo ),np.shape( local_ref )

           if ( np.sum( inlier_mask.astype(int) ) >= 10 ) & ( np.std(local_ref[undef_mask][inlier_mask]) > 0 ) :
           
              corrcoef=np.corrcoef( powerrayo[undef_mask][inlier_mask],local_ref[undef_mask][inlier_mask] )[0,1]

           else                                      :

              corrcoef=-1.0

           #if (k == 0) and (i == 348) :
           #  plt.plot(  powerrayo[undef_mask] , local_ref[undef_mask] ,'or')
           #  plt.plot(  powerrayo[undef_mask][inlier_mask],local_ref[undef_mask][inlier_mask] , 'ok')
           #  plt.show()
           #  print(corrcoef,np.sum( inlier_mask.astype(int) )/nr )
           
           #  plt.plot( r[undef_mask] , powerrayo[undef_mask] )
           #  plt.plot( r[undef_mask][inlier_mask] , local_ref[undef_mask][inlier_mask] ,'ok' )
           #  plt.plot( r[undef_mask] , local_ref[undef_mask])
           #  plt.show()


           if( ( corrcoef > corr_threshold ) & ( np.sum( inlier_mask.astype(int) )/( nr-offset ) > percent_ref_threshold ) )  :

              
              #This means that this ray is likely to be contaminated by interference.

              undef_mask = ( tmp_z[i,:,k] != undef )
              zrayo = np.power( 10.0,powerrayo  / 10.0 ) 
              z     = np.power( 10.0,ref[i,:,k] / 10.0  )

              #If the reflectivity is far from the fitted interference, and is greather than the fitted
              #Interference, then correct the power substracting the interference power.         
              tmp_mask = np.logical_and( tmp_z[i,:,k] - zrayo > 5.0 ,  undef_mask  )  
              tmp_mask = np.logical_and( z - zrayo > 0 , tmp_mask )
              ref[i, tmp_mask ,k] = 10.0*np.log10( (z - zrayo)[tmp_mask] ) 
 
              #If the reflectivity is far from the fitted interference, and is smaller than the fitted interference
              #then set that pixel as an undef pixel.
              tmp_mask = np.logical_and( tmp_z[i,:,k] - zrayo <= 5.0 , undef_mask )
              ref[i, tmp_mask ,k] = undef
              tmp_index[i, tmp_mask ,k] = 1.0 

              #if (k == 0) and (i == 348) :
              #   plt.plot( r[undef_mask] , tmp_z[i,:,k][undef_mask] )
              #   plt.plot( r[undef_mask] , z[undef_mask]      ,'ok' )
              #   plt.plot( r[undef_mask] , zrayo[undef_mask])
              #   plt.show()

           #else :
              #print(i,k)

   #Additional filter for the remaining echoes
   #consider cyclic boundary conditions in azimuth.
   for ifilter in range( npass_filter )  :

      for k in range(ne) :

         for i in range(na) :
            if  AzimuthFilter         :  #DETECT ISOLATED PIXELS IN AZIMUTH

               if ( i > 1 ) & ( i < na-2 ) :
 
                  #If we have reflectivity in only one ray but not in the neighbors this suggest an interference pattern.
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i+1,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[i+2,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0
 
               elif  i==na-1   :
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef, ref[0,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[1,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==na-3   :
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[0,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==0      :

                  tmp_mask = np.logical_and( ref[na-1,:,k] == undef , ref[i+1,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[na-2,:,k] == undef , ref[i+2,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==1      :

                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i+1,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[na-1,:,k] == undef , ref[i+2,:,k] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0


            if ElevationFilter         :   #DETECT ISOLATED PIXELS IN ELEVATION

               if ( k > 0 ) & ( k < ne-1 ) :

                  tmp_mask = np.logical_and( ref[i,:,k-1] == undef , ref[i,:,k+1] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               if ( k == 0 )                :

                  tmp_mask = np.logical_and( ref[i,:,k+2] == undef , ref[i,:,k+1] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               if ( k == ne-1 )            :

                  tmp_mask = np.logical_and( ref[i,:,k-2] == undef , ref[i,:,k-1] == undef )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

   return tmp_index 


def smooth_fft( x )   :


    import numpy as np
    import scipy.fftpack

    N = np.size( x )

    w = scipy.fftpack.rfft( x )
    #f = scipy.fftpack.rfftfreq( N , 1.0 )
    spectrum = w**2

    cutoff_idx = spectrum > (spectrum.max()/10.0)
    #w2 = w.copy()
    w[cutoff_idx] = 0

    x_smooth = scipy.fftpack.irfft(w)


    return x_smooth


#Read topography data
def read_topo( my_file )           :
    import numpy as np

    f=open(my_file,'r')

    #First read the header containing nr and na
    tmpdata=np.fromfile(f,dtype='f4',count=2)
    na=int(tmpdata[0])
    nr=int(tmpdata[1])

    my_topo=dict()
    my_topo['mean']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['max']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['min']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['number']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['range']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['azimuth']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['latitude']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['longitude']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )


    return my_topo

#Write topography data
def write_topo( my_topo , my_file )           :
    import numpy as np

    f=open(my_file,'w')

    na=my_topo['mean'].shape[0]
    nr=my_topo['mean'].shape[1]

    #First write nr and na
    np.array(na).astype('f4').tofile(f)
    np.array(nr).astype('f4').tofile(f)

    #Write the components of the my_topo dictionary.
    np.reshape( my_topo['mean'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['max'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['min'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['number'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( (my_topo['range'].data).astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['azimuth'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['latitude'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['longitude'].astype('f4') , (nr*na) ).tofile(f)

    return my_topo

#Read a raster file an return lon,lat and data.
def read_raster_data(inputfile)  :
    import osgeo.gdal as gdal
    import numpy as np

    my_raster = gdal.Open(inputfile)

    nx=my_raster.RasterXSize
    ny=my_raster.RasterYSize

    nb=my_raster.RasterCount

    my_raster_data=np.zeros((nx,ny,nb))

    #Read the data and store it into a numpy array

    for ib in range( 0 , nb )  :

        my_raster_data[:,:,ib] = my_raster.ReadAsArray(ib)


    #Get the lat and lons.
    [lon,lat]=get_lat_lon(inputfile)

    return lon , lat , my_raster_data

#Read a raster file and get lat lon structure.
def get_lat_lon(inputfile)  :
    import osgeo.gdal as gdal
    import numpy as np

    my_raster = gdal.Open(inputfile)

    gt = my_raster.GetGeoTransform()
    proj = my_raster.GetProjection()

    xres = gt[1]
    yres = gt[5]

    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * my_raster.RasterXSize) - xres * 0.5
    ymin = gt[3] + (yres * my_raster.RasterYSize) + yres * 0.5
    ymax = gt[3] - yres * 0.5

    lon=np.zeros( my_raster.RasterXSize )
    lat=np.zeros( my_raster.RasterYSize )
    for ii in range( 0 , my_raster.RasterXSize )  :
       lon[ii] = xmin + ii * xres
    for ii in range( 0 , my_raster.RasterYSize )  :
       lat[ii]= ymax + ii*yres

    my_raster = None

    # create a grid of xy coordinates in the original projection
    [lon,lat] = np.meshgrid(lon,lat)

    return lon , lat

#Download topography data and generate binary topography file.
def generate_topo_file( rlon , rlat , rrange , razimuth , raster_path , topo_file )    :
    import elevation
    import numpy as np
    import os
    from common_functions import common_functions as cf

    product='SRTM1'  #SRTM1 - 30 m res , SRTM3 - 90 m res.
         
    my_topo=dict()

    [ my_topo['range'] ,  my_topo['azimuth'] ] = np.meshgrid( rrange , razimuth )
    nr=np.size( rrange )
    na=np.size( razimuth )

    #Get lat and lon corresponding to the grid points in the polar coordinate.

    [my_topo['longitude'],my_topo['latitude']]=cf.com_ra_to_ll(cen_lon=rlon,cen_lat=rlat
                                                              ,r=my_topo['range'],a=my_topo['azimuth'],nr=nr,na=na)

    #Define grid limits
    maxlon = int( np.ceil( np.max(my_topo['longitude']) )  )
    minlon = int( np.floor( np.min(my_topo['longitude']) ) )
    maxlat = int( np.ceil( np.max(my_topo['latitude']) )  )
    minlat = int( np.floor( np.min(my_topo['latitude']) ) )

    #Download the data to a raster file.
    #west, south, east, north 

    my_topo['mean'] = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )  
    my_topo['min']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )
    my_topo['max']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )
    my_topo['number']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.int32 ) 

    #Data is downloaded and processed into 1 deg patches.
    for ilon in range( minlon , maxlon )  :
      for ilat in range( minlat , maxlat )  : 

         print('Downloading data for tile Lon=',ilon,' Lat=',ilat) 

         #If raster file is not present, then download it from the internet.
         raster_file = raster_path + '/' + product + str(ilon) + '_' + str(ilat) + '.tif'

         if ( not os.path.isfile( raster_file ) )  :
             elevation.clip(bounds=(ilon,ilat,ilon+1,ilat+1),output=raster_file,product=product)

         #Read data from raster file.
         [raster_lon,raster_lat,raster_data]=read_raster_data(raster_file)
         raster_nx=raster_data.shape[0]
         raster_ny=raster_data.shape[1]

         #Convert raster lat lon to range and azimuth.
         [raster_r,raster_a]=cf.com_ll_to_ra(cen_lon=rlon,cen_lat=rlat
                                            ,lon=raster_lon,lat=raster_lat
                                            ,nx=raster_nx,ny=raster_ny)

         #Interpolate raster data to polar coordinates surrounding the radar.
         rmin=np.min( rrange )
         dr  =rrange[1] - rrange[0]
         amin=np.min( razimuth )
         da  =razimuth[1]  - razimuth[0]

         print('Interpolating patch Lon=',ilon,' Lat=',ilat) 
         cf.com_interp_boxavereg(xini=amin,dx=da,nx=na
                                ,yini=rmin,dy=dr,ny=nr
                                ,xin=raster_a.reshape( raster_nx * raster_ny )
                                ,yin=raster_r.reshape( raster_nx * raster_ny )
                                ,datain=raster_data.reshape( raster_nx * raster_ny )
                                ,nin=raster_nx*raster_ny
                                ,data_sum=my_topo['mean']
                                ,data_max=my_topo['max']
                                ,data_min=my_topo['min']
                                ,data_n=my_topo['number'],undef=-999)

          
    #Compute the mean topography.
    mask = my_topo['number'] > 0
    my_topo['mean'][ mask ] = my_topo['mean'][ mask ] / my_topo['number'][ mask ]

    #Complete missing values using neighborhood values.
    mask = my_topo['number'] == 0  #Identify missing data points
    my_topo['mean']=cf.com_complete_missing_2d(field=my_topo['mean'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)
    my_topo['max']=cf.com_complete_missing_2d(field=my_topo['max'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)
    my_topo['min']=cf.com_complete_missing_2d(field=my_topo['min'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)


    #Write topo file in polar coordintes into binary format.
    #This data will be used for fast interpolation to the corresponding volume elevation angles.
    write_topo( my_topo , topo_file )

    return my_topo








