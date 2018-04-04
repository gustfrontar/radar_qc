#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause


def main_qc( filename , options ) :

   import sys
   import time
   sys.path.append('../fortran')

   import numpy as np
   import numpy.ma as ma
   import matplotlib.pyplot as plt
   import pyart
   from common_qc_tools  import qc  #Fortran code routines.
   from common_qc_tools  import qc_const
   from scipy.interpolate import interp2d
   import common_gdem_functions as cgf
   import netCDF4
   import os

   import matplotlib.pyplot as plt

   qc_const.undef = options['undef']   #This set the undef value for the fortran routines.
   undef          = options['undef']   #This set the undef value for the rest of the script.

   output=dict() #Initialize output dictionary.

   computed_etfilter=False   #Flag to indicate if echo top has been computed already.

#  Constant parameters

   #Codigos que permitan identificar el control que actuo en cada pixel.

   #Para la reflectividad
   QCCODE_ATTENUATION = 10
   QCCODE_SPECKLE     = 11
   QCCODE_TEXTURE     = 12
   QCCODE_RHOFILTER   = 13
   QCCODE_SIGN        = 14
   QCCODE_BLOCKING    = 15
   QCCODE_ECHOTOP     = 16
   QCCODE_ECHODEPTH   = 17

   #Para la velocidad radial
   QCCODE_DEALIAS     = 30

   #El codigo de los datos buenos para reflectividad y velocidad radial.
   QCCODE_GOOD        = 0

   #Shortcut to variable names
   name_v=options['v_name']
   name_ref=options['ref_name']
   name_rho=options['rho_name']

   radar = pyart.io.read(filename)

   #Get the nyquist velocity
   nyquistv=radar.get_nyquist_vel(0,check_uniform=True)

   #===================================================
   # RESHAPE VARIABLES
   #===================================================
   #From time,range -> azimuth,range,elevation

   startt=time.time()


   if name_ref in radar.fields :
        start=time.time()

        [ output['ref'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact'] ]=order_variable( radar , name_ref , undef )
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cref'] = np.zeros(output['ref'].shape) 

        output['cref'][:] = output['ref']                 #Initialize the corrected reflectivity array.

        output['cref'][ output['cref'] == undef ]=options['norainrefval']
        output['ref'] [ output['ref']  == undef ]=options['norainrefval']

        output['qcref'] = np.zeros(output['cref'].shape)  #Set the qc flag array to 0.

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("ref -> az,r,el",end-start) )
 
   if name_v in radar.fields  : 

        start=time.time()

        [ output['v'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , name_v , undef  )
 
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne
 
        output['cv'] = np.zeros(output['v'].shape) 

        output['cv'][:] = output['v']                     #Initialize the corrected doppler velocity array

        output['qcv'] = np.zeros(output['v'].shape)       #Set the qc flag array to 0.

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("v -> az,r,el",end-start) )

   #===================================================
   # GEOREFERENCE RADAR DATA
   #===================================================

   #Use pyart rutines to compute x, y and z at each grid point
   #Reorder data
   start=time.time()

   #dm is a dummy variable

   [ output['altitude'] , dm , dm , dm , dm , dm ] = order_variable( radar , 'altitude' , undef )
   [ output['x']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'x' , undef ) 
   [ output['y']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'y' , undef )

   #Compute distance to radar for the first azimuth (we assume that distance to radar will be the same.
   #for all other azimuths.
   output['distance']=np.power( np.power(output['x'][0,:,:],2)+np.power(output['y'][0,:,:],2) , 0.5 )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("x,y,z -> az,r,el",end-start) )

   #===================================================
   # READ AND INTERPOLATE TOPOGRAPHY DATA
   #===================================================

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

   output['topo'] = np.zeros((na,nr,ne))

   #Loop over the vertical levels.
   for ie in range( 0 , ne )   :

      ppi_r = output['distance'][:,ie]
      ppi_a = output['az']

      output['topo'][:,:,ie] =  interpolator( ppi_r, ppi_a )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("topography interpolation",end-start) )

   #===================================================
   # PRE-DEALIASING FILTER FOR RV
   #===================================================
 
   #TODO

   #===================================================
   # DEALIASING 
   #===================================================

   if ( options['ifdealias'] and (name_v in radar.fields) ) : 
     
     start=time.time()

     #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
     winddealias=pyart.correct.region_dealias.dealias_region_based(radar, interval_splits=options['interval_splits'],interval_limits=None, 
                 skip_between_rays=options['skip_between_rays'], skip_along_ray=options['skip_along_ray'],centered=True,nyquist_vel=None,
                 check_nyquist_uniform=True,gatefilter=False,rays_wrap_around=True,keep_original=False,set_limits=True,
                 vel_field=name_v,corr_vel_field=None)

     #Add wind dealias to the radar objetc.
     radar.fields['Vda']                  = winddealias
     radar.fields['Vda']['coordinates']   = radar.fields[name_v]['coordinates']
     radar.fields['Vda']['units']         = radar.fields[name_v]['units']
     radar.fields['Vda']['long_name']     = radar.fields[name_v]['long_name']
     radar.fields['Vda']['standard_name'] = radar.fields[name_v]['standard_name']

     tmp_v=np.copy( output['cv'] )
     #Re-order dealiased wind data.
     [ output['cv'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , 'Vda' , undef  )

     output['qcv'][ output['cv'] != tmp_v ]=QCCODE_DEALIAS

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("dealiasing",end-start) )

   #===================================================
   # DEALIASING BORDER FILTER
   #===================================================

   #TODO


   #===================================================
   # ECHO TOP FILTER  
   #===================================================


   if ( options['ifetfilter'] ) :
  
     start=time.time()

     if ( not computed_etfilter )     :

        tmp_z=np.zeros((nr,ne))
        tmp_d=np.zeros((nr,ne))

        tmp_z[:]=output['altitude'][0,:,:]    #Store only one RHI section of the altitutde (we will assume that the altitude of a pixel is independent
                                              #of the azimuth).
        tmp_d[:]=output['distance']
                                   
        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]
 
        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['cref'],heigth=tmp_z,rrange=tmp_d,na=na,nr=nr,ne=ne
                                                ,nx=options['etfilternx'],ny=options['etfilterny'],nz=options['etfilternz'])
        computed_etfilter = True  #In case we need any of the other variables computed in this routine.

        tmp_ref=0 #Unset tmp_ref to free some memory.
        

     #Mask the pixels with echo top values under the threshold.
     #Do not mask tose pixels where the echo top threshold in below the maximum radar height (i.e. the heigth of the maximum elevation)
     output['cref'][ np.logical_and( tmp_max_z > options['etfiltertr'] , tmp_data_3d[:,:,:,0] < options['etfiltertr'] ) ]   = undef
     output['qcref'][ np.logical_and( tmp_max_z > options['etfiltertr'] , tmp_data_3d[:,:,:,0] < options['etfiltertr'] )  ] = QCCODE_ECHOTOP

     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( options['etfilter_save'] )     : 
        output['echo_top'] = tmp_data_3d[:,:,:,0] 


     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("echo-top filter",end-start) )

   #===================================================
   # RHO HV FILTER
   #===================================================

   if options['ifrhofilter']   :

      start=time.time()

      nx=options['rhofilternx']
      ny=options['rhofilterny']
      nz=options['rhofilternz']

      if name_rho in radar.fields  :

         [ output['rho'] , dm , dm , dm , dm , dm  ]=order_variable( radar , name_rho , undef  )

         output['rho_smooth']=qc.box_functions_2d(datain=output['rho'],na=na,nr=nr,ne=ne,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

         output['cref'][ output['rho_smooth'] < options['rhofiltertr'] ] = undef
         output['qcref'][output['rho_smooth'] < options['rhofiltertr'] ] = QCCODE_RHOFILTER

      else   :
         display('Warning: could not perform RHO-HV filter because rho was not found on this file')

      if [ not options['rhofilter_save'] ] :
          output.pop('rho_smooth')
          output.pop('rho')

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format("rho filter",end-start) )

   #===================================================
   # ECHO DEPTH FILTER 
   #===================================================

   if ( options['ifedfilter']  ) :

     start=time.time()

     if ( not computed_etfilter )     :

        tmp_z=np.zeros((nr,ne))
        tmp_d=np.zeros((nr,ne))

        tmp_z[:]=output['altitude'][0,:,:]    #Store only one RHI section of the altitutde (we will assume that the altitude of a pixel is independent
                                              #of the azimuth).
        tmp_d[:]=output['distance']
        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]

        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['cref'],heigth=tmp_z,rrange=tmp_d,na=na,nr=nr,ne=ne
                                                ,nx=options['edfilternx'],ny=options['edfilterny'],nz=options['edfilternz'])
        computed_etfilter = True  #In case we need any of the other variables computed in this routine.

     #Mask the pixels with echo depth values under the threshold.
     #Do not mask tose pixels where the echo depth threshold in below the maximum radar height (i.e. the heigth of the maximum elevation)
     output['cref'] [ np.logical_and( tmp_max_z > options['edfiltertr'] , tmp_data_3d[:,:,:,2] < options['edfiltertr'] ) ] = undef
     output['qcref'][ np.logical_and( tmp_max_z > options['edfiltertr'] , tmp_data_3d[:,:,:,2] < options['edfiltertr'] ) ] = QCCODE_ECHODEPTH

     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( options['edfilter_save'] )     :
        output['echo_depth'] = tmp_data_3d[:,:,:,2] 

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format("echo-depth filter",end-start) )

 
   #===================================================
   # SPECKLE FILTER
   #===================================================

   if options['ifspfilter']  :

       start=time.time()

       if options['spfilter_ref']  :

          #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
          nx=options['spfilternx']
          ny=options['spfilterny']
          nz=options['spfilternz']
          tr=options['spfilterreftr']
          output['speckle_ref']=qc.box_functions_2d(datain=output['cref'].data,na=na,nr=nr,ne=ne,boxx=nx,boxy=ny,boxz=nz,operation='COUN',threshold=tr) 

          #Set the pixels with values below the threshold as undef. 
          output['cref'][ output['speckle_ref'] < options['spfiltertr']  ] = undef
       
          output['qcref'][ output['speckle_ref'] < options['spfiltertr'] ] = QCCODE_SPECKLE

          #If the field is not included in the output then set it to 0.
          if [ not options['spfilter_save'] ] : 
             output.pop('speckle_ref')   


       if options['spfilter_v']  :

          #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
          nx=options['spfilternx']
          ny=options['spfilterny']
          nz=options['spfilternz']
          tr=options['spfilterreftr']
          output['speckle_v']=qc.box_functions_2d(datain=output['cv'].data,na=na,nr=nr,ne=ne,boxx=nx,boxy=ny,boxz=nz,operation='COUN',threshold=tr)

          #Set the pixels with values below the threshold as undef. 
          output['cv'][ output['speckle_v'] < options['spfiltertr']  ] = undef
       
          output['qcv'][ output['speckle_v'] < options['spfiltertr'] ] = QCCODE_SPECKLE

          #If the field is not included in the output then set it to 0.
          if [ not options['spfilter_save'] ] :
             output.pop('speckle_v')


       end=time.time()

       print("The elapsed time in {:s} is {:2f}".format("speckle filter",end-start) )

   #===================================================
   # ATTENUATION FILTER
   #===================================================

   if options['ifattfilter']   :
      
      start=time.time()

      beaml=radar.range['data'][1]-radar.range['data'][0] #Get beam length

      output['attenuation']=qc.get_attenuation(var=output['cref'],na=na,nr=nr,ne=ne,beaml=beaml,cal_error=options['attcalerror'])

      #Set the pixels with values below the threshold as undef. 
      output['cref'][ output['attenuation'] < options['attfiltertr']  ] = undef

      output['qcref'][ output['attenuation'] < options['attfiltertr'] ] = QCCODE_ATTENUATION

      if  not options['attfilter_save']  :
          output.pop('attenuation')


      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format("attenuation filter",end-start) )

   #===================================================
   # TOPOGRAPHY BLOCKING FILTER
   #===================================================

   if options['ifblfilter']   :

      start=time.time()

      output['blocking']=qc.compute_blocking( radarz=output['altitude'] , topo=output['topo'] , na=na , nr=nr , ne=ne      ,
                                              radar_beam_width_v=radar.instrument_parameters['radar_beam_width_v']['data'] , 
                                              beam_length=radar.range['meters_between_gates']                              , 
                                              radarrange=radar.range['data'] , radarelev=np.unique(radar.elevation['data'])  ) 


      #Compute correction 
      if options['blocking_correction']  :
         #Correct partially blocked precipitation echoes.
         mask=np.logical_and( output['blocking'] > 0.1 , output['blocking'] <= 0.3 )
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 1.0

         mask=np.logical_and( output['blocking'] > 0.3 , output['blocking'] <= 0.4 )
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 2.0

         mask= output['blocking'] > 0.4 
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 3.0

      #Set the pixels with values below the threshold as undef. 
      output['cref'][ output['blocking'] > options['blocking_threshold']  ] = undef
      output['cv'][ output['blocking']   > options['blocking_threshold']  ] = undef

      output['qcref'][ output['blocking'] > options['blocking_threshold'] ] = QCCODE_BLOCKING

      if  not options['blocking_save']  :
          output.pop('blocking')  

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format("blocking filter",end-start) )


   #===================================================
   # DOPPLER NOISE FILTER
   #===================================================

   #TODO

   #===================================================
   # INTERFERENCE FILTER
   #===================================================

   #TODO

   #===================================================
   # LOW DOPPLER VOLOCITY FILTER
   #===================================================

   #TODO

   #===================================================
   # ADD CORRECTED DATA TO RADAR OBJECT
   #===================================================

   #TODO

   #===================================================
   # END
   #===================================================


   endt=time.time()

   print("The elapsed time in {:s} is {:2f}".format("the entire QC",endt-startt) )

   return radar , output

#===========================================================================================================
# OTRAS FUNCIONES CONTENIDAS EN ESTE MODULO
#===========================================================================================================   

#From Pyart order to ARE (Azmimuth , range , elevation )

def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   import numpy.ma as ma

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
      #Find the azimuths corresponding to the current elevation.
      azlev=azimuth[ levmask ]
      timelev=time[ levmask ]
      #Get variabile values corresponding to the current elevation
      varlev=var[ levmask , : ]

      #For the first azimuth which is a special case because it contains zero.
      az_index=np.logical_or( azlev <= ray_angle_res/2.0 , azlev >= 360 - ray_angle_res/2.0 )
     
      if ( np.sum(az_index) > 0 ) : 
         order_var[0,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
         order_time[0,ilev] = np.nanmean( timelev[ az_index ] )
         azimuth_exact[0,ilev] = np.nanmean( azlev[ az_index ] )

      #Para los que vienen despues.
      for iaz in range(1,naz) :
         #Search for all the rays that are close to the current azimuth and level.
         az_index=np.logical_and( azlev <= order_azimuth[iaz] + ray_angle_res/2.0 , azlev >= order_azimuth[iaz] - ray_angle_res/2.0 )
         if( np.sum( az_index ) > 0 ) :
            order_var[iaz,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
            order_time[iaz,ilev] = np.nanmean( timelev[ az_index ] )
            azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )
            azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )
            order_index[iaz,ilev] = az_index[0]  #If multiple levels corresponds to a single azimuth / elevation chose the first one.

   order_var[ np.isnan( order_var ) ]= undef
   order_index[ np.isnan( order_index ) ]=undef

   return order_var , order_azimuth , levels , order_time , order_index , azimuth_exact

#From ARE order to Pyart order 
def order_variable_inv (  radar , var , order_index  )  :
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

    output_var=np.ones((na,nb)) * options['undef']
       

    for ia in range(0,na)  :

       for ie in range(0,ne)  :

          if ( not order_index[ia,ie] == undef )  :
     
             output_var[ia,order_index[ia,ie]]=var[ia,:,ie] 

    return output_var


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
    np.reshape( my_topo['range'].astype('f4') , (nr*na) ).tofile(f)
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




