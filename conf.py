#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.

#======================================
# ENVIRONMENTAL VARIABLES 
#======================================

import os
import numpy as np

os.environ['OMP_STACKSIZE'] = "1G" # visible in this process + all children

#======================================
# QC PARAMETERS
#======================================
options = {}

#General

options['ref_name']='dBZ'            #Reflectivity
options['cref_name']='CdBZ'          #Corrected reflectivity (qc output)
options['v_name']='V'                #Dopper velocity
options['cv_name']='CV'              #Corrected wind (qc ouput)
options['rho_name']='RhoHV'          #Rho HV

options['norainrefval']=-0.1
options['undef']=-9.99e9

#Fuzzy logic parameter

options['w_tr']=0.5                  #When total normalized weight is greather than this value
                                     #the pixel is flagged as missing data.
#Each filter has associated an ify and ifx array that descrives the importance function.
#Each filter has also a w which determines its relative importance among the other filters.


#The weight is computed for each filter based on the importance function evaluated at the filter value.
#The relative importance of each filter is given by the importance_w.

#Dealiasing parameters (pyart)
options['ifdealias']=False

options['da_interval_split']=3
options['da_skip_between_ray']=10
options['da_skip_along_ray']=10

#Rho filter parameters  ===================================================================== 

filter_name='RhoFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter

options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([-2,0.5,0.8,2])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=10
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.7                     #Threshold for force


#Echo top filter parameters ===================================================================

filter_name='EchoTopFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter

options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,2000])    #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=11
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
                                                            #is lower than this threshold
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=3000                    #Threshold for force

#Echo depth filter parameters ===================================================================

filter_name='EchoDepthFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter

options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,2000])    #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=12
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=3000                    #Threshold for force

#Reflectivity speckle filter  parameters ==========================================================

filter_name='RefSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                           #Enable / disable filter

options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=13
options[filter_name]['reftr']=5                             #Reflectivity threshold
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force


#Doppler speckle filter  parameters ==============================================================

filter_name='DopplerSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                           #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=14
options[filter_name]['dvtr']=0.2                            #Reflectivity threshold
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force

#Attenuation parameters           ==============================================================
filter_name='Attenuation'
options[filter_name]=dict()
options[filter_name]['flag']=False    #Attenuation filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,15,25,1])           #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=15
options[filter_name]['attcalerror']=1.0                     #Calibration error
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=20                      #Threshold for force

#Blocking parameters               ==============================================================
#This filter is not part of the Fuzzy logic algorithm.
filter_name='BlFilter'
options[filter_name]=dict()
options[filter_name]=False                        #Blocking filter
options[filter_name]['blocking_correction']=True  #Wether correction will be applied for partially blocked beams.
options[filter_name]['blocking_threshold']=0.5    #Beams with blocking above this threshold will be eliminated.
options[filter_name]['save']=True                 #Save blocking factor into qc_output dictionary.


#Low elevation angles filter parameters ==============================================================
options['iflefilter']=False                       #Low elevation angles filter.

options['flfilter_minangle']=2.0     #Reflectivity with echo tops lower than this angle will be eliminated.

#Dealiasing border filter 
options['ifdabfilter']=False         #Dealiasing border filter.
options['dabfilter_boxx']=3          #Edge "expansion" in azimuth
options['dabfilter_boxy']=3          #Edge "expansion" in range
options['dabfilter_boxz']=0          #Edge "expansion" in elevation

#Low doppler velocity filter
options['ifldvfilter']=False          #Low doppler velocity filter.
options['ldvfilter_vtr']=0.2          #Velocity threshold.
options['ldvfilter_htr']=2000         #Height threshold.
options['ldvfilter_use_terrain']=True #Wether AGL height is used. 

#Detect missing parameters
options['ifmissfilter']=False   #Missing values filter


#Topography parameters

options['toporawdatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/raw/"
options['toporadardatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/radar/"



