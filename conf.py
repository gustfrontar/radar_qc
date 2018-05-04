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
options['da_texture_filter']=True    #Wether a texture filter will be applied before performing dealiasing.
options['da_texture_thr']=1          #Texture filter threshold.
options['da_texture_nx']=3
options['da_texture_ny']=3
options['da_texture_nz']=0
options['da_texture_code']=44

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
options[filter_name]['force_value']=0.5                     #Threshold for force

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
options[filter_name]['force_value']=0.5                     #Threshold for force

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
options[filter_name]['force_value']=0.5                     #Threshold for force

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

#Doppler texture filter  parameters ==============================================================

filter_name='DopplerTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=3                                #NX
options[filter_name]['ny']=3                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,50,100,200])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=15
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force

#Reflectivity texture filter  parameters ==============================================================

filter_name='ReflectivityTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=3                                #NX
options[filter_name]['ny']=3                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,50,100,200])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=16
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
options[filter_name]['code']=17
options[filter_name]['attcalerror']=1.0                     #Calibration error
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force

#Blocking parameters               ==============================================================
#This filter is not part of the Fuzzy logic algorithm.
filter_name='BlockingFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                  #Blocking filter
options[filter_name]['blocking_correction']=True    #Wether correction will be applied for partially blocked beams.
options[filter_name]['blocking_threshold']=0.5      #Beams with blocking above this threshold will be eliminated.
options[filter_name]['save']=True                   #Save blocking factor into qc_output dictionary.


#Low elevation angles filter parameters ==============================================================
filter_name='LowElevFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                           #Low Elevation angle filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['min_angle']=2.0                       #Minimun angle
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=18
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force

#Low doppler velocity filter            ==============================================================
filter_name='LowDopplerFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                              #Low Doppler Velocity filter
options[filter_name]['nx']=0                                   #NX
options[filter_name]['ny']=0                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1,0,0])            #Importance function y
options[filter_name]['ifx']=np.array([-200,-1,-0.5,0.5,1,200]) #Importance function x
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 19
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['height_thr']=2000                        #Height threshold.
options[filter_name]['use_terrain']=True                       #Wether AGL height is used.

#Interference filter            ==============================================================
#This filter is not included in the Fuzzy-logic approach.
filter_name='InterferenceFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Interference filter
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['code']= 20
options[filter_name]['offset']=100                  #Number of ranges that will be discarded in the interference fit process.
options[filter_name]['att']=0.01                    #Atmospheric gases attenuation constant.
options[filter_name]['npass_filter']=3              #Number of passes of the azimuthal continuity filter.
options[filter_name]['percent_valid_threshold']=0.4 #Rays with valid pixels over this percentaje will be examinated.
options[filter_name]['corr_threshold']=0.6          #Rays that correlates well with the interference pattern will be flagged as 
                                                    #contaminated.
options[filter_name]['dbz_threshold']=5.0           #Reflectivity threshold to count pixels which are close to the interference pattern.
options[filter_name]['percent_ref_threshold']=0.6   #If more than this percent of the ray correlates well with the interference pattern, then
                                                    #the ray is flagged as contaminated by interference.

#Dealiasing border filter            ==============================================================
filter_name='DealiasingBorderFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                              #Low Doppler Velocity filter
options[filter_name]['nx']=3                                   #NX
options[filter_name]['ny']=3                                   #NY
options[filter_name]['nz']=3                                   #NZ
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force

#Doppler Noise filter            ==============================================================
filter_name='DopplerNoiseFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                              #Low Doppler Velocity filter
options[filter_name]['nx']=1                                   #NX
options[filter_name]['ny']=1                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['nx2']=10                                 #NX
options[filter_name]['ny2']=10                                 #NY
options[filter_name]['nz2']=0                                  #NZ
options[filter_name]['threshold_1']=4
options[filter_name]['threshold_2']=15       
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['ify_1']=np.array([0,0,1,1])              #Importance function y
options[filter_name]['ifx_1']=np.array([0,2,6,10])             #Importance function x
options[filter_name]['ify_2']=np.array([0,0,1,1])              #Importance function y
options[filter_name]['ifx_2']=np.array([0,10,15,20])           #Importance function x
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force

#Detect missing parameters
options['ifmissfilter']=False   #Missing values filter

#Topography parameters

options['toporawdatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/raw/"
options['toporadardatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/radar/"

