#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.

#======================================
# ENVIRONMENTAL VARIABLES 
#======================================

import os

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

#Dealiasing parameters (pyart)
options['ifdealias']=False

options['interval_split']=3
options['skip_between_rays']=10
options['skip_along_rays']=10

#Rho filter parameters
options['ifrhofilter']=False    #Rhohv filter

options['rhofilternx']=2
options['rhofilterny']=2
options['rhofilternz']=0
options['rhofiltertr']=0.5
options['rhofilter_save']=False

#Echo top parameters           #Filter layeers with echo top below a certain threshold.
options['ifetfilter']=False      #Echo top filter

options['etfilternx']=2        #Smooth parameter in x
options['etfilterny']=2        #Smooth parameter in y
options['etfilternz']=0        #Smooth parameter in z (dummy)
options['etfiltertr']=3000     #Echo top threshold.
options['etfilter_save']=True  #Wether echo top will be included in the output

#Echo depth parameters          #Filters layers with depths below a certain threshold.
options['ifedfilter']=False     #Echo depth filter

options['edfilternx']=2         #Smooth parameter in x
options['edfilterny']=2         #Smooth parameter in y
options['edfilternz']=0         #Smooth parameter in z (dummy)
options['edfiltertr']=3000      #Echo top threshold.
options['edfilter_save']=True   #Wether echo top will be included in the output

#Speckle parameters
options['ifspfilter']=False     #Speckle filter

options['spfilternx']=2           #Box size in X NX=(2*spfilternx + 1)
options['spfilterny']=2           #Box size in Y
options['spfilternz']=0           #Box size in Z
options['spfilterreftr']=5        #Reflectivity threshold
options['spfiltertr']=0.3         #Count threshold
options['spfilter_save']=True     #Save filter fields.
options['spfilter_ref']=True      #Wether the filter will be applied to Ref.
options['spfilter_v']=True        #Wether the filter will be applied to Vr

#Attenuation parameters
options['ifattfilter']=False    #Attenuation filter

options['attfiltertr']=20.0       #Attenuation threshold in dBZ
options['attcalerror']=1.0        #Calibration error
options['attfilter_save']=True    #Save filter fields

#Blocking parameters
options['ifblfilter']=True      #Blocking filter

options['blocking_correction']=True  #Wether correction will be applied for partially blocked beams.
options['blocking_threshold']=0.5    #Beams with blocking above this threshold will be eliminated.
options['blocking_save']=True        #Save blocking factor into qc_output dictionary.

#Detect missing parameters
options['ifmissfilter']=False   #Missing values filter


#Topography parameters

options['toporawdatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/raw/"
options['toporadardatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/radar/"


#QC CODES

#Reflectivity
options['qccodes']=dict()
options['qccodes']['QCCODE_ATTENUATION']=10
options['qccodes']['QCCODE_SPECKLE']=11
options['qccodes']['QCCODE_TEXTURE']=12
options['qccodes']['QCCODE_RHOFILTER']=13
options['qccodes']['QCCODE_SIGN']=14
options['qccodes']['QCCODE_BLOCKING']=15
options['qccodes']['QCCODE_ECHOTOP']=16
options['qccodes']['QCCODE_ECHODEPTH']=17

#Doppler velocity
options['qccodes']['QCCODE_DEALIAS']=30

#El codigo de los datos buenos para reflectividad y velocidad radial.
options['qccodes']['QCCODE_GOOD']=0



