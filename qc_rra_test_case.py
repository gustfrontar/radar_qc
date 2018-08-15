#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause
# History:  Created 10-2017

import numpy as np
import os 
import glob
import sys
import matplotlib.pyplot as plt
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import radar_qc_module as rqc
import conf_rma
import conf_gematronik 

#La idea es hacer un loop sobre un grupo de carpetas / subcarpetas y archivos.
#En este caso asumo que hay 2 niveles de carpetas (por eso el /*/*/)
 
datapath='/media/jruiz/534AA1022E493735/DATOS_RADAR/RRA_radar/ANG'
 

#Generate a list with all the files contained in the folders and sub-folders.
file_list=[]
 
for (dirpath, dirnames, filenames) in os.walk(datapath):


    for filename in filenames:
        f = '/'.join([dirpath,filename])
        file_list.append(f)


#Proceed to perform qc for each file.
for ifile in file_list  :

   #Test if we have an RMA or GEMATRONIK RADAR
   if ifile.find('RMA') > 0  :
      options = conf_rma.options 
   else                   :
      options = conf_gematronik.options

   options['filename'] = ifile
   #Generate output file name
   options['filename_out'] = filename[ 0:filename.find('.nc') ] + 'corr' + '.nc' 

   options['plot']['Path'] = os.path.dirname(ifile) + '/plots/'

   os.makedirs( options['plot']['Path'],exist_ok=True)


   print(options['filename'])
   print(options['filename_out'])
   print(options['plot']['Path'])
   print(options['is_rma'])
#Performs QC operations based on options
#[radar , qc_output] = rqc.main_qc( conf.options )




