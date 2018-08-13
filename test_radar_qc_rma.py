#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause
# History:  Created 10-2017

#Add path to additional modules.
import sys
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import radar_qc_module as rqc
import conf_rma as conf

#=======================================

# read in the file, create a RadarMapDisplay object
#filename='./cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3'    #Ejemplo Parana
filename='./cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3'  #Ejemplo Cordoba

#Performs QC operations based on options
[radar , qc_output] = rqc.main_qc( filename , conf.options )

print('End of QC')




