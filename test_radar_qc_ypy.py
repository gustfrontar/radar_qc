#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.
import sys
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import numpy as np
import matplotlib.pyplot as plt
import pyart
import netCDF4
import radar_qc_module as rqc
import os

import qc_plot_tools as qcpt

import conf

#=======================================


# read in the file, create a RadarMapDisplay object
filename='./cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc'

#Performs QC operations based on options
[radar , qc_output] = rqc.main_qc( filename , conf.options )

print('End of QC')


elev=0

if conf.options['Dealiasing']['flag']  :
   qcpt.plot_dealiasing( qc_output , conf.options , figname='DealiasingTest.png' , elev=elev)
if conf.options['RhoFilter']['flag']   :
   qcpt.plot_rhofilter( qc_output , conf.options , figname='RhoFilterTest.png' , elev=elev )
if conf.options['EchoTopFilter']['flag'] :
   qcpt.plot_echotopfilter( qc_output , conf.options , figname='EchoTopFilter.png' , elev=elev )
if conf.options['EchoDepthFilter']['flag'] :
   qcpt.plot_echodepthfilter( qc_output , conf.options , figname='EchoDepthFilter.png' , elev=elev )
if conf.options['RefSpeckleFilter']['flag'] :
   qcpt.plot_refspecklefilter( qc_output , conf.options , figname='RefSpeckleFilter.png' , elev=elev )









