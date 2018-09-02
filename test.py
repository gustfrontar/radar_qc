
import os

datapath='/home/jruiz/test_qc/'

for (dirpath, dirnames, filenames) in os.walk( datapath ):

    print(dirnames,filenames)

