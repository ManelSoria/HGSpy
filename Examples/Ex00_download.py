#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
#  Download HGS database
from __future__ import print_function

import HGSpy as HGS


# Call the HGS download script to create a new HGS database
DOWNLOAD = False # Whether to rebuild or just download the database
INFO     = True  # Print information when downloading
db = HGS.data.download(download=DOWNLOAD,info=INFO)


# Now we have a new database object, let us add some interesting
# mixtures
db.add('air7921',['N2','O2'],[79,21])                              # Dry air with 79% N2 and 21% O2
db.add('air',['N2','O2','Ar','CO2'],[78.084,20.946,0.9340,0.036]) # Dry air
db.add('RP1',['C10H22','C10H18'],[61.25, 38.75])                   # Rocket propellant


# Finally we store the database inside the HGS package
db.save(HGS.HGSDATA)