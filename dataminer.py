import sys
from astropy.io import fits
import numpy as np

INST_DICT = {"NAME":['CZTI', 'LAXPC1', 'LAXPC2', 'LAXPC3'], "TABLES": [4, 1, 1, 1]}


class DATA_EXTRACTOR:
	def __init__(self, filename):
		self.FNAME = filename
		
	def GET_HEADER_INFO(self):
		HDULIST = fits.open(self.FNAME, memmap=True)
		INSTRUMENT = HDULIST[0].header['INSTRUME']
		HDULIST.close()
		if INSTRUMENT in INST_DICT["NAME"]:
			return INSTRUMENT
		else:
			print ('THIS CODE CANNOT PROCESS THE DATA FROM ', INSTRUMENT)
			return None
			sys.exit()



	def EXTRACT_DATA(self, INSTRUMENT):
	
		if INSTRUMENT == None:
			print ('CANNOT IDENTIFY THE BACKEND')
			sys.exit()
		else:
			N_HDULIST = INST_DICT["TABLES"][np.where(np.array(INST_DICT["NAME"]) == INSTRUMENT)[0][0]]
			HDULIST = fits.open(self.FNAME, memmap=True)
			TIME_ARRAY = np.array([])

			for i in range(N_HDULIST):
				TIME = HDULIST[i+1].data.field('Time')
				TIME_ARRAY = np.append(TIME_ARRAY, TIME) 
				del TIME
			HDULIST.close()
		return TIME_ARRAY


