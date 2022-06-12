import os
import argparse
import time
import h5py
#from pathlib import Path
import numpy as np
from astropy.io import fits
from fast_histogram import histogram1d
from dataminer import *
from astropy.io import ascii
#from astropy.table import Table

from readpar import *



## This part takes the input file.
Parser = argparse.ArgumentParser(description='Code takes the barycentered fits file obtained from ASTROSAT-CZTI and ASTROSAT-LAXPC to prouduce the folded pulse profile. \n Requirements:\n (A) Fits file \n (B) A parameter file. \n (C) Name of the output file.')
Parser.add_argument('-i', '--input', help='FITS file (must be barycentered)', required=True)
Parser.add_argument('-p', '--parameter', help='Parameter file', required=True)
Parser.add_argument('-o', '--output', help='Name of the output file', required=True)
Parser.add_argument('-b', '--nbin', help='Number of bins in the profile', required=True)
#parser.add_argument('-T', '--tscrunch', help='(Y/N) Scrunch all time subints into one.', required =False, default='N')
args = Parser.parse_args()

Eventfile = args.input # Reading the name of the fits file from user input 
Parameter_file = args.parameter # Reading the name of the parameter file from user input
outfile = args.output # Reading the name of the output file from user input
nbin = int(args.nbin) # Reading the number of bins provided by the user.

t1 = time.time() # recording the time stamps in order to calculate how much time it takes to compute.
num_of_sub = 8.0
###################################################################################
## Reading the parameter file, converting it to hdf5 format
## Reading the FITS file and converting it to the hdf5 format , calculating the length of the time series and deciding 
## the number of block to read from hd5 file to fold.
####################################################################################
rot_par = read_parfile(Parameter_file)
rot_par.append(nbin)

t_0 = rot_par[0]
nu = rot_par[1]
nu_dot = rot_par[2]
nu_ddot = rot_par[3]
nbin = rot_par[4]


f_par = h5py.File("parameter.hdf5", "w")
dset_par = f_par.create_dataset("para", data=rot_par, dtype='float64')

INSTRUMENT = DATA_EXTRACTOR(Eventfile).GET_HEADER_INFO()

t_fits = DATA_EXTRACTOR(Eventfile).EXTRACT_DATA(INSTRUMENT)

#t_fits = scidata.field('TIME')
num_lines = len(t_fits)
n_pulse = (np.max(t_fits) - np.min(t_fits))*nu
print ("Number of pulses in the data: ", int(n_pulse))
print ("Number of photons:", num_lines)

block_len = num_lines/num_of_sub
block_len_i = int(np.floor(block_len))

#print ('Looping will be done over', num_of_sub, 'blocks each having', block_len_i,'size')
block_len_f = block_len - np.floor(block_len)
rem_lines = int (block_len_f*num_of_sub)
#print ('Last block will have', rem_lines, 'events.')
print ('FITS file closed')
#print ('Started making the hd5 file for time-stamps')
f_time = h5py.File("time_ticks.hdf5", "w")

t_hd5_st= time.time()
#print ('Started making the hdf dataset type of the times-stamps')
dset_t = f_time.create_dataset("time_ticks", data=t_fits, dtype='float64')
t_hd5_sp = time.time()
#print ('Created the hdf5 file for time-stamps and dataset type produced it took', (t_hd5_sp-t_hd5_st)/60., 'minutes to do')
block_index = np.arange(num_of_sub)


vals_use=[86400.0, 55197.0, 0.50000, 0.16666666, 1.0]
vals_use = np.array(vals_use)


t_min = np.min(dset_t)
t0 = (dset_par[0]-vals_use[1])*vals_use[0]
del_t_min = (t_min -t0)
phase_min = dset_par[1]*del_t_min + dset_par[2]*del_t_min*del_t_min*vals_use[2] + dset_par[3]*del_t_min*del_t_min*del_t_min*vals_use[3]

profile_total = np.zeros(int (nbin))

print ('***************** FOLDING STARTED ***************')

t_fold_st = time.time()

for i in range(len(block_index)):


	dset_t_l= dset_t[i*block_len_i:(i+1)*block_len_i]
	start_index = i*block_len_i
	end_index = (i+1)*block_len_i
	

	t_cal = dset_t_l

	del_t = (t_cal -t0)

	phase_tot = dset_par[1]*del_t + dset_par[2]*del_t*del_t*vals_use[2] + dset_par[3]*del_t*del_t*del_t*vals_use[3]
	phase_tot = phase_tot - phase_min
	phase_frac = phase_tot - np.floor(phase_tot)

	nbins = np.arange(dset_par[4])

	bins = np.floor(dset_par[4]*phase_frac);
	bins_phase =	((bins +0.5)/dset_par[4]);

	#bn = np.linspace(0.0,1.0, 129)
	#print len(bn)
	#prof, bin_edges = np.histogram(phase2, bins = bn)
	
	prof = histogram1d(phase_frac,range=[0,1],bins=int (nbin))
	#prof = np.histogram(phase_frac, bins = int (nbin))
	
	profile_total = profile_total + prof
	

	#print ('Block', i, 'has been folded')



##############################################
## Now it will start folding the remaining lines
## to get the total profile we put a condition if 
## number of events are less than the number of bins 
## over which it need to be folded we 
## throw away those events else we fold.
##############################################

if (rem_lines > nbin):

	dset_t_l= dset_t[num_of_sub*block_len_i:num_lines]

	t_cal = dset_t_l

	

	del_t = (t_cal -t0)


	phase_tot = dset_par[1]*del_t + dset_par[2]*del_t*del_t*vals_use[2] + dset_par[3]*del_t*del_t*del_t*vals_use[3]
	phase_tot = phase_tot - phase_min
	phase_frac = phase_tot - np.floor(phase_tot)



	nbins = np.arange(dset_par[4])


	bins = np.floor(dset_par[4]*phase_frac);
	bins_phase =	((bins +0.5)/dset_par[4]);

	
	prof = histogram1d(phase_frac,range=[0,1],bins=int (nbin))
	
	profile_total = profile_total + prof

	#print ('folding the last block')
#else:
#	print ('The last block has lesser the number of events than that of number of phase bins you want hence we skip those events')

print ('\n*****************FINISHED FOLDING****************\n')

print ('Total number of events:', np.sum(profile_total))
	

t_fold_stop = time.time()

print ('Folding took', (t_fold_stop - t_fold_st) , 'secs.' )

t2 = time.time()

data = [np.arange(nbin)+1, profile_total ]
ascii.write(data, outfile, names=['#BINS', 'counts'], overwrite=True)


os.remove("time_ticks.hdf5")
os.remove("parameter.hdf5")


