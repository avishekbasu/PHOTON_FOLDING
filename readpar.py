import sys
from detect_delimiter import detect
	



def read_parfile(fil):
	try:
		with open(fil) as f:
			lines = f.readlines()
			for line in lines:
				delim = detect(line)
				assert delim=='\t'
				if line.split('\t')[0] == 'PEPOCH':
					t_0 = (line.split('\t')[1])
				if line.split('\t')[0] == "F0":
					nu = line.split('\t')[1]
				if line.split('\t')[0] == "F1":
					nu_dot = line.split('\t')[1]
				if line.split('\t')[0] == "F2":
					nu_ddot = line.split('\t')[1]
		return [float(t_0), float(nu), float(nu_dot), float(nu_ddot)]
					
	except AssertionError:
		print ('Parameter file has wrong format, the parameter key word and values should be tab separated')
		sys.exit()
	
