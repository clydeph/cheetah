#
#   Crawler routines for PAL/HDF5
#   Tested using Anaconda / Python 3.4
#

import glob
import lib.cfel_filetools as cfel_file



def scan_data(data_dir):

	debug = False
		
	pattern = data_dir + '/[0-9][0-9][0-9][0-9][0-9][0-9][0-9]/[0-9][0-9][0-9][0-9][0-9][0-9][0-9].h5'
	files = glob.glob(pattern)

	pattern = data_dir + '/[r,R][0-9][0-9][0-9][0-9]/[r,R][0-9][0-9][0-9][0-9].h5'
	files = files + glob.glob(pattern)

	# Create sorted file list (glob seems to return files in random order)
	files.sort()

	if debug:
		print(files)

	# Extract the run bit from HDF5 file name
	# Turn the 000000X string into an integer (for compatibility with old crawler files)
	out = []
	for filename in files:
		#thisrun = filename.split('.')[0]
		#thisrun = int(thisrun[-7::])
		#out.append(thisrun)
		thisrun = filename.split('/')[-1].split('.')[0]
		if 7 == len(thisrun) : out.append(int(thisrun))
		elif 5 == len(thisrun) : out.append(int(thisrun[1:]))


	# Find unique run values (due to multiple XTC files per run)
	run_list = list(sorted(set(out)))
	nruns = len(run_list)


	# Default status for each is ready
	status = ['Ready']*nruns

	# Loop through file names checking for '.inprogress' suffix
	for filename in files:
		if filename.endswith('.inprogress'): # FIXME: find out how PAL handles filenames when copying
			thisrun = filename.split('-')[1]
			thisrun = thisrun[1:5]
			if thisrun in run_list:
				run_indx = run_list.index(thisrun)
				status[run_indx] = 'Copying'

		if filename.endswith('.fromtape'):
			thisrun = filename.split('-')[1]
			thisrun = thisrun[1:5]
			if thisrun in run_list:
				run_indx = run_list.index(thisrun)
				status[run_indx] = 'Restoring'

	# Create the result
	result = {
		'run': run_list,
		'status' : status
	}

	if debug:
		print(result['run'])

	# Write dict to CSV file
	keys_to_save = ['run','status']
	cfel_file.dict_to_csv('data_status.csv', result, keys_to_save)

#end scan_data
