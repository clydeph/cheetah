//
//  hdf5-reader.h
//
//  Created by Anton Barty on 7/02/2014.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//  Modified by Chun Hong Yoon on 5/15/2018.

#ifndef __hdf5_reader__
#define __hdf5_reader__

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>

/*
 *  Structure for holding SACLA HDF5 metadata
 */
typedef struct {
	
	char filename[1024]; // Filename
	char version[1024]; // PAL hdf5 schema version

	char experimentID[1024]; // Run information

	hid_t   file_id;
	int32_t run_number;

	uint32_t    ndetectors;	
	uint64_t    nevents; // Detector event tags within a run

	int		*pumpLaserOn;
	int		*pumpLaserCode;
	float	*pumpLaserDelay;
	float	*photon_energy_keV;

	// Device objects within a run
	char	**detector_name;
	float	*detectorZ;
	int		*data_raw_is_float;


} h5_info_t;


/*
 *  Prototypes for functions written to read HDF5 data
 */
int HDF5_ReadHeader(const char*, h5_info_t*);
//int HDF5_Read2dDetectorFields(h5_info_t*, long);
//int HDF5_ReadEventTags(h5_info_t*, long);
//int HDF5_ReadRunInfo(h5_info_t*, long);
int HDF5_ReadImageRaw(h5_info_t*, long, void*);
int HDF5_cleanup(h5_info_t*);



#endif /* defined(__hdf5_reader__) */
