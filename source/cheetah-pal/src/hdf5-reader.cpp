//
//  hdf5-reader.cpp
//
//

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <cstddef>

#include "../include/hdf5-reader.h"

/*unsigned getRunNumber(const std::string& str)
{
  std::size_t found = str.find_last_of("/\\");
  std::size_t places = 7;
  char runNum[1024];
  strncpy( runNum, str.substr(found+1).c_str(), places );
  return atoi(runNum);
}*/


unsigned getRunNumber(const std::string& str)
{
	std::size_t found1 = str.find_last_of("/\\");
	std::size_t found2 = str.find_last_of(".");
	if (std::string::npos == found1) found1 = 0;
	else found1++;
	if (std::string::npos == found2) return -1;
	if (found2 <= found1) return -1;

	std::string fn = str.substr(found1, found2 -found1);
	if (5 == fn.length())
	{
		// /xfel/ffs/dat/.../raw_data/R1234.h5
		return atoi(fn.substr(1,4).c_str());
	}
	else if (7 == fn.length())
	{
		// /xfel/ffs/dat/.../raw_data/1234567.h5
		return atoi(fn.c_str());
	}

	return -1;
}


/*
 *  Function for parsing PAL HDF5 metadata
 */
int HDF5_ReadHeader(const char *filename, h5_info_t *result) {
    
    char h5field[1024];

    herr_t status;
    hid_t dataset;
    hid_t datatype;
    //H5T_class_t t_class;
    //H5T_order_t order;
    //size_t datatype_size;
    hsize_t dims_out[5];
    hid_t dataspace;

    int rank;
    int status_n;

    // Get run number
    result->run_number = getRunNumber(filename);
	if(result->run_number < 0)
	{
		printf("ERROR: Invalid RunNumber!\n");
		return -1;
	}
	printf("RunNumber : %d\n", result->run_number);

    // Does the file exist?
    FILE *fp = fopen(filename, "r");
    if(fp)
	{
		fclose(fp);
	}
	else
	{
		printf("ERROR: File does not exist %s\n",filename);
		exit(1);
	}

    // Open PAL HDF5 file and read in header information
    hid_t   file_id;
    file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
	if(file_id < 0)
	{
		printf("ERROR: Could not open HDF5 file %s\n",filename);
		return -1;
	}
    result->file_id = file_id;
    strcpy(result->filename, filename);

	// Get File Format version
	memset(result->version, 0, 1024);
	sprintf(h5field, "/R%04d/header/file_version", result->run_number);
	status = H5LTread_dataset_string(file_id, h5field, result->version);
	if (status < 0)
	{
		printf("ERROR: Could not read file_version\n");
		return -1;
	}
	printf("Data File Format version: %s\n", result->version);

    // Get experiment ID
	memset(result->experimentID, 0, 1024);
    sprintf(h5field, "/R%04d/header/exp_id", result->run_number);
    status = H5LTread_dataset_string(file_id, h5field, result->experimentID);
	if (status < 0)
	{
		printf("ERROR: Could not read exp_id\n");
		return -1;
	}
	printf("Exp ID: %s\n", result->experimentID);


	// Get number of events from length of N
	sprintf(h5field, "/R%04d/scan_dat/N", result->run_number);
	dataset = H5Dopen2(file_id, h5field, H5P_DEFAULT);
	if (dataset < 0)
	{
		printf("ERROR: Could not read %s\n", h5field);
		return -1;
	}
	dataspace = H5Dget_space(dataset);
	rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	status = H5Sclose(dataspace);
	status = H5Dclose(dataset);
	printf("Number of Events: %lu\n", (unsigned long)(dims_out[0]));
    result->nevents = dims_out[0];

    // Get list of photon energies
    result->photon_energy_keV = (float*) calloc( result->nevents, sizeof(float));
    sprintf(h5field, "/R%04d/scan_dat/photon_energy", result->run_number);
    status = H5LTread_dataset_float(file_id, h5field, result->photon_energy_keV);
	if (status < 0)
	{
		printf("ERROR: Could not read %s\n", h5field);
		return -1;
	}

    // Get list of pump laser code
    result->pumpLaserOn = (int*) calloc( result->nevents, sizeof(int) );
    sprintf(h5field, "/R%04d/scan_dat/laser_on", result->run_number);
    status = H5LTread_dataset_int(file_id, h5field, result->pumpLaserOn);
	if (status < 0)
	{
		printf("ERROR: Could not read %s\n", h5field);
		return -1;
	}

    // Get list of pump laser delay
    result->pumpLaserDelay = (float*) calloc( result->nevents, sizeof(float) );
    sprintf(h5field, "/R%04d/scan_dat/laser_delay_time", result->run_number);
    status = H5LTread_dataset_float(file_id, h5field, result->pumpLaserDelay);
	if (status < 0)
	{
		printf("ERROR: Could not read %s\n", h5field);
		return -1;
	}

    // Get number of detectors
    sprintf(h5field, "/R%04d/header/detector_num", result->run_number);
    status = H5LTread_dataset(file_id, h5field, H5T_NATIVE_UINT, &(result->ndetectors));
	if (status < 0)
	{
		printf("ERROR: Could not read %s\n", h5field);
		return -1;
	}
	printf("Number of Detectors: %d\n", result->ndetectors);

    // Get detector names & distances
    result->detector_name = (char**) calloc(result->ndetectors, sizeof(char*));
	result->detectorZ = (float*)calloc(result->ndetectors, sizeof(float));
	result->data_raw_is_float = (int*)calloc(result->ndetectors, sizeof(int));
    for(int i=0; i<result->ndetectors; i++)
	{
		sprintf(h5field, "/R%04d/header/detector_%d_name", result->run_number, i);
		result->detector_name[i] = (char*) calloc(1024, sizeof(char));
		status = H5LTread_dataset_string(file_id, h5field, result->detector_name[i]);
		if (status < 0)
		{
			printf("ERROR: Could not read %s\n", h5field);
			return -1;
		}
		printf("Detector_%d name: %s\n", i, result->detector_name[i]);

		sprintf(h5field, "/R%04d/header/detector_%d_distance", result->run_number, i);
		status = H5LTread_dataset(file_id, h5field, H5T_NATIVE_FLOAT, &result->detectorZ[i]);
		if (status < 0)
		{
			printf("ERROR: Could not read %s\n", h5field);
			return -1;
		}
		printf("Detector_%d Distance: %f\n", i, result->detectorZ[i]);

		sprintf(h5field, "/R%04d/scan_dat/%s_data", result->run_number, result->detector_name[i]);
		dataset = H5Dopen2(file_id, h5field, H5P_DEFAULT);
		if (dataset < 0)
		{
			printf("ERROR: Could not read %s\n", h5field);
			return -1;
		}
		datatype = H5Dget_type(dataset);
		//t_class = H5Tget_class(datatype);
		//order = H5Tget_order(datatype);
		//datatype_size  = H5Tget_size(datatype);

		dataspace = H5Dget_space(dataset);    /* dataspace handle */
		rank      = H5Sget_simple_extent_ndims(dataspace);
		status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

		printf("Detector_%d Data Shape: RANK %d, (%ld", i, rank, (unsigned long)(dims_out[0]));
		for(int j=1; j<rank; j++) printf(" x %ld", dims_out[j]);
		printf(")\n");

		char dtype_text[1024];
		size_t dtype_text_size=1024;
		H5LTdtype_to_text(datatype, dtype_text, H5LT_DDL, &dtype_text_size);
		printf("Detector_%d Data Type: %s\n", i, dtype_text);

		status = H5Tclose(datatype);
		status = H5Sclose(dataspace);
		status = H5Dclose(dataset);

		if ( 0 != strcmp(dtype_text, "H5T_STD_U16LE") 
			&& 0 != strcmp(dtype_text, "H5T_NATIVE_USHORT") 
			&& 0 != strcmp(dtype_text, "H5T_IEEE_F32LE") 
			&& 0 != strcmp(dtype_text, "H5T_NATIVE_FLOAT"))
		{
			printf("Unsupported Data Type!!\n");
			return -1;
		}

		if ( 0 == strcmp(dtype_text, "H5T_IEEE_F32LE") 
			|| 0 == strcmp(dtype_text, "H5T_NATIVE_FLOAT") ) result->data_raw_is_float[i] = true;
		else result->data_raw_is_float[i] = false;
	}


    return 0;
}

/*
 *  Function for reading all <n> 2D detectors into one massive 2D array (for passing to Cheetah or CrystFEL)
 */
int HDF5_ReadImageRaw(h5_info_t *header, long eventID, void *pBuffer) {

    hid_t	dataset;         /* handles */
    hid_t	datatype, dataspace;
    hid_t	memspace;

    hsize_t     dimsm[3];              /* memory space dimensions */
    hsize_t     dims_out[3];           /* dataset dimensions */

    hsize_t	 count[3];              /* size of the hyperslab in the file */
    hsize_t	 offset[3];             /* hyperslab offset in the file */
    hsize_t	 count_out[3];          /* size of the hyperslab in memory */
    hsize_t	 offset_out[3];         /* hyperslab offset in memory */
    int          status, status_n, rank;
    int RANK_OUT = 3;

    // FIXME: how to figure out which detector gets used for hit finding
    char    h5field[1024];
    sprintf(h5field, "/R%04d/scan_dat/%s_data", header->run_number, header->detector_name[0]);
    dataset = H5Dopen2( header->file_id, h5field, H5P_DEFAULT);
	if (0 > dataset) return -1;

    /*
     * Get datatype and dataspace handles and then query
     * dataset class, order, size, rank and dimensions.
     */
    datatype  = H5Dget_type(dataset);     /* datatype handle */
	if (0 > datatype) return -1;
    dataspace = H5Dget_space(dataset);    /* dataspace handle */
	if (0 > dataspace) return -1;
    rank      = H5Sget_simple_extent_ndims(dataspace);
	if (0 > rank) return -1;
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	if (0 > status_n) return -1;

	// Define hyperslab in the dataset.
    offset[0] = eventID;
    offset[1] = 0;
    offset[2] = 0;
    count[0]  = 1;
    count[1]  = dims_out[1];
    count[2]  = dims_out[2];
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	if (0 > status) return -1;

	// Define the memory dataspace.
    dimsm[0] = dims_out[0];
    dimsm[1] = dims_out[1];
    dimsm[2] = dims_out[2];
    memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);
	if (0 > memspace) return -1;

	// Define memory hyperslab.
    offset_out[0] = 0;
    offset_out[1] = 0;
    offset_out[2] = 0;
    count_out[0]  = 1;
    count_out[1]  = dims_out[1];
    count_out[2]  = dims_out[2];
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	if (0 > status) return -1;

	// Read data from hyperslab in the file into the hyperslab in memory and display.
    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, pBuffer);
	if (0 > status) return -1;

	// Close/release resources.
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    return 1;
}


/*
 *  Cleanup stale HDF5 references and close the file
 */
int HDF5_cleanup(h5_info_t *header) {

    std::cout << "Cleaning up HDF5 links\n";
	hid_t ids[256];
	long n_ids = H5Fget_obj_ids(header->file_id, H5F_OBJ_ALL, 256, ids);
	for (long i=0; i<n_ids; i++ ) {

		hid_t id;
		H5I_type_t type;

		id = ids[i];
		type = H5Iget_type(id);

		if ( type == H5I_GROUP )
			H5Gclose(id);
		if ( type == H5I_DATASET )
			H5Dclose(id);
		if ( type == H5I_DATASPACE )
			H5Sclose(id);
		if ( type == H5I_DATATYPE )
			H5Dclose(id);
	}

	H5Fclose(header->file_id);


	free(header->photon_energy_keV);
	free(header->pumpLaserOn);
	free(header->pumpLaserDelay);

	for(int i = 0; i < header->ndetectors; i++)
	{
		free(header->detector_name[i]);
	}
	free(header->detector_name);
	free(header->detectorZ);
	free(header->data_raw_is_float);

    return 1;
}

