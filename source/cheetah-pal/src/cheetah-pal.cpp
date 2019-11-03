//
//  cheetah-pal.cpp
//
//  Originally created by Anton Barty on 20/1/14.
//  Copyright (c) 2014 Anton Barty. All rights reserved.
//
//  Created by Chun Hong Yoon on 21/3/18.
//  Copyright (c) 2018 Chun Hong Yoon. All rights reserved.
//
//  Modified by Sang-Youn Park on 29/10/19.
//  Copyright (c) 2019 Sang-Youn Park. All rights reserved.

#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cheetah.h"
#include "hdf5-reader.h"

int main(int argc, const char * argv[])
{
	static time_t startT = 0;

	// PAL hdf5 input file and Cheetah configuration file
	char fname[1024];
	char cheetahini[1024];

	static cGlobal cheetahGlobal;
	h5_info_t *header = NULL;
	void *pBuffer = NULL;

	time_t	tnow;
	double	dtime, datarate;
	long    detID = 0;

	if (argc != 3) {
		printf("Usage: %s pal.h5 cheetah.ini\n", argv[0]);
		return -1;
	}

	time(&startT);

	// Take configuration from command line arguments
	strcpy(fname, argv[1]);
	strcpy(cheetahini, argv[2]);


	//Initialise Cheetah
	printf("Setting up Cheetah...\n");
	strcpy(cheetahGlobal.configFile, cheetahini);
	sprintf(cheetahGlobal.facility, "PAL");
	if (cheetahInit(&cheetahGlobal))
	{
		printf("Fail cheetahInit()\n");
		return -1;
	}
	printf("Done cheetahInit\n");


	// Open PAL HDF5 file
	// Read file header and information
	printf("\nStart HDF5_ReadHeader: %s\n", fname);
	header = (h5_info_t*)malloc(sizeof(h5_info_t));
	if (HDF5_ReadHeader(fname, header))
	{
		printf("Fail HDF5_ReadHeader()\n");
		return -1;
	}
	// Gather detector fields and event tags for this run
	printf("Done HDF5_ReadHeader\n\n");


	//------------------------------------------
	// Copy header information to cheetahGlobal
	//------------------------------------------
	sprintf(cheetahGlobal.experimentID, header->experimentID);
	cheetahGlobal.runNumber = header->run_number;


	// Create a buffer for holding the detector image data
	if (header->data_raw_is_float[detID])
		pBuffer = calloc(cheetahGlobal.detector[detID].pix_nn, sizeof(float));
	else
		pBuffer = calloc(cheetahGlobal.detector[detID].pix_nn, sizeof(uint16_t));
	

	// Loop through all events found in this run
	for(uint64_t eventID=0; eventID < header->nevents; eventID++)
	{
		// Cheetah: Calculate time beteeen processing of data frames
		time(&tnow);
		dtime = difftime(tnow, cheetahGlobal.tlast);
		if(dtime > 1.)
		{
			datarate = (eventID - cheetahGlobal.lastTimingFrame)/dtime;
			cheetahGlobal.lastTimingFrame = eventID;
			time(&cheetahGlobal.tlast);
			cheetahGlobal.datarate = datarate;
		}

		// Read next image
		if (HDF5_ReadImageRaw(header, eventID, pBuffer) < 0)
		{
			printf("Fail - ReadImageRaw %lud\n", eventID);
			break;
		}


		// Cheetah: Create a new eventData structure in which to place all information
		cEventData	*eventData;
		eventData = cheetahNewEvent(&cheetahGlobal);


		// Cheetah: Populate event structure with meta-data
		sprintf(eventData->filename, fname);
		eventData->frameNumber = eventID;
		eventData->runNumber = header->run_number;
		eventData->nPeaks = 0;
		eventData->pumpLaserOn = header->pumpLaserOn[eventID];
		//eventData->pumpLaserCode = header->pumpLaserCode[eventID];
		eventData->pumpLaserDelay = (double) header->pumpLaserDelay[eventID];
		eventData->photonEnergyeV = header->photon_energy_keV[eventID]*1000;
		eventData->wavelengthA = 12398.42 / eventData->photonEnergyeV; // 4.1357E-15 * 2.9979E8 * 1E10 / eV (A)
		eventData->pGlobal = &cheetahGlobal;
		eventData->fiducial = eventID; // must be unique
		sprintf(eventData->eventname, "%07lud",eventID);
		//eventData->detectorDistance = header->detectorPosition;
		eventData->detector[detID].detectorZ = header->detectorZ[detID];

		if (header->data_raw_is_float[detID])
		{
			// Jungfrau
			eventData->detector[detID].data_raw_is_float = true;
			memcpy(eventData->detector[detID].data_raw, pBuffer, cheetahGlobal.detector[detID].pix_nn * sizeof(float));
		}
		else
		{
			// Rayonix
			eventData->detector[detID].data_raw_is_float = false;
			memcpy(eventData->detector[detID].data_raw16, pBuffer, cheetahGlobal.detector[detID].pix_nn * sizeof(uint16_t));
		}

		// Cheetah: Process this event
		cheetahProcessEventMultithreaded(&cheetahGlobal, eventData);
	}

	// Clean up stale IDs and exit
	HDF5_cleanup(header);

	time_t endT;
	time(&endT);
	double diff = difftime(endT,startT);
	std::cout << "time taken: " << diff << " second(s)\n";

	// Cheetah: Cleanly exit by closing all files, releasing memory, etc.
	cheetahExit(&cheetahGlobal);

	if (pBuffer) free(pBuffer);
	if (header) free(header);

	return 0;
}

