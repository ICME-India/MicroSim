/*
 * Contains all functions required for reading user input, and
 * sharing read input with all processes.
 */

#ifndef INPUTREADER_H_
#define INPUTREADER_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#include "structures.h"
#include "utilityFunctions.h"

/*
 * Reads input from specified file and broadcasts it to all ranks
 * If the length of any input array is not known at compile time,
 * it will have to be broadcast after input is read.
 * This is done by allocating memory for all ranks other than the
 * master, and then broadcasting.
 * A log containing the read input is written to a file with the
 * same name as the input file, and is stored in the .out format
 * in the DATA folder.
 */
void readInput_MPI(domainInfo *simDomain, controls *simControls,
                   simParameters *simParams, int numRanks,
                   char *argv[]);

/*
 * Reads filling parameters from the specified file.
 * The filling routines can be found in initialise.h and its
 * related files.
 */
void readFill(fillParameters *fill, char *argv[], int rank);

#endif
