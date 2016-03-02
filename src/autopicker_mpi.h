/*
 * autopicker_mpi.h
 *
 *  Created on: Sep 18, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef AUTOPICKER_MPI_H_
#define AUTOPICKER_MPI_H_

#include "src/mpi.h"
#include "src/autopicker.h"
#include "src/parallel.h"

class AutoPickerMpi: public AutoPicker
{
private:
	MpiNode *node;

public:
	/** Destructor, calls MPI_Finalize */
    ~AutoPickerMpi()
    {
        delete node;
    }

    /** Read
     * This could take care of mpi-parallelisation-dependent variables
     */
    void read(int argc, char **argv);

    // Parallelized run function
    void run();

};

#endif /* AUTOPICKER_MPI_H_ */
