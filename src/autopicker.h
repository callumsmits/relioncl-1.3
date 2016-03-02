/*
 * autopicker.h
 *
 *  Created on: Sep 18, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef AUTOPICKER_H_
#define AUTOPICKER_H_
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/projector.h"
#include "src/ctf.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"

struct Peak
{
	int x;
	int y;
	int ref;
	double psi;
	double fom;
	double relative_fom;
};

class AutoPicker
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_in, fn_ref, fns_autopick, fn_out;

	// Pixel size (for low-pass filter and particle diameter)
	double angpix;

	// Metadata of the micrographs
	MetaDataTable MDmic;

	// Particle diameter (in Angstroms)
	double particle_diameter;
	int particle_radius2;

	// Low pass filetr cutoff (in Angstroms)
	double lowpass;

	// Original size of the reference images
	int particle_size;

	// Dimension of the filtered image
	int current_size;

	// Vector with all original reference images
	std::vector<MultidimArray<double> > Mrefs;

	// FTs of the reference images (either for autopicking or for feature calculation)
	std::vector<Projector > PPref;

	///// Autopicking stuff

	// Re-read precalculated best_localCCF and SPI arrays from disc
	bool do_read_fom_maps;

	// Write precalculated best_localCCF and SPI arrays to disc
	bool do_write_fom_maps;

	// All micrographs to autopick from
	std::vector<FileName> fn_micrographs;

	// Original size of the micrographs
	int micrograph_size, micrograph_xsize, micrograph_ysize;

	// Is density in micrograph inverted wrt templates?
	bool do_invert;

	// Correct the references for CTF effects?
	bool do_ctf;

	// Keep the CTFs unchanged until the first peak?
	bool intact_ctf_first_peak;

	// Apart from keeping particle_size/2 away from the sides, should we exclude more? E.g. to get rid of Polara bar code?
	int autopick_skip_side;

	// In-plane rotational sampling (in degrees)
	double psi_sampling;

	// Fraction of expected probability ratio to consider as peaks
	double min_fraction_expected_Pratio;

	// Number of Angstroms any 2 particle peaks need to be apart
	double min_particle_distance;

	// Size of the downsize micrographs for autopicking
	int downsize_mic;

	// Number of non-zero pixels in the circular mask, and of its inverse (for background normalisation in do_diff2)
	int nr_pixels_circular_mask, nr_pixels_circular_invmask;

	// Array with Fourier-transform of the (circular) mask, and of its inverse
	MultidimArray<Complex > Fmsk, Finvmsk;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some general stuff after reading
	void initialise();

	// General function to decide what to do
	void run();

	void autoPickOneMicrograph(FileName &fn_mic);

private:

	// Uses Roseman2003 formulae to calculate stddev under the mask through FFTs
	// The FFTs of the micrograph (Fmic), micrograph-squared (Fmic2) and the mask (Fmsk) need to be provided at downsize_mic
	// The putput (Mstddev) will be at (binned) micrograph_size
	void calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &Fmic, const MultidimArray<Complex > &Fmic2,
			MultidimArray<Complex > &Fmsk, int nr_nonzero_pixels_mask, MultidimArray<double> &Mstddev, MultidimArray<double> &Mmean);

	// Peak search for all pixels above a given threshold in the map
	void peakSearch(const MultidimArray<double> &Mccf, const MultidimArray<double> &Mpsi, int iref, int skip_side, std::vector<Peak> &peaks);

	// Now prune the coordinates: within min_particle_distance: all peaks are the same cluster
	// From each cluster, take the single peaks with the highest ccf
	// If then, there is another peaks at a distance of at least min_particle_distance: take that one as well, and so forth...
	void prunePeakClusters(std::vector<Peak> &peaks, int min_distance);


	// Only keep those peaks that are at the given distance apart from each other
	void removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance);


};


#endif /* AUTOPICKER_H_ */
