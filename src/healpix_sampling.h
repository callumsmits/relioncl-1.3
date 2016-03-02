/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef _HEALPIX_SAMPLING_HH
#define _HEALPIX_SAMPLING_HH

#include "src/Healpix_2.15a/healpix_base.h"
#include "src/metadata_table.h"
#include "src/macros.h"
#include "src/multidim_array.h"
#include "src/symmetries.h"
#include "src/euler.h"

// For the angular searches
#define NOPRIOR 0
#define PRIOR_ROTTILT_PSI 1


class HealpixSampling
{

public:
	/** Healpix sampling object */
    Healpix_Base healpix_base;

    /** Random perturbation */
    double random_perturbation;

    /** Amount of random perturbation */
    double perturbation_factor;

    /** In-plane (psi-angle) sampling rate
     */
    double psi_step;

    /** Healpix order */
    int healpix_order;

    /** Mode for orientational prior distribution
     * Note this option is not written to the STAR file, as it really belongs to mlmodel.
     * It is included here for convenience, and always needs to be set in initialise
     */
    int orientational_prior_mode;

    /* Translational search range and sampling rate
     */
    double offset_range, offset_step;

    /** Flag whether this is a real 3D sampling */
    bool is_3D;

    /** Flag whether the translations are 3D (for volume refinement) */
    bool is_3d_trans;

    /** Name of the Symmetry group */
    FileName fn_sym;

    /** List of symmetry operators */
    std::vector <Matrix2D<double> > R_repository, L_repository;

    /** Two numbers that describe the symmetry group */
    int pgGroup;
    int pgOrder;

    /** Limited tilt angle range */
    double limit_tilt;

    /** vector with the original pixel number in the healpix object */
    std::vector<int> directions_ipix;

    /** vector with sampling points described by angles */
    std::vector<Matrix1D<double> > directions_angles;

    /** vector with the psi-samples */
    std::vector<double> psi_angles;

    /** vector with the X,Y(,Z)-translations */
    std::vector<Matrix1D<double> > translations;

    /** vector with pointers to the (rot,tilt) pairs (directions) that have non-zero prior probability */
    std::vector<int> pointer_dir_nonzeroprior;

    /** vector with pointers to the psi-samples that have non-zero prior probability */
    std::vector<int> pointer_psi_nonzeroprior;

//TMP DEBUGGING: normally protected!
public:
    /** vector with the prior probabilities for those directions that have non-zero prior probability*/
    std::vector<double> directions_prior;

    /** vector with the prior probabilities for those psi-angles that have non-zero prior probability*/
    std::vector<double> psi_prior;



public:

    // Empty constructor
    HealpixSampling() {}

    // Destructor
    ~HealpixSampling()
    {
    	directions_ipix.clear();
    	directions_angles.clear();
    	psi_angles.clear();
    	pointer_dir_nonzeroprior.clear();
    	pointer_psi_nonzeroprior.clear();
    	directions_prior.clear();
    	psi_prior.clear();
    	translations.clear();

    }

    // Start from all empty vectors and meaningless parameters
    void clear();

    /** Set up the entire sampling object
     *
     * The in-plane (psi-angle) sampling is linear,
     * input_psi_sampling is modified to contain an integer number of equally-sized sampling points
     * For the 3D-case, a negative input_psi_sampling will result in a psi-sampling similar to the sqrt of a HealPix pixel area.
	 *
     * The HEALPix sampling is implemented as described by Gorski et al (2005), The Astrophysical Journal, 622:759-771
     * The order defines the number of sampling points, and thereby the angular sampling rate
     * From this paper is the following table:
     *
     * order	Npix	Theta-sampling
     * 0		12		58.6
     * 1		48		29.3
     * 2		192		14.7
     * 3		768		7.33
     * 4		3072	3.66
     * 5		12288	1.83
     * 6		49152	0.55
     * 7		196608	0.28
     * 8		786432	0.14
     * etc...
     *
     * */
    void initialise(int prior_mode, int ref_dim = -1, bool do_3d_trans = false);

    // Reset the random perturbation
    void resetRandomlyPerturbedSampling();

    // Read in all information from the command line to build the sampling object
    void read(IOParser &parser, int ori_size, int ref_dim);

    // Read CL options after a -continue statement.
    void readContinue(int argc, char **argv, bool &directions_have_changed);

    // Read in all information from a STAR file (for restarting)
    void read(FileName fn_in);

    // Write the sampling information to a STAR file
    void write(FileName fn_out);

    /* Set the non-oversampled list of translations */
    void setTranslations(double offset_step = -1., double offset_range = -1.);

    /* Set only a single translation */
    void setOneTranslation(Matrix1D<double> offset);

    /* Set the non-oversampled lists of directions and in-plane rotations */
    void setOrientations(int _order = -1, double _psi_step = -1.);

    /* Set only a single orientation */
    void setOneOrientation(double rot, double tilt, double psi);


    /* Write all orientations as a sphere in a bild file
     * Mainly useful for debugging */
    void writeAllOrientationsToBild(FileName fn_bild, std::string rgb = "1 0 0", double size = 0.025);
    void writeNonZeroPriorOrientationsToBild(FileName fn_bild, double rot_prior, double tilt_prior, std::string rgb = "0 0 1", double size = 0.025);

    /* Select all orientations with zero prior probabilities
     * store all these in the vectors pointer_dir_nonzeroprior and pointer_psi_nonzeroprior
     * Also precalculate their prior probabilities and store in directions_prior and psi_prior
     */
    void selectOrientationsWithNonZeroPriorProbability(
    		double prior_rot, double prior_tilt, double prior_psi,
    		double sigma_rot, double sigma_tilt, double sigma_psi,
    		double sigma_cutoff = 3.);

    /* Randomly reject part of te non-zero prior probabilities, so that the optimization no longer follows the steepest downward gradient
     * This procedure was inspired by Hans Elmlund's PRIME algorithm.
     */
    void randomSelectionNonZeroPriorProbability(double fraction_to_keep);


    /** Get the symmetry group of this sampling object
     */
    FileName symmetryGroup();

    /* Get the original HEALPix index for this direction
     * Note that because of symmetry-equivalence removal idir no longer corresponds to the HEALPix pixel number
     *
     */
    long int getHealPixIndex(long int idir);

    /** The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
     */
    void checkDirection(double &rot, double &tilt);

    /* Get the rot and tilt angles in the center of the ipix'th HEALPix sampling pixel
     * This involves calculations in the HEALPix library
     */
    void getDirectionFromHealPix(long int ipix, double &rot, double &tilt);

    /* Get the translational sampling step in pixels */
    double getTranslationalSampling(int adaptive_oversampling = 0);

    /* Get approximate angular sampling in degrees for any adaptive oversampling
     */
    double getAngularSampling(int adaptive_oversampling = 0);

    /* Get the number of symmetry-unique sampling points
     * Note that because of symmetry-equivalence removal this number is not the number of original HEALPix pixels
     * In the case of orientational priors, the number of directions with non-zero prior probability is returned
     */
    long int NrDirections(int oversampling_order = 0, bool include_zeroprior = false);

    /* Get the number of in-plane (psi-angle) sampling points
     */
    long int NrPsiSamplings(int oversampling_order = 0, bool include_zeroprior = false);

    /* Get the number of in-plane translational sampling points
     */
    long int NrTranslationalSamplings(int oversampling_order = 0);

    /* Get the total number of (oversampled) sampling points, i.e. all (rot, tilt, psi, xoff, yoff) quintets
    */
    long int NrSamplingPoints(int oversampling_order = 0, bool include_zeroprior = false);

    /* How often is each orientation oversampled? */
    int oversamplingFactorOrientations(int oversampling_order);

    /* How often is each translation oversampled? */
    int oversamplingFactorTranslations(int oversampling_order);

    /* Get the rot and tilt angles from the precalculated sampling_points_angles vector
     * This does not involve calculations in the HEALPix library
     * Note that because of symmetry-equivalence removal idir no longer corresponds to the HEALPix pixel number
     */
    void getDirection(long int idir, double &rot, double &tilt);

    /* Get the value for the ipsi'th precalculated psi angle
     */
    void getPsiAngle(long int ipsi, double &psi);

    /* Get the value for the itrans'th precalculated translations
     */
    void getTranslation(long int itrans, Matrix1D<double> &trans);

    /* Get the position of this sampling point in the original array */
    long int getPositionSamplingPoint(int iclass, long int idir, long int ipsi, long int itrans);

    /* Get the position of this sampling point in the oversampled array */
    long int getPositionOversampledSamplingPoint(long int ipos, int oversampling_order, int iover_rot, int iover_trans);

    /* Get the vectors of (xx, yy) for a more finely (oversampled) translational sampling
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same (xx, yy) pair as the original itrans.
     * An oversampling_order == 1 will give rise to 2*2 new (rot, tilt) pairs.
     * An oversampling_order == 2 will give rise to 4*4 new (rot, tilt) pairs.
     * etc.
     */
    void getTranslations(long int itrans, int oversampling_order,
								    std::vector<Matrix1D<double> > &my_translations);

    /* Get the vectors of (rot, tilt, psi) angle triplets for a more finely (oversampled) sampling
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same (rot, tilt, psi) triplet as the original ipix.
     * An oversampling_order == 1 will give rise to 2*2*2 new (rot, tilt, psi) triplets.
     * An oversampling_order == 2 will give rise to 4*4*4 new (rot, tilt, psi) triplets.
     * etc.
     *
     * If only_nonzero_prior is true, then only the orientations with non-zero prior probabilities will be returned
     * This is for local angular searches
     */
    void getOrientations(long int idir, long int ipsi, int oversampling_order, std::vector<Matrix1D<double> > &my_orientations);

    /** Get the prior probability for this orientation
     */
    double getPriorProbability(long int idir, long int ipsi);

    /** Get the number of the original direction from the ones with also non-zero prior probability
     */
    long int getDirectionNumberAlsoZeroPrior(long int idir);

    /** Get the number of the original direction from the ones with also non-zero prior probability
     */
    long int getPsiNumberAlsoZeroPrior(long int ipsi);

    /* Gets the vector of psi angles for a more finely (oversampled) sampling and
     * pushes each instance back into the oversampled_orientations vector with the given rot and tilt
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same psi angle as the original ipsi.
     * An oversampling_order == 1 will give rise to 2 new psi angles
     * An oversampling_order == 2 will give rise to 4 new psi angles
     * etc.
     */
    void pushbackOversampledPsiAngles(long int ipsi, int oversampling_order,
    		double rot, double tilt, std::vector<Matrix1D<double> > &oversampled_orientations);

    /* Calculate an angular distance between two sets of Euler angles */
    double calculateAngularDistance(double rot1, double tilt1, double psi1,
    		double rot2, double tilt2, double psi2);

    /* Write a BILD file describing the angular distribution
     *  R determines the radius of the sphere on which cylinders will be placed
     *  Rmax_frac determines the length of the longest cylinder (relative to R, 0.2 + +20%)
     *  width_frac determines how broad each cylinder is. frac=1 means they touch each other
     * */
    void writeBildFileOrientationalDistribution(MultidimArray<double> &pdf_direction,
    		FileName &fn_bild, double R, double offset = 0., double Rmax_frac = 0.3, double width_frac = 0.5);

private:

    /* Eliminate points from the sampling_points_vector and sampling_points_angles vectors
     * that are outside the allowed tilt range.
     * Let tilt angles range from -90 to 90, then:
     * if (limit_tilt > 0) then top views (that is with ABS(tilt) > limit_tilt) are removed and side views are kept
     * if (limit_tilt < 0) then side views (that is with ABS(tilt) < limit_tilt) are removed and top views are kept
    */
    void removePointsOutsideLimitedTiltAngles();

    /* Eliminate symmetry-equivalent points from the sampling_points_vector and sampling_points_angles vectors
        This function first calls removeSymmetryEquivalentPointsGeometric,
        and then checks each point versus all others to calculate an angular distance
        If this distance is less than 0.8 times the angular sampling, the point is deleted
        This cares care of sampling points near the edge of the geometrical considerations
    */
    void removeSymmetryEquivalentPoints(double max_ang);

    /* eliminate symmetry-related points based on simple geometrical considerations,
        symmetry group, symmetry order */
    void removeSymmetryEquivalentPointsGeometric(const int symmetry, int sym_order,
												 std::vector <Matrix1D<double> >  &sampling_points_vector);



};
//@}
#endif
