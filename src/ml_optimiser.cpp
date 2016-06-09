/***************************************************************************
 *
 * Author: "Callum Smits"
 * Victor Chang Cardiac Research Institute
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
#include "src/ml_optimiser.h"
#include <sys/stat.h>
//#define DEBUG
//#define DEBUG_CHECKSIZES
//#define CHECKSIZES
//Some global threads management variables
Mutex global_mutex, global_mutex2;
Barrier * global_barrier;
ThreadManager * global_ThreadManager;


// Global functions to work with threads
void globalGetFourierTransformsAndCtfs(ThreadArgument &thArg)
{
	((MlOptimiser*)thArg.workClass)->doThreadGetFourierTransformsAndCtfs(thArg.thread_id);
}

void globalThreadPrecalculateShiftedImagesCtfsAndInvSigma2s(ThreadArgument &thArg)
{
	((MlOptimiser*)thArg.workClass)->doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s(thArg.thread_id);
}

void globalThreadGetSquaredDifferencesAllOrientations(ThreadArgument &thArg)
{
	((MlOptimiser*)thArg.workClass)->doThreadGetSquaredDifferencesAllOrientations(thArg.thread_id);
}

void globalThreadConvertSquaredDifferencesToWeightsAllOrientations(ThreadArgument &thArg)
{
	((MlOptimiser*)thArg.workClass)->doThreadConvertSquaredDifferencesToWeightsAllOrientations(thArg.thread_id);
}

void globalThreadStoreWeightedSumsAllOrientations(ThreadArgument &thArg)
{
	((MlOptimiser*)thArg.workClass)->doThreadStoreWeightedSumsAllOrientations(thArg.thread_id);
}


/** ========================== I/O operations  =========================== */


void MlOptimiser::usage()
{

	parser.writeUsage(std::cerr);
}

void MlOptimiser::read(int argc, char **argv, int rank)
{
//#define DEBUG_READ

	parser.setCommandLine(argc, argv);

	if (checkParameter(argc, argv, "--continue"))
	{
		parser.addSection("Continue options");
		FileName fn_in = parser.getOption("--continue", "_optimiser.star file of the iteration after which to continue");
		// Read in previously calculated parameters
		if (fn_in != "")
			read(fn_in, rank);
		// And look for additional command-line options...
		parseContinue(argc, argv);
	}
	else
	{
		// Start a new run from scratch
		parseInitial(argc, argv);
	}

}

void MlOptimiser::parseContinue(int argc, char **argv)
{
#ifdef DEBUG
	std::cerr << "Entering parseContinue" << std::endl;
#endif

	int general_section = parser.addSection("General options");
	// Not all parameters are accessible here...
	FileName fn_out_new = parser.getOption("--o", "Output rootname", "OLD_ctX");
	if (fn_out_new == "OLD_ctX" || fn_out_new == fn_out )
		fn_out += "_ct" + integerToString(iter);
	else
		fn_out = fn_out_new;

	std::string fnt;
	fnt = parser.getOption("--iter", "Maximum number of iterations to perform", "OLD");
	if (fnt != "OLD")
		nr_iter = textToInteger(fnt);

	fnt = parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "OLD");
	if (fnt != "OLD")
		mymodel.tau2_fudge_factor = textToFloat(fnt);

	// Solvent flattening
	if (parser.checkOption("--flatten_solvent", "Switch on masking on the references?", "OLD"))
		do_solvent = true;

	// Check whether the mask has changed
	fnt = parser.getOption("--solvent_mask", "User-provided mask for the references", "OLD");
	if (fnt != "OLD")
		fn_mask = fnt;

	// Check whether the secondary mask has changed
	fnt = parser.getOption("--solvent_mask2", "User-provided secondary mask", "OLD");
	if (fnt != "OLD")
		fn_mask2 = fnt;

	// Check whether tau2-spectrum has changed
	fnt = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "OLD");
	if (fnt != "OLD")
		fn_tau = fnt;

	// Check whether particle diameter has changed
	fnt = parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "OLD");
	if (fnt != "OLD")
		particle_diameter = textToFloat(fnt);

	// Check whether to join the random halves again
	do_join_random_halves = parser.checkOption("--join_random_halves", "Join previously split random halves again (typically to perform a final reconstruction).");

	// Re-align movie frames
	int movie_section = parser.addSection("Re-align movie frames");

	fn_data_movie = parser.getOption("--realign_movie_frames", "Input STAR file with the movie frames", "");

	// TODO: add this to EMDL_OPTIMISER and read/write of optimiser.star
	nr_frames_per_prior = textToInteger(parser.getOption("--nr_frames_prior", "Number of movie frames to calculate running-average priors", "5"));

	// (integer-) divide running average width by 2 to have the side only
	// TODO: add this to EMDL_OPTIMISER and read/write of optimiser.star
	movie_frame_running_avg_side = textToInteger(parser.getOption("--movie_frames_running_avg", "Number of movie frames in each running average", "3")) / 2;

	// ORIENTATIONS
	int orientations_section = parser.addSection("Orientations");

	fnt = parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "OLD");
	if (fnt != "OLD")
		adaptive_oversampling = textToInteger(fnt);

	// Check whether angular sampling has changed
	// Do not do this for auto_refine, but make sure to do this when realigning movies!
	if (!do_auto_refine || fn_data_movie != "")
	{
		directions_have_changed = false;
		fnt = parser.getOption("--healpix_order", "Healpix order for the angular sampling rate on the sphere (before oversampling): hp2=15deg, hp3=7.5deg, etc", "OLD");
		if (fnt != "OLD")
		{
			int _order = textToInteger(fnt);
			if (_order != sampling.healpix_order)
			{
				directions_have_changed = true;
				sampling.healpix_order = _order;
			}
		}

		fnt = parser.getOption("--psi_step", "Angular sampling (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "OLD");
		if (fnt != "OLD")
			sampling.psi_step = textToFloat(fnt);

		fnt = parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "OLD");
		if (fnt != "OLD")
			sampling.offset_range = textToFloat(fnt);

		fnt = parser.getOption("--offset_step", "Sampling rate for origin offsets (in pixels)", "OLD");
		if (fnt != "OLD")
			sampling.offset_step = textToFloat(fnt);
	}

	fnt = parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which auto-refine procedure will use local searches", "OLD");
	if (fnt != "OLD")
		autosampling_hporder_local_searches = textToInteger(fnt);

	// Check whether the prior mode changes
	double _sigma_rot, _sigma_tilt, _sigma_psi, _sigma_off;
	int _mode;
	fnt = parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_tilt", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_tilt = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_off", "Stddev. on the translations", "OLD");
	if (fnt != "OLD")
	{
		mymodel.sigma2_offset = textToFloat(fnt) * textToFloat(fnt);
	}

	if (parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?"))
		do_skip_align = true;

	if (parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?"))
		do_skip_rotate = true;
	else
		do_skip_rotate = false; // do_skip_rotate should normally be false...

	do_skip_maximization = parser.checkOption("--skip_maximize", "Skip maximization step (only write out data.star file)?");

	int corrections_section = parser.addSection("Corrections");

	// Can only switch the following option ON, not OFF
	if (parser.checkOption("--scale", "Switch on intensity-scale corrections on image groups", "OLD"))
		do_scale_correction = true;

	// Can only switch the following option ON, not OFF
	if (parser.checkOption("--norm", "Switch on normalisation-error correction","OLD"))
		do_norm_correction = true;

	int computation_section = parser.addSection("Computation");

	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));

	fnt = parser.getOption("--pool", "Number of images to be processed together", "OLD");
	if (fnt != "OLD")
		max_nr_pool = textToInteger(fnt);

	combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");

	verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));

	int expert_section = parser.addSection("Expert options");

	fnt = parser.getOption("--strict_highres_exp", "Resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "OLD");
	if (fnt != "OLD")
		strict_highres_exp = textToFloat(fnt);

	// Debugging/analysis/hidden stuff
	do_map = !checkParameter(argc, argv, "--no_map");
	minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
	debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0."));
	debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0."));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
	// Read in initial sigmaNoise spectrum
	fn_sigma = getParameter(argc, argv, "--sigma","");
	sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1."));
	do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
	do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
	do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
	do_use_all_data = checkParameter(argc, argv, "--use_all_data");
	do_always_cc  = checkParameter(argc, argv, "--always_cc");

	do_print_metadata_labels = false;
	do_print_symmetry_ops = false;
#ifdef DEBUG
	std::cerr << "Leaving parseContinue" << std::endl;
#endif

}

void MlOptimiser::parseInitial(int argc, char **argv)
{
#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Entering "<<std::endl;
#endif

	// Read/initialise mymodel and sampling from a STAR file
    FileName fn_model = getParameter(argc, argv, "--model", "None");
	if (fn_model != "None")
	{
		mymodel.read(fn_model);
	}
	// Read in the sampling information from a _sampling.star file
    FileName fn_sampling = getParameter(argc, argv, "--sampling", "None");
	if (fn_sampling != "None")
	{
		sampling.read(fn_sampling);
	}

	// General optimiser I/O stuff
    int general_section = parser.addSection("General options");
    fn_data = parser.getOption("--i", "Input images (in a star-file or a stack)");
    fn_out = parser.getOption("--o", "Output rootname");
    nr_iter = textToInteger(parser.getOption("--iter", "Maximum number of iterations to perform", "50"));
	mymodel.pixel_size = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)"));
	mymodel.tau2_fudge_factor = textToFloat(parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "1"));
	mymodel.nr_classes = textToInteger(parser.getOption("--K", "Number of references to be refined", "1"));
    particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "-1"));
	do_zero_mask = parser.checkOption("--zero_mask","Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)");
	do_solvent = parser.checkOption("--flatten_solvent", "Perform masking on the references as well?");
	fn_mask = parser.getOption("--solvent_mask", "User-provided mask for the references (default is to use spherical mask with particle_diameter)", "None");
	fn_mask2 = parser.getOption("--solvent_mask2", "User-provided secondary mask (with its own average density)", "None");
	fn_tau = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "None");
	do_split_random_halves = parser.checkOption("--split_random_halves", "Refine two random halves of the data completely separately");
	low_resol_join_halves = textToFloat(parser.getOption("--low_resol_join_halves", "Resolution (in Angstrom) up to which the two random half-reconstructions will not be independent to prevent diverging orientations","-1"));

	// Initialisation
	int init_section = parser.addSection("Initialisation");
	fn_ref = parser.getOption("--ref", "Image, stack or star-file with the reference(s). (Compulsory for 3D refinement!)", "None");
	mymodel.sigma2_offset = textToFloat(parser.getOption("--offset", "Initial estimated stddev for the origin offsets", "3"));
	mymodel.sigma2_offset *= mymodel.sigma2_offset;

	// Perform cross-product comparison at first iteration
	do_firstiter_cc = parser.checkOption("--firstiter_cc", "Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)");
	ini_high = textToFloat(parser.getOption("--ini_high", "Resolution (in Angstroms) to which to limit refinement in the first iteration ", "-1"));

	// Set the orientations
    int orientations_section = parser.addSection("Orientations");
	// Move these to sampling
	adaptive_oversampling = textToInteger(parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "1"));
	sampling.healpix_order = textToInteger(parser.getOption("--healpix_order", "Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc", "2"));
	sampling.psi_step = textToFloat(parser.getOption("--psi_step", "Sampling rate (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "-1"));
	sampling.limit_tilt = textToFloat(parser.getOption("--limit_tilt", "Limited tilt angle: positive for keeping side views, negative for keeping top views", "-91"));
	sampling.fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	sampling.offset_range = textToFloat(parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "6"));
	sampling.offset_step = textToFloat(parser.getOption("--offset_step", "Sampling rate (before oversampling) for origin offsets (in pixels)", "2"));
	sampling.perturbation_factor = textToFloat(parser.getOption("--perturb", "Perturbation factor for the angular sampling (0=no perturb; 0.5=perturb)", "0.5"));
	do_auto_refine = parser.checkOption("--auto_refine", "Perform 3D auto-refine procedure?");
	autosampling_hporder_local_searches = textToInteger(parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which autosampling procedure will use local searches", "4"));
	parser.setSection(orientations_section);
	double _sigma_ang = textToFloat(parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "-1"));
	double _sigma_rot = textToFloat(parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
	double _sigma_tilt = textToFloat(parser.getOption("--sigma_tilt", "Stddev on the second Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
	double _sigma_psi = textToFloat(parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "-1"));
	if (_sigma_ang > 0.)
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		// the sigma-values for the orientational prior are in model (and not in sampling) because one might like to estimate them
		// from the data by calculating weighted sums of all angular differences: therefore it needs to be in wsum_model and thus in mymodel.
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = _sigma_ang * _sigma_ang;
	}
	else if (_sigma_rot > 0. || _sigma_tilt > 0. || _sigma_psi > 0.)
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot  = (_sigma_rot > 0. ) ? _sigma_rot * _sigma_rot   : 0.;
		mymodel.sigma2_tilt = (_sigma_tilt > 0.) ? _sigma_tilt * _sigma_tilt : 0.;
		mymodel.sigma2_psi  = (_sigma_psi > 0. ) ? _sigma_psi * _sigma_psi   : 0.;
	}
	else
	{
		//default
		mymodel.orientational_prior_mode = NOPRIOR;
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 0.;
	}
	do_skip_align = parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?");
	do_skip_rotate = parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?");
	do_skip_maximization = false;

	// CTF, norm, scale, bfactor correction etc.
	int corrections_section = parser.addSection("Corrections");
	do_ctf_correction = parser.checkOption("--ctf", "Perform CTF correction?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
	refs_are_ctf_corrected = parser.checkOption("--ctf_corrected_ref", "Have the input references been CTF-amplitude corrected?");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Have the data been CTF phase-flipped?");
	only_flip_phases = parser.checkOption("--only_flip_phases", "Only perform CTF phase-flipping? (default is full amplitude-correction)");
	do_norm_correction = parser.checkOption("--norm", "Perform normalisation-error correction?");
	do_scale_correction = parser.checkOption("--scale", "Perform intensity-scale corrections on image groups?");

	// Computation stuff
	// The number of threads is always read from the command line
	int computation_section = parser.addSection("Computation");
	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));
	available_memory = textToFloat(parser.getOption("--memory_per_thread", "Available RAM (in Gb) for each thread", "2"));
	max_nr_pool = textToInteger(parser.getOption("--pool", "Number of images to be processed together", "8"));
	combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");

	// Expert options
	int expert_section = parser.addSection("Expert options");
	mymodel.padding_factor = textToInteger(parser.getOption("--pad", "Oversampling factor for the Fourier transforms of the references", "2"));
	mymodel.interpolator = (parser.checkOption("--NN", "Perform nearest-neighbour instead of linear Fourier-space interpolation?")) ? NEAREST_NEIGHBOUR : TRILINEAR;
	mymodel.r_min_nn = textToInteger(parser.getOption("--r_min_nn", "Minimum number of Fourier shells to perform linear Fourier-space interpolation", "10"));
	verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));
	random_seed = textToInteger(parser.getOption("--random_seed", "Number for the random seed generator", "-1"));
	max_coarse_size = textToInteger(parser.getOption("--coarse_size", "Maximum image size for the first pass of the adaptive sampling approach", "-1"));
	adaptive_fraction = textToFloat(parser.getOption("--adaptive_fraction", "Fraction of the weights to be considered in the first pass of adaptive oversampling ", "0.999"));
	width_mask_edge = textToInteger(parser.getOption("--maskedge", "Width of the soft edge of the spherical mask (in pixels)", "5"));
	fix_sigma_noise = parser.checkOption("--fix_sigma_noise", "Fix the experimental noise spectra?");
	fix_sigma_offset = parser.checkOption("--fix_sigma_offset", "Fix the stddev in the origin offsets?");
	incr_size = textToInteger(parser.getOption("--incr_size", "Number of Fourier shells beyond the current resolution to be included in refinement", "10"));
	do_print_metadata_labels = parser.checkOption("--print_metadata_labels", "Print a table with definitions of all metadata labels, and exit");
	do_print_symmetry_ops = parser.checkOption("--print_symmetry_ops", "Print all symmetry transformation matrices, and exit");
	strict_highres_exp = textToFloat(parser.getOption("--strict_highres_exp", "Resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "-1"));
	dont_raise_norm_error = parser.checkOption("--dont_check_norm", "Skip the check whether the images are normalised correctly");



	// TODO: read/write do_always_cc in optmiser.star file!!!
	// SA-stuff
	do_sim_anneal = parser.checkOption("--sim_anneal", "Perform simulated-annealing to improve overall convergence of random starting models?");
	temp_ini = textToFloat(parser.getOption("--temp_ini", "Initial temperature (K) for simulated annealing", "1000"));
	temp_fin = textToFloat(parser.getOption("--temp_fin", "Initial temperature (K) for simulated annealing", "1"));
	do_always_cc  = parser.checkOption("--always_cc", "Perform CC-calculation in all iterations (useful for faster denovo model generation?)");

	///////////////// Special stuff for first iteration (only accessible via CL, not through readSTAR ////////////////////

	// When reading from the CL: always start at iteration 1
	iter = 0;
    // When starting from CL: always calculate initial sigma_noise
    do_calculate_initial_sigma_noise = true;
    // Start average norm correction at 1!
    mymodel.avg_norm_correction = 1.;
    // Always initialise the PDF of the directions
    directions_have_changed = true;

    // Only reconstruct and join random halves are only available when continuing an old run
    do_join_random_halves = false;

    // For auto-sampling and convergence check
    nr_iter_wo_resol_gain = 0;
    nr_iter_wo_large_hidden_variable_changes = 0;
    current_changes_optimal_classes = 9999999;
    current_changes_optimal_offsets = 999.;
    current_changes_optimal_orientations = 999.;
    smallest_changes_optimal_classes = 9999999;
    smallest_changes_optimal_offsets = 999.;
    smallest_changes_optimal_orientations = 999.;
    acc_rot = acc_trans = 999.;

    best_resol_thus_far = 1./999.;
    has_converged = false;
    has_high_fsc_at_limit = false;
    has_large_incr_size_iter_ago = 0;

    // Never realign movies from the start
    do_realign_movies = false;

    // Debugging/analysis/hidden stuff
	do_map = !checkParameter(argc, argv, "--no_map");
	minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
	debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0"));
	debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0"));
	// Read in initial sigmaNoise spectrum
	fn_sigma = getParameter(argc, argv, "--sigma","");
	do_calculate_initial_sigma_noise = (fn_sigma == "") ? true : false;
	sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1"));
	do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
	do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
	do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
	do_use_all_data = checkParameter(argc, argv, "--use_all_data");

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Done"<<std::endl;
#endif

}


void MlOptimiser::read(FileName fn_in, int rank)
{

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar entering ..."<<std::endl;
#endif

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MlOptimiser::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MD;

    // Read general stuff
    FileName fn_model, fn_model2, fn_sampling;
    MD.readStar(in, "optimiser_general");
    in.close();

    if (!MD.getValue(EMDL_OPTIMISER_OUTPUT_ROOTNAME, fn_out) ||
        !MD.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model) ||
		!MD.getValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data) ||
		!MD.getValue(EMDL_OPTIMISER_SAMPLING_STARFILE, fn_sampling) ||
        !MD.getValue(EMDL_OPTIMISER_ITERATION_NO, iter) ||
        !MD.getValue(EMDL_OPTIMISER_NR_ITERATIONS, nr_iter) ||
        !MD.getValue(EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES, do_split_random_halves) ||
        !MD.getValue(EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, low_resol_join_halves) ||
        !MD.getValue(EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING, adaptive_oversampling) ||
		!MD.getValue(EMDL_OPTIMISER_ADAPTIVE_FRACTION, adaptive_fraction) ||
		!MD.getValue(EMDL_OPTIMISER_RANDOM_SEED, random_seed) ||
		!MD.getValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, particle_diameter) ||
		!MD.getValue(EMDL_OPTIMISER_WIDTH_MASK_EDGE, width_mask_edge) ||
		!MD.getValue(EMDL_OPTIMISER_DO_ZERO_MASK, do_zero_mask) ||
		!MD.getValue(EMDL_OPTIMISER_DO_SOLVENT_FLATTEN, do_solvent) ||
		!MD.getValue(EMDL_OPTIMISER_SOLVENT_MASK_NAME, fn_mask) ||
		!MD.getValue(EMDL_OPTIMISER_SOLVENT_MASK2_NAME, fn_mask2) ||
		!MD.getValue(EMDL_OPTIMISER_TAU_SPECTRUM_NAME, fn_tau) ||
		!MD.getValue(EMDL_OPTIMISER_COARSE_SIZE, coarse_size) ||
		!MD.getValue(EMDL_OPTIMISER_MAX_COARSE_SIZE, max_coarse_size) ||
		!MD.getValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, strict_highres_exp) ||
		!MD.getValue(EMDL_OPTIMISER_INCR_SIZE, incr_size) ||
		!MD.getValue(EMDL_OPTIMISER_DO_MAP, do_map) ||
		!MD.getValue(EMDL_OPTIMISER_DO_AUTO_REFINE, do_auto_refine) ||
		!MD.getValue(EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER, autosampling_hporder_local_searches) ||
	    !MD.getValue(EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN, nr_iter_wo_resol_gain) ||
	    !MD.getValue(EMDL_OPTIMISER_BEST_RESOL_THUS_FAR, best_resol_thus_far) ||
	    !MD.getValue(EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, nr_iter_wo_large_hidden_variable_changes) ||
		!MD.getValue(EMDL_OPTIMISER_DO_SKIP_ALIGN, do_skip_align) ||
		//!MD.getValue(EMDL_OPTIMISER_DO_SKIP_ROTATE, do_skip_rotate) ||
	    !MD.getValue(EMDL_OPTIMISER_ACCURACY_ROT, acc_rot) ||
	    !MD.getValue(EMDL_OPTIMISER_ACCURACY_TRANS, acc_trans) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, current_changes_optimal_orientations) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, current_changes_optimal_offsets) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES, current_changes_optimal_classes) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, smallest_changes_optimal_orientations) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, smallest_changes_optimal_offsets) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, smallest_changes_optimal_classes) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_CONVERGED, has_converged) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, has_high_fsc_at_limit) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, has_large_incr_size_iter_ago) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_NORM, do_norm_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_SCALE, do_scale_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_CTF, do_ctf_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_REALIGN_MOVIES, do_realign_movies) ||
		!MD.getValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak) ||
		!MD.getValue(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped) ||
		!MD.getValue(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases) ||
		!MD.getValue(EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED, refs_are_ctf_corrected) ||
		!MD.getValue(EMDL_OPTIMISER_FIX_SIGMA_NOISE, fix_sigma_noise) ||
		!MD.getValue(EMDL_OPTIMISER_FIX_SIGMA_OFFSET, fix_sigma_offset) ||
		!MD.getValue(EMDL_OPTIMISER_MAX_NR_POOL, max_nr_pool) ||
		!MD.getValue(EMDL_OPTIMISER_AVAILABLE_MEMORY, available_memory))
    	REPORT_ERROR("MlOptimiser::readStar: incorrect optimiser_general table");

    if (do_split_random_halves &&
    		!MD.getValue(EMDL_OPTIMISER_MODEL_STARFILE2, fn_model2))
    	REPORT_ERROR("MlOptimiser::readStar: splitting data into two random halves, but rlnModelStarFile2 not found in optimiser_general table");

    // Initialise some stuff for first-iteration only (not relevant here...)
    do_calculate_initial_sigma_noise = false;
    do_average_unaligned = false;
    do_generate_seeds = false;
    do_firstiter_cc = false;
    ini_high = 0;

    // Initialise some of the other, hidden or debugging stuff
    minres_map = 5;
    do_bfactor = false;
    gridding_nr_iter = 10;
    debug1 = debug2 = 0.;

    // Then read in sampling, mydata and mymodel stuff
    mydata.read(fn_data);
    if (do_split_random_halves)
    {
		if (rank % 2 == 1)
			mymodel.read(fn_model);
		else
			mymodel.read(fn_model2);
    }
    else
    {
    	mymodel.read(fn_model);
    }
	sampling.read(fn_sampling);

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar done."<<std::endl;
#endif

}


void MlOptimiser::write(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model, int random_subset)
{

	FileName fn_root, fn_tmp, fn_model, fn_model2, fn_data, fn_sampling;
	std::ofstream  fh;
	if (iter > -1)
		fn_root.compose(fn_out+"_it", iter, "", 3);
	else
		fn_root = fn_out;

	// First write "main" STAR file with all information from this run
	// Do this for random_subset==0 and random_subset==1
	if (do_write_optimiser && random_subset < 2)
	{
		fn_tmp = fn_root+"_optimiser.star";
		fh.open((fn_tmp).c_str(), std::ios::out);
		if (!fh)
			REPORT_ERROR( (std::string)"MlOptimiser::write: Cannot write file: " + fn_tmp);

		// Write the command line as a comment in the header
		fh << "# RELION optimiser"<<std::endl;
		fh << "# ";
		parser.writeCommandLine(fh);

		if (do_split_random_halves && !do_join_random_halves)
		{
			fn_model  = fn_root + "_half1_model.star";
			fn_model2 = fn_root + "_half2_model.star";
		}
		else
		{
			fn_model = fn_root + "_model.star";
		}
		fn_data = fn_root + "_data.star";
		fn_sampling = fn_root + "_sampling.star";

		MetaDataTable MD;
		MD.setIsList(true);
		MD.setName("optimiser_general");
		MD.addObject();
		MD.setValue(EMDL_OPTIMISER_OUTPUT_ROOTNAME, fn_out);
		if (do_split_random_halves)
		{
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE2, fn_model2);
		}
		else
		{
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
		}
		MD.setValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data);
		MD.setValue(EMDL_OPTIMISER_SAMPLING_STARFILE, fn_sampling);
		MD.setValue(EMDL_OPTIMISER_ITERATION_NO, iter);
		MD.setValue(EMDL_OPTIMISER_NR_ITERATIONS, nr_iter);
		MD.setValue(EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES, do_split_random_halves);
		MD.setValue(EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, low_resol_join_halves);
		MD.setValue(EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING, adaptive_oversampling);
		MD.setValue(EMDL_OPTIMISER_ADAPTIVE_FRACTION, adaptive_fraction);
		MD.setValue(EMDL_OPTIMISER_RANDOM_SEED, random_seed);
		MD.setValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, particle_diameter);
		MD.setValue(EMDL_OPTIMISER_WIDTH_MASK_EDGE, width_mask_edge);
		MD.setValue(EMDL_OPTIMISER_DO_ZERO_MASK, do_zero_mask);
		MD.setValue(EMDL_OPTIMISER_DO_SOLVENT_FLATTEN, do_solvent);
		MD.setValue(EMDL_OPTIMISER_SOLVENT_MASK_NAME, fn_mask);
		MD.setValue(EMDL_OPTIMISER_SOLVENT_MASK2_NAME, fn_mask2);
		MD.setValue(EMDL_OPTIMISER_TAU_SPECTRUM_NAME, fn_tau);
		MD.setValue(EMDL_OPTIMISER_COARSE_SIZE, coarse_size);
		MD.setValue(EMDL_OPTIMISER_MAX_COARSE_SIZE, max_coarse_size);
		MD.setValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, strict_highres_exp);
		MD.setValue(EMDL_OPTIMISER_INCR_SIZE, incr_size);
		MD.setValue(EMDL_OPTIMISER_DO_MAP, do_map);
		MD.setValue(EMDL_OPTIMISER_DO_AUTO_REFINE, do_auto_refine);
		MD.setValue(EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER, autosampling_hporder_local_searches);
	    MD.setValue(EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN, nr_iter_wo_resol_gain);
	    MD.setValue(EMDL_OPTIMISER_BEST_RESOL_THUS_FAR,best_resol_thus_far);
	    MD.setValue(EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, nr_iter_wo_large_hidden_variable_changes);
		MD.setValue(EMDL_OPTIMISER_DO_SKIP_ALIGN, do_skip_align);
		MD.setValue(EMDL_OPTIMISER_DO_SKIP_ROTATE, do_skip_rotate);
	    MD.setValue(EMDL_OPTIMISER_ACCURACY_ROT, acc_rot);
	    MD.setValue(EMDL_OPTIMISER_ACCURACY_TRANS, acc_trans);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, current_changes_optimal_orientations);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, current_changes_optimal_offsets);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES, current_changes_optimal_classes);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, smallest_changes_optimal_orientations);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, smallest_changes_optimal_offsets);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, smallest_changes_optimal_classes);
	    MD.setValue(EMDL_OPTIMISER_HAS_CONVERGED, has_converged);
	    MD.setValue(EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, has_high_fsc_at_limit);
	    MD.setValue(EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, has_large_incr_size_iter_ago);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_NORM, do_norm_correction);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_SCALE, do_scale_correction);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_CTF, do_ctf_correction);
		MD.setValue(EMDL_OPTIMISER_DO_REALIGN_MOVIES, do_realign_movies);
		MD.setValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak);
		MD.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped);
		MD.setValue(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases);
		MD.setValue(EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED, refs_are_ctf_corrected);
		MD.setValue(EMDL_OPTIMISER_FIX_SIGMA_NOISE, fix_sigma_noise);
		MD.setValue(EMDL_OPTIMISER_FIX_SIGMA_OFFSET, fix_sigma_offset);
		MD.setValue(EMDL_OPTIMISER_MAX_NR_POOL, max_nr_pool);
		MD.setValue(EMDL_OPTIMISER_AVAILABLE_MEMORY, available_memory);

		MD.write(fh);
		fh.close();
	}

	// Then write the mymodel to file
	if (do_write_model)
	{
		if (do_split_random_halves && !do_join_random_halves)
			mymodel.write(fn_root + "_half" + integerToString(random_subset), sampling);
		else
			mymodel.write(fn_root, sampling);
	}

	// And write the mydata to file
	if (do_write_data)
		mydata.write(fn_root);

	// And write the sampling object
	if (do_write_sampling)
		sampling.write(fn_root);

}

char * MlOptimiser::loadProgramSource(const char *filename)
{
    struct stat statbuf;
    FILE        *fh;
    char        *source;
    
    fh = fopen(filename, "r");
    if (fh == 0)
        return 0;
    
    stat(filename, &statbuf);
    source = (char *) malloc(statbuf.st_size + 1);
    fread(source, statbuf.st_size, 1, fh);
    source[statbuf.st_size] = '\0';
    
    return source;
}

/** ========================== Initialisation  =========================== */

void MlOptimiser::initialise()
{
#ifdef DEBUG
    std::cerr<<"MlOptimiser::initialise Entering"<<std::endl;
#endif

    initialiseGeneral();

    initialiseWorkLoad();

	if (fn_sigma != "")
	{
		// Read in sigma_noise spetrum from file DEVELOPMENTAL!!! FOR DEBUGGING ONLY....
		MetaDataTable MDsigma;
		double val;
		int idx;
		MDsigma.read(fn_sigma);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
		{
			MDsigma.getValue(EMDL_SPECTRAL_IDX, idx);
			MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, val);
			if (idx < XSIZE(mymodel.sigma2_noise[0]))
				mymodel.sigma2_noise[0](idx) = val;
		}
		if (idx < XSIZE(mymodel.sigma2_noise[0]) - 1)
		{
			if (verb > 0) std::cout<< " WARNING: provided sigma2_noise-spectrum has fewer entries ("<<idx+1<<") than needed ("<<XSIZE(mymodel.sigma2_noise[0])<<"). Set rest to zero..."<<std::endl;
		}
		// Use the same spectrum for all classes
		for (int igroup = 0; igroup< mymodel.nr_groups; igroup++)
			mymodel.sigma2_noise[igroup] =  mymodel.sigma2_noise[0];

	}
	else if (do_calculate_initial_sigma_noise || do_average_unaligned)
	{
		MultidimArray<double> Mavg;

		// Calculate initial sigma noise model from power_class spectra of the individual images
		calculateSumOfPowerSpectraAndAverageImage(Mavg);

		// Set sigma2_noise and Iref from averaged poser spectra and Mavg
		setSigmaNoiseEstimatesAndSetAverageImage(Mavg);
	}

    // First low-pass filter the initial references
	if (iter == 0)
		initialLowPassFilterReferences();

	// Initialise the data_versus_prior ratio to get the initial current_size right
	if (iter == 0)
		mymodel.initialiseDataVersusPrior(fix_tau); // fix_tau was set in initialiseGeneral

	// Check minimum group size of 10 particles
	if (verb > 0)
	{
		bool do_warn = false;
		for (int igroup = 0; igroup< mymodel.nr_groups; igroup++)
		{
			if (mymodel.nr_particles_group[igroup] < 10)
			{
				std:: cout << "WARNING: There are only " << mymodel.nr_particles_group[igroup] << " particles in group " << igroup + 1 << std::endl;
				do_warn = true;
			}
		}
		if (do_warn)
		{
			std:: cout << "WARNING: You may want to consider joining some micrographs into larger groups to obtain more robust noise estimates. " << std::endl;
			std:: cout << "         You can do so by using the same rlnMicrographName label for particles from multiple different micrographs in the input STAR file. " << std::endl;
		}
	}

	// Write out initial mymodel
	write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);


	// Do this after writing out the model, so that still the random halves are written in separate files.
	if (do_realign_movies)
	{
		// Resolution seems to decrease again after 1 iteration. Therefore, just perform a single iteration until we figure out what exactly happens here...
		has_converged = true;
		// Then use join random halves
		do_join_random_halves = true;

		// If we skip the maximization step, then there is no use in using all data
		if (!do_skip_maximization)
		{
			// Use all data out to Nyquist because resolution gains may be substantial
			do_use_all_data = true;
		}
	}

#ifdef DEBUG
    std::cerr<<"MlOptimiser::initialise Done"<<std::endl;
#endif
}

void MlOptimiser::initialiseGeneral(int rank)
{

#ifdef DEBUG
	std::cerr << "Entering initialiseGeneral" << std::endl;
#endif

#ifdef TIMING
	//DIFFF = timer.setNew("difff");
	TIMING_EXP =           timer.setNew("expectation");
	TIMING_MAX =           timer.setNew("maximization");
	TIMING_RECONS =        timer.setNew("reconstruction");
	TIMING_ESP =           timer.setNew("expectationSomeParticles");
	TIMING_ESP_READ  =     timer.setNew(" - ESP: read");
	TIMING_ESP_DIFF1 =     timer.setNew(" - ESP: getAllSquaredDifferences1");
	TIMING_ESP_DIFF2 =     timer.setNew(" - ESP: getAllSquaredDifferences2");
	TIMING_DIFF_PROJ =     timer.setNew(" -  - ESPdiff2: project");
	TIMING_DIFF_SHIFT =    timer.setNew(" -  - ESPdiff2: shift");
	TIMING_DIFF_DIFF2 =    timer.setNew(" -  - ESPdiff2: diff2");
	TIMING_ESP_WEIGHT1 =   timer.setNew(" - ESP: convertDiff2ToWeights1");
	TIMING_ESP_WEIGHT2 =   timer.setNew(" - ESP: convertDiff2ToWeights2");
	TIMING_WEIGHT_EXP =    timer.setNew(" -  - ESPweight: exp");
	TIMING_WEIGHT_SORT =   timer.setNew(" -  - ESPweight: sort");
	TIMING_ESP_WSUM =      timer.setNew(" - ESP: storeWeightedSums");
	TIMING_WSUM_PROJ =     timer.setNew("  - - ESPwsum: project");
	TIMING_WSUM_DIFF2 =    timer.setNew(" -  - ESPwsum: diff2");
	TIMING_WSUM_SUMSHIFT = timer.setNew(" -  - ESPwsum: shift");
	TIMING_WSUM_BACKPROJ = timer.setNew(" -  - ESPwsum: backproject");
#endif

    if (do_use_opencl) {
        initialiseCL(rank);
    }
    
	if (do_print_metadata_labels)
	{
		if (verb > 0)
			EMDL::printDefinitions(std::cout);
		exit(0);
	}

	// Print symmetry operators to cout
	if (do_print_symmetry_ops)
	{
		if (verb > 0)
		{
			SymList SL;
			SL.writeDefinition(std::cout, sampling.symmetryGroup());
		}
		exit(0);
	}

	// Check for errors in the command-line option
	if (parser.checkForErrors(verb))
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	nr_threads_original = nr_threads;
	// If we are not continuing an old run, now read in the data and the reference images
	if (iter == 0)
	{

		// Read in the experimental image metadata
		mydata.read(fn_data);

		// Also get original size of the images to pass to mymodel.read()
		int ori_size = -1;
		mydata.MDexp.getValue(EMDL_IMAGE_SIZE, ori_size);
		if (ori_size%2 != 0)
			REPORT_ERROR("This program only works with even values for the image dimensions!");
    	mymodel.readImages(fn_ref, ori_size, mydata,
    			do_average_unaligned, do_generate_seeds, refs_are_ctf_corrected);

    	// Check consistency of EMDL_CTF_MAGNIFICATION and MEBL_CTF_DETECTOR_PIXEL_SIZE with mymodel.pixel_size
    	double mag, dstep, first_angpix, my_angpix;
    	bool has_magn = false;
    	if (mydata.MDimg.containsLabel(EMDL_CTF_MAGNIFICATION) && mydata.MDimg.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mydata.MDimg)
			{
    			mydata.MDimg.getValue(EMDL_CTF_MAGNIFICATION, mag);
    			mydata.MDimg.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
    			my_angpix = 10000. * dstep / mag;
    			if (!has_magn)
    			{
    				first_angpix = my_angpix;
    				has_magn = true;
    			}
    			else if (ABS(first_angpix - my_angpix) > 0.01)
    				REPORT_ERROR("MlOptimiser::initialiseGeneral: ERROR inconsistent magnification and detector pixel sizes in images in input STAR file");
			}
    	}
    	if (mydata.MDmic.containsLabel(EMDL_CTF_MAGNIFICATION) && mydata.MDmic.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mydata.MDmic)
			{
    			mydata.MDimg.getValue(EMDL_CTF_MAGNIFICATION, mag);
    			mydata.MDimg.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
    			my_angpix = 10000. * dstep / mag;
    			if (!has_magn)
    			{
    				first_angpix = my_angpix;
    				has_magn = true;
    			}
    			else if (ABS(first_angpix - my_angpix) > 0.01)
    				REPORT_ERROR("MlOptimiser::initialiseGeneral: ERROR inconsistent magnification and detector pixel sizes in micrographs in input STAR file");
			}
    	}
    	if (has_magn && ABS(first_angpix - mymodel.pixel_size) > 0.01)
    	{
    		if (verb > 0)
    			std::cout << "MlOptimiser::initialiseGeneral: WARNING modifying pixel size from " << mymodel.pixel_size <<" to "<<first_angpix << " based on magnification information in the input STAR file" << std::endl;
    		mymodel.pixel_size = first_angpix;
    	}

	}
	// Expand movies if fn_data_movie is given AND we were not doing expanded movies already
	else if (fn_data_movie != "" && !do_realign_movies)
	{

		if (verb > 0)
			std::cout << " Expanding current model for movie frames... " << std::endl;

		do_realign_movies = true;
		nr_iter_wo_resol_gain = -1;
		nr_iter_wo_large_hidden_variable_changes = 0;
		smallest_changes_optimal_offsets = 999.;
		smallest_changes_optimal_orientations = 999.;
		current_changes_optimal_orientations = 999.;
		current_changes_optimal_offsets = 999.;

		// If we're realigning movie frames, then now read in the metadata of the movie frames and combine with the metadata of the average images
		mydata.expandToMovieFrames(fn_data_movie);

		// Now also modify the model to contain many more groups....
		// each groups has to become Nframes groups (get Nframes from new mydata)
		mymodel.expandToMovieFrames(mydata, movie_frame_running_avg_side);

		// Don't do norm correction for realignment of movies.
		do_norm_correction = false;

	}

	if (mymodel.nr_classes > 1 && do_split_random_halves)
		REPORT_ERROR("ERROR: One cannot use --split_random_halves with more than 1 reference... You could first classify, and then refine each class separately using --random_halves.");

	if (do_join_random_halves && !do_split_random_halves)
		REPORT_ERROR("ERROR: cannot join random halves because they were not split in the previous run");

	if (do_always_join_random_halves)
		std::cout << " Joining half-reconstructions at each iteration: this is a developmental option to test sub-optimal FSC usage only! " << std::endl;

	// If fn_tau is provided, read in the tau spectrum
	fix_tau = false;
	if (fn_tau != "None")
	{
		fix_tau = true;
		mymodel.readTauSpectrum(fn_tau, verb);
	}


	// Initialise the sampling object (sets prior mode and fills translations and rotations inside sampling object)
	sampling.initialise(mymodel.orientational_prior_mode, mymodel.ref_dim, false);

	// Default max_coarse_size is original size
	if (max_coarse_size < 0)
		max_coarse_size = mymodel.ori_size;

	if (particle_diameter < 0.)
    	particle_diameter = (mymodel.ori_size - width_mask_edge) * mymodel.pixel_size;

	if (do_auto_refine)
	{
		nr_iter = 999;
		has_fine_enough_angular_sampling = false;
	}

    // For do_average_unaligned, always use initial low_pass filter
    if (do_average_unaligned && ini_high < 0.)
    {
    	// By default, use 0.07 dig.freq. low-pass filter
    	// See S.H.W. Scheres (2010) Meth Enzym.
    	ini_high = 1./mymodel.getResolution(ROUND(0.07 * mymodel.ori_size));
    }

    // Fill tabulated sine and cosine tables
    tab_sin.initialise(5000);
    tab_cos.initialise(5000);

	// For skipped alignments: set nr_pool to one to have each thread work on one particle (with its own unique sampling arrays of 1 orientation and translation)
	// Also do not perturb this orientation, nor do oversampling or priors
	if (do_skip_align || do_skip_rotate)
	{
		mymodel.orientational_prior_mode = NOPRIOR;
		sampling.orientational_prior_mode = NOPRIOR;
		adaptive_oversampling = 0;
		nr_pool = max_nr_pool = 1;
		sampling.perturbation_factor = 0.;
		sampling.random_perturbation = 0.;
		sampling.setOneOrientation(0.,0.,0.);
		directions_have_changed = true;
		if (do_realign_movies)
			nr_threads = 1; // use only one thread, as there are no particles/orientations to parallelise anyway...
		if (do_skip_align)
		{
			Matrix1D<double> offset(2);
			sampling.setOneTranslation(offset);
		}
	}

	// Resize the pdf_direction arrays to the correct size and fill with an even distribution
	if (directions_have_changed)
		mymodel.initialisePdfDirection(sampling.NrDirections(0, true));

	// Initialise the wsum_model according to the mymodel
	wsum_model.initialise(mymodel, sampling.symmetryGroup());

	// Check that number of pooled particles is not larger than 1 for local angular searches
	// Because for local searches, each particle has a different set of nonzeroprior orientations, and thus a differently sized Mweight
	// If larger than 1, just reset to 1
	if (mymodel.orientational_prior_mode != NOPRIOR && max_nr_pool > 1)
	{
		if (verb > 0)
			std::cout << " Performing local angular searches! Lowering max_nr_pool from "<<max_nr_pool<<" to 1!" << std::endl;
		max_nr_pool = 1;
	}

	// Initialise sums of hidden variable changes
	// In later iterations, this will be done in updateOverallChangesInHiddenVariables
	sum_changes_optimal_orientations = 0.;
	sum_changes_optimal_offsets = 0.;
	sum_changes_optimal_classes = 0.;
	sum_changes_count = 0.;

	// Skip scale correction if there are nor groups
	if (mymodel.nr_groups == 1)
		do_scale_correction = false;

	// Check for rlnReconstructImageName in the data.star file. If it is present, set do_use_reconstruct_images to true
	do_use_reconstruct_images = mydata.MDimg.containsLabel(EMDL_IMAGE_RECONSTRUCT_NAME);
	if (do_use_reconstruct_images && verb > 0)
		std::cout <<" Using rlnReconstructImageName from the input data.star file!" << std::endl;

#ifdef DEBUG
	std::cerr << "Leaving initialiseGeneral" << std::endl;
#endif

}

void MlOptimiser::initialiseCL(int rank) {
    char name[128];
    char extensionString[1024];
    int err;
    size_t imageWidth, imageHeight;
    cl_int CL_PreferredWGSMultiple;
    cl_ulong CL_localSize;
    
    err = clGetDeviceInfo(CL_device, CL_DEVICE_NAME, 128, name, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    err = clGetDeviceInfo(CL_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(CL_global_memsize), &CL_global_memsize, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    
    //For multi two cards sharing the same, so divide by 2...
//    CL_global_memsize /= 2;
    
    err = clGetDeviceInfo(CL_device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(CL_wgs), &CL_wgs, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    err = clGetDeviceInfo(CL_device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(CL_maxComputeUnits), &CL_maxComputeUnits, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    err = clGetDeviceInfo(CL_device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(CL_maxMemAlloc), &CL_maxMemAlloc, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    //This was for two cards sharing the same device...
//    CL_maxMemAlloc = (cl_ulong)((float)CL_maxMemAlloc * 0.5);
    
    err = clGetDeviceInfo(CL_device, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(imageWidth), &imageWidth, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    err = clGetDeviceInfo(CL_device, CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(imageHeight), &imageHeight, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    
    cl_device_local_mem_type memtype;
    err = clGetDeviceInfo(CL_device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(memtype), &memtype, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    
    err = clGetDeviceInfo(CL_device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(CL_localSize), &CL_localSize, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    
    err = clGetDeviceInfo(CL_device, CL_DEVICE_EXTENSIONS, sizeof(extensionString), extensionString, NULL);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error with clGetDeviceInfo" << std::endl;
        return;
    }
    //        std::cerr << "Extensions: " << extensionString << std::endl;
    //        if (strstr(extensionString, "cl_khr_int64_base_atomics") != NULL) {
    //            CL_atomicSupport = true;
    //        }
    if (strstr(extensionString, "cl_khr_fp64") == NULL) {
        std::cerr << "Slave " << rank << " not using CL device " << name << " - no double precision support" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_context = clCreateContext(0, 1, &CL_device, NULL, NULL, &err);
    if (!CL_context)
    {
        std::cerr << "Error: Failed to create a compute context" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_ComputeQueue = clCreateCommandQueue(CL_context, CL_device, CL_QUEUE_PROFILING_ENABLE, &err);
    if (!CL_ComputeQueue)
    {
        std::cerr << "Error: Failed to create compute command queue" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_CopyFromDeviceQueue = clCreateCommandQueue(CL_context, CL_device, CL_QUEUE_PROFILING_ENABLE, &err);
    if (!CL_CopyFromDeviceQueue)
    {
        std::cerr << "Error: Failed to create command queue to copy from device" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_CopyToDeviceQueue = clCreateCommandQueue(CL_context, CL_device, CL_QUEUE_PROFILING_ENABLE, &err);
    if (!CL_CopyToDeviceQueue)
    {
        std::cerr << "Error: Failed to create command queue to copy to device" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    //Detection of 64 bit atomics seems a bit flakey, so just try to build with atomics, if that fails switch to no atomics
    CL_atomicSupport = true;
    char *source;
    char clCode[] = {"#pragma\tOPENCL\tEXTENSION\tcl_khr_fp64\t:\tenable\n"
        "#pragma\tOPENCL\tEXTENSION\tcl_khr_int64_base_atomics\t:\tenable\n"
        "\n"
        "//__constant\tsampler_t\tsampler\t=\tCLK_FILTER_NEAREST;\n"
        "\n"
        "#define\tArow\t3\n"
        "#define\tASize\t9\n"
        "#define\tAELEM(m,gi,i,j)\t(m[gi*ASize+i*Arow+j])\n"
        "#define\tAE(m,i,j)\t(m[i*Arow+j])\n"
        "\n"
        "#define\tMODEL_ELEM(m,k,i,j)\t(m[k*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
        "\n"
        "#define\tMODEL_ELEM_2D(m,i,j)\t(m[i*modelDim.x+j])\n"
        "\n"
        "#define\tMODEL_ELEM_DOUBLE(m,k,i,j,c)\t(m[2*(k*modelDim.x*modelDim.y+i*modelDim.x+j)+c])\n"
        "\n"
        "#define\tMODEL_ELEM_2D_DOUBLE(m,i,j,c)\t(m[2*(i*modelDim.x+j)+c])\n"
        "\n"
        "#define\tOUTPUT_MODEL_ELEM(m,n,k,i,j)\t(m[n*modelDim.x*modelDim.y*modelDim.z+k*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
        "\n"
        "#define\tOUTPUT_MODEL_ELEM_2D(m,n,i,j)\t(m[n*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
        "\n"
        "#define\tOUTPUT_ELEM(m,gi,i,j)\t(m[gi*projDim.x*projDim.y+i*projDim.x+j])\n"
        "\n"
        "#define\tOUTPUT_ELEM_TR(m,gi,i,j)\t(m[gi+numOrients*(i*projDim.x+j)])\n"
        "\n"
        "#define\tOUTPUT_ELEM_WGS(m,num,i,j)\t(m[num*wgsPerFref*lgs+i*projDim.x+j])\n"
        "\n"
        "#define\tOUTPUT_ELEM_LINEAR(m,gi,i)\t(m[gi*projDim.x*projDim.y+i])\n"
        "\n"
        "#define\tINPUT_ELEM_LINEAR(m,ipart,num_trans,itrans,i)\t(m[ipart*num_trans*projDim.x*projDim.y+itrans*projDim.x*projDim.y+i])\n"
        "\n"
        "#define\tFLOOR(x)\t(((x)\t==\t(int)(x))\t?\t(int)(x):(((x)\t>\t0)\t?\t(int)(x)\t:\t\\\n"
        "(int)((x)\t-\t1)))\n"
        "\n"
        "#define\tCOMPLEX_MULT(a,\tb)\t(double2)(a.x\t*\tb.x\t-\ta.y\t*\tb.y,\ta.x\t*\tb.y\t+\ta.y\t*\tb.x)\n"
        "#define\tCOMPLEX_NORM(a)\t(dot(a,a))\n"
        "#define\tLIN_INTERP(a,\tl,\th)\t((l)\t+\t((h)\t-\t(l))\t*\t(a))\n"
        "\n"
        "#define\tDEG2RAD(a)\t(a\t*\t2\t*\tM_PI\t/\t360.0)\n"
        "\n"
        "#define\tWGS\t256\n"
        "//#define\tWGS_MULTIPLE\t4\n"
        "//#define\tGPU\n"
        "\n"
        "\n"
        "double\taccumulate(double\tinput,\tlocal\tdouble*\tbuffer)\t{\n"
        "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
        "//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
        "\t\t\t\tsize_t\tls\t=\tWGS;\n"
        "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\tbuffer[li]\t=\tinput;\n"
        "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\n"
        "\t\t\t\t#pragma\tunroll\t1\n"
        "\t\t\t\tfor(unsigned\tint\ts=ls/2;\ts>8;\ts>>=1)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\ts];\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
        "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
        "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tsum;\n"
        "\t\t\t\tif\t(li\t==\t0)\n"
        "\t\t\t\t\t\t\t\tsum\t=\tbuffer[li]\t+\tbuffer[li\t+\t1]\t+\tbuffer[li\t+\t2]\t+\tbuffer[li\t+\t3]\t+\tbuffer[li\t+\t4]\t+\tbuffer[li\t+\t5]\t+\tbuffer[li\t+\t6]\t+\tbuffer[li\t+\t7]\t+\tbuffer[li\t+\t8]\t+\tbuffer[li\t+\t9]\t+\tbuffer[li\t+\t10]\t+\tbuffer[li\t+\t11]\t+\tbuffer[li\t+\t12]\t+\tbuffer[li\t+\t13]\t+\tbuffer[li\t+\t14]\t+\tbuffer[li\t+\t15];\n"
        "\t\t\t\treturn\tsum;\n"
        "}\n"
        "\n"
        "double2\taccumulate2(double\tinput1,\tlocal\tdouble*\tbuffer1,\tdouble\tinput2,\tlocal\tdouble*\tbuffer2)\t{\n"
        "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
        "\t\t\t\t//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
        "\t\t\t\tsize_t\tls\t=\tWGS;\n"
        "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\tbuffer1[li]\t=\tinput1;\n"
        "\t\t\t\tbuffer2[li]\t=\tinput2;\n"
        "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t\n"
        "#pragma\tunroll\t1\n"
        "\t\t\t\tfor(unsigned\tint\ts=ls/2;\ts>8;\ts>>=1)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tbuffer1[li]\t+=\tbuffer1[li\t+\ts];\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tbuffer2[li]\t+=\tbuffer2[li\t+\ts];\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tsum;\n"
        "\t\t\t\tif\t(li\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\tsum.x\t=\tbuffer1[li]\t+\tbuffer1[li\t+\t1]\t+\tbuffer1[li\t+\t2]\t+\tbuffer1[li\t+\t3]\t+\tbuffer1[li\t+\t4]\t+\tbuffer1[li\t+\t5]\t+\tbuffer1[li\t+\t6]\t+\tbuffer1[li\t+\t7]\t+\tbuffer1[li\t+\t8]\t+\tbuffer1[li\t+\t9]\t+\tbuffer1[li\t+\t10]\t+\tbuffer1[li\t+\t11]\t+\tbuffer1[li\t+\t12]\t+\tbuffer1[li\t+\t13]\t+\tbuffer1[li\t+\t14]\t+\tbuffer1[li\t+\t15];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tsum.y\t=\tbuffer2[li]\t+\tbuffer2[li\t+\t1]\t+\tbuffer2[li\t+\t2]\t+\tbuffer2[li\t+\t3]\t+\tbuffer2[li\t+\t4]\t+\tbuffer2[li\t+\t5]\t+\tbuffer2[li\t+\t6]\t+\tbuffer2[li\t+\t7]\t+\tbuffer2[li\t+\t8]\t+\tbuffer2[li\t+\t9]\t+\tbuffer2[li\t+\t10]\t+\tbuffer2[li\t+\t11]\t+\tbuffer2[li\t+\t12]\t+\tbuffer2[li\t+\t13]\t+\tbuffer2[li\t+\t14]\t+\tbuffer2[li\t+\t15];\n"
        "\t\t\t\t}\n"
        "\t\t\t\treturn\tsum;\n"
        "}\n"
        "\n"
        "/*double\taccumulate(double\tinput,\tlocal\tdouble*\tbuffer)\t{\n"
        "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
        "//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
        "\t\t\t\t\t\t\t\tsize_t\tls\t=\tWGS;\n"
        "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\tbuffer[li]\t=\tinput;\n"
        "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(ls\t>=\t512)\t{\tif\t(li\t<\t256)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t256];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
        "\t\t\t\tif\t(ls\t>=\t256)\t{\tif\t(li\t<\t128)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t128];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
        "\t\t\t\tif\t(ls\t>=\t128)\t{\tif\t(li\t<\t\t64)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t64];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(li\t<\t32)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t64)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t32];\t}\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t32)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t16];\t}\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t16)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t8];\t}\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t8)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t4];\t}\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t4)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t2];\t}\n"
        "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t2)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t1];\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\treturn\tbuffer[li];\n"
        "}*/\n"
        "\n"
        "\n"
        "kernel\tvoid\tblank(global\tdouble\t*Fref)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tFref[gi]\t=\t0.0;\n"
        "}\n"
        "\n"
        "kernel\tvoid\tinitWithConstantDouble(global\tdouble\t*Fref,\tconst\tdouble\tvalue,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tFref[gi]\t=\tvalue;\n"
        "}\n"
        "\n"
        "kernel\tvoid\tinitWithConstant(global\tint\t*Fref,\tconst\tint\tvalue,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tFref[gi]\t=\tvalue;\n"
        "}\n"
        "\n"
        "\n"
        "kernel\tvoid\tgenerateOrientationMatrix(global\tdouble*\teulerOrientations,\tglobal\tdouble*\texp_R_mic,\tconst\tint\tIS_INV,\tconst\tdouble\tpaddingFactor,\tglobal\tdouble*\tmatrix,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\to\t=\tgi\t*\t3;\n"
        "\t\t\t\tsize_t\tm\t=\tgi\t*\t9;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tA[9],\tAI[9];\n"
        "\t\t\t\tdouble\tca,\tsa,\tcb,\tsb,\tcg,\tsg;\n"
        "\t\t\t\tdouble\tcc,\tcs,\tsc,\tss;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Convert\teuler\tangles\tto\tmatrix\n"
        "\t\t\t\tdouble\talpha,\tbeta,\tgamma;\n"
        "\t\t\t\t\n"
        "\t\t\t\talpha\t=\teulerOrientations[o];\n"
        "\t\t\t\tbeta\t=\teulerOrientations[o+1];\n"
        "\t\t\t\tgamma\t=\teulerOrientations[o+2];\n"
        "\t\t\t\t\n"
        "\t\t\t\talpha\t=\tDEG2RAD(alpha);\n"
        "\t\t\t\tbeta\t\t=\tDEG2RAD(beta);\n"
        "\t\t\t\tgamma\t=\tDEG2RAD(gamma);\n"
        "\t\t\t\t\n"
        "\t\t\t\tca\t=\tcos(alpha);\n"
        "\t\t\t\tcb\t=\tcos(beta);\n"
        "\t\t\t\tcg\t=\tcos(gamma);\n"
        "\t\t\t\tsa\t=\tsin(alpha);\n"
        "\t\t\t\tsb\t=\tsin(beta);\n"
        "\t\t\t\tsg\t=\tsin(gamma);\n"
        "\t\t\t\tcc\t=\tcb\t*\tca;\n"
        "\t\t\t\tcs\t=\tcb\t*\tsa;\n"
        "\t\t\t\tsc\t=\tsb\t*\tca;\n"
        "\t\t\t\tss\t=\tsb\t*\tsa;\n"
        "\t\t\t\t\n"
        "\t\t\t\tA[0]\t=\t\tcg\t*\tcc\t-\tsg\t*\tsa;\n"
        "\t\t\t\tA[1]\t=\t\tcg\t*\tcs\t+\tsg\t*\tca;\n"
        "\t\t\t\tA[2]\t=\t-cg\t*\tsb;\n"
        "\t\t\t\tA[3]\t=\t-sg\t*\tcc\t-\tcg\t*\tsa;\n"
        "\t\t\t\tA[4]\t=\t-sg\t*\tcs\t+\tcg\t*\tca;\n"
        "\t\t\t\tA[5]\t=\tsg\t*\tsb;\n"
        "\t\t\t\tA[6]\t=\t\tsc;\n"
        "\t\t\t\tA[7]\t=\t\tss;\n"
        "\t\t\t\tA[8]\t=\tcb;\n"
        "\t\t\t\t\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tASize;\ti++)\t{\n"
        "\t\t\t\t\t\t\t\tAI[i]\t=\t0;\n"
        "\t\t\t\t}\n"
        "\n"
        "\t\t\t\t//Mult\twith\texp_R_mic\n"
        "\t\t\t\t#pragma\tunroll\t1\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tArow;\ti++)\n"
        "\t\t\t\t\t\t\t\t#pragma\tunroll\t1\n"
        "\t\t\t\t\t\t\t\tfor\t(int\tj\t=\t0;\tj\t<\tArow;\tj++)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t#pragma\tunroll\t1\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tfor\t(int\tk\t=\t0;\tk\t<\tArow;\tk++)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAE(AI,i,j)\t+=\tAE(exp_R_mic,i,k)\t*\tAE(A,k,j);\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Invert\n"
        "\t\t\t\tA[0]\t=\t\t\tAI[8]*AI[4]-AI[7]*AI[5];\n"
        "\t\t\t\tA[1]\t=\t-(AI[8]*AI[1]-AI[7]*AI[2]);\n"
        "\t\t\t\tA[2]\t=\t\t\tAI[5]*AI[1]-AI[4]*AI[2];\n"
        "\t\t\t\tA[3]\t=\t-(AI[8]*AI[3]-AI[6]*AI[5]);\n"
        "\t\t\t\tA[4]\t=\t\t\tAI[8]*AI[0]-AI[6]*AI[2];\n"
        "\t\t\t\tA[5]\t=\t-(AI[5]*AI[0]-AI[3]*AI[2]);\n"
        "\t\t\t\tA[6]\t=\t\t\tAI[7]*AI[3]-AI[6]*AI[4];\n"
        "\t\t\t\tA[7]\t=\t-(AI[7]*AI[0]-AI[6]*AI[1]);\n"
        "\t\t\t\tA[8]\t=\t\t\tAI[4]*AI[0]-AI[3]*AI[1];\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\ttmp;\n"
        "\t\t\t\ttmp\t=\tAI[0]\t*\tA[0]\t+\tAI[3]\t*\tA[1]\t+\tAI[6]\t*\tA[2];\n"
        "\t\t\t\t\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\t9;\ti++)\t{\n"
        "\t\t\t\t\t\t\t\tA[i]\t/=\ttmp;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Now\ttranspose\tif\tnot\tA_INV\n"
        "\t\t\t\tAI[0]\t=\tIS_INV\t?\tA[0]\t:\tA[0];\n"
        "\t\t\t\tAI[1]\t=\tIS_INV\t?\tA[1]\t:\tA[3];\n"
        "\t\t\t\tAI[2]\t=\tIS_INV\t?\tA[2]\t:\tA[6];\n"
        "\t\t\t\tAI[3]\t=\tIS_INV\t?\tA[3]\t:\tA[1];\n"
        "\t\t\t\tAI[4]\t=\tIS_INV\t?\tA[4]\t:\tA[4];\n"
        "\t\t\t\tAI[5]\t=\tIS_INV\t?\tA[5]\t:\tA[7];\n"
        "\t\t\t\tAI[6]\t=\tIS_INV\t?\tA[6]\t:\tA[2];\n"
        "\t\t\t\tAI[7]\t=\tIS_INV\t?\tA[7]\t:\tA[5];\n"
        "\t\t\t\tAI[8]\t=\tIS_INV\t?\tA[8]\t:\tA[8];\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Multiply\tby\tpadding\tfactor\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\t9;\ti++)\t{\n"
        "\t\t\t\t\t\t\t\tmatrix[m+i]\t=\tAI[i]\t*\tpaddingFactor;\n"
        "\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tscaleMatrix(global\tdouble*\tmatrix,\tdouble\tscaleFactor,\tint\tlimit)\t{\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\tmatrix[gi]\t*=\tscaleFactor;\n"
        "\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcModelProjection_TR(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tint\tnumOrients,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
        "\t\t\t\tint\tis_neg_x;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\tdouble2\td000,\td001,\td010,\td011,\td100,\td101,\td110,\td111;\n"
        "\t\t\t\tdouble2\tdx00,\tdx01,\tdx10,\tdx11,\tdxy0,\tdxy1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
        "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
        "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tz0\t=\tFLOOR(zp);\n"
        "\t\t\t\tfz\t=\tzp\t-\tz0;\n"
        "\t\t\t\tz0\t-=\tSTARTINGZ;\n"
        "\t\t\t\tz1\t=\tz0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
        "\t\t\t\t\t\t\t\td000\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td001\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td010\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td011\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\td100\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td101\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td110\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td111\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t1)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\td000:\t%f\t%f\td001:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\td000.x,\td000.y,\td001.x,\td001.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
        "\t\t\t\t\t\t\t\tdx00\t=\tLIN_INTERP(fx,\td000,\td001);\n"
        "\t\t\t\t\t\t\t\tdx01\t=\tLIN_INTERP(fx,\td100,\td101);\n"
        "\t\t\t\t\t\t\t\tdx10\t=\tLIN_INTERP(fx,\td010,\td011);\n"
        "\t\t\t\t\t\t\t\tdx11\t=\tLIN_INTERP(fx,\td110,\td111);\n"
        "\t\t\t\t\t\t\t\tdxy0\t=\tLIN_INTERP(fy,\tdx00,\tdx10);\n"
        "\t\t\t\t\t\t\t\tdxy1\t=\tLIN_INTERP(fy,\tdx01,\tdx11);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
        "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fz,\tdxy0,\tdxy1)\t*\t(double2)(1,-1);\n"
        "\t\t\t\t\t\t\t\tOUTPUT_ELEM_TR(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
        "//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcModelProjection(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
        "\t\t\t\tint\tis_neg_x;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\tdouble2\td000,\td001,\td010,\td011,\td100,\td101,\td110,\td111;\n"
        "\t\t\t\tdouble2\tdx00,\tdx01,\tdx10,\tdx11,\tdxy0,\tdxy1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tuint\tprojSize\t=\tprojDim.x\t*\tprojDim.y;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
        "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
        "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tz0\t=\tFLOOR(zp);\n"
        "\t\t\t\tfz\t=\tzp\t-\tz0;\n"
        "\t\t\t\tz0\t-=\tSTARTINGZ;\n"
        "\t\t\t\tz1\t=\tz0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
        "\t\t\t\t\t\t\t\td000\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td001\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td010\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td011\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\td100\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td101\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td110\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td111\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t1)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\td000:\t%f\t%f\td001:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\td000.x,\td000.y,\td001.x,\td001.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
        "\t\t\t\t\t\t\t\tdx00\t=\tLIN_INTERP(fx,\td000,\td001);\n"
        "\t\t\t\t\t\t\t\tdx01\t=\tLIN_INTERP(fx,\td100,\td101);\n"
        "\t\t\t\t\t\t\t\tdx10\t=\tLIN_INTERP(fx,\td010,\td011);\n"
        "\t\t\t\t\t\t\t\tdx11\t=\tLIN_INTERP(fx,\td110,\td111);\n"
        "\t\t\t\t\t\t\t\tdxy0\t=\tLIN_INTERP(fy,\tdx00,\tdx10);\n"
        "\t\t\t\t\t\t\t\tdxy1\t=\tLIN_INTERP(fy,\tdx01,\tdx11);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
        "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fz,\tdxy0,\tdxy1)\t*\t(double2)(1,-1);\n"
        "\t\t\t\t\t\t\t\tOUTPUT_ELEM(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcModelProjection2D_TR(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tint\tnumOrients,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
        "\t\t\t\tint\tis_neg_x;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\tdouble2\td00,\td01,\td10,\td11,\tdx0,\tdx1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
        "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
        "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
        "\t\t\t\t\t\t\t\td00\t=\tMODEL_ELEM_2D(model,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td01\t=\tMODEL_ELEM_2D(model,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td10\t=\tMODEL_ELEM_2D(model,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td11\t=\tMODEL_ELEM_2D(model,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
        "\t\t\t\t\t\t\t\tdx0\t=\tLIN_INTERP(fx,\td00,\td01);\n"
        "\t\t\t\t\t\t\t\tdx1\t=\tLIN_INTERP(fx,\td10,\td11);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
        "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fy,\tdx0,\tdx1)\t*\t(double2)(1,-1);\n"
        "\t\t\t\t\t\t\t\tOUTPUT_ELEM_TR(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fy,\tdx0,\tdx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcModelProjection2D(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
        "\t\t\t\tint\tis_neg_x;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\tdouble2\td00,\td01,\td10,\td11,\tdx0,\tdx1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
        "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
        "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
        "\t\t\t\t\t\t\t\td00\t=\tMODEL_ELEM_2D(model,\ty0,\tx0);\n"
        "\t\t\t\t\t\t\t\td01\t=\tMODEL_ELEM_2D(model,\ty0,\tx1);\n"
        "\t\t\t\t\t\t\t\td10\t=\tMODEL_ELEM_2D(model,\ty1,\tx0);\n"
        "\t\t\t\t\t\t\t\td11\t=\tMODEL_ELEM_2D(model,\ty1,\tx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
        "\t\t\t\t\t\t\t\tdx0\t=\tLIN_INTERP(fx,\td00,\td01);\n"
        "\t\t\t\t\t\t\t\tdx1\t=\tLIN_INTERP(fx,\td10,\td11);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
        "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fy,\tdx0,\tdx1)\t*\t(double2)(1,-1);\n"
        "\t\t\t\t\t\t\t\tOUTPUT_ELEM(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fy,\tdx0,\tdx1);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tctfAndScaleDataPoint(global\tdouble2*\tFref,\tconst\tint\tdo_ctf_correction,\tglobal\tdouble*\tctf,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble2*\tFrefOut,\tconst\tint\tprojSize,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\tint\tFctfPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tref;\n"
        "\t\t\t\tref\t=\tFref[gi];\n"
        "//\t\t\t\tref\t=\t(double2)(1,1);\n"
        "\t\t\t\tref\t=\tdo_ctf_correction\t?\tref\t*\tctf[FctfPos]\t:\tref;\n"
        "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
        "\t\t\t\tFrefOut[gi]\t=\tref;\n"
        "//\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
        "//\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tctfAndScaleDataPoint_TR(global\tdouble2*\tFref,\tconst\tint\tdo_ctf_correction,\tglobal\tdouble*\tctf,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble2*\tFrefOut,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFctfPos\t=\tgi\t/\tnumOrients;\n"
        "//\t\t\t\tint\tFrefNo\t=\tgi\t-\tFctfPos\t*\tnumOrients;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tref;\n"
        "\t\t\t\tref\t=\tFref[gi];\n"
        "\t\t\t\t//\t\t\t\tref\t=\t(double2)(1,1);\n"
        "\t\t\t\tref\t=\tdo_ctf_correction\t?\tref\t*\tctf[FctfPos]\t:\tref;\n"
        "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
        "\t\t\t\tFrefOut[gi]\t=\tref;\n"
        "//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%3d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
        "//\t\t\t\t\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tscaleCTF(global\tdouble*\tFref,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble*\tFrefOut,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tref;\n"
        "\t\t\t\tref\t=\tFref[gi];\n"
        "\t\t\t\t//\t\t\t\tref\t=\t(double2)(1,1);\n"
        "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
        "\t\t\t\tFrefOut[gi]\t=\tref;\n"
        "\t\t\t\t//\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
        "\t\t\t\t//\t\t\t\t}\n"
        "}\n"
        "\n"
        "//Here\tglobal\tsize\tis\tnumOrientations\t*\tnumTrans\n"
        "kernel\tvoid\tcalculateDiff2AndSuma2CC_TR(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tdiff2,\tglobal\tdouble*\tsuma2,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tglobal\tchar\t*skip,\tconst\tint\tipart,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(skip[gi])\t{\n"
        "//\t\t\t\t\t\t\t\tdiff2[gi]\t=\t0;\n"
        "//\t\t\t\t\t\t\t\tsuma2[gi]\t=\t0;\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\tnumOrients;\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tnumOrients);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFimgBase\t=\t(trans\t*\tfSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\td2,\ts2;//,\tdiff2Sum,\tsuma2Sum;\n"
        "\t\t\t\td2\t=\t0;\n"
        "\t\t\t\ts2\t=\t0;\n"
        "\n"
        "\t\t\t\tfor\t(int\tFimgPos\t=\t0;\tFimgPos\t<\tfSize\t;\tFimgPos++)\t{\n"
        "\t\t\t\t\t\t\t\tdouble2\tref,\timg;\n"
        "\t\t\t\t\t\t\t\tref\t=\tFref[ti\t+\tFimgPos\t*\tnumOrients];\n"
        "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos\t+\tFimgBase];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\td2\t+=\tdot(ref,img);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\ts2\t+=\tdot(ref,ref);\n"
        "\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\tif\t((ti\t==\t0)\t&&\t(trans\t<\t2))\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFimgpos:\t%d\timgPosBase:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\ts2:\t%f\td2:\t%f\\n\",\tgi,\ttrans,\tti,\tFimgPos,\tFimgBase,\tref.x,\tref.y,\timg.x,\timg.y,\ts2,\td2);\n"
        "//\t\t\t\t\t\t\t\t}\n"
        "\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tdiff2[gi]\t=\td2;\n"
        "\t\t\t\tsuma2[gi]\t=\ts2;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tif\t(((FrefNo\t<2\t)\t&&\t(trans\t<2\t))\t||\t((FrefNo\t==\t4607)\t&&\t(trans\t==\t19)))\t{\n"
        "\t\t\t\t//\t\t\t\tif\t(isnan(diff2Sum)\t||\tisnan(suma2Sum))\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tindex:\t%d\ttrans:\t%d\tti:\t%d\trefno:\t%d\timgPos:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\td2:\t%f\ts2:\t%f\tdiff2S:\t%f\ts2S:\t%f\\n\",\tgi,\tindex,\ttrans,\tti,\tFrefNo,\tFimgPos,\tref.x,\tref.y,\timg.x,\timg.y,\td2,\ts2\t*\t1000000,\tdiff2Sum,\tsuma2Sum);\n"
        "\t\t\t\t//\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateDiff2AndSuma2CC(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tdiff2,\tglobal\tdouble*\tsuma2,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tref,\timg;\n"
        "\t\t\t\tref\t=\tFref[ti];\n"
        "\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\n"
        "//\t\t\t\tdiff2Sum\t=\tref.x\t*\timg.x;\n"
        "//\t\t\t\tdiff2Sum\t+=\tref.y\t*\timg.y;\n"
        "\n"
        "//\t\t\t\tdiff2[gi]\t=\tdot(ref,img);\n"
        "//\t\t\t\tsuma2[gi]\t=\tdot(ref,ref);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\td2,\ts2;//,\tdiff2Sum,\tsuma2Sum;\n"
        "\t\t\t\t\n"
        "\t\t\t\td2\t=\tdot(ref,img);\n"
        "//\t\t\t\tdiff2Sum\t=\taccumulate(d2,\tbuffer);\n"
        "\n"
        "\t\t\t\ts2\t=\tdot(ref,ref);\n"
        "//\t\t\t\tsuma2Sum\t=\taccumulate(s2,\tbuffer);\n"
        "\t\t\t\t\n"
        "#ifdef\tGPU\n"
        "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\t__local\tdouble\tbuffer1[WGS];\n"
        "\t\t\t\t__local\tdouble\tbuffer2[WGS];\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tsum;\n"
        "\t\t\t\tsum\t=\taccumulate2(d2,buffer1,\ts2,\tbuffer2);\n"
        "#endif\n"
        "\n"
        "//\t\t\t\tsize_t\tindex\t=\tFrefNo\t*\twgsPerFref\t+\tFimgPos\t/\tlgs\t+\ttrans\t*\twgsPerFref\t*\torientationsPerTrans;\n"
        "\t\t\t\tsize_t\tindex\t=\tgi\t/\tlgs;\n"
        "#ifdef\tGPU\n"
        "\t\t\t\tif\t(li\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\tdiff2[index]\t=\tsum.x;\n"
        "\t\t\t\t\t\t\t\tsuma2[index]\t=\tsum.y;\n"
        "\t\t\t\t}\n"
        "#else\n"
        "\t\t\t\tdiff2[index]\t=\td2;\n"
        "\t\t\t\tsuma2[index]\t=\ts2;\n"
        "#endif\n"
        "\t\t\t\t\n"
        "//\t\t\t\tif\t(((FrefNo\t<2\t)\t&&\t(trans\t<2\t))\t||\t((FrefNo\t==\t4607)\t&&\t(trans\t==\t19)))\t{\n"
        "//\t\t\t\tif\t(isnan(diff2Sum)\t||\tisnan(suma2Sum))\t{\n"
        "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tindex:\t%d\ttrans:\t%d\tti:\t%d\trefno:\t%d\timgPos:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\td2:\t%f\ts2:\t%f\tdiff2S:\t%f\ts2S:\t%f\\n\",\tgi,\tindex,\ttrans,\tti,\tFrefNo,\tFimgPos,\tref.x,\tref.y,\timg.x,\timg.y,\td2,\ts2\t*\t1000000,\tdiff2Sum,\tsuma2Sum);\n"
        "//\t\t\t\t}\n"
        "}\n"
        "\n"
        "inline\tvoid\tAtomicAdd(volatile\t__global\tdouble\t*source,\tconst\tdouble\toperand)\t{\n"
        "\t\t\t\tunion\t{\n"
        "\t\t\t\t\t\t\t\tunsigned\tlong\tintVal;\n"
        "\t\t\t\t\t\t\t\tdouble\tfloatVal;\n"
        "\t\t\t\t}\tnewVal;\n"
        "\t\t\t\tunion\t{\n"
        "\t\t\t\t\t\t\t\tunsigned\tlong\tintVal;\n"
        "\t\t\t\t\t\t\t\tdouble\tfloatVal;\n"
        "\t\t\t\t}\tprevVal;\n"
        "\t\t\t\tdo\t{\n"
        "\t\t\t\t\t\t\t\tprevVal.floatVal\t=\t*source;\n"
        "\t\t\t\t\t\t\t\tnewVal.floatVal\t=\tprevVal.floatVal\t+\toperand;\n"
        "\t\t\t\t}\twhile\t(atomic_cmpxchg((volatile\t__global\tunsigned\tlong\t*)source,\tprevVal.intVal,\tnewVal.intVal)\t!=\tprevVal.intVal);\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateDiff2AndSuma2_TR(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tdiff2,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tglobal\tchar\t*skip,\tconst\tint\tipart,\t\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t/*size_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\t__local\tdouble\tlocalMinSigma2[WGS\t*\tWGS_MULTIPLE];\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tWGS_MULTIPLE;\ti++)\t{\n"
        "\t\t\t\t\t\t\t\tint\tindex\t=\ti\t*\tWGS\t+\tli;\n"
        "\t\t\t\t\t\t\t\tif\t(index\t<\tfSize)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tlocalMinSigma2[index]\t=\t0.5\t*\tMinvsigma2[index];\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);*/\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(skip[gi])\t{\n"
        "//\t\t\t\t\t\t\t\tdiff2[gi]\t=\t0;\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\tnumOrients;\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tnumOrients);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFimgBase\t=\t(trans\t*\tfSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\td2;//,\tdiff2Sum,\tsuma2Sum;\n"
        "\t\t\t\td2\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tfor\t(int\tFimgPos\t=\t0;\tFimgPos\t<\tfSize\t;\tFimgPos++)\t{\n"
        "\t\t\t\t\t\t\t\tdouble2\tref,\timg,\tdiff;\n"
        "\t\t\t\t\t\t\t\tref\t=\tFref[ti\t+\tFimgPos\t*\tnumOrients];\n"
        "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos\t+\tFimgBase];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\td2\t+=\tdot(diff,diff)\t*\tMinvsigma2[FimgPos];\n"
        "//\t\t\t\t\t\t\t\td2\t+=\tdot(diff,diff)\t*\tlocalMinSigma2[FimgPos];\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t((ipart\t==\t1)\t&&\t(get_global_id(0)\t==\t0))\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFimgpos:\t%d\timgPosBase:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\tdiff:\t%f\t%f\td2:\t%f\tMin2:\t%f\\n\",\tgi,\ttrans,\tti,\tFimgPos,\tFimgBase,\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\td2,\tMinvsigma2[FimgPos]);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tdiff2[gi]\t=\td2;\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateDiff2AndSuma2(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tdiff2,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
        "\t\t\t\tref\t=\tFref[ti];\n"
        "\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\td2,\tdiff2Sum;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdiff\t=\tref\t-\timg;\n"
        "\t\t\t\t\n"
        "\t\t\t\td2\t=\tdot(diff,diff)\t*\t0.5\t*\tMinvsigma2[FimgPosMin];\n"
        "//\t\t\t\tdiff2Sum\t=\taccumulate(d2,\tbuffer);\n"
        "\n"
        "#ifdef\tGPU\n"
        "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
        "\t\t\t\t__local\tdouble\tbuffer[WGS];\n"
        "\t\t\t\tbuffer[li]\t=\td2;\n"
        "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t\n"
        "#pragma\tunroll\t1\n"
        "\t\t\t\tfor(unsigned\tint\ts=WGS/2;\ts>16;\ts>>=1)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\ts];\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
        "\t\t\t\t}\n"
        "\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
        "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "//\t\t\t\tsize_t\tindex\t=\tFrefNo\t*\twgsPerFref\t+\tFimgPos\t/\tlgs;\n"
        "\t\t\t\tsize_t\tindex\t=\t(gi\t/\tlgs);\n"
        "\t\t\t\tif\t(li\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\tAtomicAdd(&diff2[index],\tsum);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdiff2[index]\t=\tbuffer[li]\t+\tbuffer[li\t+\t1]\t+\tbuffer[li\t+\t2]\t+\tbuffer[li\t+\t3]\t+\tbuffer[li\t+\t4]\t+\tbuffer[li\t+\t5]\t+\tbuffer[li\t+\t6]\t+\tbuffer[li\t+\t7]\t+\tbuffer[li\t+\t8]\t+\tbuffer[li\t+\t9]\t+\tbuffer[li\t+\t10]\t+\tbuffer[li\t+\t11]\t+\tbuffer[li\t+\t12]\t+\tbuffer[li\t+\t13]\t+\tbuffer[li\t+\t14]\t+\tbuffer[li\t+\t15]\t+\tbuffer[li\t+\t16]\t+\tbuffer[li\t+\t17]\t+\tbuffer[li\t+\t18]\t+\tbuffer[li\t+\t19]\t+\tbuffer[li\t+\t20]\t+\tbuffer[li\t+\t21]\t+\tbuffer[li\t+\t22]\t+\tbuffer[li\t+\t23]\t+\tbuffer[li\t+\t24]\t+\tbuffer[li\t+\t25]\t+\tbuffer[li\t+\t26]\t+\tbuffer[li\t+\t27]\t+\tbuffer[li\t+\t28]\t+\tbuffer[li\t+\t29]\t+\tbuffer[li\t+\t30]\t+\tbuffer[li\t+\t31];\n"
        "\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\tdiff2[index]\t=\tbuffer[li];\n"
        "\t\t\t\t}\n"
        "#else\n"
        "\t\t\t\tsize_t\tindex\t=\t(gi\t/\tlgs);\n"
        "\t\t\t\tdiff2[index]\t=\td2;\n"
        "#endif\n"
        "\t\t\t\t\n"
        "//\t\t\t\tif\t(gi\t<\twgsPerFref\t*\tlgs\t*\t2)\t{\n"
        "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFrefNo:\t%d\timgPos:\t%d\timgPosMin:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\tdiff:\t%f\t%f\td2:\t%f\tMin2:\t%f\tdiff2S:\t%f\tindex:\t%d\\n\",\tgi,\ttrans,\tti,\tFrefNo,\tFimgPos,\tFimgPosMin,\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\td2,\tMinvsigma2[FimgPosMin],\tdiff2Sum,\tindex);\n"
        "//\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateSWSNoiseEstNormCorrection(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tsigma2_noise,\tglobal\tdouble*\tdata_vs_prior,\tglobal\tdouble*\tsumXA,\tglobal\tdouble*\tsumAA,\tconst\tint\tprojSize,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(projSize);\n"
        "\t\t\t\tint\tsaveNoiseCalculations,\tsaveNormCalculations;\n"
        "\n"
        "\t\t\t\tdouble\tsum\t=\t0,\txa\t=\t0,\taa\t=\t0;\n"
        "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(projSize);\n"
        "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\tint\tMresolFine\t=\tMresol_fine[FimgPosMin];\n"
        "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tprojSize)\t+\t(trans\t*\tprojSize);\n"
        "\t\t\t\t\n"
        "//\t\t\t\tif\t(MresolFine\t>\t-1)\t{\n"
        "\t\t\t\tsaveNoiseCalculations\t=\tMresolFine\t>\t-1\t?\t1\t:\t0;\n"
        "//\t\t\t\tsaveCalculations\t=\tMresolFine\t>\t-1\t?\t1\t:\t0;\n"
        "\t\t\t\tsaveNormCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\tsaveNoiseCalculations\t:\t0;\n"
        "//\t\t\t\t\t\t\t\tsaveCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\tdouble\tcurrentWeight\t=\tweight[weightIndex];\n"
        "\t\t\t\t\t\t\t\tif\t(currentWeight\t>\t0)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(projSize);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tprojSize);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
        "\t\t\t\t\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\tdouble\tw2;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\tw2\t=\tweight[weightIndex]\t*\tdot(diff,\tdiff);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tsum\t+=\tMresolFine\t>\t-1\t?\tweight[weightIndex]\t*\tdot(diff,\tdiff)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tsum\t+=\tsaveNoiseCalculations\t?\tcurrentWeight\t*\tdot(diff,\tdiff)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\txa\t+=\tsaveNormCalculations\t?\tcurrentWeight\t*\tdot(ref,img)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\taa\t+=\tsaveNormCalculations\t?\tcurrentWeight\t*\tdot(ref,ref)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(currentWeight\t>\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(isnan(sum))\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\tsum:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tcurrentWeight,\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff),\tsum);\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
        "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
        "\t\t\t\t}\n"
        "//\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tsigma2_noise[get_global_id(0)]\t=\tsum;\n"
        "\t\t\t\tsumXA[get_global_id(0)]\t=\txa;\n"
        "\t\t\t\tsumAA[get_global_id(0)]\t=\taa;\n"
        "\n"
        "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "}\n"
        "\n"
        "\n"
        "kernel\tvoid\tcalculateSWSNoiseEstimate(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tsigma2_noise,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tsum\t=\t0;\n"
        "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(lgs\t*\twgsPerFref);\n"
        "\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tMresolFine\t=\tMresol_fine[FimgPosMin];\n"
        "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
        "\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
        "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\tdouble\tw2;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\tw2\t=\tweight[weightIndex]\t*\tdot(diff,\tdiff);\n"
        "\t\t\t\t\t\t\t\tsum\t+=\tMresolFine\t>\t-1\t?\tweight[weightIndex]\t*\tdot(diff,\tdiff)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
        "//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
        "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
        "\t\t\t\t}\n"
        "\n"
        "\t\t\t\tsigma2_noise[get_global_id(0)]\t=\tsum;\n"
        "\n"
        "//\t\t\t\tif\t(sum\t>\t0)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
        "//\t\t\t\t\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateSWSNormCorrection(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tdata_vs_prior,\tglobal\tdouble*\tsumXA,\tglobal\tdouble*\tsumAA,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tsaveCalculations;\n"
        "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(lgs\t*\twgsPerFref);\n"
        "\n"
        "\t\t\t\tdouble\txa\t=\t0,\taa\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\n"
        "\t\t\t\tsaveCalculations\t=\tMresol_fine[FimgPosMin]\t>\t-1\t?\t1\t:\t0;\n"
        "\t\t\t\tsaveCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\tsaveCalculations\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
        "//\t\t\t\t\t\t\t\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
        "\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
        "\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble2\tref,\timg;\n"
        "\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
        "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\txa\t+=\tsaveCalculations\t?\tweight[weightIndex]\t*\tdot(ref,img)\t:\t0;\n"
        "\t\t\t\t\t\t\t\taa\t+=\tsaveCalculations\t?\tweight[weightIndex]\t*\tdot(ref,ref)\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
        "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tsumXA[get_global_id(0)]\t=\txa;\n"
        "\t\t\t\tsumAA[get_global_id(0)]\t=\taa;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "}\n"
        "\n"
        "//Global\tsize\there\tis\tno\tof\trotations\t*\twgs\t*\twgsPerFref\n"
        "//Sum\tfor\teach\tF\tfor\tall\ttranslations\tover\tone\torientation\n"
        "kernel\tvoid\tcalculateSWSSum(global\tdouble*\tctf,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tweight,\tglobal\tdouble2*\tFimgOut,\tglobal\tdouble*\tFweight,\tconst\tuint\tprojSize,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit,\tconst\tuint\tgi_limit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tgi_limit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tglobalStep\t=\tgi_limit\t/\t(projSize);\n"
        "\n"
        "\t\t\t\tdouble2\tsum;//\t=\t(double2)(0,\t0);\n"
        "\t\t\t\tsum\t=\tFimgOut[get_global_id(0)];\n"
        "\t\t\t\tdouble\tweightSum;//\t=\t0;\n"
        "\t\t\t\tweightSum\t=\t\tFweight[get_global_id(0)];\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(projSize);\n"
        "\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
        "\t\t\t\tint\tFrefNo\t=\tti\t/\t(projSize);\n"
        "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tprojSize)\t+\t(trans\t*\tprojSize);\n"
        "\n"
        "\t\t\t\tdouble\tmyctf\t=\tctf[FimgPosMin];\n"
        "\t\t\t\tdouble\tMinsigv2\t=\tMinvsigma2[FimgPosMin];\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Calculate\tone\tvalue\tof\tFimg\tand\tFweight\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0.0)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(projSize);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tprojSize);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tdouble2\timg;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tdouble\tweightxinvsigma2;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tweightxinvsigma2\t=\tweight[weightIndex]\t*\tmyctf\t*\tMinsigv2;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tsum\t+=\timg\t*\tweightxinvsigma2;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tweightSum\t+=\tweightxinvsigma2\t*\tmyctf;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t421)\t{\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\timgPOosMin:\t%d,\tctf:\t%f,\tMinsigv2:\t%f,\tFimg:\t%f\t%f\tSum:\t%f\t%f\twSum:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tFimgPosMin,\tmyctf,\tMinvsigma2[FimgPosMin],\timg.x,\timg.y,\tsum.x,\tsum.y,\tweightSum);\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\tgi\t+=\torientationsPerTrans\t*\tprojSize;\n"
        "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
        "//\t\t\t\t\t\t\t\tweightIndex\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tFimgOut[get_global_id(0)]\t=\tsum;\n"
        "\t\t\t\tFweight[get_global_id(0)]\t=\tweightSum;\n"
        "//\t\t\t\tFimgOut[get_global_id(0)]\t=\t1.0;\n"
        "//\t\t\t\tFweight[get_global_id(0)]\t=\t2.0;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "}\n"
        "\n"
        "\n"
        "//Backproject\n"
        "//Two\tchoices\thow\tto\tdo\tthis\n"
        "//A)\tGlobal\tsize\tis\trelated\tto\tcompute\tunits\n"
        "//\t\t\tLoop\tover\tthe\tinput\tF's\t-->\tmaps\teach\tF\tand\tweight\tfrom\ta\trotation\tto\ta\twork-group\titem\n"
        "//\t\t\tEach\twork\tgroup\titem\tthen\tneeds\tto\twrite\t8\tweights\tand\t8\tvalues\tto\tmodel\n"
        "//\t\t\tModels\tare\tfar\ttoo\tlarge\tfor\tlocal\tmemory\n"
        "//\t\t\tTo\tsynchronise\tthe\twrites\twill\tneed\tto\thave\tseparate\tmodel\tfor\teach\twork\titem\t-\tbetter\tthan\tdouble\tatomic\toperations\n"
        "//\t\t\tTo\tbe\tpractical\twill\thave\tto\treduce\twork\tgroup\tsize\tto\t1...\n"
        "\n"
        "//B)\tGlobal\tsize\there\tis\tthe\tmodel\tsize\t-->\tmaps\teach\tweight\tand\tvalue\tof\tmodel\tto\ta\twork\titem\n"
        "//\t\t\tSimple\tway\tto\timplement\tthis\tis\tthat\teach\twork\titem\tloops\tthrough\torientations\tand\tstores\tvalues\tthat\tmatch\tthe\tcurrent\tmodel\tposition\n"
        "\n"
        "//This\tis\toption\tA\n"
        "kernel\tvoid\tcalcBackProjectionNonAtomic(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble4*\tmodelValues,\tconst\tuint\tprojSize,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
        "\t\t\t\tbool\tis_neg_x;\n"
        "\t\t\t\tdouble\tdd000,\tdd001,\tdd010,\tdd011,\tdd100,\tdd101,\tdd110,\tdd111;\n"
        "\t\t\t\tdouble2\tmy_val;\n"
        "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tint\toutputModel\t=\tget_global_id(0);\n"
        "//\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\n"
        "//\t\t\t\tif\t(gi\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
        "//\t\t\t\t}\n"
        "\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\t\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tint\ti;\n"
        "\t\t\t\t\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tint\tx;\n"
        "\t\t\t\t\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\t\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\n"
        "\t\t\t\t\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\t\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\t\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
        "\t\t\t\t\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\n"
        "\t\t\t\t\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\t\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\t\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\t\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\t\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\t\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tz0\t=\tFLOOR(zp);\n"
        "\t\t\t\t\t\t\t\tfz\t=\tzp\t-\tz0;\n"
        "\t\t\t\t\t\t\t\tz0\t-=\tSTARTINGZ;\n"
        "\t\t\t\t\t\t\t\tz1\t=\tz0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
        "\t\t\t\t\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
        "\t\t\t\t\t\t\t\tmfz\t=\t1.\t-\tfz;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdd000\t=\tmfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\tdd001\t=\tmfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\tdd010\t=\tmfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\tdd011\t=\tmfz\t*\t\tfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\tdd100\t=\t\tfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\tdd101\t=\t\tfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\tdd110\t=\t\tfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\tdd111\t=\t\tfz\t*\t\tfy\t*\t\tfx;\n"
        "\n"
        "\t\t\t\t\t\t\t\tmy_weight\t=\tOUTPUT_ELEM(fWeight,FrefNo,i,x);\n"
        "\n"
        "\t\t\t\t\t\t\t\tvalid\t=\tmy_weight\t>\t0.0\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM(f2d,FrefNo,i,x);\n"
        "\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tom:\t%d\tXP:\t%.2f\tX0:\t%d\tX1:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\ti:\t%d\tx:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\tgi,\toutputModel,\txp,\tx0,\tx1,\typ,\ty0,\ty1,\tzp,\tz0,\tz1,\ti,\tx,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty0,\tx0)\t+=\tdd000\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty0,\tx1)\t+=\tdd001\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty1,\tx0)\t+=\tdd010\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty1,\tx1)\t+=\tdd011\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty0,\tx0)\t+=\tdd100\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty0,\tx1)\t+=\tdd101\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty1,\tx0)\t+=\tdd110\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty1,\tx1)\t+=\tdd111\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_weight;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tactual:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tstored:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\\n\",\tgi,\tdd000\t*\tmy_weight,\tdd001\t*\tmy_weight,\tdd010\t*\tmy_weight,\tdd011\t*\tmy_weight,\tdd100\t*\tmy_weight,\tdd101\t*\tmy_weight,\tdd110\t*\tmy_weight,\tdd111\t*\tmy_weight,\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0),\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1));\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
        "\t\t\t\t}\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcBackProjection(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble\t*modelValues,\tglobal\tdouble\t*modelWeights,\tconst\tuint\tprojSize,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
        "\t\t\t\tbool\tis_neg_x;\n"
        "\t\t\t\tdouble\tdd000,\tdd001,\tdd010,\tdd011,\tdd100,\tdd101,\tdd110,\tdd111;\n"
        "\t\t\t\tdouble2\tmy_val;\n"
        "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\toutputModel\t=\t0;\n"
        "\t\t\t\t//\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tif\t(gi\t==\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
        "\t\t\t\t//\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tz0\t=\tFLOOR(zp);\n"
        "\t\t\t\tfz\t=\tzp\t-\tz0;\n"
        "\t\t\t\tz0\t-=\tSTARTINGZ;\n"
        "\t\t\t\tz1\t=\tz0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
        "\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
        "\t\t\t\tmfz\t=\t1.\t-\tfz;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdd000\t=\tmfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\tdd001\t=\tmfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\tdd010\t=\tmfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\tdd011\t=\tmfz\t*\t\tfy\t*\t\tfx;\n"
        "\t\t\t\tdd100\t=\t\tfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\tdd101\t=\t\tfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\tdd110\t=\t\tfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\tdd111\t=\t\tfz\t*\t\tfy\t*\t\tfx;\n"
        "\t\t\t\t\n"
        "\t\t\t\tmy_weight\t=\tOUTPUT_ELEM(fWeight,FrefNo,i,x);\n"
        "\t\t\t\t\n"
        "\t\t\t\tvalid\t=\tmy_weight\t>\t0.0\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM(f2d,FrefNo,i,x);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tom:\t%d\tXP:\t%.2f\tX0:\t%d\tX1:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\ti:\t%d\tx:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\tgi,\toutputModel,\txp,\tx0,\tx1,\typ,\ty0,\ty1,\tzp,\tz0,\tz1,\ti,\tx,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty0,\tx0,\t0),\tdd000\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty0,\tx0,\t1),\tdd000\t*\tmy_val.y);\n"
        "\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty0,\tx1,\t0),\tdd001\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty0,\tx1,\t1),\tdd001\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty1,\tx0,\t0),\tdd010\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty1,\tx0,\t1),\tdd010\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty1,\tx1,\t0),\tdd011\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz0,\ty1,\tx1,\t1),\tdd011\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty0,\tx0,\t0),\tdd100\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty0,\tx0,\t1),\tdd100\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty0,\tx1,\t0),\tdd101\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty0,\tx1,\t1),\tdd101\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty1,\tx0,\t0),\tdd110\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty1,\tx0,\t1),\tdd110\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty1,\tx1,\t0),\tdd111\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_DOUBLE(modelValues,\tz1,\ty1,\tx1,\t1),\tdd111\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz0,\ty0,\tx0),\tdd000\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz0,\ty0,\tx1),\tdd001\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz0,\ty1,\tx0),\tdd010\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz0,\ty1,\tx1),\tdd011\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz1,\ty0,\tx0),\tdd100\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz1,\ty0,\tx1),\tdd101\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz1,\ty1,\tx0),\tdd110\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM(modelWeights,\tz1,\ty1,\tx1),\tdd111\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tactual:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tstored:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\\n\",\tgi,\tdd000\t*\tmy_weight,\tdd001\t*\tmy_weight,\tdd010\t*\tmy_weight,\tdd011\t*\tmy_weight,\tdd100\t*\tmy_weight,\tdd101\t*\tmy_weight,\tdd110\t*\tmy_weight,\tdd111\t*\tmy_weight,\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1));\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalcBackProjection2D(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble\t*modelValues,\tglobal\tdouble\t*modelWeights,\tconst\tuint\tprojSize,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
        "\t\t\t\tbool\tis_neg_x;\n"
        "\t\t\t\tdouble\tdd00,\tdd01,\tdd10,\tdd11;\n"
        "\t\t\t\tdouble2\tmy_val;\n"
        "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\toutputModel\t=\t0;\n"
        "\t\t\t\t//\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\t\t\t\tif\t(gi\t==\t0)\t{\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
        "\t\t\t\t//\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
        "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti;\n"
        "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tx;\n"
        "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
        "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
        "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\n"
        "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\t\n"
        "\t\t\t\t\n"
        "\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
        "\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
        "\t\t\t\t}\n"
        "\t\t\t\telse\n"
        "\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\tfx\t=\txp\t-\tx0;\n"
        "\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\tfy\t=\typ\t-\ty0;\n"
        "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\n"
        "\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
        "\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdd00\t=\tmfy\t*\tmfx;\n"
        "\t\t\t\tdd01\t=\tmfy\t*\t\tfx;\n"
        "\t\t\t\tdd10\t=\t\tfy\t*\tmfx;\n"
        "\t\t\t\tdd11\t=\t\tfy\t*\t\tfx;\n"
        "\t\t\t\t\n"
        "\t\t\t\tmy_weight\t=\tOUTPUT_ELEM(fWeight,FrefNo,i,x);\n"
        "\t\t\t\t\n"
        "\t\t\t\tvalid\t=\tmy_weight\t>\t0.0\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(valid)\t{\n"
        "\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM(f2d,FrefNo,i,x);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tom:\t%d\tXP:\t%.2f\tX0:\t%d\tX1:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\ti:\t%d\tx:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\tgi,\toutputModel,\txp,\tx0,\tx1,\typ,\ty0,\ty1,\tzp,\tz0,\tz1,\ti,\tx,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty0,\tx0,\t0),\tdd00\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty0,\tx0,\t1),\tdd00\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty0,\tx1,\t0),\tdd01\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty0,\tx1,\t1),\tdd01\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty1,\tx0,\t0),\tdd10\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty1,\tx0,\t1),\tdd10\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty1,\tx1,\t0),\tdd11\t*\tmy_val.x);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D_DOUBLE(modelValues,\ty1,\tx1,\t1),\tdd11\t*\tmy_val.y);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D(modelWeights,\ty0,\tx0),\tdd00\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D(modelWeights,\ty0,\tx1),\tdd01\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D(modelWeights,\ty1,\tx0),\tdd10\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\tAtomicAdd(&MODEL_ELEM_2D(modelWeights,\ty1,\tx1),\tdd11\t*\tmy_weight);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tactual:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tstored:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\\n\",\tgi,\tdd000\t*\tmy_weight,\tdd001\t*\tmy_weight,\tdd010\t*\tmy_weight,\tdd011\t*\tmy_weight,\tdd100\t*\tmy_weight,\tdd101\t*\tmy_weight,\tdd110\t*\tmy_weight,\tdd111\t*\tmy_weight,\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0),\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1));\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "\n"
        "//Global\there\tis\tthe\tmodel\tsize\n"
        "kernel\tvoid\tsumModel(global\tdouble4\t*modelValuesIn,\tglobal\tdouble2\t*modelValuesOut,\tglobal\tdouble\t*modelWeightsOut,\tconst\tint\ttotalModels,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tvalueSum\t=\tmodelValuesOut[gi];\n"
        "\t\t\t\tdouble\tweightSum\t=\tmodelWeightsOut[gi];\n"
        "\t\t\t\t\n"
        "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\ttotalModels;\ti++)\t{\n"
        "\t\t\t\t\t\t\t\tdouble4\tmodelValueWeight\t=\tmodelValuesIn[gi];\n"
        "\t\t\t\t\t\t\t\tvalueSum\t+=\tmodelValueWeight.xy;\n"
        "\t\t\t\t\t\t\t\tweightSum\t+=\tmodelValueWeight.z;\n"
        "\t\t\t\t\t\t\t\tgi\t+=\tlimit;\n"
        "\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tmodelValuesOut[get_global_id(0)]\t=\tvalueSum;\n"
        "\t\t\t\tmodelWeightsOut[get_global_id(0)]\t=\tweightSum;\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "//Option\tB\tseems\tbetter,\tif\tdoing\texcess\tcalculations...\n"
        "//Tried\tB\tfirst\t-\ttoo\tmany\tcollisions\t(or\ttoo\tmany\tloops)\tand\tcrashes\ton\tthe\tGPU\n"
        "kernel\tvoid\tcalcBackProjectionB(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble2*\tmodelValues,\tglobal\tdouble*\tmodelWeights,\tconst\tint\twgsPerFref,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
        "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
        "\t\t\t\tbool\tis_neg_x;\n"
        "\t\t\t\tdouble\tdd000,\tdd001,\tdd010,\tdd011,\tdd100,\tdd101,\tdd110,\tdd111;\n"
        "\t\t\t\tdouble2\tmy_val;\n"
        "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
        "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
        "\n"
        "\t\t\t\tdouble2\tmodelValue\t=\t(0,0);\n"
        "\t\t\t\tdouble\tmodelWeight\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\t//Calculate\tmodel\tcoordinates\n"
        "\t\t\t\tint\txysize\t=\tmodelDim.x\t*\tmodelDim.y;\n"
        "\t\t\t\tint\tzi\t=\tgi\t/\txysize;\n"
        "\t\t\t\tint\tremainder\t=\tgi\t-\tzi\t*\txysize;\n"
        "\t\t\t\tint\tyi\t=\tremainder\t/\tmodelDim.x;\n"
        "\t\t\t\tint\txi\t=\tremainder\t-\tyi\t*\tmodelDim.x;\n"
        "\t\t\t\t\n"
        "//\t\t\t\tif\t(gi\t==\t53100\t+\t840\t+\t1)\n"
        "//\t\t\t\tprintf(\"gi:\t%d\tMDim:\t%d\t%d\t%d\txi:\t%d\tyi:\t%d\tzi:\t%d\txysize:\t%d\\n\",\tgi,\tmodelDim.x,\tmodelDim.y,\tmodelDim.z,\txi,\tyi,\tzi,\txysize);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
        "\t\t\t\t\n"
        "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
        "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
        "\t\t\t\t\n"
        "//\t\t\t\tif\t(gi\t==\t0)\t{\n"
        "//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
        "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
        "//\t\t\t\t}\n"
        "\t\t\t\t\n"
        "\t\t\t\tfor\t(int\tFrefNo\t=\t0;\tFrefNo\t<\torientations;\tFrefNo++)\t{\n"
        "//\t\t\t\tint\tFrefNo\t=\t0;\n"
        "\t\t\t\t\t\t\t\tfor\t(int\ti=0;\ti\t<\tprojDim.y;\ti++)\n"
        "\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tint\tyValid\t=\t1;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t//\tDont\tsearch\tbeyond\tsquare\twith\tside\tmax_r\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfirst_x\t=\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tx==0\tplane\tis\tstored\ttwice\tin\tthe\tFFTW\tformat.\tDont\tset\tit\ttwice\tin\tBACKPROJECTION!\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfirst_x\t=\t1;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyValid\t=\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\ty2\t=\ty\t*\ty;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tfor\t(int\tx=first_x;\tx\t<=\tmy_r_max;\tx++)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tint\tvalid\t=\tyValid;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tOnly\tinclude\tpoints\twith\tradius\t<\tmax_r\t(exclude\tpoints\toutside\tcircle\tin\tsquare)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tr2\t=\tx\t*\tx\t+\ty2;\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\tr2\t>\tmax_r2\t?\t0\t:\tvalid;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(r2\t>\tmax_r2)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcontinue;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tthe\trelevant\tvalue\tin\tthe\tinput\timage\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM_WGS(f2d,FrefNo,i,x);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tDIRECT_A2D_ELEM(f2d,\ti,\tx);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tthe\tweight\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(Mweight\t!=\tNULL)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_weight\t=\tDIRECT_A2D_ELEM(*Mweight,\ti,\tx);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_weight\t=\tOUTPUT_ELEM_WGS(fWeight,FrefNo,i,x);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\telse:\tmy_weight\twas\talready\tinitialised\tto\t1.\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\tmy_weight\t>\t0.\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(my_weight\t>\t0.)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"Got\ta\tpositive\tweight\");\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
        "\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(xp\t<\t0)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\telse\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tx0\t=\tFLOOR(xp);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfx\t=\t1.0\t-\tfabs(xp\t-\t(double)xi);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tx1\t=\tx0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty0\t=\tFLOOR(yp);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t-=\tSTARTINGY;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfy\t=\t1.0\t-\tfabs(yp\t-\t(double)yi);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty1\t=\ty0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz0\t=\tFLOOR(zp);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t-=\tSTARTINGZ;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfz\t=\t1.0\t-\tfabs(zp\t-\t(double)zi);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz0\t-=\tSTARTINGZ;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz1\t=\tz0\t+\t1;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((xi\t==\tx0)\t||\t(xi\t==\tx1))\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((yi\t==\ty0)\t||\t(yi\t==\ty1))\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((zi\t==\tz0)\t||\t(zi\t==\tz1))\t?\tvalid\t:\t0;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(valid)\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"XP:\t%.2f\tX0:\t%d\tX1:\t%d\tXi:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tYi:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\tZi:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\txp,\tx0,\tx1,\txi,\typ,\ty0,\ty1,\tyi,\tzp,\tz0,\tz1,\tzi,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(valid)\tprintf(\"Valid\");\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmodelValue\t=\tvalid\t?\tmodelValue\t+\tfx\t*\tfy\t*\tfz\t*\tmy_val\t:\tmodelValue;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmodelWeight\t=\tvalid\t?\tmodelWeight\t+\tfx\t*\tfy\t*\tfz\t*\tmy_weight\t:\tmodelWeight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "/*\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfz\t=\t1.\t-\tfz;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd000\t=\tmfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd001\t=\tmfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd010\t=\tmfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd011\t=\tmfz\t*\t\tfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd100\t=\t\tfz\t*\tmfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd101\t=\t\tfz\t*\tmfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd110\t=\t\tfz\t*\t\tfy\t*\tmfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd111\t=\t\tfz\t*\t\tfy\t*\t\tfx;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(is_neg_x)\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tconj(my_val);\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_val;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_weight;\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_weight;\n"
        "*/\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t//\tendif\tweight>0.\n"
        "\t\t\t\t\t\t\t\t\t\t\t\t}\t//\tendif\tx-loop\n"
        "\t\t\t\t\t\t\t\t}\t//\tendif\ty-loop\n"
        "\t\t\t\t}\t//\tendif\tFrefNo\n"
        "\t\t\t\t\n"
        "\t\t\t\tmodelValues[gi]\t=\tmodelValue;\n"
        "\t\t\t\tmodelWeights[gi]\t=\tmodelWeight;\n"
        "\t\t\t\t\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculateSigma2Offset(global\tdouble\t*weights,\tconst\tdouble2\tpriorOffset,\tglobal\tdouble2\t*translations,\tglobal\tdouble\t*sigma2OffsetReturn,\tconst\tint\ttotalOrients,\tconst\tint\tgi_limit,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tgi_limit)\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tsigma2Offset\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\t0;\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\tdouble2\tmyTranslation\t=\ttranslations[trans];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble2\tdiff;\n"
        "\t\t\t\t\t\t\t\tdiff\t=\tpriorOffset\t-\tmyTranslation;\n"
        "\t\t\t\t\t\t\t\tdouble\tsum;\n"
        "\t\t\t\t\t\t\t\tsum\t=\tdot(diff,\tdiff);\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble\tweight\t=\tweights[gi];\n"
        "\t\t\t\t\t\t\t\tif\t(weight\t>\t0.0)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\tsigma2Offset\t+=\tweight\t*\tsum;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\ttrans++;\n"
        "\t\t\t\t\t\t\t\tgi\t+=\ttotalOrients;\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%f\t%f\tprior:\t%f\t%f\tdiff:\t%f\t%f\tweight\t%f\ts2offset:\t%f\\n\",\tgi,\tmyTranslation.x,\tmyTranslation.y,\tpriorOffset.x,\tpriorOffset.y,\tdiff.x,\tdiff.y,\tweight,\tsigma2Offset);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\tsigma2OffsetReturn[get_global_id(0)]\t+=\tsigma2Offset;\n"
        "}\n"
        "\n"
        "kernel\tvoid\tcalculate2DPriorOffset(global\tdouble\t*weights,\tconst\tdouble2\tpriorOffset,\tglobal\tdouble2\t*translations,\tglobal\tdouble2\t*priorOffsetReturn,\tconst\tint\ttotalOrients,\tconst\tint\tgi_limit,\tconst\tint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tgi_limit)\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\toffsetSum\t=\t0;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ttrans\t=\t0;\n"
        "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
        "\t\t\t\t\t\t\t\tdouble2\tmyTranslation\t=\ttranslations[trans];\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble2\tsum;\n"
        "\t\t\t\t\t\t\t\tsum\t=\tpriorOffset\t+\tmyTranslation;\n"
        "\t\t\t\t\t\t\t\t\n"
        "\t\t\t\t\t\t\t\tdouble\tweight\t=\tweights[gi];\n"
        "\t\t\t\t\t\t\t\tif\t(weight\t>\t0.0)\t{\n"
        "\t\t\t\t\t\t\t\t\t\t\t\toffsetSum\t+=\tweight\t*\tsum;\n"
        "\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t\t\t\t\ttrans++;\n"
        "\t\t\t\t\t\t\t\tgi\t+=\ttotalOrients;\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%f\t%f\tprior:\t%f\t%f\tdiff:\t%f\t%f\tweight\t%f\ts2offset:\t%f\\n\",\tgi,\tmyTranslation.x,\tmyTranslation.y,\tpriorOffset.x,\tpriorOffset.y,\tdiff.x,\tdiff.y,\tweight,\tsigma2Offset);\n"
        "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
        "\t\t\t\t}\n"
        "\t\t\t\tpriorOffsetReturn[get_global_id(0)]\t+=\toffsetSum;\n"
        "}\n"
        "\n"
        "kernel\tvoid\tshiftImageInFourierTransform(const\tint2\tFDim,\tconst\tdouble\toridim,\tglobal\tdouble2\t*shifts,\tglobal\tdouble2\t*Fin,\tglobal\tdouble2\t*Fout,\tconst\tuint\tlimit)\t{\n"
        "\t\t\t\t\n"
        "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
        "\t\t\t\t\n"
        "\t\t\t\tif\t(gi\t>=\tlimit)\n"
        "\t\t\t\t\t\t\t\treturn;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFSize\t=\tFDim.x\t*\tFDim.y;\n"
        "\t\t\t\tint\ttrans\t=\tgi\t/\tFSize;\n"
        "\t\t\t\tdouble2\tshift\t=\tshifts[trans];\n"
        "\t\t\t\t\n"
        "\t\t\t\tshift\t/=\t-oridim;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\tFIndex\t=\tgi\t-\ttrans\t*\tFSize;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble2\tin\t=\tFin[FIndex];\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ti\t=\tFIndex\t/\tFDim.x;\n"
        "\t\t\t\tint\tx\t=\tFIndex\t-\ti\t*\tFDim.x;\n"
        "\t\t\t\t\n"
        "\t\t\t\tint\ty\t=\t(i\t<\tFDim.x)\t?\ti\t:\ti\t-\tFDim.y;\n"
        "\t\t\t\t\n"
        "\t\t\t\tdouble\tdotp\t=\t2\t*\tM_PI_F\t*\t(x\t*\tshift.x\t+\ty\t*\tshift.y);\n"
        "\t\t\t\tdouble\ta\t=\tcos(dotp);\n"
        "\t\t\t\tdouble\tb\t=\tsin(dotp);\n"
        "\t\t\t\tdouble\tc\t=\tin.x;\n"
        "\t\t\t\tdouble\td\t=\tin.y;\n"
        "\t\t\t\tdouble\tac\t=\ta\t*\tc;\n"
        "\t\t\t\tdouble\tbd\t=\tb\t*\td;\n"
        "\t\t\t\tdouble\tab_cd\t=\t(a\t+\tb)\t*\t(c\t+\td);\n"
        "\t\t\t\tFout[gi]\t=\t(double2)(ac\t-\tbd,\tab_cd\t-\tac\t-\tbd);\n"
        "\t\t\t\t//\t\t\t\tif\t(trans\t==\t1)\n"
        "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ti:\t%d\ty:\t%d\tx:\t%d\ttrans:\t%d\tshift:\t%f\t%f\tin:\t%f\t%f\tdotp:\t%f\ta\tb\tc\td:\t%f\t%f\t%f\t%f\tout:\t%f\t%f\\n\",\tgi,\ti,\ty,\tx,\ttrans,\tshift.x,\tshift.y,\tin.x,\tin.y,\tdotp,\ta,\tb,\tc,\td,\tac\t-\tbd,\tab_cd\t-\tac\t-\tbd);\n"
        "}\n"
    };
    source = clCode;
//    source = loadProgramSource("/em/Applications/relion/relion-1.3.mod/src/ml_optimiser_exp_par_AtomicBP.cl");
//                        source = loadProgramSource("/home/calsmi/relion/relion-1.3.mod/src/ml_optimiser_exp_par_AtomicBP.cl");
    
    
    CL_program = clCreateProgramWithSource(CL_context, 1, (const char **) & source, NULL, &err);
    if (!CL_program || err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create program" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    err = clBuildProgram(CL_program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[32768];
        
        clGetProgramBuildInfo(CL_program, CL_device, CL_PROGRAM_BUILD_LOG, 32768, buffer, &len);
        if (strstr(buffer, "cmpxchg") == NULL) {
            std::cerr << "Error: Failed to build OpenCL executable - build log: " << buffer << std::endl;
            do_use_opencl = false;
            return;
        } else {
            CL_atomicSupport = false;
        }
    }
    
    if (CL_atomicSupport == false) {
        char clCode[] = {"#pragma\tOPENCL\tEXTENSION\tcl_khr_fp64\t:\tenable\n"
            "//#pragma\tOPENCL\tEXTENSION\tcl_khr_int64_base_atomics\t:\tenable\n"
            "\n"
            "//__constant\tsampler_t\tsampler\t=\tCLK_FILTER_NEAREST;\n"
            "\n"
            "#define\tArow\t3\n"
            "#define\tASize\t9\n"
            "#define\tAELEM(m,gi,i,j)\t(m[gi*ASize+i*Arow+j])\n"
            "#define\tAE(m,i,j)\t(m[i*Arow+j])\n"
            "\n"
            "#define\tMODEL_ELEM(m,k,i,j)\t(m[k*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
            "\n"
            "#define\tMODEL_ELEM_2D(m,i,j)\t(m[i*modelDim.x+j])\n"
            "\n"
            "#define\tOUTPUT_MODEL_ELEM(m,n,k,i,j)\t(m[n*modelDim.x*modelDim.y*modelDim.z+k*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
            "\n"
            "#define\tOUTPUT_MODEL_ELEM_2D(m,n,i,j)\t(m[n*modelDim.x*modelDim.y+i*modelDim.x+j])\n"
            "\n"
            "#define\tOUTPUT_ELEM(m,gi,i,j)\t(m[gi*projDim.x*projDim.y+i*projDim.x+j])\n"
            "\n"
            "#define\tOUTPUT_ELEM_TR(m,gi,i,j)\t(m[gi+numOrients*(i*projDim.x+j)])\n"
            "\n"
            "#define\tOUTPUT_ELEM_WGS(m,num,i,j)\t(m[num*wgsPerFref*lgs+i*projDim.x+j])\n"
            "\n"
            "#define\tOUTPUT_ELEM_LINEAR(m,gi,i)\t(m[gi*projDim.x*projDim.y+i])\n"
            "\n"
            "#define\tINPUT_ELEM_LINEAR(m,ipart,num_trans,itrans,i)\t(m[ipart*num_trans*projDim.x*projDim.y+itrans*projDim.x*projDim.y+i])\n"
            "\n"
            "#define\tFLOOR(x)\t(((x)\t==\t(int)(x))\t?\t(int)(x):(((x)\t>\t0)\t?\t(int)(x)\t:\t\\\n"
            "(int)((x)\t-\t1)))\n"
            "\n"
            "#define\tCOMPLEX_MULT(a,\tb)\t(double2)(a.x\t*\tb.x\t-\ta.y\t*\tb.y,\ta.x\t*\tb.y\t+\ta.y\t*\tb.x)\n"
            "#define\tCOMPLEX_NORM(a)\t(dot(a,a))\n"
            "#define\tLIN_INTERP(a,\tl,\th)\t((l)\t+\t((h)\t-\t(l))\t*\t(a))\n"
            "\n"
            "#define\tDEG2RAD(a)\t(a\t*\t2\t*\tM_PI\t/\t360.0)\n"
            "\n"
            "#define\tWGS\t256\n"
            "//#define\tWGS_MULTIPLE\t4\n"
            "//#define\tGPU\n"
            "\n"
            "\n"
            "double\taccumulate(double\tinput,\tlocal\tdouble*\tbuffer)\t{\n"
            "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
            "//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
            "\t\t\t\tsize_t\tls\t=\tWGS;\n"
            "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\tbuffer[li]\t=\tinput;\n"
            "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\n"
            "\t\t\t\t#pragma\tunroll\t1\n"
            "\t\t\t\tfor(unsigned\tint\ts=ls/2;\ts>8;\ts>>=1)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\ts];\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
            "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
            "//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tsum;\n"
            "\t\t\t\tif\t(li\t==\t0)\n"
            "\t\t\t\t\t\t\t\tsum\t=\tbuffer[li]\t+\tbuffer[li\t+\t1]\t+\tbuffer[li\t+\t2]\t+\tbuffer[li\t+\t3]\t+\tbuffer[li\t+\t4]\t+\tbuffer[li\t+\t5]\t+\tbuffer[li\t+\t6]\t+\tbuffer[li\t+\t7]\t+\tbuffer[li\t+\t8]\t+\tbuffer[li\t+\t9]\t+\tbuffer[li\t+\t10]\t+\tbuffer[li\t+\t11]\t+\tbuffer[li\t+\t12]\t+\tbuffer[li\t+\t13]\t+\tbuffer[li\t+\t14]\t+\tbuffer[li\t+\t15];\n"
            "\t\t\t\treturn\tsum;\n"
            "}\n"
            "\n"
            "double2\taccumulate2(double\tinput1,\tlocal\tdouble*\tbuffer1,\tdouble\tinput2,\tlocal\tdouble*\tbuffer2)\t{\n"
            "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
            "\t\t\t\t//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
            "\t\t\t\tsize_t\tls\t=\tWGS;\n"
            "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\tbuffer1[li]\t=\tinput1;\n"
            "\t\t\t\tbuffer2[li]\t=\tinput2;\n"
            "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t\n"
            "#pragma\tunroll\t1\n"
            "\t\t\t\tfor(unsigned\tint\ts=ls/2;\ts>8;\ts>>=1)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tbuffer1[li]\t+=\tbuffer1[li\t+\ts];\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tbuffer2[li]\t+=\tbuffer2[li\t+\ts];\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tsum;\n"
            "\t\t\t\tif\t(li\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\tsum.x\t=\tbuffer1[li]\t+\tbuffer1[li\t+\t1]\t+\tbuffer1[li\t+\t2]\t+\tbuffer1[li\t+\t3]\t+\tbuffer1[li\t+\t4]\t+\tbuffer1[li\t+\t5]\t+\tbuffer1[li\t+\t6]\t+\tbuffer1[li\t+\t7]\t+\tbuffer1[li\t+\t8]\t+\tbuffer1[li\t+\t9]\t+\tbuffer1[li\t+\t10]\t+\tbuffer1[li\t+\t11]\t+\tbuffer1[li\t+\t12]\t+\tbuffer1[li\t+\t13]\t+\tbuffer1[li\t+\t14]\t+\tbuffer1[li\t+\t15];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tsum.y\t=\tbuffer2[li]\t+\tbuffer2[li\t+\t1]\t+\tbuffer2[li\t+\t2]\t+\tbuffer2[li\t+\t3]\t+\tbuffer2[li\t+\t4]\t+\tbuffer2[li\t+\t5]\t+\tbuffer2[li\t+\t6]\t+\tbuffer2[li\t+\t7]\t+\tbuffer2[li\t+\t8]\t+\tbuffer2[li\t+\t9]\t+\tbuffer2[li\t+\t10]\t+\tbuffer2[li\t+\t11]\t+\tbuffer2[li\t+\t12]\t+\tbuffer2[li\t+\t13]\t+\tbuffer2[li\t+\t14]\t+\tbuffer2[li\t+\t15];\n"
            "\t\t\t\t}\n"
            "\t\t\t\treturn\tsum;\n"
            "}\n"
            "\n"
            "/*double\taccumulate(double\tinput,\tlocal\tdouble*\tbuffer)\t{\n"
            "\t\t\t\t//Assume\tlocal\tsize\tis\ta\tpower\tof\ttwo!\n"
            "//\t\t\t\tsize_t\tls\t=\tget_local_size(0);\n"
            "\t\t\t\t\t\t\t\tsize_t\tls\t=\tWGS;\n"
            "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\tbuffer[li]\t=\tinput;\n"
            "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(ls\t>=\t512)\t{\tif\t(li\t<\t256)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t256];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
            "\t\t\t\tif\t(ls\t>=\t256)\t{\tif\t(li\t<\t128)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t128];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
            "\t\t\t\tif\t(ls\t>=\t128)\t{\tif\t(li\t<\t\t64)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t64];\t}\tbarrier(CLK_LOCAL_MEM_FENCE);\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(li\t<\t32)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t64)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t32];\t}\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t32)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t16];\t}\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t16)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t8];\t}\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t8)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t4];\t}\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t4)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t2];\t}\n"
            "\t\t\t\t\t\t\t\tif\t(ls\t>=\t\t\t2)\t{\tbuffer[li]\t+=\tbuffer[li\t+\t\t1];\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\treturn\tbuffer[li];\n"
            "}*/\n"
            "\n"
            "\n"
            "kernel\tvoid\tblank(global\tdouble\t*Fref)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tFref[gi]\t=\t0.0;\n"
            "}\n"
            "\n"
            "kernel\tvoid\tinitWithConstantDouble(global\tdouble\t*Fref,\tconst\tdouble\tvalue,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tFref[gi]\t=\tvalue;\n"
            "}\n"
            "\n"
            "kernel\tvoid\tinitWithConstant(global\tint\t*Fref,\tconst\tint\tvalue,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tFref[gi]\t=\tvalue;\n"
            "}\n"
            "\n"
            "\n"
            "kernel\tvoid\tgenerateOrientationMatrix(global\tdouble*\teulerOrientations,\tglobal\tdouble*\texp_R_mic,\tconst\tint\tIS_INV,\tconst\tdouble\tpaddingFactor,\tglobal\tdouble*\tmatrix,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\to\t=\tgi\t*\t3;\n"
            "\t\t\t\tsize_t\tm\t=\tgi\t*\t9;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tA[9],\tAI[9];\n"
            "\t\t\t\tdouble\tca,\tsa,\tcb,\tsb,\tcg,\tsg;\n"
            "\t\t\t\tdouble\tcc,\tcs,\tsc,\tss;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Convert\teuler\tangles\tto\tmatrix\n"
            "\t\t\t\tdouble\talpha,\tbeta,\tgamma;\n"
            "\t\t\t\t\n"
            "\t\t\t\talpha\t=\teulerOrientations[o];\n"
            "\t\t\t\tbeta\t=\teulerOrientations[o+1];\n"
            "\t\t\t\tgamma\t=\teulerOrientations[o+2];\n"
            "\t\t\t\t\n"
            "\t\t\t\talpha\t=\tDEG2RAD(alpha);\n"
            "\t\t\t\tbeta\t\t=\tDEG2RAD(beta);\n"
            "\t\t\t\tgamma\t=\tDEG2RAD(gamma);\n"
            "\t\t\t\t\n"
            "\t\t\t\tca\t=\tcos(alpha);\n"
            "\t\t\t\tcb\t=\tcos(beta);\n"
            "\t\t\t\tcg\t=\tcos(gamma);\n"
            "\t\t\t\tsa\t=\tsin(alpha);\n"
            "\t\t\t\tsb\t=\tsin(beta);\n"
            "\t\t\t\tsg\t=\tsin(gamma);\n"
            "\t\t\t\tcc\t=\tcb\t*\tca;\n"
            "\t\t\t\tcs\t=\tcb\t*\tsa;\n"
            "\t\t\t\tsc\t=\tsb\t*\tca;\n"
            "\t\t\t\tss\t=\tsb\t*\tsa;\n"
            "\t\t\t\t\n"
            "\t\t\t\tA[0]\t=\t\tcg\t*\tcc\t-\tsg\t*\tsa;\n"
            "\t\t\t\tA[1]\t=\t\tcg\t*\tcs\t+\tsg\t*\tca;\n"
            "\t\t\t\tA[2]\t=\t-cg\t*\tsb;\n"
            "\t\t\t\tA[3]\t=\t-sg\t*\tcc\t-\tcg\t*\tsa;\n"
            "\t\t\t\tA[4]\t=\t-sg\t*\tcs\t+\tcg\t*\tca;\n"
            "\t\t\t\tA[5]\t=\tsg\t*\tsb;\n"
            "\t\t\t\tA[6]\t=\t\tsc;\n"
            "\t\t\t\tA[7]\t=\t\tss;\n"
            "\t\t\t\tA[8]\t=\tcb;\n"
            "\t\t\t\t\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tASize;\ti++)\t{\n"
            "\t\t\t\t\t\t\t\tAI[i]\t=\t0;\n"
            "\t\t\t\t}\n"
            "\n"
            "\t\t\t\t//Mult\twith\texp_R_mic\n"
            "\t\t\t\t#pragma\tunroll\t1\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tArow;\ti++)\n"
            "\t\t\t\t\t\t\t\t#pragma\tunroll\t1\n"
            "\t\t\t\t\t\t\t\tfor\t(int\tj\t=\t0;\tj\t<\tArow;\tj++)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t#pragma\tunroll\t1\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tfor\t(int\tk\t=\t0;\tk\t<\tArow;\tk++)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAE(AI,i,j)\t+=\tAE(exp_R_mic,i,k)\t*\tAE(A,k,j);\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Invert\n"
            "\t\t\t\tA[0]\t=\t\t\tAI[8]*AI[4]-AI[7]*AI[5];\n"
            "\t\t\t\tA[1]\t=\t-(AI[8]*AI[1]-AI[7]*AI[2]);\n"
            "\t\t\t\tA[2]\t=\t\t\tAI[5]*AI[1]-AI[4]*AI[2];\n"
            "\t\t\t\tA[3]\t=\t-(AI[8]*AI[3]-AI[6]*AI[5]);\n"
            "\t\t\t\tA[4]\t=\t\t\tAI[8]*AI[0]-AI[6]*AI[2];\n"
            "\t\t\t\tA[5]\t=\t-(AI[5]*AI[0]-AI[3]*AI[2]);\n"
            "\t\t\t\tA[6]\t=\t\t\tAI[7]*AI[3]-AI[6]*AI[4];\n"
            "\t\t\t\tA[7]\t=\t-(AI[7]*AI[0]-AI[6]*AI[1]);\n"
            "\t\t\t\tA[8]\t=\t\t\tAI[4]*AI[0]-AI[3]*AI[1];\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\ttmp;\n"
            "\t\t\t\ttmp\t=\tAI[0]\t*\tA[0]\t+\tAI[3]\t*\tA[1]\t+\tAI[6]\t*\tA[2];\n"
            "\t\t\t\t\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\t9;\ti++)\t{\n"
            "\t\t\t\t\t\t\t\tA[i]\t/=\ttmp;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Now\ttranspose\tif\tnot\tA_INV\n"
            "\t\t\t\tAI[0]\t=\tIS_INV\t?\tA[0]\t:\tA[0];\n"
            "\t\t\t\tAI[1]\t=\tIS_INV\t?\tA[1]\t:\tA[3];\n"
            "\t\t\t\tAI[2]\t=\tIS_INV\t?\tA[2]\t:\tA[6];\n"
            "\t\t\t\tAI[3]\t=\tIS_INV\t?\tA[3]\t:\tA[1];\n"
            "\t\t\t\tAI[4]\t=\tIS_INV\t?\tA[4]\t:\tA[4];\n"
            "\t\t\t\tAI[5]\t=\tIS_INV\t?\tA[5]\t:\tA[7];\n"
            "\t\t\t\tAI[6]\t=\tIS_INV\t?\tA[6]\t:\tA[2];\n"
            "\t\t\t\tAI[7]\t=\tIS_INV\t?\tA[7]\t:\tA[5];\n"
            "\t\t\t\tAI[8]\t=\tIS_INV\t?\tA[8]\t:\tA[8];\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Multiply\tby\tpadding\tfactor\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\t9;\ti++)\t{\n"
            "\t\t\t\t\t\t\t\tmatrix[m+i]\t=\tAI[i]\t*\tpaddingFactor;\n"
            "\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tscaleMatrix(global\tdouble*\tmatrix,\tdouble\tscaleFactor,\tint\tlimit)\t{\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\tmatrix[gi]\t*=\tscaleFactor;\n"
            "\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalcModelProjection_TR(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tint\tnumOrients,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
            "\t\t\t\tint\tis_neg_x;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\t\t\t\tdouble2\td000,\td001,\td010,\td011,\td100,\td101,\td110,\td111;\n"
            "\t\t\t\tdouble2\tdx00,\tdx01,\tdx10,\tdx11,\tdxy0,\tdxy1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
            "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ti;\n"
            "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tx;\n"
            "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
            "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
            "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
            "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
            "\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tz0\t=\tFLOOR(zp);\n"
            "\t\t\t\tfz\t=\tzp\t-\tz0;\n"
            "\t\t\t\tz0\t-=\tSTARTINGZ;\n"
            "\t\t\t\tz1\t=\tz0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
            "\t\t\t\t\t\t\t\td000\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td001\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td010\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td011\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\td100\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td101\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td110\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td111\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t1)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\td000:\t%f\t%f\td001:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\td000.x,\td000.y,\td001.x,\td001.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
            "\t\t\t\t\t\t\t\tdx00\t=\tLIN_INTERP(fx,\td000,\td001);\n"
            "\t\t\t\t\t\t\t\tdx01\t=\tLIN_INTERP(fx,\td100,\td101);\n"
            "\t\t\t\t\t\t\t\tdx10\t=\tLIN_INTERP(fx,\td010,\td011);\n"
            "\t\t\t\t\t\t\t\tdx11\t=\tLIN_INTERP(fx,\td110,\td111);\n"
            "\t\t\t\t\t\t\t\tdxy0\t=\tLIN_INTERP(fy,\tdx00,\tdx10);\n"
            "\t\t\t\t\t\t\t\tdxy1\t=\tLIN_INTERP(fy,\tdx01,\tdx11);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
            "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fz,\tdxy0,\tdxy1)\t*\t(double2)(1,-1);\n"
            "\t\t\t\t\t\t\t\tOUTPUT_ELEM_TR(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalcModelProjection(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
            "\t\t\t\tint\tis_neg_x;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\t\t\t\tdouble2\td000,\td001,\td010,\td011,\td100,\td101,\td110,\td111;\n"
            "\t\t\t\tdouble2\tdx00,\tdx01,\tdx10,\tdx11,\tdxy0,\tdxy1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tuint\tprojSize\t=\tprojDim.x\t*\tprojDim.y;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tvalid\t=\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ti;\n"
            "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tx;\n"
            "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
            "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
            "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
            "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
            "\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tz0\t=\tFLOOR(zp);\n"
            "\t\t\t\tfz\t=\tzp\t-\tz0;\n"
            "\t\t\t\tz0\t-=\tSTARTINGZ;\n"
            "\t\t\t\tz1\t=\tz0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
            "\t\t\t\t\t\t\t\td000\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td001\t=\tMODEL_ELEM(model,\tz0,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td010\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td011\t=\tMODEL_ELEM(model,\tz0,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\td100\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td101\t=\tMODEL_ELEM(model,\tz1,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td110\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td111\t=\tMODEL_ELEM(model,\tz1,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t1)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\td000:\t%f\t%f\td001:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\td000.x,\td000.y,\td001.x,\td001.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
            "\t\t\t\t\t\t\t\tdx00\t=\tLIN_INTERP(fx,\td000,\td001);\n"
            "\t\t\t\t\t\t\t\tdx01\t=\tLIN_INTERP(fx,\td100,\td101);\n"
            "\t\t\t\t\t\t\t\tdx10\t=\tLIN_INTERP(fx,\td010,\td011);\n"
            "\t\t\t\t\t\t\t\tdx11\t=\tLIN_INTERP(fx,\td110,\td111);\n"
            "\t\t\t\t\t\t\t\tdxy0\t=\tLIN_INTERP(fy,\tdx00,\tdx10);\n"
            "\t\t\t\t\t\t\t\tdxy1\t=\tLIN_INTERP(fy,\tdx01,\tdx11);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
            "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fz,\tdxy0,\tdxy1)\t*\t(double2)(1,-1);\n"
            "\t\t\t\t\t\t\t\tOUTPUT_ELEM(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalcModelProjection2D_TR(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tint\tnumOrients,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
            "\t\t\t\tint\tis_neg_x;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\t\t\t\tdouble2\td00,\td01,\td10,\td11,\tdx0,\tdx1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
            "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ti;\n"
            "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tx;\n"
            "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
            "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
            "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
            "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
            "\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
            "\t\t\t\t\t\t\t\td00\t=\tMODEL_ELEM_2D(model,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td01\t=\tMODEL_ELEM_2D(model,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td10\t=\tMODEL_ELEM_2D(model,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td11\t=\tMODEL_ELEM_2D(model,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
            "\t\t\t\t\t\t\t\tdx0\t=\tLIN_INTERP(fx,\td00,\td01);\n"
            "\t\t\t\t\t\t\t\tdx1\t=\tLIN_INTERP(fx,\td10,\td11);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
            "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fy,\tdx0,\tdx1)\t*\t(double2)(1,-1);\n"
            "\t\t\t\t\t\t\t\tOUTPUT_ELEM_TR(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fy,\tdx0,\tdx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalcModelProjection2D(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tmodel,\tglobal\tdouble2*\toutput,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty;\n"
            "\t\t\t\tint\tis_neg_x;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\t\t\t\tdouble2\td00,\td01,\td10,\td11,\tdx0,\tdx1;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projDim.x\t*\tprojDim.y);\n"
            "\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojDim.x\t*\tprojDim.y);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ti;\n"
            "\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tx;\n"
            "\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
            "\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t1;\n"
            "\t\t\t\t}\n"
            "\t\t\t\telse\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tis_neg_x\t=\t0;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
            "\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
            "\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
            "\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t//\tMatrix\taccess\tcan\tbe\taccelerated\tthrough\tpre-calculation\tof\tz0*xydim\tetc.\n"
            "\t\t\t\t\t\t\t\td00\t=\tMODEL_ELEM_2D(model,\ty0,\tx0);\n"
            "\t\t\t\t\t\t\t\td01\t=\tMODEL_ELEM_2D(model,\ty0,\tx1);\n"
            "\t\t\t\t\t\t\t\td10\t=\tMODEL_ELEM_2D(model,\ty1,\tx0);\n"
            "\t\t\t\t\t\t\t\td11\t=\tMODEL_ELEM_2D(model,\ty1,\tx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tSet\tthe\tinterpolated\tvalue\tin\tthe\t2D\toutput\tarray\n"
            "\t\t\t\t\t\t\t\tdx0\t=\tLIN_INTERP(fx,\td00,\td01);\n"
            "\t\t\t\t\t\t\t\tdx1\t=\tLIN_INTERP(fx,\td10,\td11);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tTake\tcomplex\tconjugated\tfor\thalf\twith\tnegative\tx\n"
            "\t\t\t\t\t\t\t\tdouble2\tconj\t=\tLIN_INTERP(fy,\tdx0,\tdx1)\t*\t(double2)(1,-1);\n"
            "\t\t\t\t\t\t\t\tOUTPUT_ELEM(output,FrefNo,i,x)\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fy,\tdx0,\tdx1);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\tdx00:\t%f\t%f\tdx01:\t%f\t%f\tdx10:\t%f\t%f\tdx11:\t%f\t%f\tdxy0:\t%f\t%f\tdxy1:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\tdx00.x,\tdx00.y,\tdx01.x,\tdx01.y,\tdx10.x,\tdx10.y,\tdx11.x,\tdx11.y,\tdxy0.x,\tdxy0.y,\tdxy1.x,\tdxy1.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"CL\tx0:\t%d\ty0:\t%d\tz0:\t%d\ti:\t%d\tx:\t%d\tresult:\t%f\t%f\\n\",\tx0,\ty0,\tz0,\ti,\tx,\tLIN_INTERP(fz,\tdxy0,\tdxy1).x,\tLIN_INTERP(fz,\tdxy0,\tdxy1).y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tv\t=\tis_neg_x\t?\tconj\t:\tLIN_INTERP(fz,\tdxy0,\tdxy1);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"pos:\t%3d\ti:\t%d\tx:\t%d\tFref:\t%f\t%f\\n\",\ti*projDim.x+x,\ti,\tx,\tv.x,\tv.y);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "\n"
            "kernel\tvoid\tctfAndScaleDataPoint(global\tdouble2*\tFref,\tconst\tint\tdo_ctf_correction,\tglobal\tdouble*\tctf,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble2*\tFrefOut,\tconst\tint\tprojSize,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\tint\tFctfPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tref;\n"
            "\t\t\t\tref\t=\tFref[gi];\n"
            "//\t\t\t\tref\t=\t(double2)(1,1);\n"
            "\t\t\t\tref\t=\tdo_ctf_correction\t?\tref\t*\tctf[FctfPos]\t:\tref;\n"
            "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
            "\t\t\t\tFrefOut[gi]\t=\tref;\n"
            "//\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
            "//\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tctfAndScaleDataPoint_TR(global\tdouble2*\tFref,\tconst\tint\tdo_ctf_correction,\tglobal\tdouble*\tctf,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble2*\tFrefOut,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFctfPos\t=\tgi\t/\tnumOrients;\n"
            "//\t\t\t\tint\tFrefNo\t=\tgi\t-\tFctfPos\t*\tnumOrients;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tref;\n"
            "\t\t\t\tref\t=\tFref[gi];\n"
            "\t\t\t\t//\t\t\t\tref\t=\t(double2)(1,1);\n"
            "\t\t\t\tref\t=\tdo_ctf_correction\t?\tref\t*\tctf[FctfPos]\t:\tref;\n"
            "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
            "\t\t\t\tFrefOut[gi]\t=\tref;\n"
            "//\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%3d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tscaleCTF(global\tdouble*\tFref,\tconst\tint\tdo_scale_correction,\tconst\tdouble\tscale,\tglobal\tdouble*\tFrefOut,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tref;\n"
            "\t\t\t\tref\t=\tFref[gi];\n"
            "\t\t\t\t//\t\t\t\tref\t=\t(double2)(1,1);\n"
            "\t\t\t\tref\t=\tdo_scale_correction\t?\tref\t*\tscale\t:\tref;\n"
            "\t\t\t\tFrefOut[gi]\t=\tref;\n"
            "\t\t\t\t//\t\t\t\tif\t(FrefNo\t==\t0)\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tctfPos:\t%d\tctf:\t%f\tFref:\t%f\t%f\\n\",\tgi,\tFctfPos,\tctf[FctfPos],\tref.x,\tref.y);\n"
            "\t\t\t\t//\t\t\t\t}\n"
            "}\n"
            "\n"
            "//Here\tglobal\tsize\tis\tnumOrientations\t*\tnumTrans\n"
            "kernel\tvoid\tcalculateDiff2AndSuma2CC_TR(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tdiff2,\tglobal\tdouble*\tsuma2,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tglobal\tchar\t*skip,\tconst\tint\tipart,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(skip[gi])\t{\n"
            "//\t\t\t\t\t\t\t\tdiff2[gi]\t=\t0;\n"
            "//\t\t\t\t\t\t\t\tsuma2[gi]\t=\t0;\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\tnumOrients;\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tnumOrients);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFimgBase\t=\t(trans\t*\tfSize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\td2,\ts2;//,\tdiff2Sum,\tsuma2Sum;\n"
            "\t\t\t\td2\t=\t0;\n"
            "\t\t\t\ts2\t=\t0;\n"
            "\n"
            "\t\t\t\tfor\t(int\tFimgPos\t=\t0;\tFimgPos\t<\tfSize\t;\tFimgPos++)\t{\n"
            "\t\t\t\t\t\t\t\tdouble2\tref,\timg;\n"
            "\t\t\t\t\t\t\t\tref\t=\tFref[ti\t+\tFimgPos\t*\tnumOrients];\n"
            "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos\t+\tFimgBase];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\td2\t+=\tdot(ref,img);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\ts2\t+=\tdot(ref,ref);\n"
            "\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\tif\t((ti\t==\t0)\t&&\t(trans\t<\t2))\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFimgpos:\t%d\timgPosBase:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\ts2:\t%f\td2:\t%f\\n\",\tgi,\ttrans,\tti,\tFimgPos,\tFimgBase,\tref.x,\tref.y,\timg.x,\timg.y,\ts2,\td2);\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tdiff2[gi]\t=\td2;\n"
            "\t\t\t\tsuma2[gi]\t=\ts2;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\t\t\t\tif\t(((FrefNo\t<2\t)\t&&\t(trans\t<2\t))\t||\t((FrefNo\t==\t4607)\t&&\t(trans\t==\t19)))\t{\n"
            "\t\t\t\t//\t\t\t\tif\t(isnan(diff2Sum)\t||\tisnan(suma2Sum))\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tindex:\t%d\ttrans:\t%d\tti:\t%d\trefno:\t%d\timgPos:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\td2:\t%f\ts2:\t%f\tdiff2S:\t%f\ts2S:\t%f\\n\",\tgi,\tindex,\ttrans,\tti,\tFrefNo,\tFimgPos,\tref.x,\tref.y,\timg.x,\timg.y,\td2,\ts2\t*\t1000000,\tdiff2Sum,\tsuma2Sum);\n"
            "\t\t\t\t//\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculateDiff2AndSuma2CC(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tdiff2,\tglobal\tdouble*\tsuma2,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tref,\timg;\n"
            "\t\t\t\tref\t=\tFref[ti];\n"
            "\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\n"
            "//\t\t\t\tdiff2Sum\t=\tref.x\t*\timg.x;\n"
            "//\t\t\t\tdiff2Sum\t+=\tref.y\t*\timg.y;\n"
            "\n"
            "//\t\t\t\tdiff2[gi]\t=\tdot(ref,img);\n"
            "//\t\t\t\tsuma2[gi]\t=\tdot(ref,ref);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\td2,\ts2;//,\tdiff2Sum,\tsuma2Sum;\n"
            "\t\t\t\t\n"
            "\t\t\t\td2\t=\tdot(ref,img);\n"
            "//\t\t\t\tdiff2Sum\t=\taccumulate(d2,\tbuffer);\n"
            "\n"
            "\t\t\t\ts2\t=\tdot(ref,ref);\n"
            "//\t\t\t\tsuma2Sum\t=\taccumulate(s2,\tbuffer);\n"
            "\t\t\t\t\n"
            "#ifdef\tGPU\n"
            "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\t__local\tdouble\tbuffer1[WGS];\n"
            "\t\t\t\t__local\tdouble\tbuffer2[WGS];\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tsum;\n"
            "\t\t\t\tsum\t=\taccumulate2(d2,buffer1,\ts2,\tbuffer2);\n"
            "#endif\n"
            "\n"
            "//\t\t\t\tsize_t\tindex\t=\tFrefNo\t*\twgsPerFref\t+\tFimgPos\t/\tlgs\t+\ttrans\t*\twgsPerFref\t*\torientationsPerTrans;\n"
            "\t\t\t\tsize_t\tindex\t=\tgi\t/\tlgs;\n"
            "#ifdef\tGPU\n"
            "\t\t\t\tif\t(li\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\tdiff2[index]\t=\tsum.x;\n"
            "\t\t\t\t\t\t\t\tsuma2[index]\t=\tsum.y;\n"
            "\t\t\t\t}\n"
            "#else\n"
            "\t\t\t\tdiff2[index]\t=\td2;\n"
            "\t\t\t\tsuma2[index]\t=\ts2;\n"
            "#endif\n"
            "\t\t\t\t\n"
            "//\t\t\t\tif\t(((FrefNo\t<2\t)\t&&\t(trans\t<2\t))\t||\t((FrefNo\t==\t4607)\t&&\t(trans\t==\t19)))\t{\n"
            "//\t\t\t\tif\t(isnan(diff2Sum)\t||\tisnan(suma2Sum))\t{\n"
            "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tindex:\t%d\ttrans:\t%d\tti:\t%d\trefno:\t%d\timgPos:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\td2:\t%f\ts2:\t%f\tdiff2S:\t%f\ts2S:\t%f\\n\",\tgi,\tindex,\ttrans,\tti,\tFrefNo,\tFimgPos,\tref.x,\tref.y,\timg.x,\timg.y,\td2,\ts2\t*\t1000000,\tdiff2Sum,\tsuma2Sum);\n"
            "//\t\t\t\t}\n"
            "}\n"
            "\n"
            "/*inline\tvoid\tAtomicAdd(volatile\t__global\tdouble\t*source,\tconst\tdouble\toperand)\t{\n"
            "\t\t\t\tunion\t{\n"
            "\t\t\t\t\t\t\t\tunsigned\tlong\tintVal;\n"
            "\t\t\t\t\t\t\t\tdouble\tfloatVal;\n"
            "\t\t\t\t}\tnewVal;\n"
            "\t\t\t\tunion\t{\n"
            "\t\t\t\t\t\t\t\tunsigned\tlong\tintVal;\n"
            "\t\t\t\t\t\t\t\tdouble\tfloatVal;\n"
            "\t\t\t\t}\tprevVal;\n"
            "\t\t\t\tdo\t{\n"
            "\t\t\t\t\t\t\t\tprevVal.floatVal\t=\t*source;\n"
            "\t\t\t\t\t\t\t\tnewVal.floatVal\t=\tprevVal.floatVal\t+\toperand;\n"
            "\t\t\t\t}\twhile\t(atomic_cmpxchg((volatile\t__global\tunsigned\tlong\t*)source,\tprevVal.intVal,\tnewVal.intVal)\t!=\tprevVal.intVal);\n"
            "}*/\n"
            "\n"
            "kernel\tvoid\tcalculateDiff2AndSuma2_TR(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tdiff2,\tconst\tint\tfSize,\tconst\tuint\tnumOrients,\tglobal\tchar\t*skip,\tconst\tint\tipart,\t\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t/*size_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\t__local\tdouble\tlocalMinSigma2[WGS\t*\tWGS_MULTIPLE];\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\tWGS_MULTIPLE;\ti++)\t{\n"
            "\t\t\t\t\t\t\t\tint\tindex\t=\ti\t*\tWGS\t+\tli;\n"
            "\t\t\t\t\t\t\t\tif\t(index\t<\tfSize)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tlocalMinSigma2[index]\t=\t0.5\t*\tMinvsigma2[index];\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);*/\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(skip[gi])\t{\n"
            "//\t\t\t\t\t\t\t\tdiff2[gi]\t=\t0;\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\tnumOrients;\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tnumOrients);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFimgBase\t=\t(trans\t*\tfSize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\td2;//,\tdiff2Sum,\tsuma2Sum;\n"
            "\t\t\t\td2\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tfor\t(int\tFimgPos\t=\t0;\tFimgPos\t<\tfSize;\tFimgPos++)\t{\n"
            "\t\t\t\t\t\t\t\tdouble2\tref,\timg,\tdiff;\n"
            "\t\t\t\t\t\t\t\tref\t=\tFref[ti\t+\tFimgPos\t*\tnumOrients];\n"
            "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos\t+\tFimgBase];\n"
            "\n"
            "\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\td2\t+=\tdot(diff,diff)\t*\tMinvsigma2[FimgPos];\n"
            "\n"
            "//\t\t\t\t\t\t\t\td2\t+=\tdot(diff,diff)\t*\tlocalMinSigma2[FimgPos];\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t((ipart\t==\t1)\t&&\t(get_global_id(0)\t==\t0))\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFimgpos:\t%d\timgPosBase:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\tdiff:\t%f\t%f\td2:\t%f\tMin2:\t%f\\n\",\tgi,\ttrans,\tti,\tFimgPos,\tFimgBase,\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\td2,\tMinvsigma2[FimgPos]);\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tdiff2[gi]\t=\td2;\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculateDiff2AndSuma2(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tdiff2,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
            "\t\t\t\tref\t=\tFref[ti];\n"
            "\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\td2,\tdiff2Sum;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdiff\t=\tref\t-\timg;\n"
            "\t\t\t\t\n"
            "\t\t\t\td2\t=\tdot(diff,diff)\t*\t0.5\t*\tMinvsigma2[FimgPosMin];\n"
            "//\t\t\t\tdiff2Sum\t=\taccumulate(d2,\tbuffer);\n"
            "\n"
            "#ifdef\tGPU\n"
            "\t\t\t\tsize_t\tli\t=\tget_local_id(0);\n"
            "\t\t\t\t__local\tdouble\tbuffer[WGS];\n"
            "\t\t\t\tbuffer[li]\t=\td2;\n"
            "\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t\n"
            "#pragma\tunroll\t1\n"
            "\t\t\t\tfor(unsigned\tint\ts=WGS/2;\ts>16;\ts>>=1)\n"
            "\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\tif\t(li\t<\ts)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\ts];\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\tbarrier(CLK_LOCAL_MEM_FENCE);\n"
            "\t\t\t\t}\n"
            "\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t1];\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t2];\n"
            "\t\t\t\t//\t\t\t\tbuffer[li]\t+=\tbuffer[li\t+\t3];\n"
            "\t\t\t\t\n"
            "\t\t\t\t\n"
            "//\t\t\t\tsize_t\tindex\t=\tFrefNo\t*\twgsPerFref\t+\tFimgPos\t/\tlgs;\n"
            "\t\t\t\tsize_t\tindex\t=\t(gi\t/\tlgs);\n"
            "\t\t\t\tif\t(li\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\tAtomicAdd(&diff2[index],\tsum);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdiff2[index]\t=\tbuffer[li]\t+\tbuffer[li\t+\t1]\t+\tbuffer[li\t+\t2]\t+\tbuffer[li\t+\t3]\t+\tbuffer[li\t+\t4]\t+\tbuffer[li\t+\t5]\t+\tbuffer[li\t+\t6]\t+\tbuffer[li\t+\t7]\t+\tbuffer[li\t+\t8]\t+\tbuffer[li\t+\t9]\t+\tbuffer[li\t+\t10]\t+\tbuffer[li\t+\t11]\t+\tbuffer[li\t+\t12]\t+\tbuffer[li\t+\t13]\t+\tbuffer[li\t+\t14]\t+\tbuffer[li\t+\t15]\t+\tbuffer[li\t+\t16]\t+\tbuffer[li\t+\t17]\t+\tbuffer[li\t+\t18]\t+\tbuffer[li\t+\t19]\t+\tbuffer[li\t+\t20]\t+\tbuffer[li\t+\t21]\t+\tbuffer[li\t+\t22]\t+\tbuffer[li\t+\t23]\t+\tbuffer[li\t+\t24]\t+\tbuffer[li\t+\t25]\t+\tbuffer[li\t+\t26]\t+\tbuffer[li\t+\t27]\t+\tbuffer[li\t+\t28]\t+\tbuffer[li\t+\t29]\t+\tbuffer[li\t+\t30]\t+\tbuffer[li\t+\t31];\n"
            "\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\tdiff2[index]\t=\tbuffer[li];\n"
            "\t\t\t\t}\n"
            "#else\n"
            "\t\t\t\tsize_t\tindex\t=\t(gi\t/\tlgs);\n"
            "\t\t\t\tdiff2[index]\t=\td2;\n"
            "#endif\n"
            "\t\t\t\t\n"
            "//\t\t\t\tif\t(gi\t<\twgsPerFref\t*\tlgs\t*\t2)\t{\n"
            "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%d\tti:\t%d\tFrefNo:\t%d\timgPos:\t%d\timgPosMin:\t%d\tref:\t%f\t%f\timg:\t%f\t%f\tdiff:\t%f\t%f\td2:\t%f\tMin2:\t%f\tdiff2S:\t%f\tindex:\t%d\\n\",\tgi,\ttrans,\tti,\tFrefNo,\tFimgPos,\tFimgPosMin,\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\td2,\tMinvsigma2[FimgPosMin],\tdiff2Sum,\tindex);\n"
            "//\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculateSWSNoiseEstNormCorrection(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tsigma2_noise,\tglobal\tdouble*\tdata_vs_prior,\tglobal\tdouble*\tsumXA,\tglobal\tdouble*\tsumAA,\tconst\tint\tprojSize,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(projSize);\n"
            "\t\t\t\tint\tsaveNoiseCalculations,\tsaveNormCalculations;\n"
            "\n"
            "\t\t\t\tdouble\tsum\t=\t0,\txa\t=\t0,\taa\t=\t0;\n"
            "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(projSize);\n"
            "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\tint\tMresolFine\t=\tMresol_fine[FimgPosMin];\n"
            "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tprojSize)\t+\t(trans\t*\tprojSize);\n"
            "\t\t\t\t\n"
            "//\t\t\t\tif\t(MresolFine\t>\t-1)\t{\n"
            "\t\t\t\tsaveNoiseCalculations\t=\tMresolFine\t>\t-1\t?\t1\t:\t0;\n"
            "//\t\t\t\tsaveCalculations\t=\tMresolFine\t>\t-1\t?\t1\t:\t0;\n"
            "\t\t\t\tsaveNormCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\tsaveNoiseCalculations\t:\t0;\n"
            "//\t\t\t\t\t\t\t\tsaveCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\t1\t:\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\tdouble\tcurrentWeight\t=\tweight[weightIndex];\n"
            "\t\t\t\t\t\t\t\tif\t(currentWeight\t>\t0)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(projSize);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tprojSize);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
            "\t\t\t\t\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\tdouble\tw2;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\tw2\t=\tweight[weightIndex]\t*\tdot(diff,\tdiff);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tsum\t+=\tMresolFine\t>\t-1\t?\tweight[weightIndex]\t*\tdot(diff,\tdiff)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tsum\t+=\tsaveNoiseCalculations\t?\tcurrentWeight\t*\tdot(diff,\tdiff)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\txa\t+=\tsaveNormCalculations\t?\tcurrentWeight\t*\tdot(ref,img)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\taa\t+=\tsaveNormCalculations\t?\tcurrentWeight\t*\tdot(ref,ref)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(currentWeight\t>\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(isnan(sum))\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\tsum:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tcurrentWeight,\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff),\tsum);\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
            "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
            "\t\t\t\t}\n"
            "//\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tsigma2_noise[get_global_id(0)]\t=\tsum;\n"
            "\t\t\t\tsumXA[get_global_id(0)]\t=\txa;\n"
            "\t\t\t\tsumAA[get_global_id(0)]\t=\taa;\n"
            "\n"
            "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "}\n"
            "\n"
            "\n"
            "kernel\tvoid\tcalculateSWSNoiseEstimate(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tsigma2_noise,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tsum\t=\t0;\n"
            "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(lgs\t*\twgsPerFref);\n"
            "\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tMresolFine\t=\tMresol_fine[FimgPosMin];\n"
            "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble2\tdiff,\tref,\timg;\n"
            "\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
            "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\tdouble\tw2;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdiff\t=\tref\t-\timg;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\tw2\t=\tweight[weightIndex]\t*\tdot(diff,\tdiff);\n"
            "\t\t\t\t\t\t\t\tsum\t+=\tMresolFine\t>\t-1\t?\tweight[weightIndex]\t*\tdot(diff,\tdiff)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
            "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
            "\t\t\t\t}\n"
            "\n"
            "\t\t\t\tsigma2_noise[get_global_id(0)]\t=\tsum;\n"
            "\n"
            "//\t\t\t\tif\t(sum\t>\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculateSWSNormCorrection(global\tdouble2*\tFref,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tweight,\tglobal\tint*\tMresol_fine,\tglobal\tdouble*\tdata_vs_prior,\tglobal\tdouble*\tsumXA,\tglobal\tdouble*\tsumAA,\tconst\tint\twgsPerFref,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\tint\tglobalStep\t=\tget_global_size(0)\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tsaveCalculations;\n"
            "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(lgs\t*\twgsPerFref);\n"
            "\n"
            "\t\t\t\tdouble\txa\t=\t0,\taa\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref)\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\n"
            "\t\t\t\tsaveCalculations\t=\tMresol_fine[FimgPosMin]\t>\t-1\t?\t1\t:\t0;\n"
            "\t\t\t\tsaveCalculations\t=\tdata_vs_prior[FimgPosMin]\t>\t3.0\t?\tsaveCalculations\t:\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(lgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tlgs\t*\twgsPerFref\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(lgs\t*\twgsPerFref);\n"
            "//\t\t\t\t\t\t\t\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tlgs\t*\twgsPerFref);\n"
            "\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tlgs\t*\twgsPerFref);\n"
            "\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble2\tref,\timg;\n"
            "\t\t\t\t\t\t\t\tref\t=\tFref[ti];\n"
            "\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\txa\t+=\tsaveCalculations\t?\tweight[weightIndex]\t*\tdot(ref,img)\t:\t0;\n"
            "\t\t\t\t\t\t\t\taa\t+=\tsaveCalculations\t?\tweight[weightIndex]\t*\tdot(ref,ref)\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
            "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tsumXA[get_global_id(0)]\t=\txa;\n"
            "\t\t\t\tsumAA[get_global_id(0)]\t=\taa;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "}\n"
            "\n"
            "//Global\tsize\there\tis\tno\tof\trotations\t*\twgs\t*\twgsPerFref\n"
            "//Sum\tfor\teach\tF\tfor\tall\ttranslations\tover\tone\torientation\n"
            "kernel\tvoid\tcalculateSWSSum(global\tdouble*\tctf,\tglobal\tdouble2*\tFimg,\tglobal\tdouble*\tMinvsigma2,\tglobal\tdouble*\tweight,\tglobal\tdouble2*\tFimgOut,\tglobal\tdouble*\tFweight,\tconst\tuint\tprojSize,\tconst\tint\torientationsPerTrans,\tconst\tint\tlimit,\tconst\tuint\tgi_limit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tgi_limit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tglobalStep\t=\tgi_limit\t/\t(projSize);\n"
            "\n"
            "\t\t\t\tdouble2\tsum;//\t=\t(double2)(0,\t0);\n"
            "\t\t\t\tsum\t=\tFimgOut[get_global_id(0)];\n"
            "\t\t\t\tdouble\tweightSum;//\t=\t0;\n"
            "\t\t\t\tweightSum\t=\t\tFweight[get_global_id(0)];\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tweightIndex\t=\tgi\t/\t(projSize);\n"
            "\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
            "\t\t\t\tint\tFrefNo\t=\tti\t/\t(projSize);\n"
            "\t\t\t\tint\tFimgPosMin\t=\tti\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\tint\tFimgPos\t=\tti\t-\t(FrefNo\t*\tprojSize)\t+\t(trans\t*\tprojSize);\n"
            "\n"
            "\t\t\t\tdouble\tmyctf\t=\tctf[FimgPosMin];\n"
            "\t\t\t\tdouble\tMinsigv2\t=\tMinvsigma2[FimgPosMin];\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Calculate\tone\tvalue\tof\tFimg\tand\tFweight\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0.0)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ttrans\t=\tgi\t/\t(projSize\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tti\t=\tgi\t-\t(trans\t*\tprojSize\t*\torientationsPerTrans);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tFrefNo\t=\tti\t/\t(projSize);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tFimgPos\t=\tFimgPosMin\t+\t(trans\t*\tprojSize);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tdouble2\timg;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\timg\t=\tFimg[FimgPos];\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tdouble\tweightxinvsigma2;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tweightxinvsigma2\t=\tweight[weightIndex]\t*\tmyctf\t*\tMinsigv2;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tsum\t+=\timg\t*\tweightxinvsigma2;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tweightSum\t+=\tweightxinvsigma2\t*\tmyctf;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t421)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\timgPOosMin:\t%d,\tctf:\t%f,\tMinsigv2:\t%f,\tFimg:\t%f\t%f\tSum:\t%f\t%f\twSum:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tFimgPosMin,\tmyctf,\tMinvsigma2[FimgPosMin],\timg.x,\timg.y,\tsum.x,\tsum.y,\tweightSum);\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(weight[weightIndex]\t>\t0)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d,\tOrient:\t%d,\tTrans:\t%d,\tweightIndex:\t%d,\tweight:\t%f,\tMr:\t%d,\tFref:\t%f\t%f\tFimg:\t%f\t%f\tdiff:\t%f\t%f\twdiff2:\t%f\\n\",\tgi,\tFrefNo,\ttrans,\tweightIndex,\tweight[weightIndex],\tMresol_fine[FimgPosMin],\tref.x,\tref.y,\timg.x,\timg.y,\tdiff.x,\tdiff.y,\tdot(diff,\tdiff));\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\tgi\t+=\torientationsPerTrans\t*\tprojSize;\n"
            "\t\t\t\t\t\t\t\tweightIndex\t+=\tglobalStep;\n"
            "//\t\t\t\t\t\t\t\tweightIndex\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tFimgOut[get_global_id(0)]\t=\tsum;\n"
            "\t\t\t\tFweight[get_global_id(0)]\t=\tweightSum;\n"
            "//\t\t\t\tFimgOut[get_global_id(0)]\t=\t1.0;\n"
            "//\t\t\t\tFweight[get_global_id(0)]\t=\t2.0;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\t\t\t\tif\t(sum\t>\t0)\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\tsum:\t%f\\n\",\tget_global_id(0),\tsum);\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "}\n"
            "\n"
            "\n"
            "//Backproject\n"
            "//Two\tchoices\thow\tto\tdo\tthis\n"
            "//A)\tGlobal\tsize\tis\trelated\tto\tcompute\tunits\n"
            "//\t\t\tLoop\tover\tthe\tinput\tF's\t-->\tmaps\teach\tF\tand\tweight\tfrom\ta\trotation\tto\ta\twork-group\titem\n"
            "//\t\t\tEach\twork\tgroup\titem\tthen\tneeds\tto\twrite\t8\tweights\tand\t8\tvalues\tto\tmodel\n"
            "//\t\t\tModels\tare\tfar\ttoo\tlarge\tfor\tlocal\tmemory\n"
            "//\t\t\tTo\tsynchronise\tthe\twrites\twill\tneed\tto\thave\tseparate\tmodel\tfor\teach\twork\titem\t-\tbetter\tthan\tdouble\tatomic\toperations\n"
            "//\t\t\tTo\tbe\tpractical\twill\thave\tto\treduce\twork\tgroup\tsize\tto\t1...\n"
            "\n"
            "//B)\tGlobal\tsize\there\tis\tthe\tmodel\tsize\t-->\tmaps\teach\tweight\tand\tvalue\tof\tmodel\tto\ta\twork\titem\n"
            "//\t\t\tSimple\tway\tto\timplement\tthis\tis\tthat\teach\twork\titem\tloops\tthrough\torientations\tand\tstores\tvalues\tthat\tmatch\tthe\tcurrent\tmodel\tposition\n"
            "\n"
            "//This\tis\toption\tA\n"
            "kernel\tvoid\tcalcBackProjection(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble4*\tmodelValues,\tconst\tuint\tprojSize,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
            "\t\t\t\tbool\tis_neg_x;\n"
            "\t\t\t\tdouble\tdd000,\tdd001,\tdd010,\tdd011,\tdd100,\tdd101,\tdd110,\tdd111;\n"
            "\t\t\t\tdouble2\tmy_val;\n"
            "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tint\toutputModel\t=\tget_global_id(0);\n"
            "//\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\n"
            "//\t\t\t\tif\t(gi\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
            "//\t\t\t\t}\n"
            "\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\t\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\ti;\n"
            "\t\t\t\t\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tx;\n"
            "\t\t\t\t\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\n"
            "\t\t\t\t\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
            "\t\t\t\t\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\n"
            "\t\t\t\t\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\t\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\t\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\t\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\t\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\t\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tz0\t=\tFLOOR(zp);\n"
            "\t\t\t\t\t\t\t\tfz\t=\tzp\t-\tz0;\n"
            "\t\t\t\t\t\t\t\tz0\t-=\tSTARTINGZ;\n"
            "\t\t\t\t\t\t\t\tz1\t=\tz0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
            "\t\t\t\t\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
            "\t\t\t\t\t\t\t\tmfz\t=\t1.\t-\tfz;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdd000\t=\tmfz\t*\tmfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd001\t=\tmfz\t*\tmfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\tdd010\t=\tmfz\t*\t\tfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd011\t=\tmfz\t*\t\tfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\tdd100\t=\t\tfz\t*\tmfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd101\t=\t\tfz\t*\tmfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\tdd110\t=\t\tfz\t*\t\tfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd111\t=\t\tfz\t*\t\tfy\t*\t\tfx;\n"
            "\n"
            "\t\t\t\t\t\t\t\tmy_weight\t=\tOUTPUT_ELEM(fWeight,FrefNo,i,x);\n"
            "\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tmy_weight\t>\t0.0\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM(f2d,FrefNo,i,x);\n"
            "\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tom:\t%d\tXP:\t%.2f\tX0:\t%d\tX1:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\ti:\t%d\tx:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\tgi,\toutputModel,\txp,\tx0,\tx1,\typ,\ty0,\ty1,\tzp,\tz0,\tz1,\ti,\tx,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty0,\tx0)\t+=\tdd000\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty0,\tx1)\t+=\tdd001\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty1,\tx0)\t+=\tdd010\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz0,\ty1,\tx1)\t+=\tdd011\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty0,\tx0)\t+=\tdd100\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty0,\tx1)\t+=\tdd101\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty1,\tx0)\t+=\tdd110\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelValues,outputModel,\tz1,\ty1,\tx1)\t+=\tdd111\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_weight;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tactual:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tstored:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\\n\",\tgi,\tdd000\t*\tmy_weight,\tdd001\t*\tmy_weight,\tdd010\t*\tmy_weight,\tdd011\t*\tmy_weight,\tdd100\t*\tmy_weight,\tdd101\t*\tmy_weight,\tdd110\t*\tmy_weight,\tdd111\t*\tmy_weight,\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx0),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty0,\tx1),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx0),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz0,\ty1,\tx1),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx0),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty0,\tx1),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx0),\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM(modelWeights,outputModel,\tz1,\ty1,\tx1));\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
            "\t\t\t\t}\n"
            "}\n"
            "\n"
            "//Global\there\tis\tthe\tmodel\tsize\n"
            "kernel\tvoid\tsumModel(global\tdouble4\t*modelValuesIn,\tglobal\tdouble2\t*modelValuesOut,\tglobal\tdouble\t*modelWeightsOut,\tconst\tint\ttotalModels,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tvalueSum\t=\tmodelValuesOut[gi];\n"
            "\t\t\t\tdouble\tweightSum\t=\tmodelWeightsOut[gi];\n"
            "\t\t\t\t\n"
            "\t\t\t\tfor\t(int\ti\t=\t0;\ti\t<\ttotalModels;\ti++)\t{\n"
            "\t\t\t\t\t\t\t\tdouble4\tmodelValueWeight\t=\tmodelValuesIn[gi];\n"
            "\t\t\t\t\t\t\t\tvalueSum\t+=\tmodelValueWeight.xy;\n"
            "\t\t\t\t\t\t\t\tweightSum\t+=\tmodelValueWeight.z;\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tlimit;\n"
            "\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tmodelValuesOut[get_global_id(0)]\t=\tvalueSum;\n"
            "\t\t\t\tmodelWeightsOut[get_global_id(0)]\t=\tweightSum;\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalcBackProjection2D(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble4*\tmodelValues,\tconst\tuint\tprojSize,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
            "\t\t\t\tbool\tis_neg_x;\n"
            "\t\t\t\tdouble\tdd00,\tdd01,\tdd10,\tdd11;\n"
            "\t\t\t\tdouble2\tmy_val;\n"
            "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tint\toutputModel\t=\tget_global_id(0);\n"
            "\t\t\t\t//\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//\t\t\t\tif\t(gi\t==\t0)\t{\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
            "\t\t\t\t//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
            "\t\t\t\t//\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tFrefNo\t=\tgi\t/\t(projSize);\n"
            "\t\t\t\t\t\t\t\tint\tFrefPos\t=\tgi\t-\t(FrefNo\t*\tprojSize);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tvalid\t=\tFrefPos\t<\tprojDim.x\t*\tprojDim.y\t?\t1\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\ti;\n"
            "\t\t\t\t\t\t\t\ti\t=\tFrefPos\t/\tprojDim.x;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tint\tx;\n"
            "\t\t\t\t\t\t\t\tx\t=\tFrefPos\t-\t(i\t*\tprojDim.x);\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tx\t<=\tmy_r_max\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tx\t*\tx\t+\ty\t*\ty\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
            "\t\t\t\t\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\t\t\t\t\tfx\t=\txp\t-\tx0;\n"
            "\t\t\t\t\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\t\t\t\t\tfy\t=\typ\t-\ty0;\n"
            "\t\t\t\t\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\t\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
            "\t\t\t\t\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdd00\t=\tmfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd01\t=\tmfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\tdd10\t=\t\tfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\tdd11\t=\t\tfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tmy_weight\t=\tOUTPUT_ELEM(fWeight,FrefNo,i,x);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tvalid\t=\tmy_weight\t>\t0.0\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tif\t(valid)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM(f2d,FrefNo,i,x);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tom:\t%d\tXP:\t%.2f\tX0:\t%d\tX1:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\ti:\t%d\tx:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\tgi,\toutputModel,\txp,\tx0,\tx1,\typ,\ty0,\ty1,\tzp,\tz0,\tz1,\ti,\tx,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM_2D(modelValues,outputModel,\ty0,\tx0)\t+=\tdd00\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM_2D(modelValues,outputModel,\ty0,\tx1)\t+=\tdd01\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM_2D(modelValues,outputModel,\ty1,\tx0)\t+=\tdd10\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tOUTPUT_MODEL_ELEM_2D(modelValues,outputModel,\ty1,\tx1)\t+=\tdd11\t*\t(double4)(my_val.x,\tmy_val.y,\tmy_weight,\t0);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(FrefNo\t==\t0)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"i:\t%3d\tactual:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\tstored:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\\n\",\tgi,\tdd000\t*\tmy_weight,\tdd001\t*\tmy_weight,\tdd010\t*\tmy_weight,\tdd011\t*\tmy_weight,\tdd100\t*\tmy_weight,\tdd101\t*\tmy_weight,\tdd110\t*\tmy_weight,\tdd111\t*\tmy_weight,\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tgi\t+=\tget_global_size(0);\n"
            "\t\t\t\t}\n"
            "}\n"
            "\n"
            "//Option\tB\tseems\tbetter,\tif\tdoing\texcess\tcalculations...\n"
            "//Tried\tB\tfirst\t-\ttoo\tmany\tcollisions\t(or\ttoo\tmany\tloops)\tand\tcrashes\ton\tthe\tGPU\n"
            "kernel\tvoid\tcalcBackProjectionB(const\tuint4\tmodelDim,\tconst\tint\tmodel_r_max,\tconst\tint2\tprojDim,\tconst\tuint\torientations,\tglobal\tdouble*\tA,\tglobal\tdouble2*\tf2d,\tglobal\tdouble*\tfWeight,\tglobal\tdouble2*\tmodelValues,\tglobal\tdouble*\tmodelWeights,\tconst\tint\twgsPerFref,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tfx,\tfy,\tfz,\tmfx,\tmfy,\tmfz,\txp,\typ,\tzp;\n"
            "\t\t\t\tint\tfirst_x,\tx0,\tx1,\ty0,\ty1,\tz0,\tz1,\ty,\ty2,\tr2;\n"
            "\t\t\t\tbool\tis_neg_x;\n"
            "\t\t\t\tdouble\tdd000,\tdd001,\tdd010,\tdd011,\tdd100,\tdd101,\tdd110,\tdd111;\n"
            "\t\t\t\tdouble2\tmy_val;\n"
            "\t\t\t\tdouble\tmy_weight\t=\t1.;\n"
            "\t\t\t\tint\tSTARTINGY,\tSTARTINGZ;\n"
            "\n"
            "\t\t\t\tdouble2\tmodelValue\t=\t(0,0);\n"
            "\t\t\t\tdouble\tmodelWeight\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\tsize_t\tlgs\t=\tget_local_size(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\t//Calculate\tmodel\tcoordinates\n"
            "\t\t\t\tint\txysize\t=\tmodelDim.x\t*\tmodelDim.y;\n"
            "\t\t\t\tint\tzi\t=\tgi\t/\txysize;\n"
            "\t\t\t\tint\tremainder\t=\tgi\t-\tzi\t*\txysize;\n"
            "\t\t\t\tint\tyi\t=\tremainder\t/\tmodelDim.x;\n"
            "\t\t\t\tint\txi\t=\tremainder\t-\tyi\t*\tmodelDim.x;\n"
            "\t\t\t\t\n"
            "//\t\t\t\tif\t(gi\t==\t53100\t+\t840\t+\t1)\n"
            "//\t\t\t\tprintf(\"gi:\t%d\tMDim:\t%d\t%d\t%d\txi:\t%d\tyi:\t%d\tzi:\t%d\txysize:\t%d\\n\",\tgi,\tmodelDim.x,\tmodelDim.y,\tmodelDim.z,\txi,\tyi,\tzi,\txysize);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmy_r_max\t=\tmin(model_r_max,\tprojDim.x\t-\t1);\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tmax_r2\t=\tmy_r_max\t*\tmy_r_max;\n"
            "\t\t\t\t\n"
            "\t\t\t\tSTARTINGY\t=\t-1\t*\t((int)modelDim.y\t-\t1)\t/\t2;\n"
            "\t\t\t\tSTARTINGZ\t=\t-1\t*\t((int)modelDim.z\t-\t1)\t/\t2;\n"
            "\t\t\t\t\n"
            "//\t\t\t\tif\t(gi\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\tprintf(\"r_max:\t%d\tSTARTINGY:\t%d\tSTARTINGZ:\t%d\\nAinv\t=\\n\",\tmy_r_max,\tSTARTINGY,\tSTARTINGZ);\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,0,0),AELEM(A,0,0,1),AELEM(A,0,0,2));\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,1,0),AELEM(A,0,1,1),AELEM(A,0,1,2));\n"
            "//\t\t\t\t\t\t\t\tprintf(\"%.6f\\t%.6f\\t%.6f\\n\",\tAELEM(A,0,2,0),AELEM(A,0,2,1),AELEM(A,0,2,2));\n"
            "//\t\t\t\t}\n"
            "\t\t\t\t\n"
            "\t\t\t\tfor\t(int\tFrefNo\t=\t0;\tFrefNo\t<\torientations;\tFrefNo++)\t{\n"
            "//\t\t\t\tint\tFrefNo\t=\t0;\n"
            "\t\t\t\t\t\t\t\tfor\t(int\ti=0;\ti\t<\tprojDim.y;\ti++)\n"
            "\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tint\tyValid\t=\t1;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t//\tDont\tsearch\tbeyond\tsquare\twith\tside\tmax_r\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tif\t(i\t<=\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfirst_x\t=\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\t\t\t\telse\tif\t(i\t>=\tprojDim.y\t-\tmy_r_max)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty\t=\ti\t-\tprojDim.y;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tx==0\tplane\tis\tstored\ttwice\tin\tthe\tFFTW\tformat.\tDont\tset\tit\ttwice\tin\tBACKPROJECTION!\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfirst_x\t=\t1;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tyValid\t=\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\ty2\t=\ty\t*\ty;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tfor\t(int\tx=first_x;\tx\t<=\tmy_r_max;\tx++)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tint\tvalid\t=\tyValid;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tOnly\tinclude\tpoints\twith\tradius\t<\tmax_r\t(exclude\tpoints\toutside\tcircle\tin\tsquare)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tr2\t=\tx\t*\tx\t+\ty2;\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\tr2\t>\tmax_r2\t?\t0\t:\tvalid;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(r2\t>\tmax_r2)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tcontinue;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tthe\trelevant\tvalue\tin\tthe\tinput\timage\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tOUTPUT_ELEM_WGS(f2d,FrefNo,i,x);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tDIRECT_A2D_ELEM(f2d,\ti,\tx);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tthe\tweight\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tif\t(Mweight\t!=\tNULL)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_weight\t=\tDIRECT_A2D_ELEM(*Mweight,\ti,\tx);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_weight\t=\tOUTPUT_ELEM_WGS(fWeight,FrefNo,i,x);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\telse:\tmy_weight\twas\talready\tinitialised\tto\t1.\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\tmy_weight\t>\t0.\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(my_weight\t>\t0.)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"Got\ta\tpositive\tweight\");\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tlogical\tcoordinates\tin\tthe\t3D\tmap\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txp\t=\tAELEM(A,FrefNo,0,0)\t*\tx\t+\tAELEM(A,FrefNo,0,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t=\tAELEM(A,FrefNo,1,0)\t*\tx\t+\tAELEM(A,FrefNo,1,1)\t*\ty;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\tAELEM(A,FrefNo,2,0)\t*\tx\t+\tAELEM(A,FrefNo,2,1)\t*\ty;\n"
            "\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tOnly\tasymmetric\thalf\tis\tstored\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(xp\t<\t0)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tGet\tcomplex\tconjugated\thermitian\tsymmetry\tpair\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\txp\t=\t-xp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t=\t-yp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t=\t-zp;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\ttrue;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\telse\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tis_neg_x\t=\tfalse;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tTrilinear\tinterpolation\t(with\tphysical\tcoords)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tSubtract\tSTARTINGY\tand\tSTARTINGZ\tto\taccelerate\taccess\tto\tdata\t(STARTINGX=0)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tIn\tthat\tway\tuse\tDIRECT_A3D_ELEM,\trather\tthan\tA3D_ELEM\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tx0\t=\tFLOOR(xp);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfx\t=\t1.0\t-\tfabs(xp\t-\t(double)xi);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tx1\t=\tx0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty0\t=\tFLOOR(yp);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\typ\t-=\tSTARTINGY;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfy\t=\t1.0\t-\tfabs(yp\t-\t(double)yi);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty0\t-=\t\tSTARTINGY;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\ty1\t=\ty0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz0\t=\tFLOOR(zp);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tzp\t-=\tSTARTINGZ;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tfz\t=\t1.0\t-\tfabs(zp\t-\t(double)zi);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz0\t-=\tSTARTINGZ;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tz1\t=\tz0\t+\t1;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((xi\t==\tx0)\t||\t(xi\t==\tx1))\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((yi\t==\ty0)\t||\t(yi\t==\ty1))\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tvalid\t=\t((zi\t==\tz0)\t||\t(zi\t==\tz1))\t?\tvalid\t:\t0;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(valid)\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"XP:\t%.2f\tX0:\t%d\tX1:\t%d\tXi:\t%d\tYP:\t%.2f\tY0:\t%d\tY1:\t%d\tYi:\t%d\tZP:\t%.2f\tZ0:\t%d\tZ1:\t%d\tZi:\t%d\tvalid:\t%d\tf2d:\t%.4f\t%.4f\tfw:\t%.4f\tfx:\t%.4f\tfy:\t%.4f\tfz:\t%.4f\\n\",\txp,\tx0,\tx1,\txi,\typ,\ty0,\ty1,\tyi,\tzp,\tz0,\tz1,\tzi,\tvalid,\tmy_val.x,\tmy_val.y,\tmy_weight,\tfx,\tfy,\tfz);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(valid)\tprintf(\"Valid\");\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tis_neg_x\t?\tmy_val\t*\t(double2)(1,-1)\t:\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmodelValue\t=\tvalid\t?\tmodelValue\t+\tfx\t*\tfy\t*\tfz\t*\tmy_val\t:\tmodelValue;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmodelWeight\t=\tvalid\t?\tmodelWeight\t+\tfx\t*\tfy\t*\tfz\t*\tmy_weight\t:\tmodelWeight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "/*\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfx\t=\t1.\t-\tfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfy\t=\t1.\t-\tfy;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmfz\t=\t1.\t-\tfz;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd000\t=\tmfz\t*\tmfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd001\t=\tmfz\t*\tmfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd010\t=\tmfz\t*\t\tfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd011\t=\tmfz\t*\t\tfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd100\t=\t\tfz\t*\tmfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd101\t=\t\tfz\t*\tmfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd110\t=\t\tfz\t*\t\tfy\t*\tmfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tdd111\t=\t\tfz\t*\t\tfy\t*\t\tfx;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tif\t(is_neg_x)\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmy_val\t=\tconj(my_val);\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tslice\tin\t3D\tweighted\tsum\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(data,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_val;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t//\tStore\tcorresponding\tweights\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty0,\tx0)\t+=\tdd000\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty0,\tx1)\t+=\tdd001\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty1,\tx0)\t+=\tdd010\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz0,\ty1,\tx1)\t+=\tdd011\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty0,\tx0)\t+=\tdd100\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty0,\tx1)\t+=\tdd101\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty1,\tx0)\t+=\tdd110\t*\tmy_weight;\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDIRECT_A3D_ELEM(weight,\tz1,\ty1,\tx1)\t+=\tdd111\t*\tmy_weight;\n"
            "*/\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t}\t//\tendif\tweight>0.\n"
            "\t\t\t\t\t\t\t\t\t\t\t\t}\t//\tendif\tx-loop\n"
            "\t\t\t\t\t\t\t\t}\t//\tendif\ty-loop\n"
            "\t\t\t\t}\t//\tendif\tFrefNo\n"
            "\t\t\t\t\n"
            "\t\t\t\tmodelValues[gi]\t=\tmodelValue;\n"
            "\t\t\t\tmodelWeights[gi]\t=\tmodelWeight;\n"
            "\t\t\t\t\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculateSigma2Offset(global\tdouble\t*weights,\tconst\tdouble2\tpriorOffset,\tglobal\tdouble2\t*translations,\tglobal\tdouble\t*sigma2OffsetReturn,\tconst\tint\ttotalOrients,\tconst\tint\tgi_limit,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tgi_limit)\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tsigma2Offset\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\t0;\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\tdouble2\tmyTranslation\t=\ttranslations[trans];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble2\tdiff;\n"
            "\t\t\t\t\t\t\t\tdiff\t=\tpriorOffset\t-\tmyTranslation;\n"
            "\t\t\t\t\t\t\t\tdouble\tsum;\n"
            "\t\t\t\t\t\t\t\tsum\t=\tdot(diff,\tdiff);\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble\tweight\t=\tweights[gi];\n"
            "\t\t\t\t\t\t\t\tif\t(weight\t>\t0.0)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\tsigma2Offset\t+=\tweight\t*\tsum;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\ttrans++;\n"
            "\t\t\t\t\t\t\t\tgi\t+=\ttotalOrients;\n"
            "//\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
            "//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%f\t%f\tprior:\t%f\t%f\tdiff:\t%f\t%f\tweight\t%f\ts2offset:\t%f\\n\",\tgi,\tmyTranslation.x,\tmyTranslation.y,\tpriorOffset.x,\tpriorOffset.y,\tdiff.x,\tdiff.y,\tweight,\tsigma2Offset);\n"
            "//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\tsigma2OffsetReturn[get_global_id(0)]\t+=\tsigma2Offset;\n"
            "}\n"
            "\n"
            "kernel\tvoid\tcalculate2DPriorOffset(global\tdouble\t*weights,\tconst\tdouble2\tpriorOffset,\tglobal\tdouble2\t*translations,\tglobal\tdouble2\t*priorOffsetReturn,\tconst\tint\ttotalOrients,\tconst\tint\tgi_limit,\tconst\tint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tgi_limit)\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\toffsetSum\t=\t0;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ttrans\t=\t0;\n"
            "\t\t\t\twhile\t(gi\t<\tlimit)\t{\n"
            "\t\t\t\t\t\t\t\tdouble2\tmyTranslation\t=\ttranslations[trans];\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble2\tsum;\n"
            "\t\t\t\t\t\t\t\tsum\t=\tpriorOffset\t+\tmyTranslation;\n"
            "\t\t\t\t\t\t\t\t\n"
            "\t\t\t\t\t\t\t\tdouble\tweight\t=\tweights[gi];\n"
            "\t\t\t\t\t\t\t\tif\t(weight\t>\t0.0)\t{\n"
            "\t\t\t\t\t\t\t\t\t\t\t\toffsetSum\t+=\tweight\t*\tsum;\n"
            "\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t\t\t\t\ttrans++;\n"
            "\t\t\t\t\t\t\t\tgi\t+=\ttotalOrients;\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\tif\t(get_global_id(0)\t==\t0)\t{\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ttrans:\t%f\t%f\tprior:\t%f\t%f\tdiff:\t%f\t%f\tweight\t%f\ts2offset:\t%f\\n\",\tgi,\tmyTranslation.x,\tmyTranslation.y,\tpriorOffset.x,\tpriorOffset.y,\tdiff.x,\tdiff.y,\tweight,\tsigma2Offset);\n"
            "\t\t\t\t\t\t\t\t//\t\t\t\t\t\t\t\t}\n"
            "\t\t\t\t}\n"
            "\t\t\t\tpriorOffsetReturn[get_global_id(0)]\t+=\toffsetSum;\n"
            "}\n"
            "\n"
            "kernel\tvoid\tshiftImageInFourierTransform(const\tint2\tFDim,\tconst\tdouble\toridim,\tglobal\tdouble2\t*shifts,\tglobal\tdouble2\t*Fin,\tglobal\tdouble2\t*Fout,\tconst\tuint\tlimit)\t{\n"
            "\t\t\t\t\n"
            "\t\t\t\tsize_t\tgi\t=\tget_global_id(0);\n"
            "\t\t\t\t\n"
            "\t\t\t\tif\t(gi\t>=\tlimit)\n"
            "\t\t\t\t\t\t\t\treturn;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFSize\t=\tFDim.x\t*\tFDim.y;\n"
            "\t\t\t\tint\ttrans\t=\tgi\t/\tFSize;\n"
            "\t\t\t\tdouble2\tshift\t=\tshifts[trans];\n"
            "\t\t\t\t\n"
            "\t\t\t\tshift\t/=\t-oridim;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\tFIndex\t=\tgi\t-\ttrans\t*\tFSize;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble2\tin\t=\tFin[FIndex];\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ti\t=\tFIndex\t/\tFDim.x;\n"
            "\t\t\t\tint\tx\t=\tFIndex\t-\ti\t*\tFDim.x;\n"
            "\t\t\t\t\n"
            "\t\t\t\tint\ty\t=\t(i\t<\tFDim.x)\t?\ti\t:\ti\t-\tFDim.y;\n"
            "\t\t\t\t\n"
            "\t\t\t\tdouble\tdotp\t=\t2\t*\tM_PI_F\t*\t(x\t*\tshift.x\t+\ty\t*\tshift.y);\n"
            "\t\t\t\tdouble\ta\t=\tcos(dotp);\n"
            "\t\t\t\tdouble\tb\t=\tsin(dotp);\n"
            "\t\t\t\tdouble\tc\t=\tin.x;\n"
            "\t\t\t\tdouble\td\t=\tin.y;\n"
            "\t\t\t\tdouble\tac\t=\ta\t*\tc;\n"
            "\t\t\t\tdouble\tbd\t=\tb\t*\td;\n"
            "\t\t\t\tdouble\tab_cd\t=\t(a\t+\tb)\t*\t(c\t+\td);\n"
            "\t\t\t\tFout[gi]\t=\t(double2)(ac\t-\tbd,\tab_cd\t-\tac\t-\tbd);\n"
            "//\t\t\t\tif\t(trans\t==\t1)\n"
            "//\t\t\t\t\t\t\t\tprintf(\"gi:\t%d\ti:\t%d\ty:\t%d\tx:\t%d\ttrans:\t%d\tshift:\t%f\t%f\tin:\t%f\t%f\tdotp:\t%f\ta\tb\tc\td:\t%f\t%f\t%f\t%f\tout:\t%f\t%f\\n\",\tgi,\ti,\ty,\tx,\ttrans,\tshift.x,\tshift.y,\tin.x,\tin.y,\tdotp,\ta,\tb,\tc,\td,\tac\t-\tbd,\tab_cd\t-\tac\t-\tbd);\n"
            "}\n"
        };
        source = clCode;
//        source = loadProgramSource("/em/Applications/relion/relion-1.3.mod/src/ml_optimiser_exp_par.cl");
//                            source = loadProgramSource("/home/calsmi/relion/relion-1.3.mod/src/ml_optimiser_exp_par.cl");
        
        CL_program = clCreateProgramWithSource(CL_context, 1, (const char **) & source, NULL, &err);
        if (!CL_program || err != CL_SUCCESS)
        {
            std::cerr << "Error: Failed to create program" << std::endl;
            do_use_opencl = false;
            return;
        }
        
        err = clBuildProgram(CL_program, 0, NULL, NULL, NULL, NULL);
        if (err != CL_SUCCESS)
        {
            size_t len;
            char buffer[32768];
            
            clGetProgramBuildInfo(CL_program, CL_device, CL_PROGRAM_BUILD_LOG, 32768, buffer, &len);
            std::cerr << "Error: Failed to build OpenCL executable - build log: " << buffer << std::endl;
            do_use_opencl = false;
            return;
        }
    }
    
    if (CL_atomicSupport) {
        std::cout << "Slave " << rank << " using OpenCL device " << name << " mem: " << CL_global_memsize << " wgs: " << CL_wgs[0] << " 2d: " << CL_wgs[1] << " 3d: " << CL_wgs[2] << " Max compute units: " << CL_maxComputeUnits << " Max mem alloc: " << CL_maxMemAlloc << " Local size/type: " << CL_localSize << "/" << memtype << " with Atomic Back projection" << std::endl;
    } else {
        std::cout << "Slave " << rank << " using OpenCL device " << name << " mem: " << CL_global_memsize << " wgs: " << CL_wgs[0] << " 2d: " << CL_wgs[1] << " 3d: " << CL_wgs[2] << " Max compute units: " << CL_maxComputeUnits << " Max mem alloc: " << CL_maxMemAlloc << " Local size/type: " << CL_localSize << "/" << memtype << std::endl;
    }
    
    CL_blank = clCreateKernel(CL_program, "blank", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - blank" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_ctfAndScaleDataPoint = clCreateKernel(CL_program, "ctfAndScaleDataPoint", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - ctfAndScaleDataPoint" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_ctfAndScaleDataPoint_TR = clCreateKernel(CL_program, "ctfAndScaleDataPoint_TR", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - ctfAndScaleDataPoint_TR" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateDiff2AndSuma2CC_TR = clCreateKernel(CL_program, "calculateDiff2AndSuma2CC_TR", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateDiff2AndSuma2CC_TR" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateDiff2AndSuma2_TR = clCreateKernel(CL_program, "calculateDiff2AndSuma2_TR", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateDiff2AndSuma2_TR" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateDiff2AndSuma2CC = clCreateKernel(CL_program, "calculateDiff2AndSuma2CC", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateDiff2AndSuma2CC" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateDiff2AndSuma2 = clCreateKernel(CL_program, "calculateDiff2AndSuma2", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateDiff2AndSuma2" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateProjections_TR = clCreateKernel(CL_program, "calcModelProjection_TR", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calc Model projections_TR" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateProjections = clCreateKernel(CL_program, "calcModelProjection", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calc Model projections" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calcModelProjection2D_TR = clCreateKernel(CL_program, "calcModelProjection2D_TR", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calc Model projections2D_TR" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calcModelProjection2D = clCreateKernel(CL_program, "calcModelProjection2D", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calc Model projections 2D" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_generateOrientationMatrix = clCreateKernel(CL_program, "generateOrientationMatrix", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - generate Orientation Matrix" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_initWithConstant = clCreateKernel(CL_program, "initWithConstant", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - initWithConstant" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_initWithConstantDouble = clCreateKernel(CL_program, "initWithConstantDouble", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - initWithConstantDouble" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateSWSNoiseEstimate = clCreateKernel(CL_program, "calculateSWSNoiseEstimate", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateSWSNoiseEstimate" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateSWSNormCorrection = clCreateKernel(CL_program, "calculateSWSNormCorrection", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateSWSNormCorrection" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateSWSSum = clCreateKernel(CL_program, "calculateSWSSum", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateSWSSum" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_scaleCTF = clCreateKernel(CL_program, "scaleCTF", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - scaleCTF" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_scaleMatrix = clCreateKernel(CL_program, "scaleMatrix", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - scaleMatrix" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calcBackProjection = clCreateKernel(CL_program, "calcBackProjection", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calcBackProjection" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calcBackProjection2D = clCreateKernel(CL_program, "calcBackProjection2D", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calcBackProjection2D" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_sumModel = clCreateKernel(CL_program, "sumModel", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - sumModel" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateSWSNoiseEstNormCorrection = clCreateKernel(CL_program, "calculateSWSNoiseEstNormCorrection", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateSWSNoiseEstNormCorrection" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculateSigma2Offset = clCreateKernel(CL_program, "calculateSigma2Offset", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculateSigma2Offset" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_calculate2DPriorOffset = clCreateKernel(CL_program, "calculate2DPriorOffset", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - calculate2DPriorOffset" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_shiftImageInFourierTransform = clCreateKernel(CL_program, "shiftImageInFourierTransform", &err);
    if (err != CL_SUCCESS)
    {
        std::cerr << "Error: Failed to create OpenCL kernel - shiftImageInFourierTransform" << std::endl;
        do_use_opencl = false;
        return;
    }
    
    CL_GSD_memoryCalculatedCoarseForIter = 0;
    CL_GSD_memoryCalculatedFineForIter = 0;
}

void MlOptimiser::initialiseWorkLoad()
{

	// Note, this function is overloaded in ml_optimiser_mpi...

	// Randomise the order of the particles
	if (random_seed == -1) random_seed = time(NULL);
    // This is for the division into random classes
	mydata.randomiseOriginalParticlesOrder(random_seed);
    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

    divide_equally(mydata.numberOfOriginalParticles(), 1, 0, my_first_ori_particle_id, my_last_ori_particle_id);

}

void MlOptimiser::calculateSumOfPowerSpectraAndAverageImage(MultidimArray<double> &Mavg, bool myverb)
{

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Entering"<<std::endl;
#endif

    int barstep, my_nr_ori_particles = my_last_ori_particle_id - my_first_ori_particle_id + 1;
	if (myverb > 0)
	{
		std::cout << " Estimating initial noise spectra " << std::endl;
		init_progress_bar(my_nr_ori_particles);
		barstep = XMIPP_MAX(1, my_nr_ori_particles / 60);
	}

	// Note the loop over the particles (part_id) is MPI-parallelized
	int nr_ori_particles_done = 0;
	Image<double> img;
	FileName fn_img;
	MultidimArray<double> ind_spectrum, sum_spectrum, count;
	// For spectrum calculation: recycle the transformer (so do not call getSpectrum all the time)
	MultidimArray<Complex > Faux;
    Matrix1D<double> f(3);
    FourierTransformer transformer;
	MetaDataTable MDimg;

	for (long int ori_part_id = my_first_ori_particle_id; ori_part_id <= my_last_ori_particle_id; ori_part_id++, nr_ori_particles_done++)
	{

		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++)
			{

				long int group_id = mydata.getGroupId(part_id, iseries);
				// TMP test for debuging
				if (group_id < 0 || group_id >= mymodel.nr_groups)
				{
					std::cerr << " group_id= " << group_id << std::endl;
					REPORT_ERROR("MlOptimiser::calculateSumOfPowerSpectraAndAverageImage: bad group_id");
				}

				// Extract the relevant MetaDataTable row from MDimg
				MDimg = mydata.getMetaDataImage(part_id, iseries);

				// Get the image filename
				MDimg.getValue(EMDL_IMAGE_NAME, fn_img);

				// Read image from disc
				img.read(fn_img);
				img().setXmippOrigin();

				// Check that the average in the noise area is approximately zero and the stddev is one
				if (!dont_raise_norm_error)
				{
					int bg_radius2 = ROUND(particle_diameter / (2. * mymodel.pixel_size));
					bg_radius2 *= bg_radius2;
					double sum = 0.;
					double sum2 = 0.;
					double nn = 0.;
					FOR_ALL_ELEMENTS_IN_ARRAY3D(img())
					{
						if (k*k+i*i+j*j > bg_radius2)
						{
							sum += A3D_ELEM(img(), k, i, j);
							sum2 += A3D_ELEM(img(), k, i, j) * A3D_ELEM(img(), k, i, j);
							nn += 1.;
						}
					}
					// stddev
					sum2 -= sum*sum/nn;
					sum2 = sqrt(sum2/nn);
					//average
					sum /= nn;

					// Average should be close to zero, i.e. max +/-50% of stddev...
					// Stddev should be close to one, i.e. larger than 0.5 and smaller than 2)
					if (ABS(sum/sum2) > 0.5 || sum2 < 0.5 || sum2 > 2.0)
					{
						std::cerr << " fn_img= " << fn_img << " bg_avg= " << sum << " bg_stddev= " << sum2 << std::endl;
						REPORT_ERROR("ERROR: It appears that these images have not been normalised to an average background value of 0 and a stddev value of 1. \n \
								Note that the average and stddev values for the background are calculated outside a circle with the particle diameter \n \
								You can use the relion_preprocess program to normalise your images \n \
								If you are sure you have normalised the images correctly (also see the RELION Wiki), you can switch off this error message using the --dont_check_norm command line option");
					}
				}

				// Apply a similar softMask as below (assume zero translations)
				if (do_zero_mask)
					softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), width_mask_edge);

				// Calculate this image's power spectrum in: ind_spectrum
				ind_spectrum.initZeros(XSIZE(img()));
			    count.initZeros(XSIZE(img()));
			    // recycle the same transformer for all images
			    transformer.FourierTransform(img(), Faux, false);
			    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			    {
			    	long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
			    	ind_spectrum(idx) += norm(dAkij(Faux, k, i, j));
			        count(idx) += 1.;
			    }
			    ind_spectrum /= count;

				// Resize the power_class spectrum to the correct size and keep sum
				ind_spectrum.resize(wsum_model.sigma2_noise[0]); // Store sum of all groups in group 0
				wsum_model.sigma2_noise[0] += ind_spectrum;
				wsum_model.sumw_group[0] += 1.;
				mymodel.nr_particles_group[group_id] += 1;


				// Also calculate average image
				if (part_id == mydata.ori_particles[my_first_ori_particle_id].particles_id[0])
					Mavg = img();
				else
					Mavg += img();

			} // end loop iseries
		} // end loop part_id (i)

		if (myverb > 0 && nr_ori_particles_done % barstep == 0)
			progress_bar(nr_ori_particles_done);

	} // end loop ori_part_id


	// Clean up the fftw object completely
	// This is something that needs to be done manually, as among multiple threads only one of them may actually do this
	transformer.cleanup();

	if (myverb > 0)
		progress_bar(my_nr_ori_particles);

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Leaving"<<std::endl;
#endif

}

void MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage(MultidimArray<double> &Mavg)
{

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Entering"<<std::endl;
#endif

	// First calculate average image
	Mavg /= wsum_model.sumw_group[0];

	// for 2D refinements set 2D average to all references
	if (do_average_unaligned)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			mymodel.Iref[iclass] = Mavg;
	}

	// Calculate sigma2_noise estimates as average of power class spectra, and subtract power spectrum of the average image from that
	if (do_calculate_initial_sigma_noise)
	{
		// Factor 2 because of 2-dimensionality of the complex plane
		mymodel.sigma2_noise[0] = wsum_model.sigma2_noise[0] / ( 2. * wsum_model.sumw_group[0] );

		// Calculate power spectrum of the average image
		MultidimArray<double> spect;
		getSpectrum(Mavg, spect, POWER_SPECTRUM);
		spect /= 2.; // because of 2-dimensionality of the complex plane

		// Now subtract power spectrum of the average image from the average power spectrum of the individual images
		spect.resize(mymodel.sigma2_noise[0]);
		mymodel.sigma2_noise[0] -= spect;

		// Set the same spectrum for all groups
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
			mymodel.sigma2_noise[igroup] = mymodel.sigma2_noise[0];
	}


#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Leaving"<<std::endl;
#endif

}

void MlOptimiser::initialLowPassFilterReferences()
{
	if (ini_high > 0.)
	{

		// Make a soft (raised cosine) filter in Fourier space to prevent artefacts in real-space
		// The raised cosine goes through 0.5 at the filter frequency and has a width of width_mask_edge fourier pixels
		double radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
		radius -= WIDTH_FMASK_EDGE / 2.;
		double radius_p = radius + WIDTH_FMASK_EDGE;
		FourierTransformer transformer;
		MultidimArray<Complex > Faux;
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			transformer.FourierTransform(mymodel.Iref[iclass], Faux);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			{
				double r = sqrt((double)(kp*kp + ip*ip + jp*jp));
				if (r < radius)
					continue;
				else if (r > radius_p)
					DIRECT_A3D_ELEM(Faux, k, i, j) = 0.;
				else
				{
					DIRECT_A3D_ELEM(Faux, k, i, j) *= 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
				}
			}
			transformer.inverseFourierTransform(Faux, mymodel.Iref[iclass]);
		}

	}

}

/** ========================== EM-Iteration  ================================= */

void MlOptimiser::iterateSetup()
{

	// Make a barrier where all working threads wait
	global_barrier = new Barrier(nr_threads - 1);

    // Create threads to start working
	global_ThreadManager = new ThreadManager(nr_threads, this);

	// Set up the thread task distributors for the particles and the orientations (will be resized later on)
	exp_ipart_ThreadTaskDistributor = new ThreadTaskDistributor(1, 1);
	exp_iorient_ThreadTaskDistributor = new ThreadTaskDistributor(1, 1);

}
void MlOptimiser::iterateWrapUp()
{

	// delete barrier, threads and task distributors
    delete global_barrier;
	delete global_ThreadManager;
    delete exp_iorient_ThreadTaskDistributor;
    delete exp_ipart_ThreadTaskDistributor;

}

void MlOptimiser::iterate()
{

	if (do_split_random_halves)
		REPORT_ERROR("ERROR: Cannot split data into random halves without using MPI!");


	// launch threads etc
	iterateSetup();

	// Update the current resolution and image sizes, and precalculate resolution pointers
	// The rest of the time this will be done after maximization and before writing output files,
	// so that current resolution is in the output files of the current iteration
	updateCurrentResolution();

	bool has_already_reached_convergence = false;
	for (iter = iter + 1; iter <= nr_iter; iter++)
    {

#ifdef TIMING
		timer.tic(TIMING_EXP);
#endif

		// SA-stuff
		if (do_sim_anneal)
		{
			double tau = -nr_iter / (std::log(temp_fin/temp_ini));
			temperature = temp_ini * exp(-iter/tau);
			std::cout << " temperature= " << temperature << std::endl;
		}

		if (do_auto_refine)
			printConvergenceStats();

        if (opencl_reset_on_next_iter) {
            if (!do_use_opencl) {
                std::cerr << "Reactivating GPU" << std::endl;
            }
            do_use_opencl = true;
        }
        
		expectation();

#ifdef TIMING
		timer.toc(TIMING_EXP);
		timer.tic(TIMING_MAX);
#endif

		if (do_skip_maximization)
		{
			// Only write data.star file and break from the iteration loop
			write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, 0);
			break;
		}

		maximization();

#ifdef TIMING
		timer.toc(TIMING_MAX);
#endif

		// Apply masks to the reference images
		// At the last iteration, do not mask the map for validation purposes
		if (do_solvent && !has_converged)
			solventFlatten();

		// Re-calculate the current resolution, do this before writing to get the correct values in the output files
		updateCurrentResolution();

		// Write output files
		write(DO_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);

		if (do_auto_refine && has_converged)
		{
			if (verb > 0)
			{
				std::cout << " Auto-refine: Refinement has converged, stopping now... " << std::endl;
				std::cout << " Auto-refine: + Final reconstruction from all particles is saved as: " <<  fn_out << "_class001.mrc" << std::endl;
				std::cout << " Auto-refine: + Final model parameters are stored in: " << fn_out << "_model.star" << std::endl;
				std::cout << " Auto-refine: + Final data parameters are stored in: " << fn_out << "_data.star" << std::endl;
				std::cout << " Auto-refine: + Final resolution (without masking) is: " << 1./mymodel.current_resolution << std::endl;
				if (acc_rot < 10.)
					std::cout << " Auto-refine: + But you may want to run relion_postprocess to mask the unfil.mrc maps and calculate a higher resolution FSC" << std::endl;
				else
				{
					std::cout << " Auto-refine: + WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles!" << std::endl;
					std::cout << " Auto-refine: + WARNING: This has been observed to lead to spurious FSC curves, so be VERY wary of inflated resolution estimates..." << std::endl;
					std::cout << " Auto-refine: + WARNING: You most probably do NOT want to publish these results!" << std::endl;
					std::cout << " Auto-refine: + WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
				}
			}
			break;
		}

		// Check whether we have converged by now
		// If we have, set do_join_random_halves and do_use_all_data for the next iteration
		if (do_auto_refine)
			checkConvergence();

#ifdef TIMING
    	if (verb > 0)
    		timer.printTimes(false);
#endif

    }

	// delete threads etc
	iterateWrapUp();
}

void MlOptimiser::expectation()
{

//#define DEBUG_EXP
#ifdef DEBUG_EXP
	std::cerr << "Entering expectation" << std::endl;
#endif

	// Initialise some stuff
	// A. Update current size (may have been changed to ori_size in autoAdjustAngularSampling) and resolution pointers
	updateImageSizeAndResolutionPointers();

	// B. Initialise Fouriertransform, set weights in wsum_model to zero, etc
	expectationSetup();

#ifdef DEBUG_EXP
	std::cerr << "Expectation: done setup" << std::endl;
#endif

	// C. Calculate expected minimum angular errors (only for 3D refinements)
	// And possibly update orientational sampling automatically
	// TODO: also implement estimate angular sampling for 3D refinements
	if (!((iter==1 && do_firstiter_cc) || do_always_cc) && !do_skip_align)
	{
		// Set the exp_metadata (but not the exp_imagedata which is not needed for calculateExpectedAngularErrors)
		int n_trials_acc = (mymodel.ref_dim==3) ? 100 : 10;
		n_trials_acc = XMIPP_MIN(n_trials_acc, mydata.numberOfOriginalParticles());
		getMetaAndImageDataSubset(0, n_trials_acc-1, false);
		calculateExpectedAngularErrors(0, n_trials_acc-1);
	}

	// D. Update the angular sampling (all nodes except master)
	if ( iter > 1 && (do_auto_refine) )
		updateAngularSampling();

	// E. Check whether everything fits into memory, possibly adjust nr_pool and setup thread task managers
	expectationSetupCheckMemory();

#ifdef DEBUG_EXP
	std::cerr << "Expectation: done setupCheckMemory" << std::endl;
#endif
	if (verb > 0)
	{
		std::cout << " Expectation iteration " << iter;
		if (!do_auto_refine)
			std::cout << " of " << nr_iter;
		std::cout << std::endl;
		init_progress_bar(mydata.numberOfOriginalParticles());
	}

	int barstep = XMIPP_MAX(1, mydata.numberOfOriginalParticles() / 60);
	long int prev_barstep = 0, nr_ori_particles_done = 0;

	// Now perform real expectation over all particles
	// Use local parameters here, as also done in the same overloaded function in MlOptimiserMpi
	long int my_first_ori_particle, my_last_ori_particle;
	while (nr_ori_particles_done < mydata.numberOfOriginalParticles())
	{

		my_first_ori_particle = nr_ori_particles_done;
		my_last_ori_particle = XMIPP_MIN(mydata.numberOfOriginalParticles() - 1, my_first_ori_particle + nr_pool - 1);

		// Get the metadata for these particles
		getMetaAndImageDataSubset(my_first_ori_particle, my_last_ori_particle);

		// perform the actual expectation step on several particles
		expectationSomeParticles(my_first_ori_particle, my_last_ori_particle);

		// Set the metadata for these particles
		setMetaDataSubset(my_first_ori_particle, my_last_ori_particle);

		// Also monitor the changes in the optimal orientations and classes
		monitorHiddenVariableChanges(my_first_ori_particle, my_last_ori_particle);

		nr_ori_particles_done += my_last_ori_particle - my_first_ori_particle + 1;

		if (verb > 0 && nr_ori_particles_done - prev_barstep > barstep)
		{
			prev_barstep = nr_ori_particles_done;
			progress_bar(nr_ori_particles_done);
		}
	}

	if (verb > 0)
		progress_bar(mydata.numberOfOriginalParticles());

	// Clean up some memory
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		mymodel.PPref[iclass].data.clear();
#ifdef DEBUG_EXP
	std::cerr << "Expectation: done " << std::endl;
#endif

}


void MlOptimiser::expectationSetup()
{
#ifdef DEBUG
	std::cerr << "Entering expectationSetup" << std::endl;
#endif

	// Re-initialise the random seed, because with a noisy_mask, inside the previous iteration different timings of different MPI nodes may have given rise to different number of calls to ran1
	// Use the iteration number so that each iteration has a different random seed
	init_random_generator(random_seed + iter);

	// Reset the random perturbation for this sampling
	sampling.resetRandomlyPerturbedSampling();

    // Initialise Projectors and fill vector with power_spectra for all classes
	mymodel.setFourierTransformMaps(!fix_tau, nr_threads);

	// Initialise all weighted sums to zero
	wsum_model.initZeros();

}

void MlOptimiser::expectationSetupCheckMemory(bool myverb)
{

	if (mymodel.orientational_prior_mode != NOPRIOR)
	{
		// First select one random direction and psi-angle for selectOrientationsWithNonZeroPriorProbability
		// This is to get an idea how many non-zero probabilities there will be
		double ran_rot, ran_tilt, ran_psi;
		int randir = (int)(rnd_unif() * sampling.NrDirections(0, true) );
		int ranpsi = (int)(rnd_unif() * sampling.NrPsiSamplings(0, true) );
		if (randir == sampling.NrDirections(0, true))
		{
			//TMP
			REPORT_ERROR("RANDIR WAS TOO BIG!!!!");
			randir--;
		}
		if (ranpsi == sampling.NrPsiSamplings(0, true))
		{
			//TMP
			REPORT_ERROR("RANPSI WAS TOO BIG!!!!");
			ranpsi--;
		}
		sampling.getDirection(randir, ran_rot, ran_tilt);
		sampling.getPsiAngle(ranpsi, ran_psi);
		// Calculate local searches for these angles
		sampling.selectOrientationsWithNonZeroPriorProbability(ran_rot, ran_tilt, ran_psi,
								sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi));
	}

	// Check whether things will fit into memory
	// Each double takes 8 bytes, and their are mymodel.nr_classes references, express in Gb
	double Gb = sizeof(double) / (1024. * 1024. * 1024.);
	// A. Calculate approximate size of the reference maps
	// Forward projector has complex data, backprojector has complex data and real weight
	double mem_references = Gb * mymodel.nr_classes * (2 * MULTIDIM_SIZE((mymodel.PPref[0]).data) + 3 * MULTIDIM_SIZE((wsum_model.BPref[0]).data));
	// B. Calculate size of the exp_Mweight matrices with (YSIZE=nr_pool, XSIZE=mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling)
	nr_pool = max_nr_pool;

	if (mydata.maxNumberOfImagesPerOriginalParticle() > 1)
	{
		// Make sure that all particles in the data set have the same number of images
		// with the same transformation matrices so that their exp_R_mic can be re-used for pooled particles
		// If there are some particles with different transformations, then just set nr_pool to one
		// TODO: optimize this for randomised particle order.....
		// Currently that will lead to pretty bad efficiency IF there are multiple different tilt angles....
		// Or perhaps just forget about pooling. If we're re-refining the orientations that will be screwed anyway...

		// First find a particle with the maxNumberOfImagesPerParticle
		long int ref_part;
		long int maxn = mydata.maxNumberOfImagesPerOriginalParticle();
		for (ref_part = 0; ref_part < mydata.numberOfParticles(); ref_part++)
		{
			if (mydata.getNrImagesInSeries(ref_part) == maxn)
				break;
		}

		// Then check the transformation matrices for all the other particles are all the same
		// Note that particles are allowed to have fewer images in their series...
		Matrix2D<double> first_R_mic, test_R_mic;
		bool is_ok = true;
		for (long int ipart = 0; ipart < mydata.numberOfParticles(); ipart++)
		{
			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(ipart); iseries++)
			{
				first_R_mic = mydata.getMicrographTransformationMatrix(ref_part, iseries);
				test_R_mic = mydata.getMicrographTransformationMatrix(ipart, iseries);
				if (!first_R_mic.equal(test_R_mic))
				{
					is_ok = false;
					break;
				}
			}
		}
		if (!is_ok)
		{
			// Don't pool particles to prevent trouble when re-using exp_R_mic...
			nr_pool = 1;
			if (myverb > 0)
				std::cout << " Switching off the pooling of particles because there are some series with distinct transformation matrices present in the data... ";
		}
	}

	double mem_pool = Gb * nr_pool * mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling, false);
	// Estimate the rest of the program at 0.1 Gb?
	double mem_rest = 0.1;
    
    //Estimate CL memory usage
    long int nr_orients = sampling.NrDirections(adaptive_oversampling, false) * sampling.NrPsiSamplings(adaptive_oversampling, false);
    long int nr_trans = sampling.NrTranslationalSamplings(adaptive_oversampling);
    long int maxMemUsage = sizeof(cl_double) * nr_orients * nr_trans;
    if (maxMemUsage > CL_maxMemAlloc) {
        maxMemUsage = CL_maxMemAlloc;
    }
    long int clMemUsage = maxMemUsage * 2 + maxMemUsage / 2;
    clMemUsage += sizeof(cl_double) * nr_orients * 3 * 2;
    clMemUsage += sizeof(cl_double2) * nr_trans * mymodel.current_size * (mymodel.current_size / 2 + 1) * 2;
    double mem_cl = clMemUsage / sizeof(cl_double) * Gb;
    
	double total_mem_Gb_exp = mem_references + mem_pool + mem_rest + mem_cl;
	// Each reconstruction has to store 1 extra complex array (Fconv) and 4 extra double arrays (Fweight, Fnewweight. vol_out and Mconv in convoluteBlobRealSpace),
	// in adddition to the double weight-array and the complex data-array of the BPref
	// That makes a total of 2*2 + 5 = 9 * a double array of size BPref
	double total_mem_Gb_max = Gb * 9 * MULTIDIM_SIZE((wsum_model.BPref[0]).data);

	bool exp_does_not_fit = false;
	if (total_mem_Gb_exp > available_memory * nr_threads_original)
	{
		double mem_for_pool = (available_memory * nr_threads_original) - mem_rest - mem_references;
		int suggested_nr_pool = FLOOR(mem_for_pool / (Gb * mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling, true)));
		if (suggested_nr_pool > 0)
		{
			if (myverb > 0)
			{
				std::cout << "Reducing nr_pool to "<< suggested_nr_pool<<" to still fit into memory" << std::endl;
			}
			nr_pool = suggested_nr_pool;
			mem_pool = Gb * nr_pool * mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling, false);
			total_mem_Gb_exp = mem_references + mem_pool + mem_rest;
		}
		else
		{
			exp_does_not_fit = true;
		}
	}

	if (myverb > 0)
	{
		// Calculate number of sampled hidden variables:
		int nr_ang_steps = CEIL(PI * particle_diameter * mymodel.current_resolution);
		double myresol_angstep = 360. / nr_ang_steps;
		std::cout << " CurrentResolution= " << 1./mymodel.current_resolution << " Angstroms, which requires orientationSampling of at least "<< myresol_angstep
				   <<" degrees for a particle of diameter "<< particle_diameter << " Angstroms"<< std::endl;
		for (int oversampling = 0; oversampling <= adaptive_oversampling; oversampling++)
		{
			std::cout << " Oversampling= " << oversampling << " NrHiddenVariableSamplingPoints= " << mymodel.nr_classes * sampling.NrSamplingPoints(oversampling, true) << std::endl;
			std::cout << " OrientationalSampling= " << sampling.getAngularSampling(oversampling)
				<< " NrOrientations= "<<sampling.NrDirections(oversampling, false)*sampling.NrPsiSamplings(oversampling, false)<<std::endl;
			std::cout << " TranslationalSampling= " << sampling.getTranslationalSampling(oversampling)
				<< " NrTranslations= "<<sampling.NrTranslationalSamplings(oversampling)<< std::endl;
			std::cout << "=============================" << std::endl;
		}
	}

	if (myverb > 0)
	{
		std::cout << " Estimated memory for expectation step  > " << total_mem_Gb_exp << " Gb, available memory = "<<available_memory * nr_threads_original<<" Gb."<<std::endl;
		std::cout << " Estimated memory for maximization step > " << total_mem_Gb_max << " Gb, available memory = "<<available_memory * nr_threads_original<<" Gb."<<std::endl;
        std::cout << " Estimated host memory used by OpenCL > " << mem_cl << "Gb" << std::endl;

		if (total_mem_Gb_max > available_memory * nr_threads_original || exp_does_not_fit)
		{
			if (exp_does_not_fit)
			std::cout << " WARNING!!! Expected to run out of memory during expectation step ...." << std::endl;
			if (total_mem_Gb_max > available_memory * nr_threads_original)
			std::cout << " WARNING!!! Expected to run out of memory during maximization step ...." << std::endl;
			std::cout << " WARNING!!! Did you set --memory_per_thread to reflect the number of Gb per core on your computer?" << std::endl;
			std::cout << " WARNING!!! If so, then check your processes are not swapping and consider running fewer MPI processors per node." << std::endl;
			std::cout << " + Available memory for each thread, as given by --memory_per_thread      : " << available_memory << " Gb" << std::endl;
			std::cout << " + Number of threads used per MPI process, as given by --j                : " << nr_threads_original << std::endl;
			std::cout << " + Available memory per MPI process 										: " << available_memory * nr_threads_original << " Gb" << std::endl;
		}
	}

	// Now that we also have nr_pool, resize the task manager for the particles

	/// When there are multiple particles for each ori_particle, then this ThreadTaskDistributor will again be resized somewhere below
	exp_ipart_ThreadTaskDistributor->resize(nr_pool, 1);

	// Also resize task manager for the orientations in case of NOPRIOR (otherwise resizing is done in doThreadGetFourierTransformsAndCtfs)
	if (do_skip_align || do_skip_rotate)
	{
		exp_iorient_ThreadTaskDistributor->resize(1, 1);
	}
	else if (mymodel.orientational_prior_mode == NOPRIOR)
	{
		long int nr_orients = sampling.NrDirections() * sampling.NrPsiSamplings();
		int threadBlockSize = (nr_orients > 100) ? 10 : 1;
		exp_iorient_ThreadTaskDistributor->resize(nr_orients, threadBlockSize);
	}
#ifdef DEBUG
	std::cerr << "Leaving expectationSetup" << std::endl;
#endif

}

void MlOptimiser::expectationSomeParticles(long int my_first_ori_particle, long int my_last_ori_particle)
{

#ifdef TIMING
	timer.tic(TIMING_ESP);
#endif

//#define DEBUG_EXPSINGLE
#ifdef DEBUG_EXPSINGLE
	std::cerr << "Entering expectationSomeParticles..." << std::endl;
#endif

#ifdef TIMING
    timer.tic(TIMING_ESP);
    timer.tic(TIMING_ESP_READ);
#endif

    // Use global variables for thread visibility
	exp_my_first_ori_particle = my_first_ori_particle;
    exp_my_last_ori_particle = my_last_ori_particle;
    exp_nr_ori_particles = exp_my_last_ori_particle - exp_my_first_ori_particle + 1;

    // Find out how many particles there are in these ori_particles
    exp_nr_particles = 0;
    for (long int i = my_first_ori_particle; i <= my_last_ori_particle; i++)
    	exp_nr_particles += mydata.ori_particles[i].particles_id.size();

    // If there are more than one particle in each ori_particle, then do these in parallel with threads
    if (nr_pool == 1 && exp_nr_particles/exp_nr_ori_particles > 1)
    {
    	int my_pool = exp_nr_particles/exp_nr_ori_particles;
    	exp_ipart_ThreadTaskDistributor->resize(my_pool, 1);
    }

    // TODO: MAKE SURE THAT ALL PARTICLES IN SomeParticles ARE FROM THE SAME AREA, SO THAT THE R_mic CAN BE RE_USED!!!

	// In the first iteration, multiple seeds will be generated
	// A single random class is selected for each pool of images, and one does not marginalise over the orientations
	// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
    // If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

    // Generally: use all references
    iclass_min = 0;
    iclass_max = mymodel.nr_classes - 1;
    // low-pass filter again and generate the seeds
    if (do_generate_seeds)
    {
    	if (do_firstiter_cc && iter == 1)
    	{
    		// In first (CC) iter, use a single reference (and CC)
    		iclass_min = iclass_max = 0;
    	}
    	else if ( (do_firstiter_cc && iter == 2) || (!do_firstiter_cc && iter == 1))
		{
			// In second CC iter, or first iter without CC: generate the seeds
    		// Now select a single random class
    		// exp_part_id is already in randomized order (controlled by -seed)
    		// WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
			iclass_min = iclass_max = divide_equally_which_group(mydata.numberOfOriginalParticles(), mymodel.nr_classes, exp_my_first_ori_particle);
		}
    }

	// TODO: think of a way to have the different images in a single series have DIFFERENT offsets!!!
	// Right now, they are only centered with a fixed relative translation!!!!

// Thid debug is a good one to step through the separate steps of the expectation to see where trouble lies....
//#define DEBUG_ESP_MEM
#ifdef DEBUG_ESP_MEM
	char c;
	std::cerr << "Before getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
	std::cin >> c;
#endif

	// Read all image of this series into memory, apply old origin offsets and store Fimg, Fctf, exp_old_xoff and exp_old_yoff in vectors./

	exp_ipart_ThreadTaskDistributor->reset();
	global_ThreadManager->run(globalGetFourierTransformsAndCtfs);

	if (do_realign_movies )//&& movie_frame_running_avg_side > 0)
	{
		calculateRunningAveragesOfMovieFrames();
	}

#ifdef DEBUG_ESP_MEM
	std::cerr << "After getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
	std::cin >> c;
#endif

	#ifdef TIMING
    timer.toc(TIMING_ESP_READ);
#endif

	// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
	exp_significant_weight.clear();
	exp_significant_weight.resize(exp_nr_particles);
	for (int n = 0; n < exp_nr_particles; n++)
		exp_significant_weight[n] = -1.;

	// Number of rotational and translational sampling points
	exp_nr_trans = sampling.NrTranslationalSamplings();

	exp_nr_dir = sampling.NrDirections();
	exp_nr_psi = sampling.NrPsiSamplings();
	exp_nr_rot = exp_nr_dir * exp_nr_psi;

	// Only perform a second pass when using adaptive oversampling
	int nr_sampling_passes = (adaptive_oversampling > 0) ? 2 : 1;

	// Pass twice through the sampling of the entire space of rot, tilt and psi
	// The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
	// Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
	// Only those sampling points will contribute to the weighted sums in the third loop below
	for (exp_ipass = 0; exp_ipass < nr_sampling_passes; exp_ipass++)
	{

		if (strict_highres_exp > 0.)
			// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
			exp_current_image_size = coarse_size;
		else if (adaptive_oversampling > 0)
			// Use smaller images in the first pass, larger ones in the second pass
			exp_current_image_size = (exp_ipass == 0) ? coarse_size : mymodel.current_size;
		else
			exp_current_image_size = mymodel.current_size;

		// Use coarse sampling in the first pass, oversampled one the second pass
		exp_current_oversampling = (exp_ipass == 0) ? 0 : adaptive_oversampling;
		exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
		exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);


#ifdef DEBUG_ESP_MEM

	std::cerr << "Before getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
	std::cin >> c;
#endif

		// Calculate the squared difference terms inside the Gaussian kernel for all hidden variables
		getAllSquaredDifferences();

#ifdef DEBUG_ESP_MEM
	std::cerr << "After getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
	std::cin >> c;
#endif

		// Now convert the squared difference terms to weights,
		// also calculate exp_sum_weight, and in case of adaptive oversampling also exp_significant_weight
		convertAllSquaredDifferencesToWeights();

#ifdef DEBUG_ESP_MEM
	std::cerr << "After convertAllSquaredDifferencesToWeights, press any key to continue... " << std::endl;
	std::cin >> c;
#endif

	}// end loop over 2 exp_ipass iterations


	// For the reconstruction step use mymodel.current_size!
	exp_current_image_size = mymodel.current_size;

#ifdef DEBUG_ESP_MEM
	std::cerr << "Before storeWeightedSums, press any key to continue... " << std::endl;
	std::cin >> c;
#endif
	storeWeightedSums();

	// Now calculate the optimal translation for each of the individual images in the series
	//if (mydata.maxNumberOfImagesPerOriginalParticle(my_first_ori_particle, my_last_ori_particle) > 1 && !(do_firstiter_cc && iter == 1))
	//	getOptimalOrientationsForIndividualImagesInSeries();

#ifdef DEBUG_ESP_MEM
	std::cerr << "After storeWeightedSums, press any key to continue... " << std::endl;
	std::cin >> c;
#endif
#ifdef DEBUG_EXPSINGLE
		std::cerr << "Leaving expectationSingleParticle..." << std::endl;
#endif

#ifdef TIMING
	timer.toc(TIMING_ESP);
#endif

}

void MlOptimiser::maximization()
{

	if (verb > 0)
	{
		std::cout << " Maximization ..." << std::endl;
		init_progress_bar(mymodel.nr_classes);
	}

	// First reconstruct the images for each class
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{
		if (mymodel.pdf_class[iclass] > 0.)
		{
			(wsum_model.BPref[iclass]).reconstruct(mymodel.Iref[iclass], gridding_nr_iter, do_map,
					mymodel.tau2_fudge_factor, mymodel.tau2_class[iclass], mymodel.sigma2_class[iclass],
					mymodel.data_vs_prior_class[iclass], mymodel.fsc_halves_class[iclass], wsum_model.pdf_class[iclass],
					false, false, nr_threads, minres_map);

		}
		else
		{
			mymodel.Iref[iclass].initZeros();
		}

		if (verb > 0)
			progress_bar(iclass);
	}

	// Then perform the update of all other model parameters
	maximizationOtherParameters();

	// Keep track of changes in hidden variables
	updateOverallChangesInHiddenVariables();

	// This doesn't really work, and I need the original priors for the polishing...
	//if (do_realign_movies)
	//	updatePriorsForMovieFrames();

	if (verb > 0)
		progress_bar(mymodel.nr_classes);

}

void MlOptimiser::maximizationOtherParameters()
{
	// Note that reconstructions are done elsewhere!
#ifdef DEBUG
	std::cerr << "Entering maximizationOtherParameters" << std::endl;
#endif

	// Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
	double sum_weight = 0.;
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		sum_weight += wsum_model.pdf_class[iclass];

	// Update average norm_correction
	if (do_norm_correction)
	{
		mymodel.avg_norm_correction = wsum_model.avg_norm_correction / sum_weight;
	}

	if (do_scale_correction && !(iter==1 && do_firstiter_cc) )
	{
		double avg_scale_correction = 0., nr_part = 0.;
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{

#ifdef DEVEL_BFAC
			// TMP
			if (verb>0)
			{
				for (int i=0; i<XSIZE(wsum_model.wsum_signal_product_spectra[igroup]); i++)
				{
					std::cout <<" igroup= "<<igroup<< " i= "<<i<<" "<<wsum_model.wsum_signal_product_spectra[igroup](i)<<" "<<wsum_model.wsum_reference_power_spectra[igroup](i)<<std::endl;
				}
			}
#endif

			double sumXA = wsum_model.wsum_signal_product_spectra[igroup].sum();
			double sumAA = wsum_model.wsum_reference_power_spectra[igroup].sum();
			if (sumAA > 0.)
				mymodel.scale_correction[igroup] = sumXA / sumAA;
			else
				mymodel.scale_correction[igroup] = 1.;
			avg_scale_correction += (double)(mymodel.nr_particles_group[igroup]) * mymodel.scale_correction[igroup];
			nr_part += (double)(mymodel.nr_particles_group[igroup]);

		}

		// Constrain average scale_correction to one.
		avg_scale_correction /= nr_part;
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{
			mymodel.scale_correction[igroup] /= avg_scale_correction;
//#define DEBUG_UPDATE_SCALE
#ifdef DEBUG_UPDATE_SCALE
			if (verb > 0)
			{
				std::cerr<< "Group "<<igroup+1<<": scale_correction= "<<mymodel.scale_correction[igroup]<<std::endl;
				for (int i = 0; i < XSIZE(wsum_model.wsum_reference_power_spectra[igroup]); i++)
					if (wsum_model.wsum_reference_power_spectra[igroup](i)> 0.)
						std::cerr << " i= " << i << " XA= " << wsum_model.wsum_signal_product_spectra[igroup](i)
											<< " A2= " << wsum_model.wsum_reference_power_spectra[igroup](i)
											<< " XA/A2= " << wsum_model.wsum_signal_product_spectra[igroup](i)/wsum_model.wsum_reference_power_spectra[igroup](i) << std::endl;

			}
#endif
		}

	}

	// Update model.pdf_class vector (for each k)
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{
		mymodel.pdf_class[iclass] = wsum_model.pdf_class[iclass] / sum_weight;

		// for 2D also update priors of translations for each class!
		if (mymodel.ref_dim == 2)
		{
			if (wsum_model.pdf_class[iclass] > 0.)
				mymodel.prior_offset_class[iclass] = wsum_model.prior_offset_class[iclass] / wsum_model.pdf_class[iclass];
			else
				mymodel.prior_offset_class[iclass].initZeros();
		}

		// Use sampling.NrDirections(0, true) to include all directions (also those with zero prior probability for any given image)
		for (int idir = 0; idir < sampling.NrDirections(0, true); idir++)
		{
			mymodel.pdf_direction[iclass](idir) = wsum_model.pdf_direction[iclass](idir) / sum_weight;
		}
	}

	// Update sigma2_offset
	// Factor 2 because of the 2-dimensionality of the xy-plane
	if (!fix_sigma_offset)
		mymodel.sigma2_offset = (wsum_model.sigma2_offset) / (2. * sum_weight);

	// TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!

	// Also refrain from updating sigma_noise after the first iteration with first_iter_cc!
	if (!fix_sigma_noise && !(iter == 1 && do_firstiter_cc))
	{
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{
			// Factor 2 because of the 2-dimensionality of the complex-plane
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mymodel.sigma2_noise[igroup])
			{
				DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) =
						DIRECT_MULTIDIM_ELEM(wsum_model.sigma2_noise[igroup], n ) /
							(2. * wsum_model.sumw_group[igroup] * DIRECT_MULTIDIM_ELEM(Npix_per_shell, n));
			}
		}
	}

	// After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
    	refs_are_ctf_corrected = true;

	// Some statistics to output
	mymodel.LL = 	wsum_model.LL;
	if ((iter==1 && do_firstiter_cc) || do_always_cc)
		mymodel.LL /= sum_weight; // this now stores the average ccf
	mymodel.ave_Pmax = wsum_model.ave_Pmax / sum_weight;

	// After the first, special iteration, apply low-pass filter of -ini_high again
	if (iter == 1 && do_firstiter_cc)
	{
		initialLowPassFilterReferences();
		if (ini_high > 0.)
		{
			// Adjust the tau2_class and data_vs_prior_class, because they were calculated on the unfiltered maps
			// This is merely a matter of having correct output in the model.star file (these values are not used in the calculations)
			double radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
			radius -= WIDTH_FMASK_EDGE / 2.;
			double radius_p = radius + WIDTH_FMASK_EDGE;

			for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			{
				for (int rr = 0; rr < XSIZE(mymodel.tau2_class[iclass]); rr++)
				{
					double r = (double)rr;
					if (r < radius)
						continue;
					else if (r > radius_p)
					{
						DIRECT_A1D_ELEM(mymodel.tau2_class[iclass], rr) = 0.;
						DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], rr) = 0.;
					}
					else
					{
						double raisedcos = 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
						DIRECT_A1D_ELEM(mymodel.tau2_class[iclass], rr) *= raisedcos * raisedcos;
						DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], rr) *= raisedcos * raisedcos;
					}
				}
			}
		}

		if (do_generate_seeds && mymodel.nr_classes > 1)
		{
			// In the first CC-iteration only a single reference was used
			// Now copy this one reference to all K references, for seed generation in the second iteration
			for (int iclass = 1; iclass < mymodel.nr_classes; iclass++)
			{
				mymodel.tau2_class[iclass] =  mymodel.tau2_class[0];
				mymodel.data_vs_prior_class[iclass] = mymodel.data_vs_prior_class[0];
				mymodel.pdf_class[iclass] = mymodel.pdf_class[0] / mymodel.nr_classes;
				mymodel.pdf_direction[iclass] = mymodel.pdf_direction[0];
				mymodel.Iref[iclass] = mymodel.Iref[0];
			}
			mymodel.pdf_class[0] /= mymodel.nr_classes;
		}

	}

#ifdef DEBUG
	std::cerr << "Leaving maximizationOtherParameters" << std::endl;
#endif
}


void MlOptimiser::solventFlatten()
{
#ifdef DEBUG
	std::cerr << "Entering MlOptimiser::solventFlatten" << std::endl;
#endif
	// First read solvent mask from disc, or pre-calculate it
	Image<double> Isolvent, Isolvent2;
    Isolvent().resize(mymodel.Iref[0]);
	Isolvent().setXmippOrigin();
	Isolvent().initZeros();
	if (fn_mask.contains("None"))
	{
		double radius = particle_diameter / (2. * mymodel.pixel_size);
		double radius_p = radius + width_mask_edge;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent())
		{
			double r = sqrt((double)(k*k + i*i + j*j));
			if (r < radius)
				A3D_ELEM(Isolvent(), k, i, j) = 1.;
			else if (r > radius_p)
				A3D_ELEM(Isolvent(), k, i, j) = 0.;
			else
			{
				A3D_ELEM(Isolvent(), k, i, j) = 0.5 - 0.5 * cos(PI * (radius_p - r) / width_mask_edge );
			}
		}
	}
	else
	{
		Isolvent.read(fn_mask);
		Isolvent().setXmippOrigin();

		if (Isolvent().computeMin() < 0. || Isolvent().computeMax() > 1.)
			REPORT_ERROR("MlOptimiser::solventFlatten: ERROR solvent mask should contain values between 0 and 1 only...");
	}

	// Also read a second solvent mask if necessary
	if (!fn_mask2.contains("None"))
	{
		Isolvent2.read(fn_mask2);
		Isolvent2().setXmippOrigin();
		if (!Isolvent2().sameShape(Isolvent()))
			REPORT_ERROR("MlOptimiser::solventFlatten ERROR: second solvent mask is of incorrect size.");
	}

	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{

		// Then apply the expanded solvent mask to the map
		mymodel.Iref[iclass] *= Isolvent();

		// Apply a second solvent mask if necessary
		// This may for example be useful to set the interior of icosahedral viruses to a constant density value that is higher than the solvent
		// Invert the solvent mask, so that an input mask can be given where 1 is the masked area and 0 is protein....
		if (!fn_mask2.contains("None"))
			softMaskOutsideMap(mymodel.Iref[iclass], Isolvent2(), true);

	} // end for iclass
#ifdef DEBUG
	std::cerr << "Leaving MlOptimiser::solventFlatten" << std::endl;
#endif

}

void MlOptimiser::updateCurrentResolution()
{
//#define DEBUG
#ifdef DEBUG
	std::cerr << "Entering MlOptimiser::updateCurrentResolution" << std::endl;
#endif


    int maxres = 0;
	if (do_map )
	{
		// Set current resolution
		if (ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc)))
		{
			maxres = ROUND(mymodel.ori_size * mymodel.pixel_size / ini_high);
		}
		else
		{
			// Calculate at which resolution shell the data_vs_prior drops below 1
			int ires;
			for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			{
				for (ires = 1; ires < mymodel.ori_size/2; ires++)
				{
					if (DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], ires) < 1.)
						break;
				}
				// Subtract one shell to be back on the safe side
				ires--;
				if (ires > maxres)
					maxres = ires;
			}

			// Never allow smaller maxres than minres_map
			maxres = XMIPP_MAX(maxres, minres_map);
		}
	}
	else
	{
		// If we are not doing MAP-estimation, set maxres to Nyquist
		maxres = mymodel.ori_size/2;
	}
    double newres = mymodel.getResolution(maxres);


    // Check whether resolution improved, if not increase nr_iter_wo_resol_gain
    //if (newres <= best_resol_thus_far)
    if (newres <= mymodel.current_resolution+0.0001) // Add 0.0001 to avoid problems due to rounding error
    	nr_iter_wo_resol_gain++;
    else
    	nr_iter_wo_resol_gain = 0;

    // Store best resolution thus far (but no longer do anything with it anymore...)
    if (newres > best_resol_thus_far)
    	best_resol_thus_far = newres;

    mymodel.current_resolution = newres;

}

void MlOptimiser::updateImageSizeAndResolutionPointers()
{

	// Increment the current_size
    // If we are far from convergence (in the initial stages of refinement) take steps of 25% the image size
    // Do this whenever the FSC at the current_size is larger than 0.2, but NOT when this is in combination with very low Pmax values,
    // in the latter case, over-marginalisation may lead to spuriously high FSCs (2 smoothed maps may look very similar at high-res: all zero!)
    //
    int maxres = mymodel.getPixelFromResolution(mymodel.current_resolution);
	if (mymodel.ave_Pmax > 0.1 && has_high_fsc_at_limit)
    {
		maxres += ROUND(0.25 * mymodel.ori_size / 2);
    }
	else
	{
		// If we are near our resolution limit, use incr_size (by default 10 shells)
		maxres += incr_size;
	}

    // Go back from resolution shells (i.e. radius) to image size, which are BTW always even...
	mymodel.current_size = maxres * 2;

	// If realigning movies: go all the way because resolution increase may be substantial
	if (do_use_all_data)
		mymodel.current_size = mymodel.ori_size;

	// current_size can never be larger than ori_size:
	mymodel.current_size = XMIPP_MIN(mymodel.current_size, mymodel.ori_size);
	// The current size is also used in wsum_model (in unpacking)
	wsum_model.current_size = mymodel.current_size;

	// Update coarse_size
	if (strict_highres_exp > 0.)
    {
    	// Strictly limit the coarse size to the one corresponding to strict_highres_exp
    	coarse_size = 2 * ROUND(mymodel.ori_size * mymodel.pixel_size / strict_highres_exp);
    }
    else if (adaptive_oversampling > 0.)
	{
    	// Dependency of coarse_size on the angular sampling used in the first pass
    	double rotated_distance = (sampling.getAngularSampling() / 360.) * PI * particle_diameter;
		double keepsafe_factor = (mymodel.ref_dim == 3) ? 1.2 : 1.5;
		double coarse_resolution = rotated_distance / keepsafe_factor;
		// Note coarse_size should be even-valued!
		coarse_size = 2 * CEIL(mymodel.pixel_size * mymodel.ori_size / coarse_resolution);
		// Coarse size can never be larger than max_coarse_size
		coarse_size = XMIPP_MIN(max_coarse_size, coarse_size);
	}
	else
		coarse_size = mymodel.current_size;
    // Coarse_size can never become bigger than current_size
    coarse_size = XMIPP_MIN(mymodel.current_size, coarse_size);

	/// Also update the resolution pointers here

	// Calculate number of pixels per resolution shell
	Npix_per_shell.initZeros(mymodel.ori_size / 2 + 1);
	MultidimArray<double> aux;
	aux.resize(mymodel.ori_size, mymodel.ori_size / 2 + 1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(aux)
	{
		int ires = ROUND(sqrt((double)(kp*kp + ip*ip + jp*jp)));
		// TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
		// Exclude points beyond XSIZE(Npix_per_shell), and exclude half of the x=0 column that is stored twice in FFTW
		if (ires < mymodel.ori_size / 2 + 1 && !(jp==0 && ip < 0))
			Npix_per_shell(ires) += 1;
	}

	Mresol_fine.resize(mymodel.current_size, mymodel.current_size / 2 + 1);
	Mresol_fine.initConstant(-1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mresol_fine)
	{
		int ires = ROUND(sqrt((double)(kp*kp + ip*ip + jp*jp)));
		// TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
		// Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
		if (ires < mymodel.current_size / 2 + 1  && !(jp==0 && ip < 0))
		{
			DIRECT_A3D_ELEM(Mresol_fine, k, i, j) = ires;
		}
	}

	Mresol_coarse.resize(coarse_size, coarse_size/ 2 + 1);
	Mresol_coarse.initConstant(-1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mresol_coarse)
	{
		int ires = ROUND(sqrt((double)(kp*kp + ip*ip + jp*jp)));
		// Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
		// exclude lowest-resolution points
		if (ires < coarse_size / 2 + 1 && !(jp==0 && ip < 0))
		{
			DIRECT_A3D_ELEM(Mresol_coarse, k, i, j) = ires;
		}
	}

//#define DEBUG_MRESOL
#ifdef DEBUG_MRESOL
	Image<double> img;
	img().resize(YSIZE(Mresol_fine),XSIZE(Mresol_fine));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
	{
		DIRECT_MULTIDIM_ELEM(img(), n) = (double)DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
	}
	img.write("Mresol_fine.mrc");
	img().resize(YSIZE(Mresol_coarse),XSIZE(Mresol_coarse));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
	{
		DIRECT_MULTIDIM_ELEM(img(), n) = (double)DIRECT_MULTIDIM_ELEM(Mresol_coarse, n);
	}
	img.write("Mresol_coarse.mrc");
#endif


#ifdef DEBUG
	std::cerr << " current_size= " << mymodel.current_size << " coarse_size= " << coarse_size << " current_resolution= " << mymodel.current_resolution << std::endl;
	std::cerr << "Leaving MlOptimiser::updateCurrentResolution" << std::endl;
#endif

}


double MlOptimiser::calculatePdfOffset(Matrix1D<double> offset, Matrix1D<double> prior)
{
	if (mymodel.sigma2_offset < 0.0001)
	{
		return (offset.sum2() > 0.) ? 0. : 1.;
	}
	else
	{
		return exp ( (offset-prior).sum2() / (-2. * mymodel.sigma2_offset) ) / ( 2. * PI * mymodel.sigma2_offset );
	}
}

void MlOptimiser::calculateRunningAveragesOfMovieFrames()
{
	std::vector<MultidimArray<Complex > > runavg_Fimgs;
	std::vector<int> count_runavg;
	MultidimArray<Complex > Fzero;
	Fzero.resize(exp_Fimgs[0]);
	Fzero.initZeros();

	// initialise the sums at zero
	for (int iimg = 0; iimg < exp_Fimgs.size(); iimg++)
	{
		runavg_Fimgs.push_back(Fzero);
		count_runavg.push_back(0);
	}

	// running avgs NOT for series!
	int iseries = 0;

//#define DEBUG_RUNAVG
#ifdef DEBUG_RUNAVG
	FourierTransformer transformer;
	MultidimArray< Complex > Fimg;
	Image<double> It;
	if (verb)
	{
		Fimg = exp_Fimgs[0];
		It().resize(YSIZE(Fimg),YSIZE(Fimg));
		transformer.inverseFourierTransform(Fimg, It());
		CenterFFT(It(), false);
		It.write("Fimg.spi");
		std::cerr << "Written Fimg" << std::endl;
	}
#endif

	// Calculate the running sums
	for (int iimg = 0; iimg < exp_Fimgs.size(); iimg++)
	{
		// Who are we?
		int my_ipart = exp_iimg_to_ipart[iimg];
		long int my_ori_part_id = exp_ipart_to_ori_part_id[my_ipart];
		long int my_part_id = exp_ipart_to_part_id[my_ipart];
		int my_frame = exp_ipart_to_ori_part_nframe[my_ipart];

#ifdef DEBUG_RUNAVG
		if (verb)
		{
			long int my_img_id = mydata.getImageId(my_part_id, iseries);
			FileName fntt;
			mydata.MDimg.getValue(EMDL_IMAGE_NAME, fntt, my_img_id);
			std::cerr << " my= " << fntt;
		}
#endif

		long int my_first_runavg_frame = XMIPP_MAX(0, my_frame - movie_frame_running_avg_side);
		long int my_last_runavg_frame = XMIPP_MIN(mydata.ori_particles[my_ori_part_id].particles_id.size() - 1, my_frame + movie_frame_running_avg_side);

		// Run over all images again and see which ones to sum
		for (int iimg2 = 0; iimg2 < exp_Fimgs.size(); iimg2++)
		{
			int other_ipart = exp_iimg_to_ipart[iimg2];
			long int other_ori_part_id = exp_ipart_to_ori_part_id[other_ipart];
			long int other_part_id = exp_ipart_to_part_id[other_ipart];
			int other_frame = exp_ipart_to_ori_part_nframe[other_ipart];

			if (my_ori_part_id == other_ori_part_id && other_frame >= my_first_runavg_frame && other_frame <= my_last_runavg_frame)
			{
#ifdef DEBUG_RUNAVG
				if (verb)
				{
					long int other_img_id = mydata.getImageId(other_part_id, iseries);
					FileName fnt, fnm, fnp;
					mydata.MDimg.getValue(EMDL_IMAGE_NAME, fnt, other_img_id);
					mydata.MDimg.getValue(EMDL_PARTICLE_ORI_NAME, fnp, other_img_id);
					mydata.MDimg.getValue(EMDL_MICROGRAPH_NAME, fnm, other_img_id);
					std::cerr << " = " << fnt<<" "<<fnm<<" "<<fnp;
				}
#endif

				// Add to sum
				runavg_Fimgs[iimg] += exp_Fimgs[iimg2];
				count_runavg[iimg] += 1;
			}
		}

#ifdef DEBUG_RUNAVG
		if (verb)
			std::cerr << std::endl;
#endif
	}

	// Calculate averages from sums and set back in exp_ vectors
	for (int iimg = 0; iimg < exp_Fimgs.size(); iimg++)
	{
		double sum = (double)count_runavg[iimg];
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_Fimgs[iimg])
		{
			DIRECT_MULTIDIM_ELEM(exp_Fimgs[iimg], n) = DIRECT_MULTIDIM_ELEM(runavg_Fimgs[iimg], n) / sum;
		}
		// Also lower the power of the images for the sigma2_noise and diff2 calculations beyond current_size....
		// sigma2_(a+b) = sigma2_(a) + sigma2_(b)
		// The actual values are lost, just hope the images obey statistics...
		exp_power_imgs[iimg] /= sum;
		exp_highres_Xi2_imgs[iimg] /= sum;
	}
#ifdef DEBUG_RUNAVG
	if (verb)
	{
		Fimg = exp_Fimgs[0];
		It().resize(YSIZE(Fimg),YSIZE(Fimg));
		transformer.inverseFourierTransform(Fimg, It());
		CenterFFT(It(), false);
		It.write("Frunavg.spi");
		std::cerr << "Written Frunavg.spi, sleeping 2 seconds..." << std::endl;
		sleep(2);

	}
#endif

}

void MlOptimiser::doThreadGetFourierTransformsAndCtfs(int thread_id)
{
	// Only first thread initialises
	if (thread_id == 0)
	{
		exp_starting_image_no.clear();
		exp_power_imgs.clear();
		exp_highres_Xi2_imgs.clear();
		exp_Fimgs.clear();
		exp_Fimgs_nomask.clear();
		exp_Fctfs.clear();
		exp_old_offset.clear();
		exp_prior.clear();
		exp_local_oldcc.clear();
		exp_ipart_to_part_id.clear();
		exp_ipart_to_ori_part_id.clear();
		exp_ipart_to_ori_part_nframe.clear();
		exp_iimg_to_ipart.clear();

		// Resize to the right size instead of using pushbacks
		exp_starting_image_no.resize(exp_nr_particles);

		// First check how many images there are in the series for each particle...
		// And calculate exp_nr_images
		exp_nr_images = 0;
		for (long int ori_part_id = exp_my_first_ori_particle, my_image_no = 0, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
	    {

#ifdef DEBUG_CHECKSIZES
			if (ori_part_id >= mydata.ori_particles.size())
			{
				std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
				REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
			}
#endif
			for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
			{
				long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
				int iipart = exp_ipart_to_part_id.size();
				exp_starting_image_no.at(iipart) = exp_nr_images;
				exp_nr_images += mydata.getNrImagesInSeries(part_id);
				exp_ipart_to_part_id.push_back(part_id);
				exp_ipart_to_ori_part_id.push_back(ori_part_id);
				exp_ipart_to_ori_part_nframe.push_back(i);
				for (int i = 0; i < mydata.getNrImagesInSeries(part_id); i++)
					exp_iimg_to_ipart.push_back(iipart);
			}
	    }
		// Then also resize vectors for all images
		exp_power_imgs.resize(exp_nr_images);
		exp_highres_Xi2_imgs.resize(exp_nr_images);
		exp_Fimgs.resize(exp_nr_images);
		exp_Fimgs_nomask.resize(exp_nr_images);
		exp_Fctfs.resize(exp_nr_images);
		exp_old_offset.resize(exp_nr_images);
		exp_prior.resize(exp_nr_images);
		exp_local_oldcc.resize(exp_nr_images);

	}
	global_barrier->wait();

	FourierTransformer transformer;
	size_t first_ipart = 0, last_ipart = 0;
	while (exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{

		for (long int ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			// the exp_ipart_ThreadTaskDistributor was set with nr_pool,
			// but some, e.g. the last, batch of pooled particles may be smaller
			if (ipart >= exp_nr_particles)
				break;

#ifdef DEBUG_CHECKSIZES
			if (ipart >= exp_ipart_to_part_id.size())
			{
				std::cerr<< "ipart= "<<ipart<<" exp_ipart_to_part_id.size()= "<< exp_ipart_to_part_id.size() <<std::endl;
				REPORT_ERROR("ipart >= exp_ipart_to_part_id.size()");
			}
#endif
			long int part_id = exp_ipart_to_part_id[ipart];

			// Prevent movies and series at the same time...
			if (mydata.getNrImagesInSeries(part_id) > 1 && do_realign_movies)
				REPORT_ERROR("Not ready yet for dealing with image series at the same time as realigning movie frames....");

			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++)
			{

				FileName fn_img;
				Image<double> img, rec_img;
				MultidimArray<Complex > Fimg, Faux;
				MultidimArray<double> Fctf;
				int my_image_no = exp_starting_image_no.at(ipart) + iseries;
				// Which group do I belong?
				int group_id = mydata.getGroupId(part_id, iseries);

				// Get the norm_correction
				double normcorr = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM);

				// Get the optimal origin offsets from the previous iteration
				Matrix1D<double> my_old_offset(2), my_prior(2);
				XX(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF);
				YY(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF);
				XX(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR);
				YY(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR);
				// Uninitialised priors were set to 999.
				if (XX(my_prior) > 998.99 && XX(my_prior) < 999.01)
					XX(my_prior) = 0.;
				if (YY(my_prior) > 998.99 && YY(my_prior) < 999.01)
					YY(my_prior) = 0.;

				// Get the old cross-correlations
				exp_local_oldcc.at(my_image_no) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL);

				// If we do local angular searches, get the previously assigned angles to center the prior
				// Only do this for the first image in the series, as this prior work per-particle, not per-image
				// All images in the series use the same rotational sampling, brought back to "exp_R_mic=identity"
				if (do_skip_align || do_skip_rotate)
				{
					// No need to block the threads global_mutex, as nr_pool will be set to 1 anyway for do_skip_align!
					if (do_skip_align)
					{
						// Rounded translations will be applied to the image upon reading,
						// set the unique translation in the sampling object to the fractional difference
						Matrix1D<double> rounded_offset = my_old_offset;
						rounded_offset.selfROUND();
						rounded_offset = my_old_offset - rounded_offset;
						sampling.setOneTranslation(rounded_offset);
					}

					// Also set the rotations
					double old_rot, old_tilt, old_psi;
					old_rot = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT);
					old_tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT);
					old_psi = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI);
					sampling.setOneOrientation(old_rot, old_tilt, old_psi);

				}
				else if (mymodel.orientational_prior_mode != NOPRIOR && iseries == 0)
				{
					// First try if there are some fixed prior angles
					double prior_rot = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT_PRIOR);
					double prior_tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT_PRIOR);
					double prior_psi = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR);

					// If there were no defined priors (i.e. their values were 999.), then use the "normal" angles
					if (prior_rot > 998.99 && prior_rot < 999.01)
						prior_rot = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT);
					if (prior_tilt > 998.99 && prior_tilt < 999.01)
						prior_tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT);
					if (prior_psi > 998.99 && prior_psi < 999.01)
						prior_psi = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI);

					// For tilted series: convert the angles back onto the untilted ones...
					// Calculate the angles back from the Euler matrix because for tilt series exp_R_mic may have changed them...
					Matrix2D<double> A, R_mic(3,3);
					R_mic(0,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_0);
					R_mic(0,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_1);
					R_mic(0,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_2);
					R_mic(1,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_0);
					R_mic(1,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_1);
					R_mic(1,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_2);
					R_mic(2,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_0);
					R_mic(2,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_1);
					R_mic(2,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_2);
					if (!R_mic.isIdentity())
					{
						Euler_angles2matrix(prior_rot, prior_tilt, prior_psi, A);
						A = R_mic.inv() * A;
						Euler_matrix2angles(A, prior_rot, prior_tilt, prior_psi);
					}

					global_mutex.lock();

					// Select only those orientations that have non-zero prior probability
					sampling.selectOrientationsWithNonZeroPriorProbability(prior_rot, prior_tilt, prior_psi,
							sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi));

					long int nr_orients = sampling.NrDirections() * sampling.NrPsiSamplings();
					if (nr_orients == 0)
					{
						std::cerr << " sampling.NrDirections()= " << sampling.NrDirections() << " sampling.NrPsiSamplings()= " << sampling.NrPsiSamplings() << std::endl;
						REPORT_ERROR("Zero orientations fall within the local angular search. Increase the sigma-value(s) on the orientations!");
					}
					int threadBlockSize = (nr_orients > 100) ? 10 : 1;

					exp_iorient_ThreadTaskDistributor->resize(nr_orients, threadBlockSize);

					global_mutex.unlock();

				}

				// Unpack the image from the imagedata
				img().resize(mymodel.ori_size, mymodel.ori_size);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
				{
					DIRECT_A2D_ELEM(img(), i, j) = DIRECT_A3D_ELEM(exp_imagedata, my_image_no, i, j);
				}
				img().setXmippOrigin();
				if (has_converged && do_use_reconstruct_images)
				{
					rec_img().resize(mymodel.ori_size, mymodel.ori_size);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
					{
						DIRECT_A2D_ELEM(rec_img(), i, j) = DIRECT_A3D_ELEM(exp_imagedata, exp_nr_images + my_image_no, i, j);
					}
					rec_img().setXmippOrigin();
				}
//#define DEBUG_SOFTMASK
#ifdef DEBUG_SOFTMASK
					Image<double> tt;
					tt()=img();
					tt.write("Fimg_unmasked.spi");
					std::cerr << "written Fimg_unmasked.spi; press any key to continue..." << std::endl;
					char c;
					std::cin >> c;
#endif
				// Apply the norm_correction term
				if (do_norm_correction)
				{
//#define DEBUG_NORM
#ifdef DEBUG_NORM
					if (normcorr < 0.001 || normcorr > 1000. || mymodel.avg_norm_correction < 0.001 || mymodel.avg_norm_correction > 1000.)
					{
						std::cerr << " ** normcorr= " << normcorr << std::endl;
						std::cerr << " ** mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
						std::cerr << " ** fn_img= " << fn_img << " part_id= " << part_id << std::endl;
						std::cerr << " ** iseries= " << iseries << " ipart= " << ipart << " part_id= " << part_id << std::endl;
						int group_id = mydata.getGroupId(part_id, iseries);
						std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << " group_id= " << group_id <<std::endl;
						std::cerr << " part_id= " << part_id << " iseries= " << iseries << std::endl;
						std::cerr << " img_id= " << img_id << std::endl;
						REPORT_ERROR("Very small or very big (avg) normcorr!");
					}
#endif
					img() *= mymodel.avg_norm_correction / normcorr;
				}

				// Apply (rounded) old offsets first
				my_old_offset.selfROUND();
				selfTranslate(img(), my_old_offset, DONT_WRAP);
				if (has_converged && do_use_reconstruct_images)
					selfTranslate(rec_img(), my_old_offset, DONT_WRAP);

				exp_old_offset.at(my_image_no) = my_old_offset;
				// Also store priors on translations
				exp_prior.at(my_image_no) = my_prior;

				// Always store FT of image without mask (to be used for the reconstruction)
				MultidimArray<double> img_aux;
				img_aux = (has_converged && do_use_reconstruct_images) ? rec_img() : img();
				CenterFFT(img_aux, true);
				transformer.FourierTransform(img_aux, Faux);
				windowFourierTransform(Faux, Fimg, mymodel.current_size);

				// Here apply the beamtilt correction if necessary
				// This will only be used for reconstruction, not for alignment
				// But beamtilt only affects very high-resolution components anyway...
				//
				double beamtilt_x = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_X);
				double beamtilt_y = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_Y);
				double Cs = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_CS);
				double V = 1000. * DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_VOLTAGE);
				double lambda = 12.2643247 / sqrt(V * (1. + V * 0.978466e-6));
				if (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.)
					selfApplyBeamTilt(Fimg, beamtilt_x, beamtilt_y, lambda, Cs, mymodel.pixel_size, mymodel.ori_size);

				exp_Fimgs_nomask.at(my_image_no) = Fimg;

				long int ori_part_id = exp_ipart_to_ori_part_id[ipart];

				MultidimArray<double> Mnoise;
				if (!do_zero_mask)
				{
					// Make a noisy background image with the same spectrum as the sigma2_noise

					// Different MPI-distributed subsets may otherwise have different instances of the random noise below,
					// because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
					// Have the seed based on the ipart, so that each particle has a different instant of the noise
					// Do this all inside a mutex for the threads, because they all use the same static variables inside ran1...
					// (So the mutex only goal is to make things exactly reproducible with the same random_seed.)
					global_mutex.lock();

					//init_random_generator(random_seed + ori_part_id);
					if (do_realign_movies)
						init_random_generator(random_seed + part_id);
					else
						init_random_generator(random_seed + ori_part_id);

					// If we're doing running averages, then the sigma2_noise was already adjusted for the running averages.
					// Undo this adjustment here in order to get the right noise in the individual frames
					MultidimArray<double> power_noise = sigma2_fudge * mymodel.sigma2_noise[group_id];
					if (do_realign_movies)
						power_noise *= (2. * movie_frame_running_avg_side + 1.);

					// Create noisy image for outside the mask
					MultidimArray<Complex > Fnoise;
					Mnoise.resize(img());
					transformer.setReal(Mnoise);
					transformer.getFourierAlias(Fnoise);
					// Fill Fnoise with random numbers, use power spectrum of the noise for its variance
					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fnoise)
					{
						int ires = ROUND( sqrt( (double)(kp * kp + ip * ip + jp * jp) ) );
						if (ires >= 0 && ires < XSIZE(Fnoise))
						{
							double sigma = sqrt(DIRECT_A1D_ELEM(power_noise, ires));
							DIRECT_A3D_ELEM(Fnoise, k, i, j).real = rnd_gaus(0., sigma);
							DIRECT_A3D_ELEM(Fnoise, k, i, j).imag = rnd_gaus(0., sigma);
						}
						else
						{
							DIRECT_A3D_ELEM(Fnoise, k, i, j) = 0.;
						}
					}
					// Back to real space Mnoise
					transformer.inverseFourierTransform();
					Mnoise.setXmippOrigin();

					// unlock the mutex now that all calss to random functions have finished
					global_mutex.unlock();

					softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), (double)width_mask_edge, &Mnoise);

				}
				else
				{
					softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), (double)width_mask_edge);
				}
#ifdef DEBUG_SOFTMASK
					tt()=img();
					tt.write("Fimg_masked.spi");
					std::cerr << "written Fimg_masked.spi; press any key to continue..." << std::endl;
					exit(1);
					std::cin >> c;
#endif

				// Inside Projector and Backprojector the origin of the Fourier Transform is centered!
				CenterFFT(img(), true);

				// Store the Fourier Transform of the image Fimg
				transformer.FourierTransform(img(), Faux);

				// Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and ori_size
				if (mymodel.current_size < mymodel.ori_size)
				{
					MultidimArray<double> spectrum;
					spectrum.initZeros(mymodel.ori_size/2 + 1);
					double highres_Xi2 = 0.;
					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
					{
						int ires = ROUND( sqrt( (double)(kp*kp + ip*ip + jp*jp) ) );
						// Skip Hermitian pairs in the x==0 column

						if (ires > 0 && ires < mymodel.ori_size/2 + 1 && !(jp==0 && ip < 0) )
						{
							double normFaux = norm(DIRECT_A3D_ELEM(Faux, k, i, j));
							DIRECT_A1D_ELEM(spectrum, ires) += normFaux;
							// Store sumXi2 from current_size until ori_size
							if (ires >= mymodel.current_size/2 + 1)
								highres_Xi2 += normFaux;
						}
					}

					// Let's use .at() here instead of [] to check whether we go outside the vectors bounds
					exp_power_imgs.at(my_image_no) = spectrum;
					exp_highres_Xi2_imgs.at(my_image_no) = highres_Xi2;
				}
				else
				{
					exp_highres_Xi2_imgs.at(my_image_no) = 0.;
				}

				// We never need any resolutions higher than current_size
				// So resize the Fourier transforms
				windowFourierTransform(Faux, Fimg, mymodel.current_size);

				// Also store its CTF
				Fctf.resize(Fimg);

				// Now calculate the actual CTF
				if (do_ctf_correction)
				{
					CTF ctf;
					ctf.setValues(DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_U),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_V),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_ANGLE),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_VOLTAGE),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_CS),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_Q0),
							DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_BFAC));

					ctf.getFftwImage(Fctf, mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size,
							ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
//#define DEBUG_CTF_FFTW_IMAGE
#ifdef DEBUG_CTF_FFTW_IMAGE
					Image<double> tt;
					tt()=Fctf;
					tt.write("relion_ctf.spi");
					std::cerr << "Written relion_ctf.spi, now exiting..." << std::endl;
					exit(1);
#endif
//#define DEBUG_GETCTF
#ifdef DEBUG_GETCTF
					std::cerr << " intact_ctf_first_peak= " << intact_ctf_first_peak << std::endl;
					ctf.write(std::cerr);
					Image<double> tmp;
					tmp()=Fctf;
					tmp.write("Fctf.spi");
					tmp().resize(mymodel.ori_size, mymodel.ori_size);
					ctf.getCenteredImage(tmp(), mymodel.pixel_size, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
					tmp.write("Fctf_cen.spi");
					std::cerr << "Written Fctf.spi, Fctf_cen.spi. Press any key to continue..." << std::endl;
					char c;
					std::cin >> c;
#endif
				}
				else
				{
					Fctf.initConstant(1.);
				}

				// Store Fimg and Fctf
				exp_Fimgs.at(my_image_no) = Fimg;
				exp_Fctfs.at(my_image_no) = Fctf;

			} // end loop iseries
		}// end loop ipart
	} // end while threadTaskDistributor

	// All threads clear out their transformer object when they are finished
	// This is to prevent a call from the first thread to fftw_cleanup, while there are still active plans in the transformer objects....
	// The multi-threaded code with FFTW objects is really a bit of a pain...
	if (thread_id != 0)
	{
		transformer.clear();
	}

	// Wait until all threads have finished
	global_barrier->wait();

	// Only the first thread cleans up the fftw-junk in the transformer object
	if (thread_id == 0)
	{
		transformer.cleanup();
	}

}

#define CL_WGS 256
//#define CL_MODEL_EU_OVERSAMPLING 2
//#define CL_VERIFY_ON_CPU
//#define CL_FULL_VERIFY
//#define CL_PRINT_SAMPLE_VALUES
//#define CL_PRINT_SPEED
#define CL_NUMBER_OF_SAMPLES_TO_PRINT 1000

void MlOptimiser::doOpenCLPrecalculateShiftedImagesCtfsAndInvSigma2s()
{
#ifdef TIMING
    timer.tic(TIMING_DIFF_SHIFT);
#endif
    
    struct timeval prepStartTV, startTV, endTV;
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif

    int err;
    
//    std::cerr << "Start Calculate shifts" << std::endl;
    
    //Need dimensions to allocate memory before starting the cycles
    int my_image_no = exp_starting_image_no[0] + exp_iseries;
    MultidimArray<Complex > Fimg;
    windowFourierTransform(exp_Fimgs[my_image_no], Fimg, exp_current_image_size);
    cl_int2 cl_FDim; cl_FDim.x = XSIZE(Fimg); cl_FDim.y = YSIZE(Fimg);
    
    cl_int cl_num_trans = exp_nr_trans * exp_nr_oversampled_trans;
    cl_mem cl_shifts_A = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_num_trans, NULL, 0);
    cl_mem cl_shifts_B = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_num_trans, NULL, 0);
    cl_mem cl_shiftsHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * cl_num_trans, NULL, 0);
    
    cl_mem cl_Fin_A = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fin_nomask_A = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fin_B = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fin_nomask_B = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, NULL, 0);
    
    cl_mem cl_Fout_A = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fout_nomask_A = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fout_B = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fout_nomask_B = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fout_Host = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    cl_mem cl_Fout_nomask_Host = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, NULL, 0);
    
    cl_event clComputeEvents[exp_nr_particles];
    
    for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
    {
        
        cl_mem cl_Fin = ipart % 2 ? cl_Fin_A : cl_Fin_B;
        cl_mem cl_Fin_nomask = ipart % 2 ? cl_Fin_nomask_A : cl_Fin_nomask_B;
        
        cl_mem cl_Fout = ipart % 2 ? cl_Fout_A : cl_Fout_B;
        cl_mem cl_Fout_nomask = ipart % 2 ? cl_Fout_nomask_A : cl_Fout_nomask_B;
        
        cl_mem cl_shifts = ipart % 2 ? cl_shifts_A : cl_shifts_B;
        
        long int part_id = exp_ipart_to_part_id[ipart];
        int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
        
        // Downsize Fimg and Fctf (again) to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
        MultidimArray<Complex > Fimg, Fshifted, Fimg_nomask, Fshifted_nomask;
        windowFourierTransform(exp_Fimgs[my_image_no], Fimg, exp_current_image_size);
        windowFourierTransform(exp_Fimgs_nomask[my_image_no], Fimg_nomask, exp_current_image_size);
        
//        std::cerr << "Start to copy Fs to device" << std::endl;
        clEnqueueWriteBuffer(CL_CopyToDeviceQueue, cl_Fin, false, 0, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, Fimg.data, 0, NULL, NULL);
        clEnqueueWriteBuffer(CL_CopyToDeviceQueue, cl_Fin_nomask, false, 0, sizeof(cl_double2) * cl_FDim.x * cl_FDim.y, Fimg_nomask.data, 0, NULL, NULL);
        
        // Also precalculate the sqrt of the sum of all Xi2
        // (Could exp_current_image_size ever be different from mymodel.current_size? Probhably therefore do it here rather than in getFourierTransforms
        if ((iter == 1 && do_firstiter_cc) || do_always_cc)
        {
            double sumxi2 = 0.;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
            {
                sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
            }
            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
            exp_local_sqrtXi2[my_image_no] = sqrt(sumxi2);
        }
        
//        std::cerr << "Start to copy shifts to device" << std::endl;
        cl_double2 *shiftsHost = (cl_double2 *)clEnqueueMapBuffer(CL_CopyToDeviceQueue, cl_shiftsHost, true, CL_MAP_WRITE, 0, sizeof(cl_double2) * cl_num_trans, 0, 0, 0, &err);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to map buffer for transfer of shifts" << std::endl;
        }
        // Store all translated variants of Fimg
        int my_trans_image = 0;
        for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
        {
            // First get the non-oversampled translations as defined by the sampling object
            std::vector<Matrix1D <double> > oversampled_translations;
            sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);
            
            // Then loop over all its oversampled relatives
            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
            {
                memcpy(&shiftsHost[my_trans_image], oversampled_translations[iover_trans].vdata, sizeof(cl_double2));
                my_trans_image++;
            }
        }
        clEnqueueUnmapMemObject(CL_CopyToDeviceQueue, cl_shiftsHost, shiftsHost, 0, NULL, NULL);
        clEnqueueCopyBuffer(CL_CopyToDeviceQueue, cl_shiftsHost, cl_shifts, 0, 0, sizeof(cl_double2) * cl_num_trans, 0, NULL, NULL);

        clFinish(CL_CopyToDeviceQueue);
        
//        std::cerr << "Start calculation" << std::endl;
        cl_uint cl_limit = cl_FDim.x * cl_FDim.y * cl_num_trans;
        cl_double cl_oridim = (double)mymodel.ori_size;
        
        //Now run the CL kernel to calculate shifts...
        err = clSetKernelArg(CL_shiftImageInFourierTransform, 0, sizeof(cl_int2), &cl_FDim);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 1, sizeof(cl_double), &cl_oridim);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 2, sizeof(cl_mem), &cl_shifts);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 3, sizeof(cl_mem), &cl_Fin);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 4, sizeof(cl_mem), &cl_Fout);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 5, sizeof(cl_uint), &cl_limit);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to set OpenCL kernel arguments - shift Image in Fourier Transform" << std::endl;
        }
        
        cl_int wgs = CL_WGS;
        
        size_t global[1], local[1];
        global[0] = wgs * (cl_limit / wgs + 1);
        local[0] = wgs;
        
        err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_shiftImageInFourierTransform, 1, NULL, global, local, 0, NULL, NULL);
        if (err)
        {
            std::cerr << "Error: Failed to execute OpenCL kernel - Precalculate:shift Image in Fourier Transform, err: " << err << std::endl;
        }
        
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 3, sizeof(cl_mem), &cl_Fin_nomask);
        err |= clSetKernelArg(CL_shiftImageInFourierTransform, 4, sizeof(cl_mem), &cl_Fout_nomask);

        err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_shiftImageInFourierTransform, 1, NULL, global, local, 0, NULL, &clComputeEvents[ipart]);
        if (err)
        {
            std::cerr << "Error: Failed to execute OpenCL kernel - Precalculate: shift Image in Fourier Transform, err: " << err << std::endl;
        }

//        std::cerr << "starting to copy back from device" << std::endl;
        if (ipart > 0) {
            
            clWaitForEvents(1, &clComputeEvents[ipart - 1]);
            
            cl_mem cl_Fout = (ipart - 1) % 2 ? cl_Fout_A : cl_Fout_B;
            cl_mem cl_Fout_nomask = (ipart - 1) % 2 ? cl_Fout_nomask_A : cl_Fout_nomask_B;

            // Store all translated variants of Fimg
            clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_Fout, cl_Fout_Host, 0, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, NULL, NULL);
            cl_double2 *Fout = (cl_double2 *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_Fout_Host, true, CL_MAP_READ, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of shifts" << std::endl;
                do_use_opencl = false;
            }

            clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_Fout_nomask, cl_Fout_nomask_Host, 0, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, NULL, NULL);
            cl_double2 *Fout_nomask = (cl_double2 *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_Fout_nomask_Host, true, CL_MAP_READ, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of shifts" << std::endl;
                do_use_opencl = false;
            }
            
            int my_image_no = exp_starting_image_no[ipart - 1] + exp_iseries;
            int my_trans_image = my_image_no * exp_nr_trans * exp_nr_oversampled_trans;
            int current_trans = 0;
            for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
            {
                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                {
                    // Need to resize arrays and then copy data...
                    exp_local_Fimgs_shifted[my_trans_image].resize(Fimg);
                    exp_local_Fimgs_shifted_nomask[my_trans_image].resize(Fimg_nomask);
                    
                    memcpy(exp_local_Fimgs_shifted[my_trans_image].data, &Fout[current_trans * cl_FDim.x * cl_FDim.y], sizeof(cl_double2) * cl_FDim.x * cl_FDim.y);
                    memcpy(exp_local_Fimgs_shifted_nomask[my_trans_image].data, &Fout_nomask[current_trans * cl_FDim.x * cl_FDim.y], sizeof(cl_double2) * cl_FDim.x * cl_FDim.y);
                    
                    my_trans_image++;
                    current_trans++;
                }
            }
            
            clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_Fout_Host, Fout, 0, NULL, NULL);
            clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_Fout_nomask_Host, Fout_nomask, 0, NULL, NULL);
            
        }
        
        if (ipart == exp_nr_particles - 1) {
            
            clWaitForEvents(1, &clComputeEvents[ipart]);
            
            cl_mem cl_Fout = (ipart) % 2 ? cl_Fout_A : cl_Fout_B;
            cl_mem cl_Fout_nomask = (ipart) % 2 ? cl_Fout_nomask_A : cl_Fout_nomask_B;
            
            // Store all translated variants of Fimg
            clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_Fout, cl_Fout_Host, 0, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, NULL, NULL);
            cl_double2 *Fout = (cl_double2 *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_Fout_Host, true, CL_MAP_READ, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of shifts" << std::endl;
                do_use_opencl = false;
            }
            
            clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_Fout_nomask, cl_Fout_nomask_Host, 0, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, NULL, NULL);
            cl_double2 *Fout_nomask = (cl_double2 *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_Fout_nomask_Host, true, CL_MAP_READ, 0, sizeof(cl_double2) * cl_num_trans * cl_FDim.x * cl_FDim.y, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of shifts" << std::endl;
                do_use_opencl = false;
            }
            
            int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
            int my_trans_image = my_image_no * exp_nr_trans * exp_nr_oversampled_trans;
            int current_trans = 0;
            for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
            {
                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                {
                    // Need to resize arrays and then copy data...
                    exp_local_Fimgs_shifted[my_trans_image].resize(Fimg);
                    exp_local_Fimgs_shifted_nomask[my_trans_image].resize(Fimg_nomask);
                    
                    memcpy(exp_local_Fimgs_shifted[my_trans_image].data, &Fout[current_trans * cl_FDim.x * cl_FDim.y], sizeof(cl_double2) * cl_FDim.x * cl_FDim.y);
                    memcpy(exp_local_Fimgs_shifted_nomask[my_trans_image].data, &Fout_nomask[current_trans * cl_FDim.x * cl_FDim.y], sizeof(cl_double2) * cl_FDim.x * cl_FDim.y);
                    
                    my_trans_image++;
                    current_trans++;
                }
            }
            
            clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_Fout_Host, Fout, 0, NULL, NULL);
            clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_Fout_nomask_Host, Fout_nomask, 0, NULL, NULL);
            
        }
        
    }

#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedTime = (double)(endMicros - startMicros) / 1000000;
    
    
    double gpuElapsedTime = elapsedTime;
    
    unsigned long ops;
    ops = cl_num_trans * cl_FDim.x * cl_FDim.y * 15;
    std::cerr << "Shift Time: " << elapsedTime << "s, ops: " << ops << " trans: " << cl_num_trans << " GFLOPS: " << (double)ops / elapsedTime / 1e9 << std::endl;
#endif

#ifdef CL_VERIFY_ON_CPU
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif

    for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
    {
        
        long int part_id = exp_ipart_to_part_id[ipart];
#ifdef DEBUG_CHECKSIZES
        if (ipart >= exp_starting_image_no.size())
        {
            std::cerr<< "ipart= "<<ipart<<" exp_starting_image_no.size()= "<< exp_starting_image_no.size() <<std::endl;
            REPORT_ERROR("ipart >= exp_starting_image_no.size()");
        }
#endif
        int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
        
#ifdef DEBUG_CHECKSIZES
        if (my_image_no >= exp_Fimgs.size())
        {
            std::cerr<< "my_image_no= "<<my_image_no<<" exp_Fimgs.size()= "<< exp_Fimgs.size() <<std::endl;
            std::cerr << " exp_nr_trans= " << exp_nr_trans << " exp_nr_oversampled_trans= " << exp_nr_oversampled_trans << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
            REPORT_ERROR("my_image_no >= exp_Fimgs.size()");
        }
#endif
        // Downsize Fimg and Fctf (again) to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
        MultidimArray<Complex > Fimg, Fshifted, Fimg_nomask, Fshifted_nomask;
        windowFourierTransform(exp_Fimgs[my_image_no], Fimg, exp_current_image_size);
        windowFourierTransform(exp_Fimgs_nomask[my_image_no], Fimg_nomask, exp_current_image_size);
        
        
        // Also precalculate the sqrt of the sum of all Xi2
        // (Could exp_current_image_size ever be different from mymodel.current_size? Probhably therefore do it here rather than in getFourierTransforms
        if ((iter == 1 && do_firstiter_cc) || do_always_cc)
        {
            double sumxi2 = 0.;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
            {
                sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
            }
            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
            exp_local_sqrtXi2[my_image_no] = sqrt(sumxi2);
        }
        
        // Store all translated variants of Fimg
        int my_trans_image = my_image_no * exp_nr_trans * exp_nr_oversampled_trans;
        for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
        {
            // First get the non-oversampled translations as defined by the sampling object
            std::vector<Matrix1D <double> > oversampled_translations;
            sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);
            
#ifdef DEBUG_CHECKSIZES
            if (oversampled_translations.size() != exp_nr_oversampled_trans)
            {
                std::cerr<< "oversampled_translations.size()= "<<oversampled_translations.size()<<" exp_nr_oversampled_trans= "<< exp_nr_oversampled_trans <<std::endl;
                REPORT_ERROR("oversampled_translations.size() != exp_nr_oversampled_trans");
            }
#endif
            
            // Then loop over all its oversampled relatives
            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
            {
                //#define DEBUG_SHIFTS
#ifdef DEBUG_SHIFTS
                Image<double> It;
                std::cerr << " iover_trans= " << iover_trans << " XX(oversampled_translations[iover_trans] )= " << XX(oversampled_translations[iover_trans] ) << " YY(oversampled_translations[iover_trans] )= " << YY(oversampled_translations[iover_trans] ) << std::endl;
#endif
                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
//                if (my_trans_image == 0) {
//                    shiftImageInFourierTransformDebug(Fimg, Fshifted, tab_sin, tab_cos, (double)mymodel.ori_size, oversampled_translations[iover_trans]);
//                } else {
                    shiftImageInFourierTransform(Fimg, Fshifted, tab_sin, tab_cos, (double)mymodel.ori_size, oversampled_translations[iover_trans]);
//                }
                shiftImageInFourierTransform(Fimg_nomask, Fshifted_nomask, tab_sin, tab_cos, (double)mymodel.ori_size, oversampled_translations[iover_trans]);
                
#ifdef DEBUG_SHIFTS
                FourierTransformer transformer;
                It().resize(YSIZE(Fimg), YSIZE(Fimg));
                transformer.inverseFourierTransform(Fimg, It());
                CenterFFT(It(), false);
                It.write("Fimg.spi");
                transformer.inverseFourierTransform(Fshifted, It());
                CenterFFT(It(), false);
                It.write("Fshifted.spi");
                std::cerr << "Written Fimg and Fshifted, press any key to continue..." << std::endl;
                char c;
                std::cin >> c;
#endif
                //Compare to what was calculated by the CL code...
                bool Fshifted_correct = true;
                bool Fshifted_nomask_correct = true;
                
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fshifted)
                {
                    double diff_real = (DIRECT_MULTIDIM_ELEM(Fshifted, n)).real - (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[my_trans_image], n)).real;
                    double diff_imag = (DIRECT_MULTIDIM_ELEM(Fshifted, n)).imag - (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[my_trans_image], n)).imag;
                    
                    if ((fabs(diff_real) > 0.001) || (fabs(diff_imag) > 0.001)) {
                        Fshifted_correct = false;
//                        if (my_trans_image == 148) {
//                            std::cerr << n << ":" << "Fimg: " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag << std::endl;
//                            std::cerr << "CPU: " << (DIRECT_MULTIDIM_ELEM(Fshifted, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fshifted, n)).imag << " CL: " << (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[my_trans_image], n)).real << " " << (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[my_trans_image], n)).imag << std::endl;
//                        }
                    }
                }
                
                if (Fshifted_correct) {
//                    std::cerr << "Correct shifts for trans image " << my_trans_image << " ipart: " << ipart << std::endl;
                } else {
                    std::cerr << "Wrong shifts for trans image " << my_trans_image << " ipart: " << ipart << std::endl;
//                    std::cerr << "Shift: " << XX(oversampled_translations[iover_trans]) << " " << YY(oversampled_translations[iover_trans]) << std::endl;
                }
                
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fshifted_nomask)
                {
                    double diff_real = (DIRECT_MULTIDIM_ELEM(Fshifted_nomask, n)).real - (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[my_trans_image], n)).real;
                    double diff_imag = (DIRECT_MULTIDIM_ELEM(Fshifted_nomask, n)).imag - (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[my_trans_image], n)).imag;
                    
                    if ((fabs(diff_real) > 0.001) || (fabs(diff_imag) > 0.001)) {
                        Fshifted_nomask_correct = false;
                    }
                }
                if (Fshifted_nomask_correct) {
//                    std::cerr << "Correct shifts for nomask trans image " << my_trans_image << " ipart: " << ipart << std::endl;
                } else {
                    std::cerr << "Wrong shifts for nomask trans image " << my_trans_image << " ipart: " << ipart << std::endl;
                }
                
                // Store the shifted image
                exp_local_Fimgs_shifted[my_trans_image] = Fshifted;
                exp_local_Fimgs_shifted_nomask[my_trans_image] = Fshifted_nomask;
                my_trans_image++;
            }
        }
    }
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedCPUTime = (double)(endMicros - startMicros) / 1000000;
    std::cerr << "Shift CPU Time: " << elapsedCPUTime << "s GPU time: " << gpuElapsedTime << " s, equivalent to: " << elapsedCPUTime / gpuElapsedTime << " cores" << std::endl;
#endif
#endif
    
    for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
    {
        long int part_id = exp_ipart_to_part_id[ipart];
        int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
        
        // Also store downsized Fctfs
        // In the second pass of the adaptive approach this will have no effect,
        // since then exp_current_image_size will be the same as the size of exp_Fctfs
        
        MultidimArray<double> Fctf;
        windowFourierTransform(exp_Fctfs[my_image_no], Fctf, exp_current_image_size);
        exp_local_Fctfs[my_image_no] = Fctf;
        
        // Get micrograph id (for choosing the right sigma2_noise)
        int group_id = mydata.getGroupId(part_id, exp_iseries);
        
        MultidimArray<double> Minvsigma2;
        Minvsigma2.initZeros(YSIZE(Fimg), XSIZE(Fimg));
        MultidimArray<int> * myMresol = (YSIZE(Fimg) == coarse_size) ? &Mresol_coarse : &Mresol_fine;
        
        // With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*myMresol)
        {
            int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
            // Exclude origin (ires==0) from the Probability-calculation
            // This way we are invariant to additive factors
            if (ires > 0)
            {
                DIRECT_MULTIDIM_ELEM(Minvsigma2, n) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
            }
        }
        
        exp_local_Minvsigma2s[my_image_no] = Minvsigma2;
    }

    
    for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
    {
        clReleaseEvent(clComputeEvents[ipart]);
    }
    clReleaseMemObject(cl_shifts_A);
    clReleaseMemObject(cl_shifts_B);
    clReleaseMemObject(cl_shiftsHost);
    clReleaseMemObject(cl_Fin_A);
    clReleaseMemObject(cl_Fin_B);
    clReleaseMemObject(cl_Fin_nomask_A);
    clReleaseMemObject(cl_Fin_nomask_B);
    clReleaseMemObject(cl_Fout_A);
    clReleaseMemObject(cl_Fout_B);
    clReleaseMemObject(cl_Fout_Host);
    clReleaseMemObject(cl_Fout_nomask_A);
    clReleaseMemObject(cl_Fout_nomask_B);
    clReleaseMemObject(cl_Fout_nomask_Host);
    
#ifdef TIMING
    timer.toc(TIMING_DIFF_SHIFT);
#endif
    
}

void MlOptimiser::doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s(int thread_id)
{
#ifdef TIMING
	timer.tic(TIMING_DIFF_SHIFT);
#endif


	size_t first_ipart = 0, last_ipart = 0;
	while (exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long int ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			// the exp_ipart_ThreadTaskDistributor was set with nr_pool,
			// but some, e.g. the last, batch of pooled particles may be smaller
			if (ipart >= exp_nr_particles)
				break;

			long int part_id = exp_ipart_to_part_id[ipart];
#ifdef DEBUG_CHECKSIZES
			if (ipart >= exp_starting_image_no.size())
			{
				std::cerr<< "ipart= "<<ipart<<" exp_starting_image_no.size()= "<< exp_starting_image_no.size() <<std::endl;
				REPORT_ERROR("ipart >= exp_starting_image_no.size()");
			}
#endif
			int my_image_no = exp_starting_image_no[ipart] + exp_iseries;

#ifdef DEBUG_CHECKSIZES
			if (my_image_no >= exp_Fimgs.size())
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" exp_Fimgs.size()= "<< exp_Fimgs.size() <<std::endl;
				std::cerr << " exp_nr_trans= " << exp_nr_trans << " exp_nr_oversampled_trans= " << exp_nr_oversampled_trans << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
				REPORT_ERROR("my_image_no >= exp_Fimgs.size()");
			}
#endif
			// Downsize Fimg and Fctf (again) to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
			MultidimArray<Complex > Fimg, Fshifted, Fimg_nomask, Fshifted_nomask;
			windowFourierTransform(exp_Fimgs[my_image_no], Fimg, exp_current_image_size);
			windowFourierTransform(exp_Fimgs_nomask[my_image_no], Fimg_nomask, exp_current_image_size);


			// Also precalculate the sqrt of the sum of all Xi2
			// (Could exp_current_image_size ever be different from mymodel.current_size? Probhably therefore do it here rather than in getFourierTransforms
			if ((iter == 1 && do_firstiter_cc) || do_always_cc)
			{
				double sumxi2 = 0.;
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
				{
					sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
				}
				// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
				exp_local_sqrtXi2[my_image_no] = sqrt(sumxi2);
			}

			// Store all translated variants of Fimg
			int my_trans_image = my_image_no * exp_nr_trans * exp_nr_oversampled_trans;
			for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
			{
				// First get the non-oversampled translations as defined by the sampling object
				std::vector<Matrix1D <double> > oversampled_translations;
				sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);

#ifdef DEBUG_CHECKSIZES
				if (oversampled_translations.size() != exp_nr_oversampled_trans)
				{
					std::cerr<< "oversampled_translations.size()= "<<oversampled_translations.size()<<" exp_nr_oversampled_trans= "<< exp_nr_oversampled_trans <<std::endl;
					REPORT_ERROR("oversampled_translations.size() != exp_nr_oversampled_trans");
				}
#endif

				// Then loop over all its oversampled relatives
				for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
				{
//#define DEBUG_SHIFTS
#ifdef DEBUG_SHIFTS
					Image<double> It;
					std::cerr << " iover_trans= " << iover_trans << " XX(oversampled_translations[iover_trans] )= " << XX(oversampled_translations[iover_trans] ) << " YY(oversampled_translations[iover_trans] )= " << YY(oversampled_translations[iover_trans] ) << std::endl;
#endif
					// Shift through phase-shifts in the Fourier transform
					// Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
					shiftImageInFourierTransform(Fimg, Fshifted, tab_sin, tab_cos, (double)mymodel.ori_size, oversampled_translations[iover_trans]);
					shiftImageInFourierTransform(Fimg_nomask, Fshifted_nomask, tab_sin, tab_cos, (double)mymodel.ori_size, oversampled_translations[iover_trans]);

#ifdef DEBUG_SHIFTS
					FourierTransformer transformer;
					It().resize(YSIZE(Fimg), YSIZE(Fimg));
					transformer.inverseFourierTransform(Fimg, It());
					CenterFFT(It(), false);
					It.write("Fimg.spi");
					transformer.inverseFourierTransform(Fshifted, It());
					CenterFFT(It(), false);
					It.write("Fshifted.spi");
					std::cerr << "Written Fimg and Fshifted, press any key to continue..." << std::endl;
					char c;
					std::cin >> c;
#endif
					// Store the shifted image
					exp_local_Fimgs_shifted[my_trans_image] = Fshifted;
					exp_local_Fimgs_shifted_nomask[my_trans_image] = Fshifted_nomask;
					my_trans_image++;
				}
			}


			// Also store downsized Fctfs
			// In the second pass of the adaptive approach this will have no effect,
			// since then exp_current_image_size will be the same as the size of exp_Fctfs
#ifdef DEBUG_CHECKSIZES
			if (my_image_no >= exp_Fctfs.size())
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" exp_Fctfs.size()= "<< exp_Fctfs.size() <<std::endl;
				REPORT_ERROR("my_image_no >= exp_Fctfs.size()");
			}
#endif

			MultidimArray<double> Fctf;
			windowFourierTransform(exp_Fctfs[my_image_no], Fctf, exp_current_image_size);
			exp_local_Fctfs[my_image_no] = Fctf;

			// Get micrograph id (for choosing the right sigma2_noise)
			int group_id = mydata.getGroupId(part_id, exp_iseries);

			MultidimArray<double> Minvsigma2;
			Minvsigma2.initZeros(YSIZE(Fimg), XSIZE(Fimg));
			MultidimArray<int> * myMresol = (YSIZE(Fimg) == coarse_size) ? &Mresol_coarse : &Mresol_fine;

#ifdef DEBUG_CHECKSIZES
			if (!Minvsigma2.sameShape(*myMresol))
			{
				std::cerr<< "!Minvsigma2.sameShape(*myMresol)= "<<!Minvsigma2.sameShape(*myMresol) <<std::endl;
				REPORT_ERROR("!Minvsigma2.sameShape(*myMresol)");
			}
#endif
			// With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*myMresol)
			{
				int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
				// Exclude origin (ires==0) from the Probability-calculation
				// This way we are invariant to additive factors
				if (ires > 0)
				{
					DIRECT_MULTIDIM_ELEM(Minvsigma2, n) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
				}
			}

#ifdef DEBUG_CHECKSIZES
			if (my_image_no >= exp_local_Minvsigma2s.size())
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" exp_local_Minvsigma2s.size()= "<< exp_local_Minvsigma2s.size() <<std::endl;
				REPORT_ERROR("my_image_no >= exp_local_Minvsigma2s.size()");
			}
#endif
			exp_local_Minvsigma2s[my_image_no] = Minvsigma2;
		}
	}

	// Wait until all threads are finsished
	global_barrier->wait();

#ifdef TIMING
	timer.toc(TIMING_DIFF_SHIFT);
#endif


}

bool MlOptimiser::isSignificantAnyParticleAnyTranslation(long int iorient)
{

	for (long int ipart = 0; ipart < YSIZE(exp_Mcoarse_significant); ipart++)
	{
		long int ihidden = iorient * exp_nr_trans;
		for (long int itrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
		{
#ifdef DEBUG_CHECKSIZES
			if (ihidden >= XSIZE(exp_Mcoarse_significant))
			{
				std::cerr << " ihidden= " << ihidden << " XSIZE(exp_Mcoarse_significant)= " << XSIZE(exp_Mcoarse_significant) << std::endl;
				std::cerr << " iorient= " << iorient << " itrans= " << itrans << " exp_nr_trans= " << exp_nr_trans << std::endl;
				REPORT_ERROR("ihidden > XSIZE: ");
			}
#endif
			if (DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden))
				return true;
		}
	}
	return false;

}

void MlOptimiser::doThreadGetSquaredDifferencesAllOrientations(int thread_id)
{
#ifdef DEBUG_THREAD
    std::cerr << "entering doThreadGetAllSquaredDifferences" << std::endl;
#endif
    
    // Local variables
    std::vector<double> thisthread_min_diff2;
    std::vector< Matrix1D<double> > oversampled_orientations, oversampled_translations;
    MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_shift;
    MultidimArray<double> Fctf, Minvsigma2;
    Matrix2D<double> A;
    
    // Initialise local mindiff2 for thread-safety
    thisthread_min_diff2.clear();
    thisthread_min_diff2.resize(exp_nr_particles, 99.e99);
    Fref.resize(exp_local_Fimgs_shifted[0]);
    Frefctf.resize(exp_local_Fimgs_shifted[0]);
    
    // THESE TWO FOR LOOPS WILL BE PARALLELISED USING THREADS...
    // exp_iclass loop does not always go from 0 to nr_classes!
    long int iorientclass_offset = exp_iclass * exp_nr_rot;
    
    size_t first_iorient = 0, last_iorient = 0;
    while (exp_iorient_ThreadTaskDistributor->getTasks(first_iorient, last_iorient))
    {
        for (long int iorient = first_iorient; iorient <= last_iorient; iorient++)
        {
            
            long int iorientclass = iorientclass_offset + iorient;
            long int idir = iorient / exp_nr_psi;
            long int ipsi = iorient % exp_nr_psi;
            // Get prior for this direction and skip calculation if prior==0
            double pdf_orientation;
            if (mymodel.orientational_prior_mode == NOPRIOR)
            {
#ifdef DEBUG_CHECKSIZES
                if (idir >= XSIZE(mymodel.pdf_direction[exp_iclass]))
                {
                    std::cerr<< "idir= "<<idir<<" XSIZE(mymodel.pdf_direction[exp_iclass])= "<< XSIZE(mymodel.pdf_direction[exp_iclass]) <<std::endl;
                    REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
                }
#endif
                pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
            }
            else
            {
                pdf_orientation = sampling.getPriorProbability(idir, ipsi);
            }
            
            // In the first pass, always proceed
            // In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
            // if so, proceed with projecting the reference in that direction
            bool do_proceed = (exp_ipass==0) ? true : isSignificantAnyParticleAnyTranslation(iorientclass);
            
            if (do_proceed && pdf_orientation > 0.)
            {
                // Now get the oversampled (rot, tilt, psi) triplets
                // This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
                sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_orientations);
                
#ifdef DEBUG_CHECKSIZES
                if (exp_nr_oversampled_rot != oversampled_orientations.size())
                {
                    std::cerr<< "exp_nr_oversampled_rot= "<<exp_nr_oversampled_rot<<" oversampled_orientations.size()= "<< oversampled_orientations.size() <<std::endl;
                    REPORT_ERROR("exp_nr_oversampled_rot != oversampled_orientations.size()");
                }
#endif
                // Loop over all oversampled orientations (only a single one in the first pass)
                for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
                {
                    
                    // Get the Euler matrix
                    Euler_angles2matrix(XX(oversampled_orientations[iover_rot]),
                                        YY(oversampled_orientations[iover_rot]),
                                        ZZ(oversampled_orientations[iover_rot]), A);
                    
                    // Take tilt-series into account
                    A = (exp_R_mic * A).inv();
                    
                    // Project the reference map (into Fref)
#ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (thread_id == 0)
                        timer.tic(TIMING_DIFF_PROJ);
#endif
                    (mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_INV);
#ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (thread_id == 0)
                        timer.toc(TIMING_DIFF_PROJ);
#endif
                    
                    /// Now that reference projection has been made loop over someParticles!
                    for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                    {
                        // loop over all particles inside this ori_particle
                        for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                        {
                            long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                            
                            bool is_last_image_in_series = mydata.getNrImagesInSeries(part_id) == (exp_iseries + 1);
                            // Which number was this image in the combined array of exp_iseries and part_id
                            long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                            
#ifdef DEBUG_CHECKSIZES
                            if (my_image_no >= exp_local_Minvsigma2s.size())
                            {
                                std::cerr<< "my_image_no= "<<my_image_no<<" exp_local_Minvsigma2s.size()= "<< exp_local_Minvsigma2s.size() <<std::endl;
                                REPORT_ERROR("my_image_no >= exp_local_Minvsigma2.size()");
                            }
#endif
                            Minvsigma2 = exp_local_Minvsigma2s[my_image_no];
                            
                            // Apply CTF to reference projection
                            if (do_ctf_correction && refs_are_ctf_corrected)
                            {
                                
#ifdef DEBUG_CHECKSIZES
                                if (my_image_no >= exp_local_Fctfs.size())
                                {
                                    std::cerr<< "my_image_no= "<<my_image_no<<" exp_local_Fctfs.size()= "<< exp_local_Fctfs.size() <<std::endl;
                                    REPORT_ERROR("my_image_no >= exp_local_Fctfs.size()");
                                }
                                if (MULTIDIM_SIZE(Fref) != MULTIDIM_SIZE(exp_local_Fctfs[my_image_no]))
                                {
                                    std::cerr<< "MULTIDIM_SIZE(Fref)= "<<MULTIDIM_SIZE(Fref)<<" MULTIDIM_SIZE()= "<< MULTIDIM_SIZE(exp_local_Fctfs[my_image_no]) <<std::endl;
                                    REPORT_ERROR("MULTIDIM_SIZE(Fref) != MULTIDIM_SIZE(exp_local_Fctfs[my_image_no)");
                                }
                                
#endif
                                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
                                {
                                    DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[my_image_no], n);
                                }
                            }
                            else
                                Frefctf = Fref;
                            
                            if (do_scale_correction)
                            {
                                int group_id = mydata.getGroupId(part_id, exp_iseries);
#ifdef DEBUG_CHECKSIZES
                                if (group_id >= mymodel.scale_correction.size())
                                {
                                    std::cerr<< "group_id= "<<group_id<<" mymodel.scale_correction.size()= "<< mymodel.scale_correction.size() <<std::endl;
                                    REPORT_ERROR("group_id >= mymodel.scale_correction.size()");
                                }
#endif
                                double myscale = mymodel.scale_correction[group_id];
                                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
                                {
                                    DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
                                }
                            }
                            
                            for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
                            {
                                long int ihidden = iorientclass * exp_nr_trans + itrans;
                                
#ifdef DEBUG_CHECKSIZES
                                if (exp_ipass > 0 && ihidden >= XSIZE(exp_Mcoarse_significant))
                                {
                                    std::cerr<< "ihidden= "<<ihidden<<" XSIZE(exp_Mcoarse_significant)= "<< XSIZE(exp_Mcoarse_significant) <<std::endl;
                                    REPORT_ERROR("ihidden >= XSIZE(exp_Mcoarse_significant)");
                                }
#endif
                                // In the first pass, always proceed
                                // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                                bool do_proceed = (exp_ipass == 0) ? true : DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden);

                                if (do_proceed)
                                {
                                    
                                    sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);
                                    for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                                    {
#ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (thread_id == 0)
                                            timer.tic(TIMING_DIFF_DIFF2);
#endif
                                        // Get the shifted image
                                        long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
                                        itrans * exp_nr_oversampled_trans + iover_trans;
                                        
#ifdef DEBUG_CHECKSIZES
                                        if (ishift >= exp_local_Fimgs_shifted.size())
                                        {
                                            std::cerr<< "ishift= "<<ishift<<" exp_local_Fimgs_shifted.size()= "<< exp_local_Fimgs_shifted.size() <<std::endl;
                                            std::cerr << " itrans= " << itrans << std::endl;
                                            std::cerr << " ipart= " << ipart << std::endl;
                                            std::cerr << " exp_nr_oversampled_trans= " << exp_nr_oversampled_trans << " exp_nr_trans= " << exp_nr_trans << " iover_trans= " << iover_trans << std::endl;
                                            REPORT_ERROR("ishift >= exp_local_Fimgs_shifted.size()");
                                        }
#endif
                                        
                                        Fimg_shift = exp_local_Fimgs_shifted[ishift];
                                        //#define DEBUG_GETALLDIFF2
#ifdef DEBUG_GETALLDIFF2
                                        if (verb> 0)
                                        {
                                            FourierTransformer transformer;
                                            Image<double> tt;
                                            tt().resize(exp_current_image_size, exp_current_image_size);
                                            transformer.inverseFourierTransform(Fimg_shift, tt());
                                            CenterFFT(tt(),false);
                                            tt.write("Fimg_shift.spi");
                                            transformer.inverseFourierTransform(Frefctf, tt());
                                            CenterFFT(tt(),false);
                                            tt.write("Fref.spi");
                                            tt()=Minvsigma2;
                                            tt.write("Minvsigma2.spi");
                                            std::cerr << "written Minvsigma2.spi" << std::endl;
                                            
                                            char c;
                                            std::cerr << "Written Fimg_shift.spi and Fref.spi. Press any key to continue..." << std::endl;
                                            std::cin >> c;
                                            exit(1);
                                        }
#endif
                                        
                                        double diff2;
                                        if ((iter == 1 && do_firstiter_cc) || do_always_cc)
                                        {
                                            // Do not calculate squared-differences, but signal product
                                            // Negative values because smaller is worse in this case
                                            diff2 = 0.;
                                            double suma2 = 0.;
                                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                                            {
                                                diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                                diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                                suma2 += norm(DIRECT_MULTIDIM_ELEM(Frefctf, n));
                                            }
                                            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
                                            diff2 /= sqrt(suma2) * exp_local_sqrtXi2[my_image_no];
                                        }
                                        else
                                        {
                                            
#ifdef DEBUG_CHECKSIZES
                                            if (my_image_no >= exp_highres_Xi2_imgs.size())
                                            {
                                                std::cerr<< "my_image_no= "<<my_image_no<<" exp_highres_Xi2_imgs.size()= "<< exp_highres_Xi2_imgs.size() <<std::endl;
                                                REPORT_ERROR("my_image_no >= exp_highres_Xi2_imgs.size()");
                                            }
#endif
                                            
                                            // Calculate the actual squared difference term of the Gaussian probability function
                                            // If current_size < mymodel.ori_size diff2 is initialised to the sum of
                                            // all |Xij|2 terms that lie between current_size and ori_size
                                            // Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
                                            diff2 = exp_highres_Xi2_imgs[my_image_no] / 2.;
                                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                                            {
                                                double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                                double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                                diff2 += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
                                            }
                                            
                                        }
                                        
#ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (thread_id == 0)
                                            timer.toc(TIMING_DIFF_DIFF2);
#endif
                                        
                                        // Store all diff2 in exp_Mweight
                                        long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                             iover_rot, iover_trans);
                                        
                                        //#define DEBUG_DIFF2_ISNAN
#ifdef DEBUG_DIFF2_ISNAN
                                        if (std::isnan(diff2))
                                        {
                                            global_mutex.lock();
                                            std::cerr << " ipart= " << ipart << std::endl;
                                            std::cerr << " diff2= " << diff2 << " thisthread_min_diff2[ipart]= " << thisthread_min_diff2[ipart] << " ipart= " << ipart << std::endl;
                                            std::cerr << " exp_highres_Xi2_imgs[my_image_no]= " << exp_highres_Xi2_imgs[my_image_no] << std::endl;
                                            std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
                                            std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
                                            std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
                                            std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
                                            std::cerr << " ihidden_over= " << ihidden_over << " XSIZE(Mweight)= " << XSIZE(exp_Mweight) << std::endl;
                                            int group_id = mydata.getGroupId(part_id, exp_iseries);
                                            std::cerr << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
                                            if (std::isnan(mymodel.scale_correction[group_id]))
                                            {
                                                for (int i=0; i < mymodel.scale_correction.size(); i++)
                                                    std::cerr << " i= " << i << " mymodel.scale_correction[i]= " << mymodel.scale_correction[i] << std::endl;
                                            }
                                            std::cerr << " group_id= " << group_id << std::endl;
                                            Image<double> It;
                                            It()=Minvsigma2;
                                            It.write("Minvsigma2.spi");
                                            std::cerr << "written Minvsigma2.spi" << std::endl;
                                            std::cerr << "Frefctf shape= "; Frefctf.printShape(std::cerr);
                                            std::cerr << "Fimg_shift shape= "; Fimg_shift.printShape(std::cerr);
                                            It()=exp_local_Fctfs[my_image_no];
                                            It.write("exp_local_Fctf.spi");
                                            std::cerr << "written exp_local_Fctf.spi" << std::endl;
                                            FourierTransformer transformer;
                                            Image<double> tt;
                                            tt().resize(exp_current_image_size, exp_current_image_size);
                                            transformer.inverseFourierTransform(Fimg_shift, tt());
                                            CenterFFT(tt(),false);
                                            tt.write("Fimg_shift.spi");
                                            std::cerr << "written Fimg_shift.spi" << std::endl;
                                            FourierTransformer transformer2;
                                            tt().initZeros();
                                            transformer2.inverseFourierTransform(Frefctf, tt());
                                            CenterFFT(tt(),false);
                                            tt.write("Frefctf.spi");
                                            std::cerr << "written Frefctf.spi" << std::endl;
                                            FourierTransformer transformer3;
                                            tt().initZeros();
                                            transformer3.inverseFourierTransform(Fref, tt());
                                            CenterFFT(tt(),false);
                                            tt.write("Fref.spi");
                                            std::cerr << "written Fref.spi" << std::endl;
                                            std::cerr << " A= " << A << std::endl;
                                            std::cerr << " exp_R_mic= " << exp_R_mic << std::endl;
                                            std::cerr << "written Frefctf.spi" << std::endl;
                                            REPORT_ERROR("diff2 is not a number");
                                            global_mutex.unlock();
                                        }
#endif
                                        //#define DEBUG_VERBOSE
#ifdef DEBUG_VERBOSE
                                        global_mutex.lock();
                                        if (verb > 0)
                                        {
                                            std::cout << " rot= " << XX(oversampled_orientations[iover_rot]) << " tilt= "<< YY(oversampled_orientations[iover_rot]) << " psi= " << ZZ(oversampled_orientations[iover_rot]) << std::endl;
                                            std::cout << " xoff= " <<XX(oversampled_translations[iover_trans]) <<" yoff= "<<YY(oversampled_translations[iover_trans])<<std::endl;
                                            std::cout << " ihidden_over= " << ihidden_over << " diff2= " << diff2 << " thisthread_min_diff2[ipart]= " << thisthread_min_diff2[ipart] << std::endl;
                                        }
                                        global_mutex.unlock();
#endif
#ifdef DEBUG_CHECKSIZES
                                        if (ihidden_over >= XSIZE(exp_Mweight) )
                                        {
                                            std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
                                            std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
                                            std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
                                            std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
                                            std::cerr << " ihidden_over= " << ihidden_over << " XSIZE(Mweight)= " << XSIZE(exp_Mweight) << std::endl;
                                            REPORT_ERROR("ihidden_over >= XSIZE(Mweight)");
                                        }
#endif
                                        
                                        if (exp_iseries == 0)
                                            DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;
                                        else
                                            DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += diff2;
                                        
#ifdef DEBUG_CHECKSIZES
                                        if (ipart >= thisthread_min_diff2.size())
                                        {
                                            std::cerr<< "ipart= "<<ipart<<" thisthread_min_diff2.size()= "<< thisthread_min_diff2.size() <<std::endl;
                                            REPORT_ERROR("ipart >= thisthread_min_diff2.size() ");
                                        }
#endif
                                        // Keep track of minimum of all diff2, only for the last image in this series
                                        diff2 = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                        //std::cerr << " exp_ipass= " << exp_ipass << " exp_iclass= " << exp_iclass << " diff2= " << diff2 << std::endl;
                                        if (is_last_image_in_series && diff2 < thisthread_min_diff2[ipart])
                                            thisthread_min_diff2[ipart] = diff2;
                                        
                                    } // end loop iover_trans
                                } // end if do_proceed translations
                            } // end loop itrans
                        } // end loop part_id (i)
                    } // end loop ori_part_id
                }// end loop iover_rot
            } // end if do_proceed orientations
        } // end loop iorient
    } // end while task distribution
    
    
    // Now inside a mutex set the minimum of the squared differences among all threads
#ifdef DEBUG_CHECKSIZES
    if (thisthread_min_diff2.size() != exp_min_diff2.size())
    {
        std::cerr<< "thisthread_min_diff2.size()= "<<thisthread_min_diff2.size()<<" exp_min_diff2.size()= "<< exp_min_diff2.size() <<std::endl;
        REPORT_ERROR("thisthread_min_diff2.size() != exp_min_diff2.size()");
    }
#endif
    
    global_mutex.lock();
    for (int i = 0; i < exp_min_diff2.size(); i++)
    {
        if (thisthread_min_diff2[i] < exp_min_diff2[i])
        {
            exp_min_diff2[i] = thisthread_min_diff2[i];
        }
    }
    global_mutex.unlock();
    
    // Wait until all threads have finished
    global_barrier->wait();
    
#ifdef DEBUG_THREAD
    std::cerr << "leaving doThreadGetAllSquaredDifferences" << std::endl;
#endif
    
}



void MlOptimiser::enqueueCLMemoryWrite(cl_mem destination, void *source, size_t s) {
    int err;
    err = clEnqueueWriteBuffer(CL_CopyToDeviceQueue, destination, true, 0, s, source, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        std::cerr << "Error writing data to OpenCL device, error: " << err;
    }
    clFinish(CL_CopyToDeviceQueue);
}

void MlOptimiser::blankCLMemory(cl_mem destination, size_t s) {
 
    int err;
    err = clSetKernelArg(CL_blank, 0, sizeof(cl_mem), &destination);
    
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to set OpenCL kernel arguments - blank" << std::endl;
    }
    
    //    std::cerr << "Finished loading CL data" << std::endl;
    int wgs = CL_WGS;
    size_t global[1], local[1];
    //    global[0] = wgs * ceil((float)limit / (float)wgs);;
    global[0] = s;
    local[0] = wgs;
    
    err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_blank, 1, NULL, global, local, 0, NULL, NULL);
    if (err)
    {
        std::cerr << "Error: Failed to execute OpenCL kernel - blank, error: " << err << std::endl;
    }
//    clFinish(CL_ComputeQueue);
}

void MlOptimiser::initCLMemoryWithConstant(cl_mem destination, size_t s, int value) {
    
    cl_uint limit = s;
    int err;
    err = clSetKernelArg(CL_initWithConstant, 0, sizeof(cl_mem), &destination);
    err |= clSetKernelArg(CL_initWithConstant, 1, sizeof(cl_int), &value);
    err |= clSetKernelArg(CL_initWithConstant, 2, sizeof(cl_uint), &limit);
    
    
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to set OpenCL kernel arguments - initWithConstant" << std::endl;
    }
    
    //    std::cerr << "Finished loading CL data" << std::endl;
    int wgs = CL_WGS;
    size_t global[1], local[1];
    //    global[0] = wgs * ceil((float)limit / (float)wgs);;
    global[0] = wgs * ceil((float)s / (float)wgs);
    local[0] = wgs;
    
    err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_initWithConstant, 1, NULL, global, local, 0, NULL, NULL);
    if (err)
    {
        std::cerr << "Error: Failed to execute OpenCL kernel - blank, error: " << err << std::endl;
    }
//    clFinish(CL_ComputeQueue);
}

void MlOptimiser::initCLMemoryWithConstantDouble(cl_mem destination, size_t s, double value) {
    
    cl_uint limit = s;
    int err;
    err = clSetKernelArg(CL_initWithConstantDouble, 0, sizeof(cl_mem), &destination);
    err |= clSetKernelArg(CL_initWithConstantDouble, 1, sizeof(cl_double), &value);
    err |= clSetKernelArg(CL_initWithConstantDouble, 2, sizeof(cl_uint), &limit);

    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to set OpenCL kernel arguments - initWithConstant" << std::endl;
    }
    
    //    std::cerr << "Finished loading CL data" << std::endl;
    int wgs = CL_WGS;
    size_t global[1], local[1];
    //    global[0] = wgs * ceil((float)limit / (float)wgs);;
    global[0] = wgs * ceil((float)s / (float)wgs);
    local[0] = wgs;
    
    err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_initWithConstantDouble, 1, NULL, global, local, 0, NULL, NULL);
    if (err)
    {
        std::cerr << "Error: Failed to execute OpenCL kernel - blank, error: " << err << std::endl;
    }
//    clFinish(CL_ComputeQueue);
}

void MlOptimiser::doOpenCLGetSquaredDifferencesAllOrientations()
{
#ifdef DEBUG_OPENCL
    std::cerr << "entering doOpenCLetAllSquaredDifferences" << std::endl;
#endif
    
#ifdef TIMING
    timer.tic(TIMING_DIFF_PROJ);
#endif
    


    // Local variables
    std::vector<double> thisthread_min_diff2;
    std::vector< Matrix1D<double> > oversampled_orientations, oversampled_translations;
    MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_shift;
    MultidimArray<double> Fctf, Minvsigma2;
    Matrix2D<double> A;
    struct timeval prepStartTV, startTV, endTV;
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&prepStartTV, NULL);
#endif
    
    // Initialise local mindiff2 for thread-safety
    Fref.resize(exp_local_Fimgs_shifted[0]);
    Frefctf.resize(exp_local_Fimgs_shifted[0]);
    
    
    if (((exp_ipass == 0) && (CL_GSD_memoryCalculatedCoarseForIter < iter)) ||
        ((exp_ipass == 1) && (CL_GSD_memoryCalculatedFineForIter < iter))) {
        unsigned long device_memory = CL_global_memsize;
        
        //Subtract a 256 Mb space, just in case...
        device_memory -= (1 << 28);
        
        //Subtract various small things
        //Model
        MultidimArray<Complex > modelData = (mymodel.PPref[exp_iclass]).data;
        int model_size = XSIZE(modelData) * YSIZE(modelData) * ZSIZE(modelData);
        device_memory -= sizeof(cl_double2) * model_size;
        
        //local Fctf's and Minsigmav2
        int FSize = XSIZE(Fref) * YSIZE(Fref);
        int wgsMultipleToFitFSize = CL_WGS * (FSize / CL_WGS + 1);
        device_memory -= 4 * sizeof(cl_double) * wgsMultipleToFitFSize;
        
        //exp_r_mic
        device_memory -= sizeof(cl_double) * 9;
        
        //Fimg shifts
        //worst case scenario here...
        int num_trans = exp_nr_trans * exp_nr_oversampled_trans;
        int numFsTrans = CL_WGS * ceil((float)num_trans * FSize / (float)CL_WGS);
        device_memory -= 2 * sizeof(cl_double2) * numFsTrans;
        
        //Now calculate max orients with the remaining memory
        int maxOrientsIntoGlobal = device_memory / (sizeof(cl_double) * (4 * (FSize) + 4 * num_trans + 9 + 2 * 3) + sizeof(cl_char) * (2 * num_trans));
        
//        std::cerr << "Global memsize: " << CL_global_memsize << " memory after small items: " << device_memory << " meaning maxOrients: " << maxOrientsIntoGlobal << std::endl;
        
        int largestAlloc = 2 * FSize > num_trans ? 2 * FSize : num_trans;
        int maxOrientsPerAlloc = CL_maxMemAlloc / (sizeof(cl_double) * largestAlloc);

        if (exp_ipass == 0) {
            if (maxOrientsIntoGlobal < maxOrientsPerAlloc) {
                CL_GSD_maxOrientsPerCycleCoarse = maxOrientsIntoGlobal;
            } else {
                CL_GSD_maxOrientsPerCycleCoarse = maxOrientsPerAlloc;
            }
            CL_GSD_memoryCalculatedCoarseForIter = iter;
        } else {
            if (maxOrientsIntoGlobal < maxOrientsPerAlloc) {
                CL_GSD_maxOrientsPerCycleFine = maxOrientsIntoGlobal;
            } else {
                CL_GSD_maxOrientsPerCycleFine = maxOrientsPerAlloc;
            }
            CL_GSD_memoryCalculatedFineForIter = iter;
        }
        
//        std::cerr << "Max alloc: " << CL_maxMemAlloc << " meaning maxOrients: " << maxOrientsPerAlloc << " Final coarse maxOrients: " << CL_GSD_maxOrientsPerCycleCoarse << " Final fine maxOrients: " << CL_GSD_maxOrientsPerCycleFine << std::endl;
    }
    
    int maxOrientsPerCycle;
    if (exp_ipass == 0) {
        maxOrientsPerCycle = CL_GSD_maxOrientsPerCycleCoarse;
//        maxOrientsPerCycle = 3000;
    } else {
        maxOrientsPerCycle = CL_GSD_maxOrientsPerCycleFine;
//        maxOrientsPerCycle = 3000;
    }
    
    // THESE TWO FOR LOOPS WILL BE PARALLELISED USING THREADS...
    // exp_iclass loop does not always go from 0 to nr_classes!
    

    //Count the nuumber of particles
    int partCount = 0;
    for (long int ori_part_id = exp_my_first_ori_particle; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
    {
        // loop over all particles inside this ori_particle
        //                    std::cerr << "Number of particles for ori_part_id: " << ori_part_id << " is: " << mydata.ori_particles[ori_part_id].particles_id.size() << std::endl;
        for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, partCount++)
        {
        }
    }
    
    //First get orientations in a matrix for OpenCL
    //Need to first count the number of orientations...
    long int iorientclass_offset = exp_iclass * exp_nr_rot;
    
    size_t first_iorient = 0, last_iorient = 0;
    long int nr_orients = sampling.NrDirections() * sampling.NrPsiSamplings();
    long int totalOrients = 0;

    //Get number of orients for openCL save data...
    std::vector<bool> translationUsedInAnyRotation;
    translationUsedInAnyRotation.resize(exp_nr_trans * exp_nr_oversampled_trans, false);

    std::vector<int> rotationsForParticle;
    rotationsForParticle.resize(partCount);
    
    int totalPosToCheck = 0;
    
    for (long int iorient = 0; iorient < nr_orients; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        long int idir = iorient / exp_nr_psi;
        long int ipsi = iorient % exp_nr_psi;
        // Get prior for this direction and skip calculation if prior==0
        double pdf_orientation;
        if (mymodel.orientational_prior_mode == NOPRIOR)
        {
            pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
        }
        else
        {
            pdf_orientation = sampling.getPriorProbability(idir, ipsi);
        }
        
        // In the first pass, always proceed
        // In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
        // if so, proceed with projecting the reference in that direction
        bool do_proceed = (exp_ipass==0) ? true : isSignificantAnyParticleAnyTranslation(iorientclass);
        
        if (do_proceed && pdf_orientation > 0.)
        {
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                bool anySignificantWeightTranslations = false;
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        bool particleAnySignificantWeightTranslations = false;
                        
                        for (long int itrans = 0, currentTrans = 0; itrans < exp_nr_trans; itrans++)
                        {
                            long int ihidden = iorientclass * exp_nr_trans + itrans;
                            
                            // In the first pass, always proceed
                            // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                            bool do_proceed = (exp_ipass == 0) ? true : DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden);
                            
                                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                                {
                                    
                                    if (do_proceed) {
                                        anySignificantWeightTranslations = true;
                                        particleAnySignificantWeightTranslations = true;
                                        translationUsedInAnyRotation[currentTrans] = true;
//                                        particleTranslationUsedInAnyRotation[ipart * exp_nr_trans * exp_nr_oversampled_trans + currentTrans] = true;
                                        totalPosToCheck++;
                                        //                                    std::cerr << "ipart: " << ipart << " itrans: " << itrans << " iover_trans: " << iover_trans << " current_trans: " << currentTrans << std::endl;
                                    }
                                    currentTrans++;
                                } // end loop iover_trans

                        } // end loop itrans
                        if (particleAnySignificantWeightTranslations) {
                            rotationsForParticle[ipart]++;
                        }
                    }// end loop part_id (i)
                } // end loop ori_part_id
                if (anySignificantWeightTranslations) {
                    totalOrients++;
                }
            }// end if iover_rot
        }// end loop do_proceed
    }
    
    int usedTrans = 0;
    for (int i = 0; i < exp_nr_trans * exp_nr_oversampled_trans; i++) {
        if (translationUsedInAnyRotation[i]) {
            usedTrans++;
        }
    }
    
    
//    if (CL_GSD_memoryCalculatedForIter < iter) {
    
//        std::cerr << "Total orients: " << totalOrients << " Total trans: " << usedTrans << " orients x trans x particles = " << totalOrients * usedTrans * 8 << " totalPosToCheck: " << totalPosToCheck << " = " << (float)totalPosToCheck / (float)(totalOrients * usedTrans * 8) * 100.0 << "%" << std::endl;
//        if (exp_ipass == 1) {
//            CL_GSD_memoryCalculatedForIter = iter;
//        }
//    }

    
    if (totalOrients == 0) {
//        std::cerr << "No orients!" << std::endl;
        return;
    }

    //    std::cerr << "Rotations per particle" << std::endl;
    
    //Make all allocations for orients the size of maxOrients - only need to allocate memory once then
    int maxOrients = 0;
    for (int i = 0; i < partCount; i++) {
        //        int partUsedTrans = 0;
        //        for (int j = 0; j < exp_nr_trans * exp_nr_oversampled_trans; j++) {
        //            if (particleTranslationUsedInAnyRotation[i * exp_nr_trans * exp_nr_oversampled_trans + j]) {
        //                partUsedTrans++;
        //            }
        //        }
        //        std::cerr << "Particle: " << i << " rotations: " << rotationsForParticle[i] << " trans: " << partUsedTrans << std::endl;
        if (maxOrients < rotationsForParticle[i]) {
            maxOrients = rotationsForParticle[i];
        }
    }

    if (maxOrients > maxOrientsPerCycle) {
        maxOrients = maxOrientsPerCycle;
    }
    
#pragma mark Cycles Variables
    int numCycles = 0;
    std::vector<int> particleForCycle;
    std::vector<int> orientsPerCycle;
    std::vector<int> startOrientPerCycle;
    for (int i = 0; i < partCount; i++) {
        int particleRotations = rotationsForParticle[i];
        int cyclesThisParticle = 0;
        while (particleRotations > 0) {
            numCycles++;
            if (particleRotations >= maxOrientsPerCycle) {
                particleRotations -= maxOrientsPerCycle;
                orientsPerCycle.push_back(maxOrientsPerCycle);
            } else {
                orientsPerCycle.push_back(particleRotations);
                particleRotations = 0;
            }
            particleForCycle.push_back(i);
            startOrientPerCycle.push_back(cyclesThisParticle * maxOrientsPerCycle);
            cyclesThisParticle++;
        }
    }
    
//    std::cerr << "Num cycles: " << numCycles << std::endl;
//    for (int i = 0; i < numCycles; i++) {
//        std::cerr << "Cycle " << i << " particle: " << particleForCycle[i] << " orients: " << orientsPerCycle[i] << " start orient: " << startOrientPerCycle[i] << std::endl;
//    }
    
    //Load model into OpenCL device Memory
    MultidimArray<Complex > modelData = (mymodel.PPref[exp_iclass]).data;
    int model_size = XSIZE(modelData) * YSIZE(modelData) * ZSIZE(modelData);
    cl_mem cl_modelData = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * model_size, NULL, 0);
    enqueueCLMemoryWrite(cl_modelData, modelData.data, sizeof(cl_double2) * model_size);
    
    //    std::cerr << "Model mapped, size (bytes): " << sizeof(cl_double2) * model_size << std::endl;

    //Allocate memory for the orientations
    cl_uint2 projDim; projDim.x = XSIZE(Fref); projDim.y = YSIZE(Fref);
//    if (iter > 7) {
//        std::cerr << "Total orients: " << totalOrients << " required memory: " << totalOrients * 2 * sizeof(cl_double2) * projDim.x * projDim.y << std::endl;
//    }

//    std::cerr << "Proj Dim X: " << projDim.x << " y: " << projDim.y << std::endl;
//    cl_mem cl_projOut = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * totalOrients * projDim.x * projDim.y, NULL, 0);

    int err;
    int wgs = CL_WGS;
    int numWGSPerProj = ceil((float)(projDim.x * projDim.y) / (float)wgs);
    int numFs = wgs * ceil((float)maxOrients * projDim.x * projDim.y / (float)wgs);
//    std::cerr << "ProjDim: " << projDim.x * projDim.y << " wgs: " << wgs << " numWGSPerProj: " << numWGSPerProj << std::endl;
    
//    if ((iter > CL_Proj_alloc_iter) || (CL_Proj_allocated_orients != totalOrients)) {
//        std::cerr << "Reallocating pass 0 memory" << std::endl;
//        if (CL_iter_Proj) {
//            clReleaseMemObject(CL_iter_Proj);
//        }
//        CL_iter_Proj = clCreateImage(CL_context, CL_MEM_READ_ONLY, &fmt, &desc, NULL, &err);
        cl_mem CL_iter_Proj = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * numFs, NULL, 0);
        blankCLMemory(CL_iter_Proj, 2 * numFs);
//    }

//    cl_mem cl_Proj_ipass1;
//    if (exp_ipass > 0) {
//        cl_Proj_ipass1 = clCreateImage(CL_context, CL_MEM_READ_ONLY, &fmt, &desc, NULL, &err);
//        cl_Proj_ipass1 = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * totalOrients * numWGSPerProj * wgs, NULL, 0);
//        blankCLMemory(cl_Proj_ipass1, 2 * totalOrients * numWGSPerProj * wgs);
//    }
    cl_mem cl_FrefCTF = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * numFs, NULL, 0);
    
    //Allocate memory for the elements associated with each particle
    double *scales = (double *)malloc(sizeof(double) * partCount);
    //    cl_mem cl_scale = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * partCount, NULL, 0);
    //    desc.image_height = partCount;
    //    desc.image_width = projDim.x * projDim.y / 2;
    //    cl_mem cl_local_Fctfs = clCreateImage(CL_context, CL_MEM_READ_ONLY, &fmt, &desc, NULL, &err);
    cl_mem cl_local_FctfsA = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * numWGSPerProj * wgs, NULL, 0);
    cl_mem cl_Minvsigma2A = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * numWGSPerProj * wgs, NULL, 0);
    cl_mem cl_local_FctfsB = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * numWGSPerProj * wgs, NULL, 0);
    cl_mem cl_Minvsigma2B = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * numWGSPerProj * wgs, NULL, 0);
    blankCLMemory(cl_Minvsigma2A, numWGSPerProj * wgs);
    blankCLMemory(cl_local_FctfsA, numWGSPerProj * wgs);
    blankCLMemory(cl_Minvsigma2B, numWGSPerProj * wgs);
    blankCLMemory(cl_local_FctfsB, numWGSPerProj * wgs);
    
    //    double **FctfArray = (double **)malloc(sizeof(double *) * partCount);
    long int *myImageArray = (long int *)malloc(sizeof(long int) * partCount);
    
    //Allocate memory for the Fimg_shifts
    //    int totalTrans = exp_nr_trans * exp_nr_oversampled_trans * partCount;
    int totalTrans = usedTrans * partCount;
    
    //    cl_mem cl_FimgShift = clCreateImage(CL_context, CL_MEM_READ_ONLY, &fmt, &desc, NULL, &err);
    //    clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * totalTrans * projDim.x * projDim.y, NULL, 0);
    long int *FimgShiftArray = (long int *)malloc(sizeof(long int) * totalTrans);
//    cl_mem cl_FimgShift = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * numWGSPerProj * wgs, NULL, 0);
    //    std::cerr << "Image memory allocated" << std::endl;
    
    int *lastImageInSeries = (int *)malloc(sizeof(int) * partCount);
    
    long int *iTransClassArray = (long int *)malloc(sizeof(long int) * totalTrans);
    long int *iOverTransClassArray = (long int *)malloc(sizeof(long int) * totalTrans);
    
    //Now gather the data
    int iUsedTrans = 0;
    for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
    {
        // loop over all particles inside this ori_particle
        for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
        {
            long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
            
            bool is_last_image_in_series = mydata.getNrImagesInSeries(part_id) == (exp_iseries + 1);
            lastImageInSeries[ipart] = is_last_image_in_series;
            // Which number was this image in the combined array of exp_iseries and part_id
            long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
            myImageArray[ipart] = my_image_no;
            
            // Apply CTF to reference projection
            if (do_ctf_correction && refs_are_ctf_corrected)
            {
                //                origin[1] = ipart;
                //                region[0] = projDim.x * projDim.y / 2;
                //                clEnqueueWriteImage(CL_queue, cl_local_Fctfs, true, origin, region, 0, 0, (const void *)exp_local_Fctfs[my_image_no].data, 0, NULL, NULL);
                //                clEnqueueWriteBuffer(CL_queue, cl_local_Fctfs, true, ipart * sizeof(cl_double) * numWGSPerProj * wgs, sizeof(cl_double) * projDim.x * projDim.y, exp_local_Fctfs[my_image_no].data, 0, NULL, NULL);
                //                FctfArray[ipart] = exp_local_Fctfs[my_image_no].data;
            }
            
            
            //Scale
            if (do_scale_correction)
            {
                int group_id = mydata.getGroupId(part_id, exp_iseries);
                double myscale = mymodel.scale_correction[group_id];
                scales[ipart] = myscale;
                
                //                clEnqueueWriteBuffer(CL_queue, cl_scale, true, ipart * sizeof(cl_double), sizeof(cl_double), &myscale, 0, NULL, NULL);
                
            }
            
            //Other data
            //            if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
            //            if (false) {
            //                double sqrtXi2 = exp_local_sqrtXi2[my_image_no];
            //                clEnqueueWriteBuffer(CL_queue, cl_sqrtXi2, true, ipart * sizeof(cl_double), sizeof(cl_double), &sqrtXi2, 0, NULL, NULL);
            //            } else {
            //                Minvsigma2 = exp_local_Minvsigma2s[my_image_no];
            //                clEnqueueWriteImage(CL_queue, cl_Minvsigma2, true, origin, region, 0, 0, (const void *)Minvsigma2.data, 0, NULL, NULL);
            //                clEnqueueWriteBuffer(CL_queue, cl_Minvsigma2, true, ipart * sizeof(cl_double) * projDim.x * projDim.y, sizeof(cl_double) * projDim.x * projDim.y, Minvsigma2.data, 0, NULL, NULL);
            
            //                double exp_highres_Xi2 = exp_highres_Xi2_imgs[my_image_no];
            //                clEnqueueWriteBuffer(CL_queue, cl_exp_highres_Xi2_imgs, true, ipart * sizeof(cl_double), sizeof(cl_double), &exp_highres_Xi2, 0, NULL, NULL);
            
            //            }
            
            
            //Now gather itrans for the particle
            for (long int itrans = 0, currentTrans = 0; itrans < exp_nr_trans; itrans++)
            {
                
                //                sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);
                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                {
                    // Get the shifted image
                    long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
                    itrans * exp_nr_oversampled_trans + iover_trans;
                    
                    if (translationUsedInAnyRotation[currentTrans]) {
                        //                        std::cerr << "iused trans: " << iUsedTrans <<  " itrans: " << itrans << " iover_trans: " << iover_trans << " ipart: " << ipart << std::endl;
                        iTransClassArray[iUsedTrans] = itrans;
                        iOverTransClassArray[iUsedTrans] = iover_trans;
                        FimgShiftArray[iUsedTrans] = ishift;
                        iUsedTrans++;
                    }
                    //                    Fimg_shift = exp_local_Fimgs_shifted[ishift];
                    //                    clEnqueueWriteImage(CL_queue, cl_FimgShift, true, origin, region, 0, 0, (const void *)Fimg_shift.data, 0, NULL, NULL);
                    //                    FimgArray[trans_current] = Fimg_shift.data;
                    currentTrans++;
                    
                    //                    if ((ipart == 1) && (iover_trans == 0) && (itrans == 0)) {
                    //                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                    //                        {
                    //                                double Fimgreal = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                    //                                double Fimgimag = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                    //                                std::cout << "ipart: " << ipart << " itrans: " << itrans << " transcurrent: " << trans_current << " i: " << n << " Fimg: " << Fimgreal << " " << Fimgimag << std::endl;
                    //                        }
                    //                    }
                    
                    
                }
            }
        }
    }
    clFinish(CL_ComputeQueue);
    //    if (iter > 7)
    //        std::cerr << "Image data loaded" << std::endl;
    
    cl_int cl_do_ctf_correction = (do_ctf_correction && refs_are_ctf_corrected);
    cl_int cl_do_scale_correction = do_scale_correction;
    cl_int cl_num_particles = partCount;
    cl_int cl_num_trans = usedTrans;
    cl_int cl_IS_INV = IS_INV;
    cl_double cl_paddingFactor = (double)(mymodel.PPref[exp_iclass]).padding_factor;
    
    
    //    cl_mem cl_diff2Out = clCreateImage(CL_context, CL_MEM_WRITE_ONLY, &fmt, &desc, NULL, &err);
    //    cl_mem cl_diff2Out = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * totalOrients * cl_num_particles * cl_num_trans, NULL, 0);
    cl_mem cl_diff2Tmp = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
    cl_mem cl_diff2TmpHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
    
    cl_mem cl_diff2TmpB = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
    //    cl_mem cl_diff2TmpHostB = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * totalOrients * numWGSPerProj * cl_num_trans, NULL, 0);
    
    cl_mem cl_suma2Tmp;
    cl_mem cl_suma2TmpB;
    cl_mem cl_suma2TmpHost;
    if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
        cl_suma2Tmp = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
        cl_suma2TmpB = clCreateBuffer(CL_context, CL_MEM_WRITE_ONLY, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
        cl_suma2TmpHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * maxOrients * cl_num_trans, NULL, 0);
    }
    
    
    //    if (iter > 7)
    //        std::cerr << "Result memory allocated for device" << std::endl;
    
    int numFsTrans = wgs * ceil((float)cl_num_trans * projDim.x * projDim.y / (float)wgs);
    cl_mem cl_FimgShiftAllA = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * numFsTrans, NULL, 0);
    cl_mem cl_FimgShiftAllHostA = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * numFsTrans, NULL, 0);
    cl_mem cl_FimgShiftAllB = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * numFsTrans, NULL, 0);
    cl_mem cl_FimgShiftAllHostB = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * numFsTrans, NULL, 0);

    
    //Each matrix is 3x3, as is exp_R_mic
    cl_mem cl_Ainv = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * maxOrients * 9, NULL, 0);
    cl_mem cl_exp_R_mic = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * 9, NULL, 0);
    clEnqueueWriteBuffer(CL_CopyToDeviceQueue, cl_exp_R_mic, true, 0, sizeof(cl_double) * 9, exp_R_mic.mdata, 0, NULL, NULL);
    
    cl_mem cl_eulerAnglesA = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * maxOrients * 3, NULL, 0);
    cl_mem cl_eulerAnglesHostA = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * maxOrients * 3, NULL, 0);
    cl_mem cl_eulerAnglesB = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * maxOrients * 3, NULL, 0);
    cl_mem cl_eulerAnglesHostB = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * maxOrients * 3, NULL, 0);
    
    long int *iorientClassArrayA = (long int *)malloc(sizeof(long int) * maxOrients);
    long int *iOverOrientClassArrayA = (long int *)malloc(sizeof(long int) * maxOrients);
    long int *iorientClassArrayB = (long int *)malloc(sizeof(long int) * maxOrients);
    long int *iOverOrientClassArrayB = (long int *)malloc(sizeof(long int) * maxOrients);
    

    cl_mem cl_skipCalculationA = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(char) * maxOrients * usedTrans, NULL, 0);
    cl_mem cl_skipCalculationHostA = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(char) * maxOrients * usedTrans, NULL, 0);
    cl_mem cl_skipCalculationB = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(char) * maxOrients * usedTrans, NULL, 0);
    cl_mem cl_skipCalculationHostB = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(char) * maxOrients * usedTrans, NULL, 0);
    
    
    
//    if (iter > 7) {
//        std::cerr << "Projection memory allocated" << std::endl;
//    }

    cl_event clComputeEvents[numCycles];
//    cl_event clCopyToDeviceEvents[partCount];
    
#ifdef TIMING
    timer.toc(TIMING_DIFF_PROJ);
#endif
#ifdef TIMING
    timer.tic(TIMING_DIFF_DIFF2);
#endif

#ifdef CL_PRINT_SPEED
    long int orientTranslations = 0;
#endif
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif

#pragma mark Cycle start
    for (int iCycle = 0; iCycle < numCycles; iCycle++) {
        
        int ipart = particleForCycle[iCycle];
        
        //Set totalOrients to the value for this particle
//        totalOrients = rotationsForParticle[ipart];
        totalOrients = orientsPerCycle[iCycle];
        int startOrientThisCycle = startOrientPerCycle[iCycle];
        
        numFs = wgs * ceil((float)totalOrients * projDim.x * projDim.y / (float)wgs);
        
        long int *currentIOrientClassArray = iCycle % 2 ? iorientClassArrayA : iorientClassArrayB;
        long int *currentIOverOrientClassArray = iCycle % 2 ? iOverOrientClassArrayA : iOverOrientClassArrayB;

        cl_mem cl_eulerAngles = iCycle % 2 ? cl_eulerAnglesA : cl_eulerAnglesB;
        cl_mem cl_eulerAnglesHost = iCycle % 2 ? cl_eulerAnglesHostA : cl_eulerAnglesHostB;
        
        cl_mem cl_skipCalculation = iCycle % 2 ? cl_skipCalculationA : cl_skipCalculationB;
        cl_mem cl_skipCalculationHost = iCycle % 2 ? cl_skipCalculationHostA : cl_skipCalculationHostB;
        
        cl_mem cl_FimgShiftAll = iCycle % 2 ? cl_FimgShiftAllA : cl_FimgShiftAllB;
        cl_mem cl_FimgShiftAllHost = iCycle % 2 ? cl_FimgShiftAllHostA : cl_FimgShiftAllHostB;
        
        cl_mem cl_local_Fctfs = iCycle % 2 ? cl_local_FctfsA : cl_local_FctfsB;
        cl_mem cl_Minvsigma2 = iCycle % 2 ? cl_Minvsigma2A : cl_Minvsigma2B;

        long int my_image_no = myImageArray[ipart];

        if (totalOrients > 0) {
            
        double *eulerAnglesHost = (double *)clEnqueueMapBuffer(CL_CopyToDeviceQueue, cl_eulerAnglesHost, true, CL_MAP_WRITE, 0, sizeof(cl_double) * maxOrients * 3, 0, 0, 0, &err);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to map buffer for transfer of euler angles" << std::endl;
            do_use_opencl = false;
        }
        char *skipCalculation = (char *)clEnqueueMapBuffer(CL_CopyToDeviceQueue, cl_skipCalculationHost, true, CL_MAP_WRITE, 0, sizeof(char) * maxOrients * usedTrans, 0, 0, 0, &err);
        memset(skipCalculation, 1, sizeof(char) * maxOrients * usedTrans);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to map buffer for calculation optimisation" << std::endl;
            do_use_opencl = false;
        }

        
        //Now fill the orientation array
        int currentOrients = 0;
        for (long int iorient = 0; iorient < nr_orients; iorient++)
        {
            
            long int iorientclass = iorientclass_offset + iorient;
            long int idir = iorient / exp_nr_psi;
            long int ipsi = iorient % exp_nr_psi;
            // Get prior for this direction and skip calculation if prior==0
            double pdf_orientation;
            if (mymodel.orientational_prior_mode == NOPRIOR)
            {
                pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
            }
            else
            {
                pdf_orientation = sampling.getPriorProbability(idir, ipsi);
            }
            
            // In the first pass, always proceed
            // In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
            // if so, proceed with projecting the reference in that direction
            bool do_proceed = (exp_ipass==0) ? true : isSignificantAnyParticleAnyTranslation(iorientclass);
            
            if (do_proceed && pdf_orientation > 0.)
            {
                // Now get the oversampled (rot, tilt, psi) triplets
                // This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
                sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_orientations);
                
                // Loop over all oversampled orientations (only a single one in the first pass)
                for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
                {
                    bool anySignificantWeightTranslations = false;
                    
                    for (long int itrans = 0, currentTrans = 0, iUsedTrans = 0; itrans < exp_nr_trans; itrans++)
                    {
                        long int ihidden = iorientclass * exp_nr_trans + itrans;
                        
                        // In the first pass, always proceed
                        // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                        bool do_proceed = (exp_ipass == 0) ? true : DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden);
                        
                        for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                        {
                            if (do_proceed)
                            {
                                
                                //                                    if (exp_ipass > 0) {
                                //                                        std::cerr << "ipart: " << ipart << " orient: " << iorient << " trans: " << itrans << " index: " << ipart * totalOrients * usedTrans + iUsedTrans * totalOrients + currentOrients << " current orient: " << currentOrients << " current trans: " << iUsedTrans << std::endl;
                                //                                    }
                                if ((currentOrients >= startOrientThisCycle) && (currentOrients < startOrientThisCycle + totalOrients)) {
                                    skipCalculation[iUsedTrans * totalOrients + currentOrients - startOrientThisCycle] = 0;
                                }
                                anySignificantWeightTranslations = true;
#ifdef CL_PRINT_SPEED
                                orientTranslations++;
#endif
                            }
                            if (translationUsedInAnyRotation[currentTrans]) {
                                iUsedTrans++;
                            }
                            
                            currentTrans++;
                            
                        } // end loop iover_trans
                    } // end loop itrans
                    if (anySignificantWeightTranslations) {
//                        if ((iter == 2) && (exp_ipass ==1) && (ipart == 1)) {
//                            std::cerr << "currentOrient: " << currentOrients << " iorient: " << iorientclass << " ioverrot: " << iover_rot << std::endl;
//                        }
                        if ((currentOrients >= startOrientThisCycle) && (currentOrients < startOrientThisCycle + totalOrients)) {
                            currentIOverOrientClassArray[currentOrients - startOrientThisCycle] = iover_rot;
                            currentIOrientClassArray[currentOrients - startOrientThisCycle] = iorientclass;
                            // Get the Euler matrix
                            memcpy((eulerAnglesHost+(currentOrients - startOrientThisCycle) * 3), oversampled_orientations[iover_rot].vdata, sizeof(cl_double) * 3);
                        }
                        
                        currentOrients++;
                    }
                }
            }
        }
        //        if (exp_ipass == 0) {
        //            CL_Proj_allocated_orients = currentOrients;
        //            CL_Proj_alloc_iter = iter;
        //            std::cerr << "Projections regenerated" << std::endl;
        //        }
        clEnqueueUnmapMemObject(CL_CopyToDeviceQueue, cl_eulerAnglesHost, eulerAnglesHost, 0, NULL, NULL);
        clEnqueueCopyBuffer(CL_CopyToDeviceQueue, cl_eulerAnglesHost, cl_eulerAngles, 0, 0, sizeof(cl_double) * totalOrients * 3, 0, NULL, NULL);
        
        clEnqueueUnmapMemObject(CL_CopyToDeviceQueue, cl_skipCalculationHost, skipCalculation, 0, NULL, NULL);
        clEnqueueCopyBuffer(CL_CopyToDeviceQueue, cl_skipCalculationHost, cl_skipCalculation, 0, 0, sizeof(char) * totalOrients * usedTrans, 0, NULL, NULL);

        double *ipartCTF = exp_local_Fctfs[my_image_no].data;
        MultidimArray<double > MinSigma = exp_local_Minvsigma2s[my_image_no] * 0.5;
        double *minsigv2 = MinSigma.data;
        
        clEnqueueWriteBuffer(CL_CopyToDeviceQueue, cl_local_Fctfs, true, 0, sizeof(cl_double) * projDim.x * projDim.y, ipartCTF, 0, NULL, NULL);
        clEnqueueWriteBuffer(CL_CopyToDeviceQueue, cl_Minvsigma2, true, 0, sizeof(cl_double) * projDim.x * projDim.y, minsigv2, 0, NULL, NULL);

        cl_double2 *FimgShiftAllHost = (cl_double2 *)clEnqueueMapBuffer(CL_CopyToDeviceQueue, cl_FimgShiftAllHost, true, CL_MAP_WRITE, 0, sizeof(cl_double2) * numFsTrans, 0, 0, 0, &err);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to map buffer for transfer of Fimg shifts" << std::endl;
            do_use_opencl = false;
        }
        memset(FimgShiftAllHost, 0, sizeof(cl_double2) * numFsTrans);
        for (int itrans = 0; itrans < cl_num_trans; itrans++) {
            long int ishift = FimgShiftArray[ipart * cl_num_trans + itrans];
            Fimg_shift = exp_local_Fimgs_shifted[ishift];
            memcpy(FimgShiftAllHost + itrans * projDim.x * projDim.y, Fimg_shift.data, sizeof(cl_double2) * projDim.x * projDim.y);
        }
        
        clEnqueueUnmapMemObject(CL_CopyToDeviceQueue, cl_FimgShiftAllHost, FimgShiftAllHost, 0, NULL, NULL);
        clEnqueueCopyBuffer(CL_CopyToDeviceQueue, cl_FimgShiftAllHost, cl_FimgShiftAll, 0, 0, sizeof(cl_double2) * numFsTrans, 0, NULL, NULL);
        
                clFinish(CL_CopyToDeviceQueue);
        }
        //    if (iter > 7)
        //        std::cerr << "Projection memory loaded" << std::endl;
        //    }
        
        
        //    if (iter > 7)
        //        std::cerr << "Result memory allocated on host" << std::endl;
        //    cl_double *diff2Tmp = (cl_double *)malloc(sizeof(cl_double) * totalOrients * numWGSPerProj);
        //    cl_double *suma2Tmp = (cl_double *)malloc(sizeof(cl_double) * totalOrients * numWGSPerProj);
        
        
        //    gettimeofday(&endTV, NULL);
        //    unsigned long fstartMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
        //    unsigned long fendMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
        
        //    double felapsedTime = (double)(fendMicros - fstartMicros) / 1000000;
        //    std::cerr << "Prep time: " << felapsedTime << " seconds" << std::endl;
        
        
        
        
//        clWaitForEvents(1, &clCopyToDeviceEvents[ipart]);
        if (totalOrients > 0)
            blankCLMemory(CL_iter_Proj, 2 * numFs);
        
        cl_uint cl_limit = totalOrients;
        
        err = clSetKernelArg(CL_generateOrientationMatrix, 0, sizeof(cl_mem), &cl_eulerAngles);
        err |= clSetKernelArg(CL_generateOrientationMatrix, 1, sizeof(cl_mem), &cl_exp_R_mic);
        err |= clSetKernelArg(CL_generateOrientationMatrix, 2, sizeof(cl_int), &cl_IS_INV);
        err |= clSetKernelArg(CL_generateOrientationMatrix, 3, sizeof(cl_double), &cl_paddingFactor);
        err |= clSetKernelArg(CL_generateOrientationMatrix, 4, sizeof(cl_mem), &cl_Ainv);
        err |= clSetKernelArg(CL_generateOrientationMatrix, 5, sizeof(cl_uint), &cl_limit);
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to set OpenCL kernel arguments - generate matrix" << std::endl;
        }
        
        size_t global[1], local[1];
        global[0] = wgs * ceil((float)totalOrients / (float)wgs);
        local[0] = wgs;
        
        if (totalOrients >0)
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_generateOrientationMatrix, 1, NULL, global, local, 0, NULL, NULL);
        if (err)
        {
            std::cerr << "Error: Failed to execute OpenCL kernel - GSD: matrix generation, err: " << err << std::endl;
        }
        //    clFinish(CL_queue);
        
        cl_limit = totalOrients * projDim.x * projDim.y;
        cl_uint4 modelDim; modelDim.x = XSIZE(modelData); modelDim.y = YSIZE(modelData); modelDim.z = ZSIZE(modelData);
        cl_int model_r_max = (mymodel.PPref[exp_iclass]).r_max;
        
        if (mymodel.ref_dim == 2) {
            err = clSetKernelArg(CL_calcModelProjection2D_TR, 0, sizeof(cl_uint4), &modelDim);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 1, sizeof(cl_uint), &model_r_max);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 2, sizeof(cl_uint2), &projDim);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 3, sizeof(cl_mem), &cl_Ainv);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 4, sizeof(cl_mem), &cl_modelData);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 5, sizeof(cl_mem), &CL_iter_Proj);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 6, sizeof(cl_uint), &totalOrients);
            err |= clSetKernelArg(CL_calcModelProjection2D_TR, 7, sizeof(cl_uint), &cl_limit);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments - model projection" << std::endl;
            }
            
            //            std::cerr << "Finished loading CL data - model projection" << std::endl;
            global[0] = numFs;
            local[0] = wgs;
            
            if (totalOrients >0)
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcModelProjection2D_TR, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - GSD: model projection, err: " << err << std::endl;
            }
        } else if (mymodel.ref_dim == 3) {
            err = clSetKernelArg(CL_calculateProjections_TR, 0, sizeof(cl_uint4), &modelDim);
            err |= clSetKernelArg(CL_calculateProjections_TR, 1, sizeof(cl_uint), &model_r_max);
            err |= clSetKernelArg(CL_calculateProjections_TR, 2, sizeof(cl_uint2), &projDim);
            err |= clSetKernelArg(CL_calculateProjections_TR, 3, sizeof(cl_mem), &cl_Ainv);
            err |= clSetKernelArg(CL_calculateProjections_TR, 4, sizeof(cl_mem), &cl_modelData);
            err |= clSetKernelArg(CL_calculateProjections_TR, 5, sizeof(cl_mem), &CL_iter_Proj);
            err |= clSetKernelArg(CL_calculateProjections_TR, 6, sizeof(cl_uint), &totalOrients);
            err |= clSetKernelArg(CL_calculateProjections_TR, 7, sizeof(cl_uint), &cl_limit);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments - model projection" << std::endl;
            }
            
            //            std::cerr << "Finished loading CL data - model projection" << std::endl;
            global[0] = numFs;
            local[0] = wgs;
            
            if (totalOrients >0)
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateProjections_TR, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - GSD: model projection, err: " << err << std::endl;
            }
        }
        //    clFinish(CL_queue);
        
        //    if (iter > 7) {
        //        std::cerr << "Finished calculating projections" << std::endl;
        //    }
        
        cl_double scaleVal = scales[ipart];
        cl_limit = totalOrients * projDim.x * projDim.y;
        cl_uint cl_fSize = projDim.x * projDim.y;
        
//        std::cerr << "Fsize: " << cl_fSize << " x: " << projDim.x << " y: " << projDim.y << std::endl;
        
        err = 0;
        //            if (exp_ipass == 0) {
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 0, sizeof(cl_mem), &CL_iter_Proj);
        //            } else {
        //                err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 0, sizeof(cl_mem), &cl_Proj_ipass1);
        //            }
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 1, sizeof(cl_int), &cl_do_ctf_correction);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 2, sizeof(cl_mem), &cl_local_Fctfs);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 3, sizeof(cl_int), &cl_do_scale_correction);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 4, sizeof(cl_double), &scaleVal);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 5, sizeof(cl_mem), &cl_FrefCTF);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 6, sizeof(cl_uint), &cl_fSize);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 7, sizeof(cl_uint), &totalOrients);
        err |= clSetKernelArg(CL_ctfAndScaleDataPoint_TR, 8, sizeof(cl_uint), &cl_limit);
        
        if (err != CL_SUCCESS) {
            std::cerr << "Error: Failed to set OpenCL kernel arguments" << std::endl;
        }
        
        //                std::cerr << "Finished loading CL data" << std::endl;
        //            int wgs = 1;
//        size_t global[1], local[1];
        //    global[0] = wgs * ceil((float)limit / (float)wgs);;
        global[0] = numFs;
        local[0] = wgs;
        
        if (totalOrients >0)
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_ctfAndScaleDataPoint_TR, 1, NULL, global, local, 0, NULL, NULL);
        if (err)
        {
            std::cerr << "Error: Failed to execute OpenCL kernel - GSD: ctf and scale, err: " << err << std::endl;
        }
        //            clFinish(CL_queue);
        //            std::cerr << "Finished ctf and Scale" << std::endl;

        
        if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
            //    if (false) {
            //            std::cerr << "Starting ctfs" << std::endl;
//            long int my_image_no = myImageArray[ipart];
//            double *ipartCTF = exp_local_Fctfs[my_image_no].data;
            double sqrtXi2 = exp_local_sqrtXi2[my_image_no];
            
//            clEnqueueWriteBuffer(CL_ComputeQueue, cl_local_Fctfs, true, 0, sizeof(cl_double) * projDim.x * projDim.y, ipartCTF, 0, NULL, NULL);
            //            clFinish(CL_queue);
            //                        std::cerr << "Loaded ctf" << std::endl;
            
            
            cl_mem cl_currentDiff2 = iCycle % 2 ? cl_diff2Tmp : cl_diff2TmpB;
            cl_mem cl_currentSuma2 = iCycle % 2 ? cl_suma2Tmp : cl_suma2TmpB;
            
            cl_limit = totalOrients * cl_num_trans;
            
            err = 0;
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 0, sizeof(cl_mem), &cl_FrefCTF);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 1, sizeof(cl_mem), &cl_FimgShiftAll);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 2, sizeof(cl_mem), &cl_currentDiff2);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 3, sizeof(cl_mem), &cl_currentSuma2);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 4, sizeof(cl_int), &cl_fSize);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 5, sizeof(cl_uint), &totalOrients);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 6, sizeof(cl_mem), &cl_skipCalculation);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 7, sizeof(cl_int), &ipart);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2CC_TR, 8, sizeof(cl_uint), &cl_limit);
            
            
            global[0] = wgs * ceil((float)totalOrients * projDim.x * projDim.y * cl_num_trans / (float)wgs);
            
            if (totalOrients >0)
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateDiff2AndSuma2CC_TR, 1, NULL, global, local, 0, NULL, &clComputeEvents[iCycle]);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - GSD: diff2 and suma2, err: " << err << std::endl;
            }
            //                            clFinish(CL_queue);
            
            //            std::cerr << "Finished ipart: " << ipart << std::endl;
            
            if (iCycle > 0) {
                //Asynchronously transfer and process results from previous run while executing kernel
                
                int totalOrients = orientsPerCycle[iCycle - 1];
                int ipart = particleForCycle[iCycle - 1];

                if (totalOrients > 0) {
                    clWaitForEvents(1, &clComputeEvents[iCycle - 1]);
                    
                    cl_mem cl_previousDiff2 = (iCycle - 1) % 2 ? cl_diff2Tmp : cl_diff2TmpB;
                    cl_mem cl_previousSuma2 = (iCycle - 1) % 2 ? cl_suma2Tmp : cl_suma2TmpB;
                    
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousDiff2, cl_diff2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousSuma2, cl_suma2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    
                    double *diff2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_diff2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2" << std::endl;
                        do_use_opencl = false;
                    }
                    double *suma2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_suma2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of suma2" << std::endl;
                        do_use_opencl = false;
                    }
                    cl_mem cl_skipCalculationHost = (iCycle - 1) % 2 ? cl_skipCalculationHostA : cl_skipCalculationHostB;
                    char *skipCalculation = (char *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_skipCalculationHost, true, CL_MAP_READ, 0, sizeof(char) * maxOrients * usedTrans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2 optimisations" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    clFinish(CL_CopyFromDeviceQueue);
                    
                    long int my_image_no = myImageArray[ipart];
                    double sqrtXi2 = exp_local_sqrtXi2[my_image_no];
                    
                    //                clEnqueueReadBuffer(CL_queue, cl_diff2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, diff2Tmp, 0, NULL, NULL);
                    //                clEnqueueReadBuffer(CL_queue, cl_suma2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, suma2Tmp, 0, NULL, NULL);
                    
                    for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                        for (int iorient = 0; iorient < totalOrients; iorient++) {
                            if (skipCalculation[itrans*totalOrients + iorient] == 0) {
                                
                                double d2Sum = 0, s2Sum = 0;
                                //                        for (int iWG = 0; iWG < numWGSPerProj; iWG++) {
                                //                            d2Sum += diff2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                            s2Sum += suma2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                        if ((itrans == 0) && (iorient == 0)) {
                                //                            std::cerr << "itrans: " << itrans << " iorient: " << iorient << " iwg: " << iWG << " d2Sum: " << d2Sum << " s2Sum: " << s2Sum << std::endl;
                                //                        }
                                //                        }
                                d2Sum += diff2Tmp[itrans*totalOrients + iorient];
                                s2Sum += suma2Tmp[itrans*totalOrients + iorient];
                                
                                d2Sum *= -1.0;
                                if (s2Sum == 0.0) {
                                    d2Sum = 0;
                                } else {
                                    d2Sum /= sqrt(s2Sum) * sqrtXi2;
                                }
                                
                                
                                long int *iorientClassArray = (iCycle - 1) % 2 ? iorientClassArrayA : iorientClassArrayB;
                                long int *iOverOrientClassArray = (iCycle - 1) % 2 ? iOverOrientClassArrayA : iOverOrientClassArrayB;
                                
                                long int over_rot = iOverOrientClassArray[iorient];
                                long int iorientclass = iorientClassArray[iorient];
                                long int trans = iTransClassArray[itrans];
                                long int over_trans = iOverTransClassArray[itrans];
                                
                                long int ihidden = iorientclass * exp_nr_trans + trans;
                                long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                     over_rot, over_trans);
                                
                                //                    if (exp_ipass == 1) {
                                //                                                std::cerr << "itrans: " << itrans << " iorient: " << iorient << " ritrans: " << trans << " riOverTrans: " << over_trans << " riOverRot: " << over_rot << " ipart: " << ipart - 1 << " iorientclass: " << iorientclass << " d2: " << d2Sum << std::endl;
                                //                    }
                                
                                if (std::isnan(d2Sum) && do_use_opencl) {
                                    do_use_opencl = false;
                                    std::cerr << "Issue with GPU - disabling opencl and repeating calculations..." << std::endl;
                                }
                                
                                if (exp_iseries == 0)
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = d2Sum;
                                else
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += d2Sum;
                                
                                //                        if ((ipart - 1 == 0) && (iorient == 0) && (itrans == 0)) {
                                //                            std::cerr << "D2: " << diff2Tmp[itrans*totalOrients + iorient] << " S2: " << s2Sum << " Final D2: " << d2Sum << " Stored: " << DIRECT_A2D_ELEM(exp_Mweight, ipart - 1, ihidden_over) << std::endl;
                                //                        }
                                
                                d2Sum = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                
                                
                                if ((lastImageInSeries[ipart]) && d2Sum < exp_min_diff2[ipart]) {
                                    exp_min_diff2[ipart] = d2Sum;
                                }
                            }
                        }
                    }
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_skipCalculationHost, skipCalculation, 0, NULL, NULL);

                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_diff2TmpHost, diff2Tmp, 0, NULL, NULL);
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_suma2TmpHost, suma2Tmp, 0, NULL, NULL);
                }
            }
            if (iCycle == numCycles - 1) {
                //Asynchronously transfer and process results from previous run while executing kernel
                
                if (totalOrients > 0) {
                    clWaitForEvents(1, &clComputeEvents[iCycle]);
                    
                    cl_mem cl_previousDiff2 = (iCycle) % 2 ? cl_diff2Tmp : cl_diff2TmpB;
                    cl_mem cl_previousSuma2 = (iCycle) % 2 ? cl_suma2Tmp : cl_suma2TmpB;
                    
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousDiff2, cl_diff2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousSuma2, cl_suma2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    
                    double *diff2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_diff2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2" << std::endl;
                        do_use_opencl = false;
                    }
                    double *suma2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_suma2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of suma2" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    cl_mem cl_skipCalculationHost = (iCycle) % 2 ? cl_skipCalculationHostA : cl_skipCalculationHostB;
                    char *skipCalculation = (char *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_skipCalculationHost, true, CL_MAP_READ, 0, sizeof(char) * maxOrients * usedTrans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2 optimisations" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    clFinish(CL_CopyFromDeviceQueue);
                    
                    
                    //                clEnqueueReadBuffer(CL_queue, cl_diff2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, diff2Tmp, 0, NULL, NULL);
                    //                clEnqueueReadBuffer(CL_queue, cl_suma2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, suma2Tmp, 0, NULL, NULL);
                    
                    for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                        for (int iorient = 0; iorient < totalOrients; iorient++) {
                            if (skipCalculation[itrans*totalOrients + iorient] == 0) {
                                
                                double d2Sum = 0, s2Sum = 0;
                                //                        for (int iWG = 0; iWG < numWGSPerProj; iWG++) {
                                //                            d2Sum += diff2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                            s2Sum += suma2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                        if ((itrans == 0) && (iorient == 0)) {
                                //                            std::cerr << "itrans: " << itrans << " iorient: " << iorient << " iwg: " << iWG << " d2Sum: " << d2Sum << " s2Sum: " << s2Sum << std::endl;
                                //                        }
                                //                        }
                                d2Sum += diff2Tmp[itrans*totalOrients + iorient];
                                s2Sum += suma2Tmp[itrans*totalOrients + iorient];
                                
                                d2Sum *= -1.0;
                                if (s2Sum == 0.0) {
                                    d2Sum = 0;
                                } else {
                                    d2Sum /= sqrt(s2Sum) * sqrtXi2;
                                }
                                
                                long int *iorientClassArray = (iCycle) % 2 ? iorientClassArrayA : iorientClassArrayB;
                                long int *iOverOrientClassArray = (iCycle) % 2 ? iOverOrientClassArrayA : iOverOrientClassArrayB;
                                
                                long int over_rot = iOverOrientClassArray[iorient];
                                long int iorientclass = iorientClassArray[iorient];
                                long int trans = iTransClassArray[itrans];
                                long int over_trans = iOverTransClassArray[itrans];
                                
                                long int ihidden = iorientclass * exp_nr_trans + trans;
                                long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                     over_rot, over_trans);
                                
                                //                    if (exp_ipass == 1) {
                                //                        std::cerr << "itrans: " << itrans << " iorient: " << iorient << " ritrans: " << ritrans << " riOverTrans: " << riOvertrans << " riOverRot: " << riOverRot << " ipart: " << ipart << " iorientclass: " << iorientclass << std::endl;
                                //                    }
                                
                                if (std::isnan(d2Sum) && do_use_opencl) {
                                    do_use_opencl = false;
                                    std::cerr << "Issue with GPU - disabling opencl and repeating calculations..." << std::endl;
                                }
                                
                                if (exp_iseries == 0)
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = d2Sum;
                                else
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += d2Sum;
                                
                                d2Sum = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                
                                if ((lastImageInSeries[ipart]) && d2Sum < exp_min_diff2[ipart]) {
                                    exp_min_diff2[ipart] = d2Sum;
                                }
                            }
                        }
                    }
                    
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_diff2TmpHost, diff2Tmp, 0, NULL, NULL);
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_suma2TmpHost, suma2Tmp, 0, NULL, NULL);
                }
            }
        } else {
            //        double elapsedTimeInner =0;
            //            clFinish(CL_queue);
            //            std::cerr << "Loaded ctf and Minsigmav2" << std::endl;
            
            cl_mem cl_currentDiff2 = iCycle % 2 ? cl_diff2Tmp : cl_diff2TmpB;
            cl_limit = totalOrients * cl_num_trans;
            
            err = 0;
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 0, sizeof(cl_mem), &cl_FrefCTF);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 1, sizeof(cl_mem), &cl_FimgShiftAll);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 2, sizeof(cl_mem), &cl_Minvsigma2);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 3, sizeof(cl_mem), &cl_currentDiff2);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 4, sizeof(cl_int), &cl_fSize);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 5, sizeof(cl_uint), &totalOrients);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 6, sizeof(cl_mem), &cl_skipCalculation);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 7, sizeof(cl_int), &ipart);
            err |= clSetKernelArg(CL_calculateDiff2AndSuma2_TR, 8, sizeof(cl_uint), &cl_limit);
            
            
            global[0] = wgs * ceil((float)totalOrients * projDim.x * projDim.y * cl_num_trans / (float)wgs);
            
            if (totalOrients >0)
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateDiff2AndSuma2_TR, 1, NULL, global, local, 0, NULL, &clComputeEvents[iCycle]);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - GSD: diff2, err: " << err << std::endl;
            }
            //                clFinish(CL_queue);
            
            
            //            gettimeofday(&startTV, NULL);
            if (iCycle > 0) {
                //Asynchronously transfer and process results from previous run while executing kernel

                int totalOrients = orientsPerCycle[iCycle - 1];
                int ipart = particleForCycle[iCycle - 1];

                if (totalOrients > 0) {
                    
                    clWaitForEvents(1, &clComputeEvents[iCycle - 1]);
                    
                    cl_mem cl_previousDiff2 = (iCycle - 1) % 2 ? cl_diff2Tmp : cl_diff2TmpB;
                    
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousDiff2, cl_diff2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    
                    double *diff2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_diff2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    cl_mem cl_skipCalculationHost = (iCycle - 1) % 2 ? cl_skipCalculationHostA : cl_skipCalculationHostB;
                    char *skipCalculation = (char *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_skipCalculationHost, true, CL_MAP_READ, 0, sizeof(char) * maxOrients * usedTrans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2 optimisations" << std::endl;
                        do_use_opencl = false;
                    }
                    clFinish(CL_CopyFromDeviceQueue);
                    
                    long int my_image_no = myImageArray[ipart];
                    
                    //                clEnqueueReadBuffer(CL_queue, cl_diff2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, diff2Tmp, 0, NULL, NULL);
                    //                clEnqueueReadBuffer(CL_queue, cl_suma2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, suma2Tmp, 0, NULL, NULL);
                    
                    for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                        for (int iorient = 0; iorient < totalOrients; iorient++) {
                            
                            if (skipCalculation[itrans*totalOrients + iorient] == 0) {
                                
                                double d2Sum = diff2Tmp[itrans*totalOrients + iorient];
                                //                        for (int iWG = 0; iWG < numWGSPerProj; iWG++) {
                                //                            d2Sum += diff2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                        }
                                d2Sum += exp_highres_Xi2_imgs[my_image_no] / 2.;
                                //                        diff2Return[iorient * totalTrans + ipart * cl_num_trans + itrans] = d2Sum;
                                
                                long int *iorientClassArray = (iCycle - 1) % 2 ? iorientClassArrayA : iorientClassArrayB;
                                long int *iOverOrientClassArray = (iCycle - 1) % 2 ? iOverOrientClassArrayA : iOverOrientClassArrayB;
                                
                                long int over_rot = iOverOrientClassArray[iorient];
                                long int iorientclass = iorientClassArray[iorient];
                                long int trans = iTransClassArray[itrans];
                                long int over_trans = iOverTransClassArray[itrans];
                                
                                long int ihidden = iorientclass * exp_nr_trans + trans;
                                long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                     over_rot, over_trans);
                                
                                //                    if (exp_ipass == 1) {
                                //                    }
                                
                                if (std::isnan(d2Sum) && do_use_opencl) {
                                    do_use_opencl = false;
                                    std::cerr << "Issue with GPU - disabling opencl and repeating calculations..." << std::endl;
                                }
                                
                                if (exp_iseries == 0)
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = d2Sum;
                                else
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += d2Sum;
                                
                                d2Sum = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                
                                //                        if ((exp_ipass == 1) && (ipart - 1 == 1)) {
                                
                                //                            std::cerr << "itrans: " << itrans << " iorient: " << iorient << " orient: " << iorientclass << " over_rot: " << over_rot << " trans: " << trans << " over_trans: " << over_trans << " ipart: " << ipart - 1 << " cl_result: " << d2Sum << std::endl;
                                //                        }
                                
                                if ((lastImageInSeries[ipart]) && d2Sum < exp_min_diff2[ipart]) {
                                    exp_min_diff2[ipart] = d2Sum;
                                }
                            }
                        }
                    }
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_skipCalculationHost, skipCalculation, 0, NULL, NULL);

                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_diff2TmpHost, diff2Tmp, 0, NULL, NULL);
                }
            }
            
            if (iCycle == numCycles - 1) {
                
                if (totalOrients > 0) {
                    clWaitForEvents(1, &clComputeEvents[iCycle]);
                    
                    cl_mem cl_previousDiff2 = (iCycle) % 2 ? cl_diff2Tmp : cl_diff2TmpB;
                    
                    clEnqueueCopyBuffer(CL_CopyFromDeviceQueue, cl_previousDiff2, cl_diff2TmpHost, 0, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, NULL, NULL);
                    
                    double *diff2Tmp = (double *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_diff2TmpHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    cl_mem cl_skipCalculationHost = iCycle % 2 ? cl_skipCalculationHostA : cl_skipCalculationHostB;
                    char *skipCalculation = (char *)clEnqueueMapBuffer(CL_CopyFromDeviceQueue, cl_skipCalculationHost, true, CL_MAP_READ, 0, sizeof(char) * maxOrients * usedTrans, 0, 0, 0, &err);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error: Failed to map buffer for transfer of diff2 optimisations" << std::endl;
                        do_use_opencl = false;
                    }
                    
                    clFinish(CL_CopyFromDeviceQueue);
                    
                    
                    //                clEnqueueReadBuffer(CL_queue, cl_diff2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, diff2Tmp, 0, NULL, NULL);
                    //                clEnqueueReadBuffer(CL_queue, cl_suma2Tmp, true, 0, sizeof(cl_double) * totalOrients * numWGSPerProj, suma2Tmp, 0, NULL, NULL);
                    
                    for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                        for (int iorient = 0; iorient < totalOrients; iorient++) {
                            
                            if (skipCalculation[itrans*totalOrients + iorient] == 0) {
                                
                                double d2Sum = diff2Tmp[itrans*totalOrients + iorient];
                                //                        for (int iWG = 0; iWG < numWGSPerProj; iWG++) {
                                //                            d2Sum += diff2Tmp[itrans*totalOrients*numWGSPerProj + iorient*numWGSPerProj + iWG];
                                //                        }
                                d2Sum += exp_highres_Xi2_imgs[my_image_no] / 2.;
                                //                        diff2Return[iorient * totalTrans + ipart * cl_num_trans + itrans] = d2Sum;
                                
                                long int *iorientClassArray = (iCycle) % 2 ? iorientClassArrayA : iorientClassArrayB;
                                long int *iOverOrientClassArray = (iCycle) % 2 ? iOverOrientClassArrayA : iOverOrientClassArrayB;
                                
                                long int over_rot = iOverOrientClassArray[iorient];
                                long int iorientclass = iorientClassArray[iorient];
                                long int trans = iTransClassArray[itrans];
                                long int over_trans = iOverTransClassArray[itrans];
                                
                                long int ihidden = iorientclass * exp_nr_trans + trans;
                                long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                     over_rot, over_trans);
                                
                                //                    if (exp_ipass == 1) {
                                //                        std::cerr << "itrans: " << itrans << " iorient: " << iorient << " ritrans: " << ritrans << " riOverTrans: " << riOvertrans << " riOverRot: " << riOverRot << " ipart: " << ipart << " iorientclass: " << iorientclass << std::endl;
                                //                    }
                                
                                if (std::isnan(d2Sum) && do_use_opencl) {
                                    do_use_opencl = false;
                                    std::cerr << "Issue with GPU - disabling opencl and repeating calculations..." << std::endl;
                                }
                                
                                if (exp_iseries == 0)
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = d2Sum;
                                else
                                    DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += d2Sum;
                                
                                d2Sum = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                
                                if ((lastImageInSeries[ipart]) && d2Sum < exp_min_diff2[ipart]) {
                                    exp_min_diff2[ipart] = d2Sum;
                                }
                            }
                        }
                    }
                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_skipCalculationHost, skipCalculation, 0, NULL, NULL);

                    clEnqueueUnmapMemObject(CL_CopyFromDeviceQueue, cl_diff2TmpHost, diff2Tmp, 0, NULL, NULL);
                }
            }
            //            gettimeofday(&endTV, NULL);
            //            unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            
            //            elapsedTimeInner += (double)(endMicros - startMicros) / 1000000;
        }
        
        //        std::cerr << "Inner cpu time: " << elapsedTimeInner << std::endl;
    }//end particle loop
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedTime = (double)(endMicros - startMicros) / 1000000;
    
    unsigned long startPrepMicros = 1000000*prepStartTV.tv_sec + prepStartTV.tv_usec;
    
    double gpuElapsedTime = (double)(endMicros - startPrepMicros) / 1000000;

    totalOrients = 0;
    for (int iCycle = 0; iCycle < numCycles; iCycle++) {
        totalOrients += orientsPerCycle[iCycle];
    }
    
    unsigned long ops;
    if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
        ops = orientTranslations * projDim.x * projDim.y * 8 + totalOrients * (projDim.x * projDim.y * (4 + 57) + 168);
    } else {
        ops = orientTranslations * projDim.x * projDim.y * 8 + totalOrients * (projDim.x * projDim.y * (4 + 57) + 168);
    }
    std::cerr << "GSD Time: " << elapsedTime << "s, ops: " << ops << " orients: " << totalOrients << " trans: " << totalTrans << " GFLOPS: " << (double)ops / elapsedTime / 1e9 << std::endl;
#endif
    
#ifdef TIMING
    timer.toc(TIMING_DIFF_DIFF2);
#endif
    
    //    std::cerr << "Finished executing kernel!" << std::endl;
    
    //    cl_double2 *projOut = (cl_double2 *)malloc(sizeof(cl_double2) * projDim.x * projDim.y);
    //    clEnqueueReadBuffer(CL_queue, cl_projOut, true, 0, sizeof(cl_double2) * projDim.x * projDim.y, projOut, 0, NULL, NULL);
    
#ifdef TIMING
    timer.tic(TIMING_DIFF_PROJ);
#endif
    //    cl_double2 *diff2Return = (cl_double2 *)malloc(sizeof(cl_double2) * totalOrients * cl_num_particles * cl_num_trans);
    //    origin[0] = 0; origin[1] = 0; origin[2] = 0;
    //    region[0] = cl_num_particles * cl_num_trans; region[1] = totalOrients; region[2] = 1;
    //    clEnqueueReadImage(CL_queue, cl_diff2Out, true, origin, region, 0, 0, diff2Return, 0, NULL, NULL);
    
    //    clEnqueueReadBuffer(CL_queue, cl_diff2Out, true, 0, sizeof(cl_double) * totalOrients * cl_num_particles * cl_num_trans, diff2Return, 0, NULL, NULL);
    //    clFinish(CL_queue);
    
    //    clReleaseMemObject(cl_diff2Out);
    //    if (exp_ipass > 0) {
    //        clReleaseMemObject(cl_Proj_ipass1);
    //    }
//    gettimeofday(&startTV, NULL);
    
    clFinish(CL_CopyFromDeviceQueue);
    
    for (int iCycle = 0; iCycle < numCycles; iCycle++) {
        if (orientsPerCycle[iCycle] > 0) {
            clReleaseEvent(clComputeEvents[iCycle]);
        }
    }
    
    clReleaseMemObject(CL_iter_Proj);
    clReleaseMemObject(cl_modelData);
    clReleaseMemObject(cl_Ainv);
    clReleaseMemObject(cl_eulerAnglesA);
    clReleaseMemObject(cl_eulerAnglesHostA);
    clReleaseMemObject(cl_eulerAnglesB);
    clReleaseMemObject(cl_eulerAnglesHostB);
    clReleaseMemObject(cl_exp_R_mic);
    //    clReleaseMemObject(cl_projOut);
    clReleaseMemObject(cl_local_FctfsA);
    clReleaseMemObject(cl_local_FctfsB);
    //    clReleaseMemObject(cl_scale);
//    clReleaseMemObject(cl_FimgShift);
    clReleaseMemObject(cl_FimgShiftAllA);
    clReleaseMemObject(cl_FimgShiftAllHostA);
    clReleaseMemObject(cl_FimgShiftAllB);
    clReleaseMemObject(cl_FimgShiftAllHostB);
    clReleaseMemObject(cl_FrefCTF);
    clReleaseMemObject(cl_diff2Tmp);
    clReleaseMemObject(cl_diff2TmpHost);
    clReleaseMemObject(cl_diff2TmpB);
//    clReleaseMemObject(cl_diff2TmpHostB);
    if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
        clReleaseMemObject(cl_suma2TmpHost);
        clReleaseMemObject(cl_suma2Tmp);
        clReleaseMemObject(cl_suma2TmpB);
    }
    //    if ((iter == 1 && do_firstiter_cc) || do_always_cc) {
    //    if (false) {
    //        clReleaseMemObject(cl_sqrtXi2);
    //    } else {
    //        clReleaseMemObject(cl_projCTF);
    clReleaseMemObject(cl_Minvsigma2A);
    clReleaseMemObject(cl_Minvsigma2B);
    clReleaseMemObject(cl_skipCalculationA);
    clReleaseMemObject(cl_skipCalculationHostA);
    clReleaseMemObject(cl_skipCalculationB);
    clReleaseMemObject(cl_skipCalculationHostB);
    //        clReleaseMemObject(cl_exp_highres_Xi2_imgs);
    //    }
    free(scales);
    free(myImageArray);
    free(FimgShiftArray);
    free(iorientClassArrayA);
    free(iOverOrientClassArrayA);
    free(iorientClassArrayB);
    free(iOverOrientClassArrayB);
    free(iTransClassArray);
    free(iOverTransClassArray);
    free(lastImageInSeries);
//    free(diff2Tmp);
//    free(suma2Tmp);
    
    //    for (int i = 0; i < results_to_check; i++) {
    //        std::cout << "i: " << i << " Diff2: " << diff2Return[i] << std::endl;
    //    }

#ifdef TIMING
    timer.toc(TIMING_DIFF_PROJ);
#endif
//    gettimeofday(&endTV, NULL);
//    unsigned long estartMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
//    unsigned long eendMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
//    double eelapsedTime = (double)(eendMicros - estartMicros) / 1000000;
//    std::cerr << "Clean-up time: " << eelapsedTime << " seconds" << std::endl;

//    std::cerr << "Finished saving results" << std::endl;
    
    //Do the one CTF and scale
    // Apply CTF to reference projection
/*    if (do_ctf_correction && refs_are_ctf_corrected)
    {
        
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
        {
            DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[test_my_image_no], n);
            //                        std::cerr << "Out: " << DIRECT_MULTIDIM_ELEM(Frefctf, n).real << " " << DIRECT_MULTIDIM_ELEM(Frefctf, n).imag << " input: " << DIRECT_MULTIDIM_ELEM(Fref, n).real << " " << DIRECT_MULTIDIM_ELEM(Fref, n).imag << " ctf: " << DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[test_my_image_no], n) << std::endl;
        }
    }
    else
        Frefctf = Fref;
    
    if (do_scale_correction)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
        {
            DIRECT_MULTIDIM_ELEM(Frefctf, n) *= mymodel.scale_correction[cl_group_id];
        }
    }*/
    
    
//    for (int k = 0; k < projDim.y; k++) {
//        for (int j = 0; j < projDim.x ; j++) {
//            std::cout << " Real CL: " << projOut[k * projDim.x + j].x << std::endl;
//            std::cout << " imag CL: " << projOut[k * projDim.x + j].y << std::endl;
//                                    std::cout << "X " << x << " Y " << y << " Real true: " << DIRECT_A2D_ELEM(Fref, x, y).real << " Real CL: " << projOut[x * projDim.y + y].x << std::endl;
//                                    std::cout << "X " << x << " Y " << y << " Imag true: " << DIRECT_A2D_ELEM(Fref, x, y).imag << " imag CL: " << projOut[x * projDim.y + y].y << std::endl;
//            if ((DIRECT_A2D_ELEM(Frefctf, k, j).real != projOut[k * projDim.x + j].x) || (DIRECT_A2D_ELEM(Frefctf, k, j).imag != projOut[k * projDim.x + j].y)) {
//                std::cout << "K " << k << " J " << j << "True: " << DIRECT_A2D_ELEM(Frefctf, k, j).real << " " << DIRECT_A2D_ELEM(Frefctf, k, j).imag << " Real diff: " << DIRECT_A2D_ELEM(Frefctf, k, j).real - projOut[k * projDim.x + j].x << " Imag diff: " << DIRECT_A2D_ELEM(Frefctf, k, j).imag - projOut[k * projDim.x + j].y << std::endl;
//            }
//        }
//    }

    //Cpu code...
//    if ((iter == 1 && do_firstiter_cc) || do_always_cc) {

//    } else {
#ifdef CL_VERIFY_ON_CPU
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif
#endif
#ifndef CL_VERIFY_ON_CPU
    if (!do_use_opencl) {
#endif
//    if (iter > 1) {
        bool gdsAnyWrong = false;
        unsigned long totalDiffs = 0;
        unsigned long totalCorrect = 0;
//#ifdef CL_PRINT_SAMPLE_VALUES
        int numSampleValuesPrinted = 0;
//#endif
        bool printedFirstOrient = false;
    for (long int iorient = 0; iorient < nr_orients; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        long int idir = iorient / exp_nr_psi;
        long int ipsi = iorient % exp_nr_psi;
        // Get prior for this direction and skip calculation if prior==0
        double pdf_orientation;
        if (mymodel.orientational_prior_mode == NOPRIOR)
        {
#ifdef DEBUG_CHECKSIZES
            if (idir >= XSIZE(mymodel.pdf_direction[exp_iclass]))
            {
                std::cerr<< "idir= "<<idir<<" XSIZE(mymodel.pdf_direction[exp_iclass])= "<< XSIZE(mymodel.pdf_direction[exp_iclass]) <<std::endl;
                REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
            }
#endif
            pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
        }
        else
        {
            pdf_orientation = sampling.getPriorProbability(idir, ipsi);
        }
        
        // In the first pass, always proceed
        // In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
        // if so, proceed with projecting the reference in that direction
        bool do_proceed = (exp_ipass==0) ? true : isSignificantAnyParticleAnyTranslation(iorientclass);
        
        if (do_proceed && pdf_orientation > 0.)
        {
            // Now get the oversampled (rot, tilt, psi) triplets
            // This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
            sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_orientations);
            
#ifdef DEBUG_CHECKSIZES
            if (exp_nr_oversampled_rot != oversampled_orientations.size())
            {
                std::cerr<< "exp_nr_oversampled_rot= "<<exp_nr_oversampled_rot<<" oversampled_orientations.size()= "<< oversampled_orientations.size() <<std::endl;
                REPORT_ERROR("exp_nr_oversampled_rot != oversampled_orientations.size()");
            }
#endif
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                
                // Get the Euler matrix
                Euler_angles2matrix(XX(oversampled_orientations[iover_rot]),
                                    YY(oversampled_orientations[iover_rot]),
                                    ZZ(oversampled_orientations[iover_rot]), A);
                
                // Take tilt-series into account
                A = (exp_R_mic * A).inv();
                
                // Project the reference map (into Fref)
#ifdef TIMING
                // Only time one thread, as I also only time one MPI process
                
                    timer.tic(TIMING_DIFF_PROJ);
#endif
                (mymodel.PPref[exp_iclass]).get2DFourierTransformD(Fref, A, IS_INV);
#ifdef TIMING
                // Only time one thread, as I also only time one MPI process
                
                    timer.toc(TIMING_DIFF_PROJ);
#endif
                
                
                
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
//                    std::cerr << "Number of particles for ori_part_id: " << ori_part_id << " is: " << mydata.ori_particles[ori_part_id].particles_id.size() << std::endl;
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        bool is_last_image_in_series = mydata.getNrImagesInSeries(part_id) == (exp_iseries + 1);
                        // Which number was this image in the combined array of exp_iseries and part_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        
#ifdef DEBUG_CHECKSIZES
                        if (my_image_no >= exp_local_Minvsigma2s.size())
                        {
                            std::cerr<< "my_image_no= "<<my_image_no<<" exp_local_Minvsigma2s.size()= "<< exp_local_Minvsigma2s.size() <<std::endl;
                            REPORT_ERROR("my_image_no >= exp_local_Minvsigma2.size()");
                        }
#endif
                        Minvsigma2 = exp_local_Minvsigma2s[my_image_no];
                        
                        // Apply CTF to reference projection
                        if (do_ctf_correction && refs_are_ctf_corrected)
                        {
                            
#ifdef DEBUG_CHECKSIZES
                            if (my_image_no >= exp_local_Fctfs.size())
                            {
                                std::cerr<< "my_image_no= "<<my_image_no<<" exp_local_Fctfs.size()= "<< exp_local_Fctfs.size() <<std::endl;
                                REPORT_ERROR("my_image_no >= exp_local_Fctfs.size()");
                            }
                            if (MULTIDIM_SIZE(Fref) != MULTIDIM_SIZE(exp_local_Fctfs[my_image_no]))
                            {
                                std::cerr<< "MULTIDIM_SIZE(Fref)= "<<MULTIDIM_SIZE(Fref)<<" MULTIDIM_SIZE()= "<< MULTIDIM_SIZE(exp_local_Fctfs[my_image_no]) <<std::endl;
                                REPORT_ERROR("MULTIDIM_SIZE(Fref) != MULTIDIM_SIZE(exp_local_Fctfs[my_image_no)");
                            }
                            
#endif
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
                            {
                                DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[my_image_no], n);
                            }
                        }
                        else
                            Frefctf = Fref;
                        
                        if (do_scale_correction)
                        {
                            int group_id = mydata.getGroupId(part_id, exp_iseries);
#ifdef DEBUG_CHECKSIZES
                            if (group_id >= mymodel.scale_correction.size())
                            {
                                std::cerr<< "group_id= "<<group_id<<" mymodel.scale_correction.size()= "<< mymodel.scale_correction.size() <<std::endl;
                                REPORT_ERROR("group_id >= mymodel.scale_correction.size()");
                            }
#endif
                            double myscale = mymodel.scale_correction[group_id];
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
                            {
                                DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
                            }
                        }
                        
//                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
//                        {
//                            if (!printedFirstOrient) {
//                                double Frefreal = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
//                                double Frefimag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
//                                std::cout << " i: " << n << " Fref: " << Frefreal << " " << Frefimag << std::endl;
//                            }
//                        }

                        
                        for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
                        {
                            long int ihidden = iorientclass * exp_nr_trans + itrans;
                            
                            // In the first pass, always proceed
                            // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                            bool do_proceed = (exp_ipass == 0) ? true : DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden);

                            if (do_proceed)
                            {
                                
                                sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations);
                                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                                {
#ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    
                                        timer.tic(TIMING_DIFF_DIFF2);
#endif
                                    // Get the shifted image
                                    long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
                                    itrans * exp_nr_oversampled_trans + iover_trans;
                                    
                                    
                                    Fimg_shift = exp_local_Fimgs_shifted[ishift];
                                    //#define DEBUG_GETALLDIFF2
                                    
                                    double diff2;
                                    if ((iter == 1 && do_firstiter_cc) || do_always_cc)
//                                    if (false)
                                    {
                                        // Do not calculate squared-differences, but signal product
                                        // Negative values because smaller is worse in this case
                                        diff2 = 0.;
                                        double suma2 = 0.;
                                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                                        {
                                            diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                            diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                            suma2 += norm(DIRECT_MULTIDIM_ELEM(Frefctf, n));
//                                            if (!printedFirstOrient && (ipart == 1)) {
//                                            if ((iorient == 0) && (ipart < 2) && (itrans == 0)) {
//                                                double Frefreal = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
//                                                double Frefimag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
//                                                double Fimgreal = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
//                                                double Fimgimag = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
//                                                std::cout << "iorient: " << iorient << " ipart: " << ipart << " itrans: " << itrans << " i: " << n << " Fref: " << Frefreal << " " << Frefimag << " Fimg: " << Fimgreal << " " << Fimgimag << "suma2: " << suma2 << " diff2: " << diff2 << std::endl;
//                                            }

                                        }
                                        printedFirstOrient = true;
                                        // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
                                        diff2 /= sqrt(suma2) * exp_local_sqrtXi2[my_image_no];
                                    }
                                    else
                                    {
                                        
                                        
                                        // Calculate the actual squared difference term of the Gaussian probability function
                                        // If current_size < mymodel.ori_size diff2 is initialised to the sum of
                                        // all |Xij|2 terms that lie between current_size and ori_size
                                        // Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
                                        diff2 = 0;
                                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                                        {
                                            double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                            double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                            diff2 += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
//                                            if (!printedFirstOrient && (ipart == 1)) {
//                                                double Frefreal = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
//                                                double Frefimag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
//                                                double Fimgreal = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
//                                                double Fimgimag = (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
//                                                double Minvsigma2n = DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
//                                                std::cout << "ipart: " << ipart << " itrans: " << itrans << " i: " << n << " Fref: " << Frefreal << " " << Frefimag << " Fimg: " << Fimgreal << " " << Fimgimag << " Minsigma2: " << Minvsigma2n << "diff2: " << diff2 << std::endl;
//                                            }

                                        }
                                        diff2 += exp_highres_Xi2_imgs[my_image_no] / 2.;
//                                        if (!printedFirstOrient && (ipart == 1)) {
//                                        printedFirstOrient = true;
//                                        }
                                        
                                    }
                                    
#ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process

                                        timer.toc(TIMING_DIFF_DIFF2);
#endif
                                    
                                    // Store all diff2 in exp_Mweight
                                    long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                         iover_rot, iover_trans);
                                    
                                    if (do_use_opencl) {
                                        double testVal = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                        if ((fabs(diff2 - testVal) > 0.0000000001f) || (std::isnan(DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over)))) {
                                            gdsAnyWrong = true;
                                            totalDiffs++;
                                            DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;
#ifdef CL_PRINT_SAMPLE_VALUES
                                            if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                                                std::cerr << "Wrong - Actual Diff2: " << diff2 << " OpenCL diff2: " << DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) << " iOrient: " << iorient << " itrans: " << itrans << " ipart: " << ipart << std::endl;
                                                numSampleValuesPrinted++;
                                            }
#endif
                                        } else {
                                            totalDiffs++;
                                            totalCorrect++;
#ifdef CL_PRINT_SAMPLE_VALUES
                                                if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                                                std::cerr << "Right - Actual Diff2: " << diff2 << " OpenCL diff2: " << DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) << " iOrient: " << iorient << " itrans: " << itrans << " ipart: " << ipart << std::endl;
                                                numSampleValuesPrinted++;
                                            }
#endif
                                        }
                                    } else {
                                        //                                    if (diff2 < -4.0e10) {
                                        //                                        std::cerr << "Likely issue..." << std::endl;
                                        //                                    }
                                        
                                        //#define DEBUG_DIFF2_ISNAN
                                        //#define DEBUG_VERBOSE
                                        
                                        if (exp_iseries == 0)
                                            DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;
                                        else
                                            DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) += diff2;
                                        
                                        // Keep track of minimum of all diff2, only for the last image in this series
                                        diff2 = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                        //std::cerr << " exp_ipass= " << exp_ipass << " exp_iclass= " << exp_iclass << " diff2= " << diff2 << std::endl;
                                        if (is_last_image_in_series && diff2 < exp_min_diff2[ipart])
                                            exp_min_diff2[ipart] = diff2;
                                    }

                                } // end loop iover_trans
                            } // end if do_proceed translations
#ifdef CL_FULL_VERIFY
                            else {
                                if (iter > 1) {
                                //for first value in exp_iseries these should be equal to init_constant
                                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                                {
                                    long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
                                                                                                         iover_rot, iover_trans);
                                    double diff2 = (double)-999.;
                                    double testVal = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                    double diff = diff2 - testVal;
                                    if ((fabs(diff) > 0.0000000001f) || (std::isnan(DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over)))) {
                                        gdsAnyWrong = true;
                                        totalDiffs++;
//#ifdef CL_PRINT_SAMPLE_VALUES
                                        if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                                            std::cerr << "Wrong - Actual Diff2: " << diff2 << " OpenCL diff2: " << DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) << " iOrient: " << iorient << " itrans: " << itrans << " ipart: " << ipart << " Difference: " << diff << " fabs(diff): " << fabs(diff) << std::endl;
                                            numSampleValuesPrinted++;
                                        }
//#endif
                                    } else {
                                        totalDiffs++;
                                        totalCorrect++;
#ifdef CL_PRINT_SAMPLE_VALUES
                                        if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                                            std::cerr << "Right - Actual Diff2: " << diff2 << " OpenCL diff2: " << DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) << " iOrient: " << iorient << " itrans: " << itrans << " ipart: " << ipart << std::endl;
                                            numSampleValuesPrinted++;
                                        }
#endif
                                    }
                                }
                                }
                            } // end if do_proceed translations
#endif
                        } // end loop itrans
                    } // end loop part_id (i)
                } // end loop ori_part_id
            }// end loop iover_rot
        } // end if do_proceed orientations
    } // end loop iorient
#ifndef CL_VERIFY_ON_CPU
    }
#endif
#ifdef CL_VERIFY_ON_CPU
    if (gdsAnyWrong == true) {
        std::cerr << "Wrong - incorrect diff2 values calculated, " << totalCorrect << "/" << totalDiffs << " = " << (float)totalCorrect / (float)totalDiffs * 100.0f << "% correct" << std::endl;
    } else {
//        std::cerr << "Right - all diff2 values correct" << std::endl;
    }
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedCPUTime = (double)(endMicros - startMicros) / 1000000;
    std::cerr << "Diff2 CPU Time: " << elapsedCPUTime << "s GPU time: " << gpuElapsedTime << " s, equivalent to: " << elapsedCPUTime / gpuElapsedTime << " cores" << std::endl;
#endif
#endif
    // Now inside a mutex set the minimum of the squared differences among all threads
    
    /*global_mutex.lock();
    for (int i = 0; i < exp_min_diff2.size(); i++)
    {
        if (thisthread_min_diff2[i] < exp_min_diff2[i])
        {
            exp_min_diff2[i] = thisthread_min_diff2[i];
        }
    }
    global_mutex.unlock();
    
    // Wait until all threads have finished
    global_barrier->wait();
    */
    
#ifdef DEBUG_OPENCL
    std::cerr << "leaving doThreadGetAllSquaredDifferences" << std::endl;
#endif
    
}


void MlOptimiser::getAllSquaredDifferences()
{

#ifdef TIMING
	if (exp_ipass == 0) timer.tic(TIMING_ESP_DIFF1);
	else timer.tic(TIMING_ESP_DIFF2);
#endif

//#define DEBUG_GETALLDIFF2
#ifdef DEBUG_GETALLDIFF2
	std::cerr << " ipass= " << exp_ipass << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
	std::cerr << " sampling.NrPsiSamplings(exp_current_oversampling)= " << sampling.NrPsiSamplings(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.NrTranslationalSamplings(exp_current_oversampling)= " << sampling.NrTranslationalSamplings(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.NrSamplingPoints(exp_current_oversampling)= " << sampling.NrSamplingPoints(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.oversamplingFactorOrientations(exp_current_oversampling)= "<<sampling.oversamplingFactorOrientations(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.oversamplingFactorTranslations(exp_current_oversampling)= "<<sampling.oversamplingFactorTranslations(exp_current_oversampling) << std::endl;
#endif

	// Initialise min_diff and exp_Mweight for this pass
	exp_Mweight.resize(exp_nr_particles, mymodel.nr_classes * sampling.NrSamplingPoints(exp_current_oversampling, false));
	exp_Mweight.initConstant(-999.);
	if (exp_ipass==0)
		exp_Mcoarse_significant.clear();

	exp_min_diff2.clear();
	exp_min_diff2.resize(exp_nr_particles);
	for (int n = 0; n < exp_nr_particles; n++)
		exp_min_diff2[n] = 99.e99;

	// Use pre-sized vectors instead of push_backs!!
	exp_local_Fimgs_shifted.clear();
	exp_local_Fimgs_shifted.resize(exp_nr_images * sampling.NrTranslationalSamplings(exp_current_oversampling));
	exp_local_Fimgs_shifted_nomask.clear();
	exp_local_Fimgs_shifted_nomask.resize(exp_nr_images * sampling.NrTranslationalSamplings(exp_current_oversampling));
	exp_local_Minvsigma2s.clear();
	exp_local_Minvsigma2s.resize(exp_nr_images);
	exp_local_Fctfs.clear();
	exp_local_Fctfs.resize(exp_nr_images);
	exp_local_sqrtXi2.clear();
	exp_local_sqrtXi2.resize(exp_nr_images);
    
	// TODO: MAKE SURE THAT ALL PARTICLES IN SomeParticles ARE FROM THE SAME AREA, SO THAT THE R_mic CAN BE RE_USED!!!

	//for (exp_iseries = 0; exp_iseries < mydata.getNrImagesInSeries(part_id); exp_iseries++)
	for (exp_iseries = 0; exp_iseries < mydata.getNrImagesInSeries((mydata.ori_particles[exp_my_first_ori_particle]).particles_id[0]); exp_iseries++)
	{
		// Get all shifted versions of the (downsized) images, their (downsized) CTFs and their inverted Sigma2 matrices
        if (do_use_opencl) {
            doOpenCLPrecalculateShiftedImagesCtfsAndInvSigma2s();
        } else {
            exp_ipart_ThreadTaskDistributor->reset(); // reset thread distribution tasks
            global_ThreadManager->run(globalThreadPrecalculateShiftedImagesCtfsAndInvSigma2s);
        }

        
		// Get micrograph transformation matrix. Note that for all pooled particles (exp_my_first_particle-exp_my_last_particle)
		// the same exp_R_mic will be used in order to re-use the reference projections
		// This is the reason why all pooled particles should come from the same micrograph
		// TODO: THAT STILL NEEDS TO BE CONFIRMED!!!! CURRENTLY NO CHECK ON SAME-PARTICLENAME IN EACH POOL!!!
		// WORKAROUND FOR NOW: just set --pool 1
		//exp_R_mic = mydata.getMicrographTransformationMatrix((mydata.ori_particles[exp_my_first_ori_particle]).particles_id[0], exp_iseries);
		int my_image_no = exp_starting_image_no[0] + exp_iseries;
		// Get micrograph transformation matrix
		exp_R_mic.resize(3,3);
		exp_R_mic(0,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_0);
		exp_R_mic(0,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_1);
		exp_R_mic(0,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_2);
		exp_R_mic(1,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_0);
		exp_R_mic(1,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_1);
		exp_R_mic(1,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_2);
		exp_R_mic(2,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_0);
		exp_R_mic(2,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_1);
		exp_R_mic(2,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_2);

		// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
		for (exp_iclass = iclass_min; exp_iclass <= iclass_max; exp_iclass++)
		{
			if (mymodel.pdf_class[exp_iclass] > 0.)
			{
                if (do_use_opencl) {
//                    struct timeval curTV;
//                    gettimeofday(&curTV, NULL);
//                    std::cerr << "Enter diff2 time: " << curTV.tv_sec << "." << curTV.tv_usec << std::endl;
                    doOpenCLGetSquaredDifferencesAllOrientations();
//                    gettimeofday(&curTV, NULL);
//                    std::cerr << "Exit diff2 time: " << curTV.tv_sec << "." << curTV.tv_usec << std::endl;
                } else {
                    exp_iorient_ThreadTaskDistributor->reset(); // reset thread distribution tasks
                    global_ThreadManager->run(globalThreadGetSquaredDifferencesAllOrientations);
                }

			} // end if mymodel.pdf_class[iclass] > 0.
        } // end loop iclass
	} // end loop iseries

#ifdef DEBUG_GETALLDIFF2b
	for (long int part_id = exp_my_first_particle, ipart = 0; part_id <= exp_my_last_particle; part_id++, ipart++)
	{
		if (exp_min_diff2[ipart] < 0.)
		{
			std::cerr << "Negative min_diff2...." << std::endl;
			std::cerr << " ipart= " << ipart << " part_id= "<<part_id<<std::endl;
			std::cerr << " do_firstiter_cc= " << do_firstiter_cc << std::endl;
			int group_id = mydata.getGroupId(part_id, 0);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
		}
	}
#endif

#ifdef TIMING
	if (exp_ipass == 0) timer.toc(TIMING_ESP_DIFF1);
	else timer.toc(TIMING_ESP_DIFF2);
#endif

}

void MlOptimiser::doThreadConvertSquaredDifferencesToWeightsAllOrientations(int thread_id)
{
#ifdef DEBUG_THREAD
    std::cerr << "entering doThreadConvertSquaredDifferencesToWeightsAllOrientations" << std::endl;
#endif


    // Store local sum of weights for this thread and then combined all threads at the end of this function inside a mutex.
	double thisthread_sumweight = 0.;

	// exp_iclass loop does not always go from 0 to nr_classes!
	long int iorientclass_offset = exp_iclass * exp_nr_rot;

	size_t first_iorient = 0, last_iorient = 0;
	while (exp_iorient_ThreadTaskDistributor->getTasks(first_iorient, last_iorient))
	{
		for (long int iorient = first_iorient; iorient <= last_iorient; iorient++)
		{

			double pdf_orientation;
			long int iorientclass = iorientclass_offset + iorient;
			long int idir = iorient / exp_nr_psi;
			long int ipsi = iorient % exp_nr_psi;

			// Get prior for this direction
			if (mymodel.orientational_prior_mode == NOPRIOR)
			{
#ifdef DEBUG_CHECKSIZES
				if (idir >= XSIZE(mymodel.pdf_direction[exp_iclass]))
				{
					std::cerr<< "idir= "<<idir<<" XSIZE(mymodel.pdf_direction[exp_iclass])= "<< XSIZE(mymodel.pdf_direction[exp_iclass]) <<std::endl;
					REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
				}
#endif
				pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
			}
			else
			{
				pdf_orientation = sampling.getPriorProbability(idir, ipsi);
			}

			// Loop over all translations
			for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
			{

				long int ihidden = iorientclass * exp_nr_trans + itrans;

				// To speed things up, only calculate pdf_offset at the coarse sampling.
				// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
				Matrix1D<double> my_offset, my_prior;
				sampling.getTranslation(itrans, my_offset);
				// Convert offsets back to Angstroms to calculate PDF!
				// TODO: if series, then have different exp_old_xoff for each my_image_no....
				// WHAT TO DO WITH THIS?!!!

				double pdf_offset;
				if (mymodel.ref_dim == 2)
					pdf_offset = calculatePdfOffset(exp_old_offset[exp_iimage] + my_offset, mymodel.prior_offset_class[exp_iclass]);
				else
					pdf_offset = calculatePdfOffset(exp_old_offset[exp_iimage] + my_offset, exp_prior[exp_iimage]);

				// TMP DEBUGGING
				if (mymodel.orientational_prior_mode != NOPRIOR && (pdf_offset==0. || pdf_orientation==0.))
				{
					global_mutex.lock();
					std::cerr << " pdf_offset= " << pdf_offset << " pdf_orientation= " << pdf_orientation << std::endl;
					std::cerr << " exp_ipart= " << exp_ipart << " exp_part_id= " << exp_part_id << std::endl;
					std::cerr << " iorient= " << iorient << " idir= " << idir << " ipsi= " << ipsi << std::endl;
					std::cerr << " exp_nr_psi= " << exp_nr_psi << " exp_nr_dir= " << exp_nr_dir << " exp_nr_trans= " << exp_nr_trans << std::endl;
					for (long int i = 0; i < sampling.directions_prior.size(); i++)
						std::cerr << " sampling.directions_prior["<<i<<"]= " << sampling.directions_prior[i] << std::endl;
					for (long int i = 0; i < sampling.psi_prior.size(); i++)
						std::cerr << " sampling.psi_prior["<<i<<"]= " << sampling.psi_prior[i] << std::endl;
					REPORT_ERROR("ERROR! pdf_offset==0.|| pdf_orientation==0.");
					global_mutex.unlock();
				}
				if (exp_nr_oversampled_rot == 0)
					REPORT_ERROR("exp_nr_oversampled_rot == 0");
				if (exp_nr_oversampled_trans == 0)
					REPORT_ERROR("exp_nr_oversampled_trans == 0");


#ifdef TIMING
				// Only time one thread, as I also only time one MPI process
				if (thread_id == 0)
					timer.tic(TIMING_WEIGHT_EXP);
#endif

				// Now first loop over iover_rot, because that is the order in exp_Mweight as well
				long int ihidden_over = ihidden * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
				for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
				{
					// Then loop over iover_trans
					for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, ihidden_over++)
					{

#ifdef DEBUG_CHECKSIZES
						if (ihidden_over >= XSIZE(exp_Mweight))
						{
							std::cerr<< "ihidden_over= "<<ihidden_over<<" XSIZE(Mweight)= "<< XSIZE(exp_Mweight) <<std::endl;
							REPORT_ERROR("ihidden_over >= XSIZE(exp_Mweight)");
						}
#endif

						// Only exponentiate for determined values of exp_Mweight
						// (this is always true in the first pass, but not so in the second pass)
						// Only deal with this sampling point if its weight was significant
#ifdef DEBUG_CHECKSIZES
						if (exp_iimage >= YSIZE(exp_Mweight))
						{
							std::cerr<< "exp_iimage= "<<exp_iimage<<" YSIZE(exp_Mweight)= "<< YSIZE(exp_Mweight) <<std::endl;
							std::cerr << " exp_ipart= " << exp_ipart << std::endl;
							std::cerr << " DIRECT_A2D_ELEM(exp_Mweight, exp_iimage, ihidden_over)= " << DIRECT_A2D_ELEM(exp_Mweight, exp_iimage, ihidden_over) << std::endl;
							std::cerr << " DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over)= " << DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over) << std::endl;
							REPORT_ERROR("exp_iimage >= YSIZE(exp_Mweight)");
						}
#endif

						if (DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over) < 0.)
						{
							DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over) = 0.;
						}
						else
						{
							double weight = pdf_orientation * pdf_offset;

							double diff2 = DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over) - exp_min_diff2[exp_ipart];

							// next line because of numerical precision of exp-function
							if (diff2 > 700.) weight = 0.;
							// TODO: use tabulated exp function?
							else weight *= exp(-diff2);
//#define DEBUG_PSIANGLE_PDISTRIBUTION
#ifdef DEBUG_PSIANGLE_PDISTRIBUTION
							std::cout << ipsi*360./sampling.NrPsiSamplings() << " "<< weight << std::endl;
#endif
							// Store the weight
							DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, ihidden_over) = weight;

#ifdef DEBUG_CHECKSIZES
							if (std::isnan(weight))
							{
								global_mutex.lock();
								std::cerr<< "weight= "<<weight<<" is not a number! " <<std::endl;
								std::cerr << " exp_min_diff2[exp_ipart]= " << exp_min_diff2[exp_ipart] << std::endl;
								std::cerr << " exp_ipart= " << exp_ipart << std::endl;
								std::cerr << " exp_part_id= " << exp_part_id << std::endl;
								std::cerr << " mydata.getNrImagesInSeries(exp_part_id)= " << mydata.getNrImagesInSeries(exp_part_id) << std::endl;
								long int my_image_no = exp_iimage + 0;
								std::cerr << " my_image_no= " << my_image_no << std::endl;
								std::cerr << " exp_iimage= " << exp_iimage << std::endl;
								std::cerr << " DIRECT_A2D_ELEM(exp_Mweight, my_image_no, ihidden_over)= " << DIRECT_A2D_ELEM(exp_Mweight, my_image_no, ihidden_over) << std::endl;
								REPORT_ERROR("weight is not a number");
								global_mutex.unlock();
							}
#endif

							// Keep track of sum and maximum of all weights for this particle
							// Later add all to exp_thisparticle_sumweight, but inside this loop sum to local thisthread_sumweight first
							thisthread_sumweight += weight;

						} // end if/else exp_Mweight < 0.
					} // end loop iover_trans
				}// end loop iover_rot
#ifdef TIMING
				// Only time one thread, as I also only time one MPI process
				if (thread_id == 0)
					timer.toc(TIMING_WEIGHT_EXP);
#endif
			} // end loop itrans

		} // end loop iorient
	} // end while task distributor

	// Now inside a mutex update the sum of all weights
	global_mutex.lock();
	exp_thisparticle_sumweight += thisthread_sumweight;
	global_mutex.unlock();

	// Wait until all threads have finished
	global_barrier->wait();



#ifdef DEBUG_THREAD
    std::cerr << "leaving doThreadConvertSquaredDifferencesToWeightsAllOrientations" << std::endl;
#endif
}

void MlOptimiser::convertAllSquaredDifferencesToWeights()
{

#ifdef TIMING
	if (exp_ipass == 0) timer.tic(TIMING_ESP_WEIGHT1);
	else timer.tic(TIMING_ESP_WEIGHT2);
#endif

	// Convert the squared differences into weights
	// Note there is only one weight for each part_id, because a whole series of images is treated as one particle

	// Initialising...
	exp_sum_weight.resize(exp_nr_particles);
	for (int i = 0; i < exp_nr_particles; i++)
		exp_sum_weight[i] = 0.;

//#define DEBUG_CONVERTDIFF2W
#ifdef DEBUG_CONVERTDIFF2W
	double max_weight = -1.;
	double opt_psi, opt_xoff, opt_yoff;
	int opt_iover_rot, opt_iover_trans, opt_ipsi, opt_itrans;
	long int opt_ihidden, opt_ihidden_over;
#endif

	//TMP DEBUGGING
	//DEBUGGING_COPY_exp_Mweight = exp_Mweight;

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	exp_iimage = 0;
	exp_ipart = 0;
	for (long int ori_part_id = exp_my_first_ori_particle; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
	{
		// loop over all particles inside this ori_particle
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, exp_ipart++)
		{
			exp_part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			exp_thisparticle_sumweight = 0.;

            
			if ((iter == 1 && do_firstiter_cc) || do_always_cc)
			{
				// Binarize the squared differences array to skip marginalisation
				// Note this loop is not threaded. This is not so important because it will only be executed in the 1st iteration and is fast anyway
				double mymindiff2 = 99.e10, mymaxprob = -99.e10;
				long int myminidx = -1;
				// Find the smallest element in this row of exp_Mweight
				for (long int i = 0; i < XSIZE(exp_Mweight); i++)
				{
					double cc = DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, i);
					// ignore non-determined cc
					if (cc == -999.)
						continue;

					if (do_sim_anneal && iter > 1)
					{
						// P_accept = exp ( - (CCold -CC)/temperature)
						// cc is negative value, so use "+ cc"
						double my_prob = rnd_unif() * exp(-(exp_local_oldcc[exp_ipart] + cc)/temperature);
						if (my_prob > mymaxprob)
						{
							mymaxprob = my_prob;
							mymindiff2 = cc;
							myminidx = i;
						}
					}
					else
					{
						// just search for the maximum
						if (cc < mymindiff2)
						{
							mymindiff2 = cc;
							myminidx = i;
						}
					}
				}
				// Set all except for the best hidden variable to zero and the smallest element to 1
				for (long int i = 0; i < XSIZE(exp_Mweight); i++)
					DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, i)= 0.;

				DIRECT_A2D_ELEM(exp_Mweight, exp_ipart, myminidx)= 1.;
				exp_thisparticle_sumweight += 1.;

			}
			else
			{
				for (exp_iclass = iclass_min; exp_iclass <= iclass_max; exp_iclass++)
				{

					// The loops over all orientations are parallelised using threads
					exp_iorient_ThreadTaskDistributor->reset(); // reset thread distribution tasks
					global_ThreadManager->run(globalThreadConvertSquaredDifferencesToWeightsAllOrientations);

				} // end loop iclass

			} // end else iter==1 && do_firstiter_cc

			// Keep track of number of processed images
			exp_iimage += mydata.getNrImagesInSeries(exp_part_id);

			//Store parameters for this particle
			exp_sum_weight[exp_ipart] = exp_thisparticle_sumweight;

			// Check the sum of weights is not zero
// On a Mac, the isnan function does not compile. Just uncomment the define statement, as this is merely a debugging statement
//#define MAC_OSX
#ifndef MAC_OSX
			if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
			{
				std::cerr << " exp_thisparticle_sumweight= " << exp_thisparticle_sumweight << std::endl;
				Image<double> It;
				It() = exp_Mweight;
				It.write("Mweight.spi");
				//It() = DEBUGGING_COPY_exp_Mweight;
				//It.write("Mweight_copy.spi");
				It().resize(exp_Mcoarse_significant);
				if (MULTIDIM_SIZE(It()) > 0)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(It())
					{
						if (DIRECT_MULTIDIM_ELEM(exp_Mcoarse_significant, n))
							DIRECT_MULTIDIM_ELEM(It(), n) = 1.;
						else
							DIRECT_MULTIDIM_ELEM(It(), n) = 0.;
					}
					It.write("Mcoarse_significant.spi");
				}
				std::cerr << " exp_part_id= " << exp_part_id << "exp_iimage="<<exp_iimage<<std::endl;
				int group_id = mydata.getGroupId(exp_part_id, 0);
				std::cerr << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
				std::cerr << " exp_ipass= " << exp_ipass << std::endl;
				std::cerr << " sampling.NrDirections(0, true)= " << sampling.NrDirections(0, true)
						<< " sampling.NrDirections(0, false)= " << sampling.NrDirections(0, false) << std::endl;
				std::cerr << " sampling.NrPsiSamplings(0, true)= " << sampling.NrPsiSamplings(0, true)
						<< " sampling.NrPsiSamplings(0, false)= " << sampling.NrPsiSamplings(0, false) << std::endl;
				std::cerr << " mymodel.sigma2_noise[exp_ipart]= " << mymodel.sigma2_noise[exp_ipart] << std::endl;
				std::cerr << " wsum_model.sigma2_noise[exp_ipart]= " << wsum_model.sigma2_noise[exp_ipart] << std::endl;
				if (mymodel.orientational_prior_mode == NOPRIOR)
					std::cerr << " wsum_model.pdf_direction[exp_ipart]= " << wsum_model.pdf_direction[exp_ipart] << std::endl;
				if (do_norm_correction)
				{
					std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
					std::cerr << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
				}

				std::cerr << "written out Mweight.spi" << std::endl;
				std::cerr << " exp_thisparticle_sumweight= " << exp_thisparticle_sumweight << std::endl;
				std::cerr << " exp_min_diff2[exp_ipart]= " << exp_min_diff2[exp_ipart] << std::endl;
				REPORT_ERROR("ERROR!!! zero sum of weights....");
			}
#endif

		} // end loop part_id (i)
	} // end loop ori_part_id

	// The remainder of this function is not threaded.

	// Initialise exp_Mcoarse_significant
	if (exp_ipass==0)
		exp_Mcoarse_significant.resize(exp_nr_particles, XSIZE(exp_Mweight));

	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
	exp_significant_weight.clear();
	exp_significant_weight.resize(exp_nr_particles, 0.);
	for (long int ori_part_id = exp_my_first_ori_particle, my_image_no = 0, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
	{
		// loop over all particles inside this ori_particle
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

	#ifdef TIMING
	timer.tic(TIMING_WEIGHT_SORT);
	#endif
			MultidimArray<double> sorted_weight;
			// Get the relevant row for this particle
			exp_Mweight.getRow(ipart, sorted_weight);

			// Only select non-zero probabilities to speed up sorting
			long int np = 0;
            
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
			{
				if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
				{
					DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
					np++;
				}
			}
			sorted_weight.resize(np);

			// Sort from low to high values
            sorted_weight.sort();

	#ifdef TIMING
	timer.toc(TIMING_WEIGHT_SORT);
	#endif
			double frac_weight = 0.;
			double my_significant_weight;
			long int my_nr_significant_coarse_samples = 0;
			for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
			{
				if (exp_ipass==0) my_nr_significant_coarse_samples++;
				my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
				frac_weight += my_significant_weight;
				if (frac_weight > adaptive_fraction * exp_sum_weight[ipart])
					break;
			}

	#ifdef DEBUG_SORT
			// Check sorted array is really sorted
			double prev = 0.;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
			{
				if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) < prev)
				{
					Image<double> It;
					It()=sorted_weight;
					It() *= 10000;
					It.write("sorted_weight.spi");
					std::cerr << "written sorted_weight.spi" << std::endl;
					REPORT_ERROR("Error in sorting!");
				}
				prev=DIRECT_MULTIDIM_ELEM(sorted_weight, n);
			}
	#endif

			if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
			{
				std::cerr << " ipart= " << ipart << " adaptive_fraction= " << adaptive_fraction << std::endl;
				std::cerr << " frac-weight= " << frac_weight << std::endl;
				std::cerr << " exp_sum_weight[ipart]= " << exp_sum_weight[ipart] << std::endl;
				Image<double> It;
				std::cerr << " XSIZE(exp_Mweight)= " << XSIZE(exp_Mweight) << std::endl;
				It()=exp_Mweight;
				It() *= 10000;
				It.write("Mweight2.spi");
				std::cerr << "written Mweight2.spi" << std::endl;
				std::cerr << " np= " << np << std::endl;
				It()=sorted_weight;
				It() *= 10000;
				std::cerr << " XSIZE(sorted_weight)= " << XSIZE(sorted_weight) << std::endl;
				if (XSIZE(sorted_weight) > 0)
				{
					It.write("sorted_weight.spi");
					std::cerr << "written sorted_weight.spi" << std::endl;
				}
				REPORT_ERROR("my_nr_significant_coarse_samples == 0");
			}

			if (exp_ipass==0)
			{
				// Store nr_significant_coarse_samples for all images in this series
				for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++, my_image_no++)
				{
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NR_SIGN) = (double)my_nr_significant_coarse_samples;
				}

				// Keep track of which coarse samplings were significant were significant for this particle
				for (int ihidden = 0; ihidden < XSIZE(exp_Mcoarse_significant); ihidden++)
				{
					if (DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden) >= my_significant_weight)
						DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden) = true;
					else
						DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden) = false;
				}

			}
			exp_significant_weight[ipart] = my_significant_weight;
	#ifdef DEBUG_OVERSAMPLING
			std::cerr << " sum_weight[ipart]= " << exp_sum_weight[ipart] << " my_significant_weight= " << my_significant_weight << std::endl;
			std::cerr << " my_nr_significant_coarse_samples= " << my_nr_significant_coarse_samples << std::endl;
			std::cerr << " ipass= " << exp_ipass << " Pmax="<<DIRECT_A1D_ELEM(sorted_weight,XSIZE(sorted_weight) - 1)/frac_weight
					<<" nr_sign_sam= "<<nr_significant_samples<<" sign w= "<<exp_significant_weight<< "sum_weight= "<<exp_sum_weight<<std::endl;
	#endif

		} // end loop part_id (i)
	} // end loop ori_part_id


#ifdef DEBUG_CONVERTDIFF2W
	//Image<double> tt;
	//tt()=sorted_weight;
	//tt.write("sorted_weight.spi");
	//std::cerr << "written sorted_weight.spi" << std::endl;
	std::cerr << " ipass= " << exp_ipass << " exp_part_id= " << exp_part_id << std::endl;
	std::cerr << " diff2w: opt_xoff= " << opt_xoff << " opt_yoff= " << opt_yoff << " opt_psi= " << opt_psi << std::endl;
	std::cerr << " diff2w: opt_iover_rot= " << opt_iover_rot << " opt_iover_trans= " << opt_iover_trans << " opt_ipsi= " << opt_ipsi << std::endl;
	std::cerr << " diff2w: opt_itrans= " << opt_itrans << " opt_ihidden= " << opt_ihidden << " opt_ihidden_over= " << opt_ihidden_over << std::endl;
	std::cerr << "significant_weight= " << exp_significant_weight << " max_weight= " << max_weight << std::endl;
	std::cerr << "nr_significant_coarse_samples= " << nr_significant_coarse_samples <<std::endl;
	debug2 = (double)opt_ihidden_over;
#endif

#ifdef TIMING
	if (exp_ipass == 0) timer.toc(TIMING_ESP_WEIGHT1);
	else timer.toc(TIMING_ESP_WEIGHT2);
#endif

}

void MlOptimiser::doOpenCLStoreWeightedSumsAllOrientationsCalculations(int first_ori, int last_ori) {
#ifdef DEBUG_THREAD
    std::cerr << "entering doOpenCLStoreWeightedSumsAllOrientations" << std::endl;
#endif
    
#ifdef TIMING
    timer.tic(TIMING_WSUM_PROJ);
#endif
    
    MultidimArray<double> zeroArray;
    std::vector<MultidimArray<double> > cl_wsum_sigma2_noise, cl_wsum_scale_correction_XA, cl_wsum_scale_correction_AA;
    struct timeval startTV, endTV;
    std::vector<double> cl_wsum_norm_correction;
    
    std::vector<MultidimArray<double> > cl_wsum_pdf_direction;
    std::vector<double> cl_sumw_group, cl_wsum_pdf_class, cl_wsum_prior_offsetx_class, cl_wsum_prior_offsety_class, cl_max_weight;
    double cl_wsum_sigma2_offset;
    MultidimArray<double> cl_metadata;
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif
    
    // Wsum_sigma_noise2 is a 1D-spectrum for each group
    zeroArray.initZeros(mymodel.ori_size/2 + 1);
    cl_wsum_sigma2_noise.resize(mymodel.nr_groups);
    for (int n = 0; n < mymodel.nr_groups; n++)
    {
        cl_wsum_sigma2_noise[n] = zeroArray;
    }
    
    cl_wsum_scale_correction_XA.resize(exp_nr_particles);
    cl_wsum_scale_correction_AA.resize(exp_nr_particles);
    for (int n = 0; n < exp_nr_particles; n++)
    {
        cl_wsum_scale_correction_XA[n] = zeroArray;
        cl_wsum_scale_correction_AA[n] = zeroArray;
    }
    
    // wsum_norm_correction is a double for each particle
    cl_wsum_norm_correction.resize(exp_nr_particles, 0.);
    
    // wsum_pdf_direction is a 1D-array (of length sampling.NrDirections(0, true)) for each class
    zeroArray.initZeros(sampling.NrDirections(0, true));
    cl_wsum_pdf_direction.resize(mymodel.nr_classes);
    for (int n = 0; n < mymodel.nr_classes; n++)
    {
        cl_wsum_pdf_direction[n] = zeroArray;
    }
    // sumw_group is a double for each group
    cl_sumw_group.resize(mymodel.nr_groups, 0.);
    // wsum_pdf_class is a double for each class
    cl_wsum_pdf_class.resize(mymodel.nr_classes, 0.);
    if (mymodel.ref_dim == 2)
    {
        cl_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.);
        cl_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.);
    }
    // max_weight is a double for each particle
    cl_max_weight.resize(exp_nr_particles, 0.);
    // wsum_sigma2_offset is just a double
    cl_wsum_sigma2_offset = 0.;
    // metadata is a 2D array of nr_particles x METADATA_LINE_LENGTH
    cl_metadata.initZeros(exp_metadata);
    
    
    //    std::cerr << "starting CL SWS" << std::endl;
    std::vector< Matrix1D<double> > oversampled_orientations, oversampled_translations;
    Matrix2D<double> A;
    MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_shift, Fimg_shift_nomask;
    MultidimArray<double> Minvsigma2, Mctf, Fweight;
    double rot, tilt, psi;
    bool have_warned_small_scale = false;
    
    // Initialising...
    Fref.resize(exp_Fimgs[0]);
    Frefctf.resize(exp_Fimgs[0]);
    Fweight.resize(exp_Fimgs[0]);
    
    // Initialise Mctf to all-1 for if !do_ctf_corection
    Mctf.resize(exp_Fimgs[0]);
    Mctf.initConstant(1.);
    
    // Initialise Minvsigma2 to all-1 for if !do_map
    Minvsigma2.resize(exp_Fimgs[0]);
    Minvsigma2.initConstant(1.);
    
    long int nr_orients = sampling.NrDirections() * sampling.NrPsiSamplings();
    long int iorientclass_offset = exp_iclass * exp_nr_rot;
    
    
    int thread_id = 0;
    
    int wgs = CL_WGS;
    
    //Get number of orients for openCL save data...
    std::vector<bool> translationUsedInAnyRotation;
    translationUsedInAnyRotation.resize(exp_nr_trans * exp_nr_oversampled_trans, false);
    
    long int totalOrients = 0;
    for (long int iorient = first_ori; iorient < last_ori; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        
        // Only proceed if any of the particles had any significant coarsely sampled translation
        if (isSignificantAnyParticleAnyTranslation(iorientclass))
        {
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                bool anySignificantWeightTranslations = false;
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        // Which number was this image in the combined array of iseries and part_idpart_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        int group_id = mydata.getGroupId(part_id, exp_iseries);
                        
                        
                        long int ihidden = iorientclass * exp_nr_trans;
                        for (long int itrans = 0, currentTrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
                        {
                            
                            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                            {
                                
                                // Only deal with this sampling point if its weight was significant
                                long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                iover_rot * exp_nr_oversampled_trans + iover_trans;
                                
                                double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                if (weight >= exp_significant_weight[ipart]) {
                                    anySignificantWeightTranslations = true;
                                    translationUsedInAnyRotation[currentTrans] = true;
                                }
                                
                                currentTrans++;
                                
                            } // end loop iover_trans
                        } // end loop itrans
                    }// end loop part_id (i)
                } // end loop ori_part_id
                if (anySignificantWeightTranslations) {
                    totalOrients++;
                }
            }// end if iover_rot
        }// end loop do_proceed
    }
    
    int usedTrans = 0;
    for (int i = 0; i < exp_nr_trans * exp_nr_oversampled_trans; i++) {
        if (translationUsedInAnyRotation[i]) {
            usedTrans++;
        }
    }
    
    
    //Load model into OpenCL device Memory
    cl_mem cl_modelData, cl_outputModelValues, cl_outputModelWeights, cl_modelValues, cl_modelWeights;
    cl_uint4 modelDim;
    unsigned int outputModelElements;
    int cl_modelEUOversampling;
    int model_size;
    if (!do_skip_maximization) {
        MultidimArray<Complex > modelData = (mymodel.PPref[exp_iclass]).data;
        model_size = XSIZE(modelData) * YSIZE(modelData) * ZSIZE(modelData);
        //        std::cerr << "Model Size: " << model_size << std::endl;
        //        std::cerr << "Buffer size for Back projection: " << sizeof(double) * 3 * model_size << " Temp back proj buffer size: " << sizeof(double) * 3 * model_size * CL_maxComputeUnits * 3 << std::endl;
        cl_modelData = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * model_size, NULL, 0);
        enqueueCLMemoryWrite(cl_modelData, modelData.data, sizeof(cl_double2) * model_size);
        modelDim.x = XSIZE(modelData); modelDim.y = YSIZE(modelData); modelDim.z = ZSIZE(modelData);
        
        //Also allocate memory for the weighted sum model
        
        //If atomics are not supported, create temporary arrays for for each compute unit
        if (CL_atomicSupport == false) {
            //Want this to be as large as possible for speed, so size it to fit within the maxMemAlloc for the device
            cl_modelEUOversampling = CL_maxMemAlloc / (sizeof(cl_double4) * model_size * CL_maxComputeUnits);
            
            outputModelElements = wgs * ceil((float)(model_size * CL_maxComputeUnits * cl_modelEUOversampling) / (float)wgs);
            
            //Check that rounding up to wgs hasn't exceeded the maxMemAlloc
            if (sizeof(cl_double4) * outputModelElements > CL_maxMemAlloc) {
                cl_modelEUOversampling--;
                outputModelElements = wgs * ceil((float)(model_size * CL_maxComputeUnits * cl_modelEUOversampling) / (float)wgs);
            }
            
            //            std::cerr << "Iter: " << iter << " EU oversampling: " << cl_modelEUOversampling << std::endl;
            //If there isn't enough room, don't bother to use this device...
            if (cl_modelEUOversampling <= 0) {
                do_use_opencl = false;
                clReleaseMemObject(cl_modelData);
                //Need to still do these SWS calculations
                exp_iorient_ThreadTaskDistributor->reset(); // reset thread distribution tasks
                global_ThreadManager->run(globalThreadStoreWeightedSumsAllOrientations);
                return;
            }
            
            cl_outputModelValues = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double4) * outputModelElements, NULL, 0);
            //        cl_outputModelWeights = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * outputModelElements, NULL, 0);
            blankCLMemory(cl_outputModelValues, 4 * outputModelElements);
            //        blankCLMemory(cl_outputModelWeights, outputModelElements);
        }
        
        cl_modelValues = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * model_size, NULL, 0);
        cl_modelWeights = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * model_size, NULL, 0);
        modelData = (wsum_model.BPref[exp_iclass]).data;
        MultidimArray<double > modelWeight = (wsum_model.BPref[exp_iclass]).weight;
        enqueueCLMemoryWrite(cl_modelValues, modelData.data, sizeof(cl_double2) * model_size);
        enqueueCLMemoryWrite(cl_modelWeights, modelWeight.data, sizeof(cl_double) * model_size);
        
    }
    
    cl_uint2 projDim; projDim.x = XSIZE(Fref); projDim.y = YSIZE(Fref);
    cl_uint projSize = projDim.x * projDim.y;
    int err;
    int numWGSPerProj = ceil((float)(projDim.x * projDim.y) / (float)wgs);
    
    int iterProjSize = wgs * ceil((float)(totalOrients * projSize) / (float)wgs);
    cl_mem CL_iter_Proj = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * iterProjSize, NULL, 0);
    blankCLMemory(CL_iter_Proj, 2 * iterProjSize);
    
    
    cl_mem cl_FrefCTF = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * totalOrients * projSize, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to allocate Fref memory, error: " << err << std::endl;
    }
    
    //Each matrix is 3x3, as is exp_R_mic
    cl_mem cl_Ainv = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * totalOrients * 9, NULL, 0);
    cl_mem cl_exp_R_mic = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * 9, NULL, 0);
    clEnqueueWriteBuffer(CL_ComputeQueue, cl_exp_R_mic, true, 0, sizeof(cl_double) * 9, exp_R_mic.mdata, 0, NULL, NULL);
    
    clFinish(CL_ComputeQueue);
    
    //    std::cerr << "Transferred exp_r_mic" << std::endl;
    
    cl_mem cl_eulerAngles = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * totalOrients * 3, NULL, 0);
    
    cl_mem cl_eulerAnglesHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * totalOrients * 3, NULL, 0);
    double *eulerAnglesHost = (double *)clEnqueueMapBuffer(CL_ComputeQueue, cl_eulerAnglesHost, true, CL_MAP_WRITE, 0, sizeof(cl_double) * totalOrients * 3, 0, 0, 0, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to map buffer for transfer of euler angles" << std::endl;
        do_use_opencl = false;
    }
    
    long int *iorientClassArray = (long int *)malloc(sizeof(long int) * totalOrients);
    
    int currentOrients = 0;
    
    //Count the number of particles
    int partCount = 0;
    for (long int ori_part_id = exp_my_first_ori_particle; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
    {
        // loop over all particles inside this ori_particle
        //                    std::cerr << "Number of particles for ori_part_id: " << ori_part_id << " is: " << mydata.ori_particles[ori_part_id].particles_id.size() << std::endl;
        for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, partCount++)
        {
        }
    }
    
    double *scales = (double *)malloc(sizeof(double) * partCount);
    cl_mem cl_local_Fctfs = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * projSize, NULL, 0);
    cl_mem cl_Mctfs = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * projSize, NULL, 0);
    cl_mem cl_Minvsigma2 = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * projSize, NULL, 0);
    //    blankCLMemory(cl_Minvsigma2, projSize);
    //    blankCLMemory(cl_local_Fctfs, projSize);
    //    blankCLMemory(cl_Mctfs, projSize);
    
    //    std::cerr << "Setup ctf arrays" << std::endl;
    
    long int *myImageArray = (long int *)malloc(sizeof(long int) * partCount);
    
    //    int totalTrans = exp_nr_trans * exp_nr_oversampled_trans * partCount;
    int totalTrans = usedTrans * partCount;
    
    long int *FimgShiftArray = (long int *)malloc(sizeof(long int) * totalTrans);
    
    int *groupIdArray = (int *)malloc(sizeof(int) * partCount);
    
    int iUsedTrans = 0;
    for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
    {
        // loop over all particles inside this ori_particle
        for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
        {
            long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
            
            // Which number was this image in the combined array of iseries and part_idpart_id
            long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
            myImageArray[ipart] = my_image_no;
            
            int group_id = mydata.getGroupId(part_id, exp_iseries);
            groupIdArray[ipart] = group_id;
            
            if (!do_skip_maximization)
            {
                
                if (do_scale_correction)
                {
                    // TODO: implemenent B-factor as well...
                    double myscale = mymodel.scale_correction[group_id];
                    if (myscale > 10000.)
                    {
                        std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << " my_image_no= " << my_image_no << std::endl;
                        REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
                    }
                    else if (myscale < 0.001)
                    {
                        
                        if (!have_warned_small_scale)
                        {
                            std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << myscale <<
                            "); Use larger groups for more stable scale estimates." << std::endl;
                            have_warned_small_scale = true;
                        }
                        myscale = 0.001;
                    }
                    scales[ipart] = myscale;
                }
            } // end if !do_skip_maximization
            
            for (long int itrans = 0, currentTrans = 0; itrans < exp_nr_trans; itrans++)
            {
                
                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                {
                    
                    if (!do_skip_maximization)
                    {
                        // Get the shifted image
                        long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
                        itrans * exp_nr_oversampled_trans + iover_trans;
                        if (translationUsedInAnyRotation[currentTrans]) {
                            //                            std::cerr << "Used - iTrans: " << itrans << " iOvertrans: " << iover_trans << " ishift: " << ishift << " iUsedTrans: " << iUsedTrans << std::endl;
                            FimgShiftArray[iUsedTrans] = ishift;
                            iUsedTrans++;
                        }
                        currentTrans++;
                    }
                }
            }
        }
    }
    
    //    std::cerr << "Transferred trans data" << std::endl;
    //    std::cerr << "Going to allocate weight arrays - alloc size: " << sizeof(double) * totalTrans * totalOrients << std::endl;
    
    cl_mem cl_weights = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * totalOrients * usedTrans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    cl_mem cl_weightsHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * totalOrients * usedTrans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    
    
    //    double weights[partCount][totalTrans * totalOrients];
    double *weights = (double *)malloc(sizeof(double) * totalTrans * totalOrients);
    memset(weights, 0, sizeof(double) * totalOrients * totalTrans);
    
    //    std::cerr << "Allocated weights arrays" << std::endl;
    
    // exp_iclass loop does not always go from 0 to nr_classes!
    size_t first_iorient = 0, last_iorient = 0;
    for (long int iorient = first_ori; iorient < last_ori; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        
        // Only proceed if any of the particles had any significant coarsely sampled translation
        if (isSignificantAnyParticleAnyTranslation(iorientclass))
        {
            
            long int idir = iorient / exp_nr_psi;
            long int ipsi = iorient % exp_nr_psi;
            
            // Now get the oversampled (rot, tilt, psi) triplets
            // This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
            sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_orientations);
            
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                
                bool anySignificantWeightTranslations = false;
                
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        // Which number was this image in the combined array of iseries and part_idpart_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        int group_id = mydata.getGroupId(part_id, exp_iseries);
                        
                        
                        long int ihidden = iorientclass * exp_nr_trans;
                        int trans_current = 0;
                        int iUsedTrans = 0;
                        for (long int itrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
                        {
                            
                            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                            {
                                
                                // Only deal with this sampling point if its weight was significant
                                long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                iover_rot * exp_nr_oversampled_trans + iover_trans;
                                
                                double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                if (weight >= exp_significant_weight[ipart]) {
                                    
                                    int index = iUsedTrans * totalOrients + currentOrients;
                                    //                                    std::cerr << "Current Orient: " << currentOrients << " iOrient: " << iorient << " iOverRot: " << iover_rot << " iTrans: " << itrans << " iOverTrans: " << iover_trans << " ipart: " << ipart << " current Trans: " << trans_current << " index: " << ipart * totalOrients * usedTrans + index << " iUsedTrans: " << iUsedTrans << std::endl;
                                    weight /= exp_sum_weight[ipart];
                                    //                                    weights[ipart][index] = weight;
                                    weights[ipart * totalOrients * usedTrans + index] = weight;
                                    
                                    anySignificantWeightTranslations = true;
                                }
                                if (translationUsedInAnyRotation[trans_current]) {
                                    iUsedTrans++;
                                }
                                trans_current++;
                                
                                
                            } // end loop iover_trans
                        } // end loop itrans
                    }// end loop part_id (i)
                } // end loop ori_part_id
                
                if (anySignificantWeightTranslations) {
                    iorientClassArray[currentOrients] = iorientclass;
                    memcpy((eulerAnglesHost+currentOrients * 3), oversampled_orientations[iover_rot].vdata, sizeof(cl_double) * 3);
                    
                    currentOrients++;
                }
            }// end if iover_rot
        }// end loop do_proceed
        
    } // end loop ipsi
    
    clEnqueueUnmapMemObject(CL_ComputeQueue, cl_eulerAnglesHost, eulerAnglesHost, 0, NULL, NULL);
    clEnqueueCopyBuffer(CL_ComputeQueue, cl_eulerAnglesHost, cl_eulerAngles, 0, 0, sizeof(cl_double) * totalOrients * 3, 0, NULL, NULL);
    clFinish(CL_ComputeQueue);
    
    //    std::cerr << "Transferred euler angle arrays" << std::endl;
    
    cl_int cl_do_ctf_correction = (do_ctf_correction && refs_are_ctf_corrected);
    cl_int cl_do_scale_correction = do_scale_correction;
    cl_int cl_num_particles = partCount;
    //    cl_int cl_num_trans = exp_nr_trans * exp_nr_oversampled_trans;
    cl_int cl_num_trans = usedTrans;
    cl_int cl_IS_INV = IS_INV;
    cl_double cl_paddingFactor = (double)(mymodel.PPref[exp_iclass]).padding_factor;
    cl_uint cl_limit = totalOrients;
    
    cl_mem cl_FimgShiftAll = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * projSize * cl_num_trans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    cl_mem cl_FimgShiftAllHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * projSize * cl_num_trans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    cl_mem cl_FimgShiftNoMaskAll = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * projSize * cl_num_trans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    cl_mem cl_FimgShiftNoMaskAllHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * projSize * cl_num_trans, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error allocating memory, error: " << err << std::endl;
    }
    
    
    cl_mem cl_Mresol_fine = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_int) * projSize, NULL, 0);
    initCLMemoryWithConstant(cl_Mresol_fine, projSize, -1);
    int mresolSize = XSIZE(Mresol_fine) * YSIZE(Mresol_fine);
    clEnqueueWriteBuffer(CL_ComputeQueue, cl_Mresol_fine, true, 0, sizeof(int) * mresolSize, Mresol_fine.data, 0, NULL, NULL);
    clFinish(CL_ComputeQueue);
    
    cl_mem cl_sigma2_noise = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    
    cl_mem cl_sumXA = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    cl_mem cl_sumAA = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    
    //    cl_mem cl_norm_correction = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs, NULL, 0);
    //    blankCLMemory(cl_Mresol_fine, numWGSPerProj * wgs);
    blankCLMemory(cl_sigma2_noise, wgs * projSize);
    blankCLMemory(cl_sumXA, wgs * projSize);
    blankCLMemory(cl_sumAA, wgs * projSize);
    
    cl_mem cl_sigma2_noise2 = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    
    cl_mem cl_sumXA2 = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    cl_mem cl_sumAA2 = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs * projSize, NULL, 0);
    
    //    cl_mem cl_norm_correction = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * wgs, NULL, 0);
    //    blankCLMemory(cl_Mresol_fine, numWGSPerProj * wgs);
    blankCLMemory(cl_sigma2_noise2, wgs * projSize);
    blankCLMemory(cl_sumXA2, wgs * projSize);
    blankCLMemory(cl_sumAA2, wgs * projSize);
    
    unsigned int sumSize = wgs * ceil((float)totalOrients * projSize / (float)wgs);
    cl_mem cl_FimgOut = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * sumSize, NULL, 0);
    cl_mem cl_FweightOut = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * sumSize, NULL, 0);
    blankCLMemory(cl_FimgOut, 2 * sumSize);
    blankCLMemory(cl_FweightOut, sumSize);
    
    //            blankCLMemory(cl_norm_correction, wgs);
    
    cl_mem cl_dataVsPrior = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double) * projSize, NULL, 0);
    double *dataVsPrior = (double *)malloc(sizeof(double) * projSize);
    memset(dataVsPrior, 0, sizeof(double) * projSize);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
    {
        int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
        
        if (ires > -1)
        {
            dataVsPrior[n] = DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[exp_iclass], ires);
            //            std::cerr << "n: " << n << " dataVsPrior: " << dataVsPrior[n] << std::endl;
        }
    }
    clEnqueueWriteBuffer(CL_ComputeQueue, cl_dataVsPrior, true, 0, sizeof(cl_double) * projSize, dataVsPrior, 0, NULL, NULL);
    
    std::vector< std::vector< Matrix1D<double> > > precalculated_oversampled_translations;
    precalculated_oversampled_translations.resize(exp_nr_trans);
    for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
    {
        sampling.getTranslations(itrans, adaptive_oversampling, precalculated_oversampled_translations[itrans]);
    }
    
    cl_mem cl_oversampled_translations = clCreateBuffer(CL_context, CL_MEM_READ_ONLY, sizeof(cl_double2) * cl_num_trans, NULL, 0);
    cl_mem cl_oversampled_translationsHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double2) * cl_num_trans, NULL, 0);
    cl_double2 *oversampled_translationsHost = (cl_double2 *)clEnqueueMapBuffer(CL_ComputeQueue, cl_oversampled_translationsHost, true, CL_MAP_WRITE, 0, sizeof(cl_double2) * cl_num_trans, 0, 0, 0, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to map buffer for transfer of euler angles" << std::endl;
        do_use_opencl = false;
    }
    int currentTrans = 0;
    for (int trans = 0; trans < exp_nr_trans; trans++) {
        for (int ot = 0; ot < exp_nr_oversampled_trans; ot++) {
            if (translationUsedInAnyRotation[trans * exp_nr_oversampled_trans + ot]) {
                memcpy(&oversampled_translationsHost[currentTrans], precalculated_oversampled_translations[trans][ot].vdata, sizeof(cl_double2));
                //                std::cerr << "trans: " << trans << " ot: " << ot << " index: " << trans * exp_nr_oversampled_trans + ot << " copied: " << oversampled_translationsHost[currentTrans].x << " " << oversampled_translationsHost[currentTrans].y << " orig: " << precalculated_oversampled_translations[trans][ot] << std::endl;
                currentTrans++;
            }
        }
    }
    clEnqueueUnmapMemObject(CL_ComputeQueue, cl_oversampled_translationsHost, oversampled_translationsHost, 0, NULL, NULL);
    clEnqueueCopyBuffer(CL_ComputeQueue, cl_oversampled_translationsHost, cl_oversampled_translations, 0, 0, sizeof(cl_double2) * cl_num_trans, 0, NULL, NULL);
    
    int num_orients_wgs_multiple = wgs * (totalOrients / wgs + 1);
    cl_mem cl_sigma2_offset_return = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double) * num_orients_wgs_multiple, NULL, 0);
    blankCLMemory(cl_sigma2_offset_return, num_orients_wgs_multiple);
    
    cl_mem cl_prior_offset_return;
    if (mymodel.ref_dim == 2) {
        cl_prior_offset_return = clCreateBuffer(CL_context, CL_MEM_READ_WRITE, sizeof(cl_double2) * num_orients_wgs_multiple, NULL, 0);
        blankCLMemory(cl_prior_offset_return, 2 * num_orients_wgs_multiple);
    }
    
    clFinish(CL_ComputeQueue);
    free(dataVsPrior);
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    unsigned long sstartMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    unsigned long sendMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double selapsedTime = (double)(sendMicros - sstartMicros) / 1000000;
    std::cerr << "SWS Prep time: " << selapsedTime << " seconds" << std::endl;
#endif
    
#ifdef TIMING
    timer.toc(TIMING_WSUM_PROJ);
#endif
#ifdef TIMING
    timer.tic(TIMING_WSUM_DIFF2);
#endif
    
#pragma mark CL calculations start
    //    std::cerr << "Generating orients" << std::endl;
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif
    
    err = clSetKernelArg(CL_generateOrientationMatrix, 0, sizeof(cl_mem), &cl_eulerAngles);
    err |= clSetKernelArg(CL_generateOrientationMatrix, 1, sizeof(cl_mem), &cl_exp_R_mic);
    err |= clSetKernelArg(CL_generateOrientationMatrix, 2, sizeof(cl_int), &cl_IS_INV);
    err |= clSetKernelArg(CL_generateOrientationMatrix, 3, sizeof(cl_double), &cl_paddingFactor);
    err |= clSetKernelArg(CL_generateOrientationMatrix, 4, sizeof(cl_mem), &cl_Ainv);
    err |= clSetKernelArg(CL_generateOrientationMatrix, 5, sizeof(cl_uint), &cl_limit);
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to set OpenCL kernel arguments - generate matrix" << std::endl;
    }
    
    size_t global[1], local[1];
    global[0] = wgs * ceil((float)totalOrients / (float)wgs);
    local[0] = wgs;
    
    err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_generateOrientationMatrix, 1, NULL, global, local, 0, NULL, NULL);
    if (err)
    {
        std::cerr << "Error: Failed to execute OpenCL kernel - SWS: matrix generation, err: " << err << std::endl;
    }
    clFinish(CL_ComputeQueue);
    
    //    std::cerr << "Generating projections" << std::endl;
    if (!do_skip_maximization) {
        
        cl_int model_r_max = (mymodel.PPref[exp_iclass]).r_max;
        cl_limit = totalOrients * projSize;
        
        if (mymodel.ref_dim == 2) {
            err = clSetKernelArg(CL_calcModelProjection2D, 0, sizeof(cl_uint4), &modelDim);
            err |= clSetKernelArg(CL_calcModelProjection2D, 1, sizeof(cl_uint), &model_r_max);
            err |= clSetKernelArg(CL_calcModelProjection2D, 2, sizeof(cl_uint2), &projDim);
            err |= clSetKernelArg(CL_calcModelProjection2D, 3, sizeof(cl_mem), &cl_Ainv);
            err |= clSetKernelArg(CL_calcModelProjection2D, 4, sizeof(cl_mem), &cl_modelData);
            err |= clSetKernelArg(CL_calcModelProjection2D, 5, sizeof(cl_mem), &CL_iter_Proj);
            err |= clSetKernelArg(CL_calcModelProjection2D, 6, sizeof(cl_int), &cl_limit);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments - model projection" << std::endl;
            }
            
            //                std::cerr << "Finished loading CL data - model projection" << std::endl;
            global[0] = wgs * ceil((float)totalOrients * projSize / (float)wgs);
            local[0] = wgs;
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcModelProjection2D, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model projection, err: " << err << std::endl;
            }
        } else if (mymodel.ref_dim == 3) {
            err = clSetKernelArg(CL_calculateProjections, 0, sizeof(cl_uint4), &modelDim);
            err |= clSetKernelArg(CL_calculateProjections, 1, sizeof(cl_uint), &model_r_max);
            err |= clSetKernelArg(CL_calculateProjections, 2, sizeof(cl_uint2), &projDim);
            err |= clSetKernelArg(CL_calculateProjections, 3, sizeof(cl_mem), &cl_Ainv);
            err |= clSetKernelArg(CL_calculateProjections, 4, sizeof(cl_mem), &cl_modelData);
            err |= clSetKernelArg(CL_calculateProjections, 5, sizeof(cl_mem), &CL_iter_Proj);
            err |= clSetKernelArg(CL_calculateProjections, 6, sizeof(cl_int), &cl_limit);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments - model projection" << std::endl;
            }
            
            //                std::cerr << "Finished loading CL data - model projection" << std::endl;
            global[0] = wgs * ceil((float)totalOrients * projSize / (float)wgs);
            local[0] = wgs;
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateProjections, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model projection, err: " << err << std::endl;
            }
        }
        clFinish(CL_ComputeQueue);
    }
    
    
    cl_event clEvents[partCount];
    for (int ipart = 0; ipart < partCount; ipart++) {
        //        gettimeofday(&startTV, NULL);
        
        long int my_image_no = myImageArray[ipart];
        
        if (!do_skip_maximization) {
            //            std::cerr << "Starting ctf etc, particle: " << ipart << " numWGPerProj: " << numWGSPerProj << std::endl;
            double *ipartCTF = exp_local_Fctfs[my_image_no].data;
            
            if (do_ctf_correction) {
                clEnqueueWriteBuffer(CL_ComputeQueue, cl_local_Fctfs, true, 0, sizeof(cl_double) * projDim.x * projDim.y, ipartCTF, 0, NULL, NULL);
            } else {
                initCLMemoryWithConstantDouble(cl_local_Fctfs, projDim.x * projDim.y, 1.0);
            }
            //            std::cerr << "Loaded ctf - from image number: " << my_image_no << " First value: " << ipartCTF[1] << std::endl;
            
            double *minsigv2 = exp_local_Minvsigma2s[my_image_no].data;
            clEnqueueWriteBuffer(CL_ComputeQueue, cl_Minvsigma2, true, 0, sizeof(cl_double) * projDim.x * projDim.y, minsigv2, 0, NULL, NULL);
            clFinish(CL_ComputeQueue);
            
            //            std::cerr << "Starting ctf" << std::endl;
            cl_double scaleVal = scales[ipart];
            
            cl_limit = totalOrients * projSize;
            err = 0;
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 0, sizeof(cl_mem), &CL_iter_Proj);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 1, sizeof(cl_int), &cl_do_ctf_correction);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 2, sizeof(cl_mem), &cl_local_Fctfs);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 3, sizeof(cl_int), &cl_do_scale_correction);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 4, sizeof(cl_double), &scaleVal);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 5, sizeof(cl_mem), &cl_FrefCTF);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 6, sizeof(cl_int), &projSize);
            err |= clSetKernelArg(CL_ctfAndScaleDataPoint, 7, sizeof(cl_uint), &cl_limit);
            
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments" << std::endl;
            }
            
            //                std::cerr << "Finished loading CL data" << std::endl;
            //            int wgs = 1;
            size_t global[1], local[1];
            //    global[0] = wgs * ceil((float)limit / (float)wgs);;
            global[0] = wgs * ceil((float)totalOrients * projSize / (float)wgs);
            local[0] = wgs;
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_ctfAndScaleDataPoint, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: ctf and scale, err: " << err << std::endl;
                //                std::cerr << "Orients: " << totalOrients << " projSize: " << projSize << " global: " << global[0] << std::endl;
            }
            //            clFinish(CL_ComputeQueue);
            
            //            gettimeofday(&endTV, NULL);
            //            unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            double elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS CTF done elapsed Time: " << elapsedTime << std::endl;
            
            
            //            std::cerr << "Starting to scale Fctfs..." << std::endl;
            //Also scale the ctf for later use...
            cl_limit = projSize;
            err = 0;
            err |= clSetKernelArg(CL_scaleCTF, 0, sizeof(cl_mem), &cl_local_Fctfs);
            err |= clSetKernelArg(CL_scaleCTF, 1, sizeof(cl_int), &cl_do_scale_correction);
            err |= clSetKernelArg(CL_scaleCTF, 2, sizeof(cl_double), &scaleVal);
            err |= clSetKernelArg(CL_scaleCTF, 3, sizeof(cl_mem), &cl_Mctfs);
            err |= clSetKernelArg(CL_scaleCTF, 4, sizeof(cl_uint), &cl_limit);
            
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to set OpenCL kernel arguments" << std::endl;
            }
            
            global[0] = numWGSPerProj * wgs;
            local[0] = wgs;
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_scaleCTF, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: scale ctfs, err: " << err << std::endl;
            }
            //            clFinish(CL_ComputeQueue);
            
            //            std::cerr << "Transfer of weights and Fimgshifts to device" << std::endl;
            //Transfer data to the device
            //            clEnqueueWriteBuffer(CL_queue, cl_weights, true, 0, sizeof(double) * totalOrients * cl_num_trans, weights + ipart * cl_num_trans * totalOrients, 0, NULL, NULL);
            //            clFinish(CL_queue);
            
            //            gettimeofday(&endTV, NULL);
            //            startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS CTF scale done elapsed Time: " << elapsedTime << std::endl;
            
            cl_double *weightsHost = (cl_double *)clEnqueueMapBuffer(CL_ComputeQueue, cl_weightsHost, true, CL_MAP_WRITE, 0, sizeof(double) * totalOrients * cl_num_trans, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of weights" << std::endl;
                do_use_opencl = false;
            }
            memcpy(weightsHost, weights + ipart * cl_num_trans * totalOrients, sizeof(double) * totalOrients * cl_num_trans);
            clEnqueueUnmapMemObject(CL_ComputeQueue, cl_weightsHost, weightsHost, 0, NULL, NULL);
            clEnqueueCopyBuffer(CL_ComputeQueue, cl_weightsHost, cl_weights, 0, 0, sizeof(double) * totalOrients * cl_num_trans, 0, NULL, NULL);
            //            clFinish(CL_ComputeQueue);
            //            std::cerr << "Copied Weights" << std::endl;
            
            
            cl_limit = totalOrients * cl_num_trans;
            cl_double2 prior_offset_diff;
            if (mymodel.ref_dim == 2) {
                //For 2d also calculate the prior offsets
                prior_offset_diff.x = XX(exp_old_offset[my_image_no]);
                prior_offset_diff.y = YY(exp_old_offset[my_image_no]);
                
                err = 0;
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 0, sizeof(cl_mem), &cl_weights);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 1, sizeof(cl_double2), &prior_offset_diff);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 2, sizeof(cl_mem), &cl_oversampled_translations);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 3, sizeof(cl_mem), &cl_prior_offset_return);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 4, sizeof(cl_int), &totalOrients);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 5, sizeof(cl_int), &totalOrients);
                err |= clSetKernelArg(CL_calculate2DPriorOffset, 6, sizeof(cl_int), &cl_limit);
                
                global[0] = wgs * (totalOrients / wgs + 1);
                
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculate2DPriorOffset, 1, NULL, global, local, 0, NULL, NULL);
                if (err)
                {
                    std::cerr << "Error: Failed to execute OpenCL prior offset kernel, err: " << err << std::endl;
                }
                
                
                prior_offset_diff.x = XX(mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no]);
                prior_offset_diff.y = YY(mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no]);
            } else if (mymodel.ref_dim == 3) {
                prior_offset_diff.x = XX(exp_prior[my_image_no] - exp_old_offset[my_image_no]);
                prior_offset_diff.y = YY(exp_prior[my_image_no] - exp_old_offset[my_image_no]);
            }
            
            err = 0;
            err |= clSetKernelArg(CL_calculateSigma2Offset, 0, sizeof(cl_mem), &cl_weights);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 1, sizeof(cl_double2), &prior_offset_diff);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 2, sizeof(cl_mem), &cl_oversampled_translations);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 3, sizeof(cl_mem), &cl_sigma2_offset_return);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 4, sizeof(cl_int), &totalOrients);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 5, sizeof(cl_int), &totalOrients);
            err |= clSetKernelArg(CL_calculateSigma2Offset, 6, sizeof(cl_int), &cl_limit);
            
            global[0] = wgs * (totalOrients / wgs + 1);
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSigma2Offset, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL sigma2 offset kernel, err: " << err << std::endl;
            }
            
            cl_double2 *FimgShiftAllHost = (cl_double2 *)clEnqueueMapBuffer(CL_ComputeQueue, cl_FimgShiftAllHost, true, CL_MAP_WRITE, 0, sizeof(cl_double2) * projSize * cl_num_trans, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of Fimg shifts" << std::endl;
                do_use_opencl = false;
            }
            //            memset(FimgShiftAllHost, 0, sizeof(cl_double2) * numWGSPerProj * wgs * cl_num_trans);
            for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                //                std::cerr << "Copying weight: " << itrans << std::endl;
                long int ishift = FimgShiftArray[ipart * cl_num_trans + itrans];
                Fimg_shift = exp_local_Fimgs_shifted[ishift];
                memcpy(FimgShiftAllHost + itrans * projSize, Fimg_shift.data, sizeof(cl_double2) * projDim.x * projDim.y);
            }
            
            clEnqueueUnmapMemObject(CL_ComputeQueue, cl_FimgShiftAllHost, FimgShiftAllHost, 0, NULL, NULL);
            clEnqueueCopyBuffer(CL_ComputeQueue, cl_FimgShiftAllHost, cl_FimgShiftAll, 0, 0, sizeof(cl_double2) * projSize * cl_num_trans, 0, NULL, NULL);
            //            clFinish(CL_ComputeQueue);
            //            std::cerr <<  "Transferred image shifts" << std::endl;
            
            //            gettimeofday(&endTV, NULL);
            //            startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS Fimg transfer done elapsed Time: " << elapsedTime << std::endl;
            
            cl_mem cl_current_sigma2_noise = ipart % 2 ? cl_sigma2_noise : cl_sigma2_noise2;
            
            cl_uint clLimit = totalOrients * projSize * cl_num_trans;
            
            cl_mem cl_current_sumXA = ipart % 2 ? cl_sumXA : cl_sumXA2;
            cl_mem cl_current_sumAA = ipart % 2 ? cl_sumAA : cl_sumAA2;
            
            err = 0;
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 0, sizeof(cl_mem), &cl_FrefCTF);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 1, sizeof(cl_mem), &cl_FimgShiftAll);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 2, sizeof(cl_mem), &cl_weights);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 3, sizeof(cl_mem), &cl_Mresol_fine);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 4, sizeof(cl_mem), &cl_current_sigma2_noise);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 5, sizeof(cl_mem), &cl_dataVsPrior);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 6, sizeof(cl_mem), &cl_current_sumXA);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 7, sizeof(cl_mem), &cl_current_sumAA);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 8, sizeof(cl_int), &projSize);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 9, sizeof(cl_int), &totalOrients);
            err |= clSetKernelArg(CL_calculateSWSNoiseEstNormCorrection, 10, sizeof(cl_uint), &clLimit);
            
            global[0] = wgs * projSize;
            
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSWSNoiseEstNormCorrection, 1, NULL, global, local, 0, NULL, NULL);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: noise estimation and norm correction, err: " << err << std::endl;
            }
            
            /*
             err = 0;
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 0, sizeof(cl_mem), &cl_FrefCTF);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 1, sizeof(cl_mem), &cl_FimgShiftAll);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 2, sizeof(cl_mem), &cl_weights);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 3, sizeof(cl_mem), &cl_Mresol_fine);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 4, sizeof(cl_mem), &cl_current_sigma2_noise);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 5, sizeof(cl_int), &numWGSPerProj);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 6, sizeof(cl_int), &totalOrients);
             err |= clSetKernelArg(CL_calculateSWSNoiseEstimate, 7, sizeof(cl_uint), &clLimit);
             
             
             global[0] = CL_maxComputeUnits * numWGSPerProj * wgs;
             
             err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSWSNoiseEstimate, 1, NULL, global, local, 0, NULL, NULL);
             if (err)
             {
             std::cerr << "Error: Failed to execute OpenCL kernel, err: " << err << std::endl;
             }
             clFinish(CL_ComputeQueue);
             //            std::cerr << "Noise Calculations complete" << std::endl;
             //            double norm_correction;
             //            MultidimArray<double> noiseResult;
             //            noiseResult.initZeros(mymodel.ori_size/2 + 1);
             
             //            gettimeofday(&endTV, NULL);
             //            unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
             //            unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
             //            double elapsedTime = (double)(endMicros - startMicros) / 1000000;
             //            std::cerr << "SWS Noise estimate done elapsed Time: " << elapsedTime << std::endl;
             
             if (do_scale_correction) {
             cl_mem cl_current_sumXA = ipart % 2 ? cl_sumXA : cl_sumXA2;
             cl_mem cl_current_sumAA = ipart % 2 ? cl_sumAA : cl_sumAA2;
             
             err = 0;
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 0, sizeof(cl_mem), &cl_FrefCTF);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 1, sizeof(cl_mem), &cl_FimgShiftAll);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 2, sizeof(cl_mem), &cl_weights);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 3, sizeof(cl_mem), &cl_Mresol_fine);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 4, sizeof(cl_mem), &cl_dataVsPrior);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 5, sizeof(cl_mem), &cl_current_sumXA);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 6, sizeof(cl_mem), &cl_current_sumAA);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 7, sizeof(cl_int), &numWGSPerProj);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 8, sizeof(cl_int), &totalOrients);
             err |= clSetKernelArg(CL_calculateSWSNormCorrection, 9, sizeof(cl_uint), &clLimit);
             
             err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSWSNormCorrection, 1, NULL, global, local, 0, NULL, NULL);
             if (err)
             {
             std::cerr << "Error: Failed to execute OpenCL kernel, err: " << err << std::endl;
             }
             //                clFinish(CL_ComputeQueue);
             //                std::cerr << "Norm Calculations complete" << std::endl;
             
             }*/
            
            //            gettimeofday(&endTV, NULL);
            //            startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS Norm correction done elapsed Time: " << elapsedTime << std::endl;
            
            
            //            std::cerr << "Transferring Fimgshift no mask" << std::endl;
            cl_double2 *FimgShiftNoMaskAllHost = (cl_double2 *)clEnqueueMapBuffer(CL_ComputeQueue, cl_FimgShiftNoMaskAllHost, true, CL_MAP_WRITE, 0, sizeof(cl_double2) * projSize * cl_num_trans, 0, 0, 0, &err);
            if (err != CL_SUCCESS) {
                std::cerr << "Error: Failed to map buffer for transfer of Fimg shifts" << std::endl;
                do_use_opencl = false;
            }
            //            memset(FimgShiftNoMaskAllHost, 0, sizeof(cl_double2) * numWGSPerProj * wgs * cl_num_trans);
            for (int itrans = 0; itrans < cl_num_trans; itrans++) {
                //                std::cerr << "Transferring shift: " << itrans << std::endl;
                long int ishift = FimgShiftArray[ipart * cl_num_trans + itrans];
                Fimg_shift = exp_local_Fimgs_shifted_nomask[ishift];
                memcpy(FimgShiftNoMaskAllHost + itrans * projSize, Fimg_shift.data, sizeof(cl_double2) * projDim.x * projDim.y);
            }
            
            clEnqueueUnmapMemObject(CL_ComputeQueue, cl_FimgShiftNoMaskAllHost, FimgShiftNoMaskAllHost, 0, NULL, NULL);
            clEnqueueCopyBuffer(CL_ComputeQueue, cl_FimgShiftNoMaskAllHost, cl_FimgShiftNoMaskAll, 0, 0, sizeof(cl_double2) * projSize * cl_num_trans, 0, NULL, NULL);
            //            clFinish(CL_ComputeQueue);
            
            //            gettimeofday(&endTV, NULL);
            //            startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS Fimgshift transfer done elapsed Time: " << elapsedTime << std::endl;
            
            //            std::cerr << "Starting SWS Sum" << std::endl;
            cl_uint cl_gi_limit = totalOrients * projSize;
            err = 0;
            err |= clSetKernelArg(CL_calculateSWSSum, 0, sizeof(cl_mem), &cl_Mctfs);
            err |= clSetKernelArg(CL_calculateSWSSum, 1, sizeof(cl_mem), &cl_FimgShiftNoMaskAll);
            err |= clSetKernelArg(CL_calculateSWSSum, 2, sizeof(cl_mem), &cl_Minvsigma2);
            err |= clSetKernelArg(CL_calculateSWSSum, 3, sizeof(cl_mem), &cl_weights);
            err |= clSetKernelArg(CL_calculateSWSSum, 4, sizeof(cl_mem), &cl_FimgOut);
            err |= clSetKernelArg(CL_calculateSWSSum, 5, sizeof(cl_mem), &cl_FweightOut);
            err |= clSetKernelArg(CL_calculateSWSSum, 6, sizeof(cl_int), &projSize);
            err |= clSetKernelArg(CL_calculateSWSSum, 7, sizeof(cl_int), &totalOrients);
            err |= clSetKernelArg(CL_calculateSWSSum, 8, sizeof(cl_uint), &clLimit);
            err |= clSetKernelArg(CL_calculateSWSSum, 9, sizeof(cl_uint), &cl_gi_limit);
            
            global[0] = wgs * ceil((float)totalOrients * projSize / (float)wgs);
            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSWSSum, 1, NULL, global, local, 0, NULL, &clEvents[ipart]);
            if (err)
            {
                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: calculate sws sum, err: " << err << std::endl;
            }
            //            clFinish(CL_ComputeQueue);
            
            /*
             std::cerr << "Starting SWS Sum B" << std::endl;
             cl_gi_limit = totalOrients * projSize;
             err = 0;
             err |= clSetKernelArg(CL_calculateSWSSum, 0, sizeof(cl_mem), &cl_Mctfs);
             err |= clSetKernelArg(CL_calculateSWSSum, 1, sizeof(cl_mem), &cl_FimgShiftNoMaskAll);
             err |= clSetKernelArg(CL_calculateSWSSum, 2, sizeof(cl_mem), &cl_Minvsigma2);
             err |= clSetKernelArg(CL_calculateSWSSum, 3, sizeof(cl_mem), &cl_weights);
             err |= clSetKernelArg(CL_calculateSWSSum, 4, sizeof(cl_mem), &cl_FimgOut);
             err |= clSetKernelArg(CL_calculateSWSSum, 5, sizeof(cl_mem), &cl_FweightOut);
             err |= clSetKernelArg(CL_calculateSWSSum, 6, sizeof(cl_int), &projSize);
             err |= clSetKernelArg(CL_calculateSWSSum, 7, sizeof(cl_int), &totalOrients);
             err |= clSetKernelArg(CL_calculateSWSSum, 8, sizeof(cl_uint), &clLimit);
             err |= clSetKernelArg(CL_calculateSWSSum, 9, sizeof(cl_uint), &cl_gi_limit);
             
             wgs = 1;
             global[0] = wgs * ceil((float)totalOrients * projSize / (float)wgs);
             local[0] = wgs;
             err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calculateSWSSum, 1, NULL, global, local, 0, NULL, &clEvents[ipart]);
             if (err)
             {
             std::cerr << "Error: Failed to execute OpenCL kernel, err: " << err << std::endl;
             }
             clFinish(CL_ComputeQueue);
             wgs = CL_WGS;*/
            
            //            gettimeofday(&endTV, NULL);
            //            startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
            //            endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
            //            elapsedTime = (double)(endMicros - startMicros) / 1000000;
            //            std::cerr << "SWS Sum done elapsed Time: " << elapsedTime << std::endl;
            
            if (ipart == partCount - 1) {
                if (!do_skip_maximization) {
                    if (CL_atomicSupport == true) {
                        //Calculate the weighted model - only after the final particle has been processed
                        
                        cl_int model_r_max = (mymodel.PPref[exp_iclass]).r_max;
                        
                        cl_uint limit = totalOrients * projSize;
                        
                        if (mymodel.ref_dim == 2) {
                            err = clSetKernelArg(CL_calcBackProjection2D, 0, sizeof(cl_uint4), &modelDim);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 1, sizeof(cl_uint), &model_r_max);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 2, sizeof(cl_uint2), &projDim);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 3, sizeof(cl_uint), &totalOrients);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 4, sizeof(cl_mem), &cl_Ainv);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 5, sizeof(cl_mem), &cl_FimgOut);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 6, sizeof(cl_mem), &cl_FweightOut);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 7, sizeof(cl_mem), &cl_modelValues);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 8, sizeof(cl_mem), &cl_modelWeights);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 9, sizeof(cl_uint), &projSize);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 9, sizeof(cl_int), &wgs);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 10, sizeof(cl_uint), &model_size);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 10, sizeof(cl_uint), &limit);
                            if (err != CL_SUCCESS) {
                                std::cerr << "Error: Failed to set OpenCL kernel arguments - model back projection" << std::endl;
                            }
                            
                            //        std::cerr << "Finished loading CL data - model projection" << std::endl;
                            //                    global[0] = wgs * ceil((float)(model_size) / (float)wgs);;
                            //                    local[0] = wgs;
                            
                            global[0] = wgs * ceil((float)(totalOrients * projSize) / (float)wgs);
                            //                    local[0] = 1;
                            
                            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcBackProjection2D, 1, NULL, global, local, 0, NULL, NULL);
                            if (err)
                            {
                                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model back projection, err: " << err << std::endl;
                            }
                        } else if (mymodel.ref_dim == 3) {
                            err = clSetKernelArg(CL_calcBackProjection, 0, sizeof(cl_uint4), &modelDim);
                            err |= clSetKernelArg(CL_calcBackProjection, 1, sizeof(cl_uint), &model_r_max);
                            err |= clSetKernelArg(CL_calcBackProjection, 2, sizeof(cl_uint2), &projDim);
                            err |= clSetKernelArg(CL_calcBackProjection, 3, sizeof(cl_uint), &totalOrients);
                            err |= clSetKernelArg(CL_calcBackProjection, 4, sizeof(cl_mem), &cl_Ainv);
                            err |= clSetKernelArg(CL_calcBackProjection, 5, sizeof(cl_mem), &cl_FimgOut);
                            err |= clSetKernelArg(CL_calcBackProjection, 6, sizeof(cl_mem), &cl_FweightOut);
                            err |= clSetKernelArg(CL_calcBackProjection, 7, sizeof(cl_mem), &cl_modelValues);
                            err |= clSetKernelArg(CL_calcBackProjection, 8, sizeof(cl_mem), &cl_modelWeights);
                            err |= clSetKernelArg(CL_calcBackProjection, 9, sizeof(cl_uint), &projSize);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 9, sizeof(cl_int), &wgs);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 10, sizeof(cl_uint), &model_size);
                            err |= clSetKernelArg(CL_calcBackProjection, 10, sizeof(cl_uint), &limit);
                            if (err != CL_SUCCESS) {
                                std::cerr << "Error: Failed to set OpenCL kernel arguments - model back projection" << std::endl;
                            }
                            
                            //        std::cerr << "Finished loading CL data - model projection" << std::endl;
                            //                    global[0] = wgs * ceil((float)(model_size) / (float)wgs);;
                            //                    local[0] = wgs;
                            
                            global[0] = wgs * ceil((float)(totalOrients * projSize) / (float)wgs);
                            //                    local[0] = 1;
                            
                            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcBackProjection, 1, NULL, global, local, 0, NULL, NULL);
                            if (err)
                            {
                                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model back projection, err: " << err << std::endl;
                            }
                        }
                    } else {
                        //Calculate the weighted model - only after the final particle has been processed
                        cl_int model_r_max = (mymodel.PPref[exp_iclass]).r_max;
                        
                        cl_uint limit = totalOrients * projSize;
                        
                        if (mymodel.ref_dim == 2) {
                            err = clSetKernelArg(CL_calcBackProjection2D, 0, sizeof(cl_uint4), &modelDim);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 1, sizeof(cl_uint), &model_r_max);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 2, sizeof(cl_uint2), &projDim);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 3, sizeof(cl_uint), &totalOrients);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 4, sizeof(cl_mem), &cl_Ainv);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 5, sizeof(cl_mem), &cl_FimgOut);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 6, sizeof(cl_mem), &cl_FweightOut);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 7, sizeof(cl_mem), &cl_outputModelValues);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 8, sizeof(cl_uint), &projSize);
                            err |= clSetKernelArg(CL_calcBackProjection2D, 9, sizeof(cl_uint), &limit);
                            if (err != CL_SUCCESS) {
                                std::cerr << "Error: Failed to set OpenCL kernel arguments - model back projection" << std::endl;
                            }
                            
                            global[0] = CL_maxComputeUnits * cl_modelEUOversampling;
                            local[0] = 1;
                            
                            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcBackProjection2D, 1, NULL, global, local, 0, NULL, NULL);
                            if (err)
                            {
                                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model back projection, err: " << err << std::endl;
                            }
                            
                        } else if (mymodel.ref_dim == 3) {
                            err = clSetKernelArg(CL_calcBackProjection, 0, sizeof(cl_uint4), &modelDim);
                            err |= clSetKernelArg(CL_calcBackProjection, 1, sizeof(cl_uint), &model_r_max);
                            err |= clSetKernelArg(CL_calcBackProjection, 2, sizeof(cl_uint2), &projDim);
                            err |= clSetKernelArg(CL_calcBackProjection, 3, sizeof(cl_uint), &totalOrients);
                            err |= clSetKernelArg(CL_calcBackProjection, 4, sizeof(cl_mem), &cl_Ainv);
                            err |= clSetKernelArg(CL_calcBackProjection, 5, sizeof(cl_mem), &cl_FimgOut);
                            err |= clSetKernelArg(CL_calcBackProjection, 6, sizeof(cl_mem), &cl_FweightOut);
                            err |= clSetKernelArg(CL_calcBackProjection, 7, sizeof(cl_mem), &cl_outputModelValues);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 8, sizeof(cl_mem), &cl_outputModelWeights);
                            err |= clSetKernelArg(CL_calcBackProjection, 8, sizeof(cl_uint), &projSize);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 9, sizeof(cl_int), &wgs);
                            //                    err |= clSetKernelArg(CL_calcBackProjection, 10, sizeof(cl_uint), &model_size);
                            err |= clSetKernelArg(CL_calcBackProjection, 9, sizeof(cl_uint), &limit);
                            if (err != CL_SUCCESS) {
                                std::cerr << "Error: Failed to set OpenCL kernel arguments - model back projection" << std::endl;
                            }
                            
                            //        std::cerr << "Finished loading CL data - model projection" << std::endl;
                            //                    global[0] = wgs * ceil((float)(model_size) / (float)wgs);;
                            //                    local[0] = wgs;
                            
                            global[0] = CL_maxComputeUnits * cl_modelEUOversampling;
                            local[0] = 1;
                            
                            err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_calcBackProjection, 1, NULL, global, local, 0, NULL, NULL);
                            if (err)
                            {
                                std::cerr << "Error: Failed to execute OpenCL kernel - SWS: model back projection, err: " << err << std::endl;
                            }
                        }
                        
                        limit = model_size;
                        cl_int totalModels = CL_maxComputeUnits * cl_modelEUOversampling;
                        
                        err = clSetKernelArg(CL_sumModel, 0, sizeof(cl_mem), &cl_outputModelValues);
                        //                    err |= clSetKernelArg(CL_sumModel, 1, sizeof(cl_mem), &cl_outputModelWeights);
                        err |= clSetKernelArg(CL_sumModel, 1, sizeof(cl_mem), &cl_modelValues);
                        err |= clSetKernelArg(CL_sumModel, 2, sizeof(cl_mem), &cl_modelWeights);
                        err |= clSetKernelArg(CL_sumModel, 3, sizeof(cl_int), &totalModels);
                        err |= clSetKernelArg(CL_sumModel, 4, sizeof(cl_uint), &limit);
                        
                        global[0] = wgs * ceil((float)(model_size) / (float)wgs);
                        local[0] = wgs;
                        
                        err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_sumModel, 1, NULL, global, local, 0, NULL, NULL);
                        if (err)
                        {
                            std::cerr << "Error: Failed to execute OpenCL kernel - model sum, err: " << err << std::endl;
                        }
                        
                        //                    clFinish(CL_ComputeQueue);
                    }
                }
                
            }
            
#pragma mark CL save data
            if (ipart > 0) {
                //                std::cerr << "Started to copy data back" << std::endl;
                clWaitForEvents(1, &clEvents[ipart - 1]);
                
                cl_mem cl_previous_sigma2_noise = (ipart - 1) % 2 ? cl_sigma2_noise : cl_sigma2_noise2;
                
                double *sigma2_noiseReturn = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sigma2_noise, true, 0, sizeof(cl_double) * wgs * projSize, sigma2_noiseReturn, 0, NULL, NULL);
                //            clEnqueueReadBuffer(CL_CopyFromDeviceQueue cl_sigma2_noise, true, 0, sizeof(double) * projDim.x * projDim.y, noiseResult.data, 0, NULL, NULL);
                clFinish(CL_CopyFromDeviceQueue);
                if (err != CL_SUCCESS) {
                    std::cerr << "Error retrieving sigma2 noise results, error code: " << err << std::endl;
                }
                
                int group_id = groupIdArray[ipart - 1];
                for (int i = 0; i < wgs; i++) {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                    {
                        int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                        cl_wsum_norm_correction[ipart - 1] += sigma2_noiseReturn[i * projSize + n];
                        
                        if (ires > -1) {
                            DIRECT_MULTIDIM_ELEM(cl_wsum_sigma2_noise[group_id], ires) += sigma2_noiseReturn[i * projSize + n];
                        }
                    }
                }
                free(sigma2_noiseReturn);
                
                if (do_scale_correction) {
                    cl_mem cl_previous_sumXA = (ipart - 1) % 2 ? cl_sumXA : cl_sumXA2;
                    cl_mem cl_previous_sumAA = (ipart - 1) % 2 ? cl_sumAA : cl_sumAA2;
                    
                    
                    double *xa_Return = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                    err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sumXA, true, 0, sizeof(cl_double) * wgs * projSize, xa_Return, 0, NULL, NULL);
                    
                    double *aa_Return = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                    err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sumAA, true, 0, sizeof(cl_double) * wgs * projSize, aa_Return, 0, NULL, NULL);
                    clFinish(CL_CopyFromDeviceQueue);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error retrieving norm results, error code: " << err << std::endl;
                    }
                    
                    for (int i = 0; i < wgs; i++) {
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                        {
                            int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                            
                            if (ires > -1)
                            {
                                DIRECT_A1D_ELEM(cl_wsum_scale_correction_XA[ipart - 1], ires) += xa_Return[i * projSize + n];
                                
                                DIRECT_A1D_ELEM(cl_wsum_scale_correction_AA[ipart - 1], ires) += aa_Return[i * projSize + n];
                            }
                        }
                    }
                    free(xa_Return);
                    free(aa_Return);
                    
                }
                
            }
            if (ipart == partCount - 1) {
                
                //Scale the A for all rotations back so can be used when storing data
                cl_double cl_invPadding = 1.0 / cl_paddingFactor;
                cl_limit = 9 * totalOrients;
                
                err = clSetKernelArg(CL_scaleMatrix, 0, sizeof(cl_mem), &cl_Ainv);
                err |= clSetKernelArg(CL_scaleMatrix, 1, sizeof(cl_double), &cl_invPadding);
                err |= clSetKernelArg(CL_scaleMatrix, 2, sizeof(cl_uint), &cl_limit);
                if (err != CL_SUCCESS) {
                    std::cerr << "Error: Failed to set OpenCL kernel arguments - scale matrix, err: " << err << std::endl;
                }
                
                global[0] = wgs * ceil((float)(9 * totalOrients) / (float)wgs);
                
                err = clEnqueueNDRangeKernel(CL_ComputeQueue, CL_scaleMatrix, 1, NULL, global, local, 0, NULL, NULL);
                if (err)
                {
                    std::cerr << "Error: Failed to execute OpenCL kernel - SWS: scale matrix, err: " << err << std::endl;
                }
                //                clFinish(CL_ComputeQueue);
                
                clWaitForEvents(1, &clEvents[ipart]);
                
                if (mymodel.ref_dim == 2) {
                    cl_double2 *prior_offsetReturn = (cl_double2 *)malloc(sizeof(cl_double2) * num_orients_wgs_multiple);
                    err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_prior_offset_return, true, 0, sizeof(cl_double2) * num_orients_wgs_multiple, prior_offsetReturn, 0, NULL, NULL);
                    //            clEnqueueReadBuffer(CL_queue, cl_sigma2_noise, true, 0, sizeof(double) * projDim.x * projDim.y, noiseResult.data, 0, NULL, NULL);
                    clFinish(CL_CopyFromDeviceQueue);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error retrieving prior offset x results, error code: " << err << std::endl;
                    }
                    
                    for (int i = 0; i < totalOrients; i++) {
                        cl_wsum_prior_offsetx_class[exp_iclass] += prior_offsetReturn[i].x;
                        cl_wsum_prior_offsety_class[exp_iclass] += prior_offsetReturn[i].y;
                    }
                    
                    free(prior_offsetReturn);
                }
                
                double *sigma2_offsetReturn = (double *)malloc(sizeof(cl_double) * num_orients_wgs_multiple);
                err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_sigma2_offset_return, true, 0, sizeof(cl_double) * num_orients_wgs_multiple, sigma2_offsetReturn, 0, NULL, NULL);
                //            clEnqueueReadBuffer(CL_queue, cl_sigma2_noise, true, 0, sizeof(double) * projDim.x * projDim.y, noiseResult.data, 0, NULL, NULL);
                clFinish(CL_CopyFromDeviceQueue);
                if (err != CL_SUCCESS) {
                    std::cerr << "Error retrieving sigma2 offset results, error code: " << err << std::endl;
                }
                
                for (int i = 0; i < totalOrients; i++) {
                    cl_wsum_sigma2_offset += sigma2_offsetReturn[i];
                }
                free(sigma2_offsetReturn);
                
                cl_mem cl_previous_sigma2_noise = (ipart) % 2 ? cl_sigma2_noise : cl_sigma2_noise2;
                
                double *sigma2_noiseReturn = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sigma2_noise, true, 0, sizeof(cl_double) * wgs * projSize, sigma2_noiseReturn, 0, NULL, NULL);
                //            clEnqueueReadBuffer(CL_queue, cl_sigma2_noise, true, 0, sizeof(double) * projDim.x * projDim.y, noiseResult.data, 0, NULL, NULL);
                clFinish(CL_CopyFromDeviceQueue);
                if (err != CL_SUCCESS) {
                    std::cerr << "Error retrieving sigma2 noise results, error code: " << err << std::endl;
                }
                
                int group_id = groupIdArray[ipart];
                for (int i = 0; i < wgs; i++) {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                    {
                        int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                        cl_wsum_norm_correction[ipart] += sigma2_noiseReturn[i * projSize + n];
                        
                        if (ires > -1) {
                            DIRECT_MULTIDIM_ELEM(cl_wsum_sigma2_noise[group_id], ires) += sigma2_noiseReturn[i * projSize + n];
                        }
                    }
                }
                free(sigma2_noiseReturn);
                
                if (do_scale_correction) {
                    cl_mem cl_previous_sumXA = (ipart) % 2 ? cl_sumXA : cl_sumXA2;
                    cl_mem cl_previous_sumAA = (ipart) % 2 ? cl_sumAA : cl_sumAA2;
                    
                    
                    double *xa_Return = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                    err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sumXA, true, 0, sizeof(cl_double) * wgs * projSize, xa_Return, 0, NULL, NULL);
                    
                    double *aa_Return = (double *)malloc(sizeof(cl_double) * wgs * projSize);
                    err = clEnqueueReadBuffer(CL_CopyFromDeviceQueue, cl_previous_sumAA, true, 0, sizeof(cl_double) * wgs * projSize, aa_Return, 0, NULL, NULL);
                    clFinish(CL_CopyFromDeviceQueue);
                    if (err != CL_SUCCESS) {
                        std::cerr << "Error retrieving norm results, error code: " << err << std::endl;
                    }
                    
                    for (int i = 0; i < wgs; i++) {
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                        {
                            int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                            
                            if (ires > -1)
                            {
                                DIRECT_A1D_ELEM(cl_wsum_scale_correction_XA[ipart], ires) += xa_Return[i * projSize + n];
                                
                                DIRECT_A1D_ELEM(cl_wsum_scale_correction_AA[ipart], ires) += aa_Return[i * projSize + n];
                            }
                        }
                    }
                    free(xa_Return);
                    free(aa_Return);
                    
                }
            }
            //            for (int i = 0; i < CL_maxComputeUnits * numWGSPerProj * wgs; i++) {
            //                exp_wsum_norm_correction[ipart] += sigma2_noiseReturn[i];
            //                if (ipart == 0) {
            //                    std::cerr << "i: " << i << " sum: " << exp_wsum_norm_correction[ipart] << std::endl;
            //                }
            //            }
            
            
            //            exp_wsum_norm_correction[ipart] += norm_correction;
            //            wsum_model.sigma2_noise[ipart] += noiseResult;
            
        }
    }
    
    
    cl_double2 *modelValues = (cl_double2 *)malloc(sizeof(cl_double2) * model_size);
    cl_double *modelWeights = (cl_double *)malloc(sizeof(cl_double) * model_size);
    err |= clEnqueueReadBuffer(CL_ComputeQueue, cl_modelValues, true, 0, sizeof(cl_double2) * model_size, modelValues, 0, NULL, NULL);
    err |= clEnqueueReadBuffer(CL_ComputeQueue, cl_modelWeights, true, 0, sizeof(cl_double) * model_size, modelWeights, 0, NULL, NULL);
    clFinish(CL_ComputeQueue);
    if (err != CL_SUCCESS) {
        std::cerr << "Error retrieving Fimg, Fweight, model value and weight results, error code: " << err << std::endl;
    }
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    unsigned long startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    unsigned long endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedTime = (double)(endMicros - startMicros) / 1000000;
    unsigned long ops;
    //Generate matrix
    ops = totalOrients * 168;
    //CTF and scale and projections
    ops += partCount * (totalOrients * projDim.x * projDim.y * (4 + 57) + projDim.x * projDim.y);
    //Noise estimate
    ops += partCount * (totalOrients * projDim.x * projDim.y * cl_num_trans * 7);
    //SWS Norm correction
    ops += partCount * (totalOrients * projDim.x * projDim.y * cl_num_trans * 10);
    //Sum shifts
    ops += partCount * (totalOrients * projDim.x * projDim.y * cl_num_trans * 6);
    std::cerr << "SWS Time: " << elapsedTime << "s, ops: " << ops << " orients: " << totalOrients << " trans: " << usedTrans << " GFLOPS: " << (double)ops / elapsedTime / 1e9 << std::endl;
    gettimeofday(&startTV, NULL);
#endif
#ifdef TIMING
    timer.toc(TIMING_WSUM_DIFF2);
#endif
#ifdef TIMING
    timer.tic(TIMING_WSUM_SUMSHIFT);
#endif
    
#pragma mark CL cleanup
    
    cl_mem cl_AinvHost = clCreateBuffer(CL_context, CL_MEM_ALLOC_HOST_PTR, sizeof(cl_double) * totalOrients * 9, NULL, 0);
    clEnqueueCopyBuffer(CL_ComputeQueue, cl_Ainv, cl_AinvHost, 0, 0, sizeof(cl_double) * totalOrients * 9, 0, NULL, NULL);
    
    double *AinvHost = (double *)clEnqueueMapBuffer(CL_ComputeQueue, cl_AinvHost, true, CL_MAP_READ, 0, sizeof(cl_double) * totalOrients * 9, 0, 0, 0, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Error: Failed to map buffer for transfer of Ainv" << std::endl;
        do_use_opencl = false;
    }
    clFinish(CL_ComputeQueue);
    
    currentOrients = 0;
    A.resize(3, 3);
    Fimg.initZeros(Fref);
    Fweight.initZeros(Fref);
    //    std::cerr << "Completed precalculations of translations" << std::endl;
    double test;
    for (long int iorient = first_ori; iorient < last_ori; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        
        
        // Only proceed if any of the particles had any significant coarsely sampled translation
        if (isSignificantAnyParticleAnyTranslation(iorientclass))
        {
            
            long int idir = iorient / exp_nr_psi;
            long int ipsi = iorient % exp_nr_psi;
            
            // Now get the oversampled (rot, tilt, psi) triplets
            // This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
            //            sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_orientations);
            
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                
                //                rot = XX(oversampled_orientations[iover_rot]);
                //                tilt = YY(oversampled_orientations[iover_rot]);
                //                psi = ZZ(oversampled_orientations[iover_rot]);
                // Get the Euler matrix
                //                Euler_angles2matrix(rot, tilt, psi, A);
                
                // Take tilt-series into account
                //                A = (exp_R_mic * A).inv();
                
                bool anySignificantWeightTranslations = false;
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        // Which number was this image in the combined array of iseries and part_idpart_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        int group_id = mydata.getGroupId(part_id, exp_iseries);
                        
                        
                        long int ihidden = iorientclass * exp_nr_trans;
                        int trans_current = 0;
                        for (long int itrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
                        {
                            /*
                             bool getTranslations = false;
                             for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                             {
                             
                             // Only deal with this sampling point if its weight was significant
                             long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                             iover_rot * exp_nr_oversampled_trans + iover_trans;
                             
                             double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                             if (weight >= exp_significant_weight[ipart]) {
                             getTranslations = true;
                             }
                             }
                             
                             if (getTranslations) {
                             sampling.getTranslations(itrans, adaptive_oversampling, oversampled_translations);
                             }*/
                            
                            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                            {
                                
                                // Only deal with this sampling point if its weight was significant
                                long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                iover_rot * exp_nr_oversampled_trans + iover_trans;
                                
                                double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                if (weight >= exp_significant_weight[ipart]) {
                                    
                                    weight /= exp_sum_weight[ipart];
                                    
                                    if (!do_skip_maximization) {
                                        //Only large loops using the CL device - do the rest on the CPU
                                        // Store sum of weights for this group
                                        cl_sumw_group[group_id] += weight;
                                        
                                        // Store weights for this class and orientation
                                        cl_wsum_pdf_class[exp_iclass] += weight;
                                        
                                        if (mymodel.ref_dim ==2)
                                        {
                                            // Also store weighted offset differences for prior_offsets of each class
                                            //                                            cl_wsum_prior_offsetx_class[exp_iclass] += weight * XX(exp_old_offset[my_image_no] + precalculated_oversampled_translations[itrans][iover_trans]);
                                            //                                            cl_wsum_prior_offsety_class[exp_iclass] += weight * YY(exp_old_offset[my_image_no] + precalculated_oversampled_translations[itrans][iover_trans]);
                                            
                                            // Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
                                            //                                                                                                          reduces to double
                                            //double                 double     std::vector<Matrix1D<double> >           std::vector<Matrix1D<double> > std::vector< std::vector< Matrix1D<double> > >
                                            //                                            cl_wsum_sigma2_offset += weight * ((mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no] - precalculated_oversampled_translations[itrans][iover_trans]).sum2());
                                            
                                            //                                            test += weight * ((mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no] - precalculated_oversampled_translations[itrans][iover_trans]).sum2());
                                            
                                            //                                            if (currentOrients < 1) {
                                            //                                                std::cerr << "i: " << trans_current << " weight: " << weight << " prior_offset: " << mymodel.prior_offset_class[exp_iclass] << " old_offset: " << exp_old_offset[my_image_no] << " trans: " << precalculated_oversampled_translations[itrans][iover_trans] << " sum: " << (mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no] - precalculated_oversampled_translations[itrans][iover_trans]).sum2() << " s2: " << test << std::endl;
                                            //                                            }
                                        }
                                        else
                                        {
                                            // Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
                                            //                                            cl_wsum_sigma2_offset += weight * ((exp_prior[my_image_no] - exp_old_offset[my_image_no] - precalculated_oversampled_translations[itrans][iover_trans]).sum2());
                                            
                                        }
                                        
                                        
                                        // Store weight for this direction of this class
                                        if (mymodel.orientational_prior_mode == NOPRIOR)
                                        {
                                            DIRECT_MULTIDIM_ELEM(cl_wsum_pdf_direction[exp_iclass], idir) += weight;
                                        }
                                        else
                                        {
                                            // In the case of orientational priors, get the original number of the direction back
                                            long int mydir = sampling.getDirectionNumberAlsoZeroPrior(idir);
                                            DIRECT_MULTIDIM_ELEM(cl_wsum_pdf_direction[exp_iclass], mydir) += weight;
                                        }
                                    } //end do_skip_maximisation
                                    
                                    if (weight > cl_max_weight[ipart])
                                    {
                                        
                                        // Store optimal image parameters
                                        cl_max_weight[ipart] = weight;
                                        
                                        //retrieve angles calculated in OpenCL
                                        memcpy(A.mdata, AinvHost + currentOrients * 9, sizeof(double) * 9);
                                        //                                        A /= (double)(mymodel.PPref[exp_iclass]).padding_factor;
                                        
                                        // Calculate the angles back from the Euler matrix because for tilt series exp_R_mic may have changed them...
                                        //std::cerr << " ORI rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
                                        //Don't need to invert because is already inverted from the CL calculations
                                        Euler_matrix2angles(A.inv(), rot, tilt, psi);
                                        //std::cerr << " BACK rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
                                        
                                        
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_ROT) = rot;
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_TILT) = tilt;
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PSI) = psi;
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_XOFF) = XX(exp_old_offset[my_image_no]) + XX(precalculated_oversampled_translations[itrans][iover_trans]);
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_YOFF) = YY(exp_old_offset[my_image_no]) + YY(precalculated_oversampled_translations[itrans][iover_trans]);
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_CLASS) = (double)exp_iclass + 1;
                                        DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PMAX) = cl_max_weight[ipart];
                                    }
                                    
                                    anySignificantWeightTranslations = true;
                                }
                                trans_current++;
                                
                            } // end loop iover_trans
                        } // end loop itrans
                    }// end loop part_id (i)
                } // end loop ori_part_id
                
                if (anySignificantWeightTranslations) {
                    //                    MultidimArray<Complex > FimgTest;
                    //                    MultidimArray<double> FweightTest;
                    
                    // Initialising...
                    //                    FweightTest.resize(exp_Fimgs[0]);
                    //                    FimgTest.initZeros(Fref);
                    //                    FweightTest.initZeros(Fref);
                    
                    if (!do_skip_maximization) {
                        
                        //                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
                        //                        {
                        //                            double myweight = FweightOut[currentOrients * wgs * numWGSPerProj + n];
                        //                            (DIRECT_MULTIDIM_ELEM(Fimg, n)).real += FimgOut[currentOrients * wgs * numWGSPerProj + n].x;
                        //                            (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += FimgOut[currentOrients * wgs * numWGSPerProj + n].y;
                        //                            (DIRECT_MULTIDIM_ELEM(Fweight, n)) = 1.0;
                        //
                        //                        }
                        //                        std::cerr << "Current orient: " << iorient << std::endl;
                        //                        std::cerr << "Orig A: " << A << std::endl;
                        //                        memcpy(A.mdata, AinvHost + currentOrients * 9, sizeof(double) * 9);
                        //                        A /= (double)(mymodel.PPref[exp_iclass]).padding_factor;
                        //                        std::cerr << "New  A: " << A << std::endl;
                        //                        memcpy(Fimg.data, &FimgOut[currentOrients * wgs * numWGSPerProj], sizeof(cl_double2) * projDim.x * projDim.y);
                        //                        memcpy(Fweight.data, &FweightOut[currentOrients * wgs * numWGSPerProj], sizeof(cl_double) * projDim.x * projDim.y);
                        //                        std::cerr << "Copied data for orient: " << currentOrients << std::endl;
                        //                        (wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_INV, &Fweight);
                    }
                    
                    currentOrients++;
                }
            }// end if iover_rot
        }// end loop do_proceed
        
    } // end loop ipsi
    
    clEnqueueUnmapMemObject(CL_ComputeQueue, cl_AinvHost, AinvHost, 0, NULL, NULL);
    clReleaseMemObject(cl_AinvHost);
#ifdef TIMING
    timer.toc(TIMING_WSUM_SUMSHIFT);
#endif
#ifdef TIMING
    timer.tic(TIMING_WSUM_BACKPROJ);
#endif
    
    
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    unsigned long fstartMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    unsigned long fendMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double felapsedTime = (double)(fendMicros - fstartMicros) / 1000000;
    std::cerr << "SWS Clean up time: " << felapsedTime << " seconds" << std::endl;
#endif
    
#pragma mark CPU verify host calculations
#ifdef CL_VERIFY_ON_CPU
#ifdef CL_PRINT_SPEED
    gettimeofday(&startTV, NULL);
#endif
    // Make local copies of weighted sums (excepts BPrefs, which are too big)
    // so that there are not too many mutex locks below
    std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_scale_correction_XA, thr_wsum_scale_correction_AA, thr_wsum_pdf_direction;
    std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class, thr_max_weight;
    double thr_wsum_sigma2_offset;
    MultidimArray<double> thr_metadata;
    
    // Wsum_sigma_noise2 is a 1D-spectrum for each group
    zeroArray.initZeros(mymodel.ori_size/2 + 1);
    thr_wsum_sigma2_noise.resize(mymodel.nr_groups);
    for (int n = 0; n < mymodel.nr_groups; n++)
    {
        thr_wsum_sigma2_noise[n] = zeroArray;
    }
    // scale-correction terms are a spectrum for each particle
    thr_wsum_scale_correction_XA.resize(exp_nr_particles);
    thr_wsum_scale_correction_AA.resize(exp_nr_particles);
    for (int n = 0; n < exp_nr_particles; n++)
    {
        thr_wsum_scale_correction_XA[n] = zeroArray;
        thr_wsum_scale_correction_AA[n] = zeroArray;
    }
    // wsum_pdf_direction is a 1D-array (of length sampling.NrDirections(0, true)) for each class
    zeroArray.initZeros(sampling.NrDirections(0, true));
    thr_wsum_pdf_direction.resize(mymodel.nr_classes);
    for (int n = 0; n < mymodel.nr_classes; n++)
    {
        thr_wsum_pdf_direction[n] = zeroArray;
    }
    // wsum_norm_correction is a double for each particle
    thr_wsum_norm_correction.resize(exp_nr_particles, 0.);
    // sumw_group is a double for each group
    thr_sumw_group.resize(mymodel.nr_groups, 0.);
    // wsum_pdf_class is a double for each class
    thr_wsum_pdf_class.resize(mymodel.nr_classes, 0.);
    if (mymodel.ref_dim == 2)
    {
        thr_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.);
        thr_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.);
    }
    // max_weight is a double for each particle
    thr_max_weight.resize(exp_nr_particles, 0.);
    // wsum_sigma2_offset is just a double
    thr_wsum_sigma2_offset = 0.;
    // metadata is a 2D array of nr_particles x METADATA_LINE_LENGTH
    thr_metadata.initZeros(exp_metadata);
    
    //Copy weighted image sums from CL device for verification
    cl_double2 *FimgOut = (cl_double2 *)malloc(sizeof(cl_double2) * totalOrients * projSize);
    double *FweightOut = (double *)malloc(sizeof(double) * totalOrients * projSize);
    err = clEnqueueReadBuffer(CL_ComputeQueue, cl_FimgOut, true, 0, sizeof(cl_double2) * totalOrients * projSize, FimgOut, 0, NULL, NULL);
    err |= clEnqueueReadBuffer(CL_ComputeQueue, cl_FweightOut, true, 0, sizeof(cl_double) * totalOrients * projSize, FweightOut, 0, NULL, NULL);
    
    
    currentOrients = 0;
    
    
    // exp_iclass loop does not always go from 0 to nr_classes!
    //    long int iorientclass_offset = exp_iclass * exp_nr_rot;
    for (long int iorient = first_ori; iorient < last_ori; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        
        // Only proceed if any of the particles had any significant coarsely sampled translation
        //        std::cerr << "Orient: " << iorient << std::endl;
        
        if (isSignificantAnyParticleAnyTranslation(iorientclass))
        {
            
            long int idir = iorient / exp_nr_psi;
            long int ipsi = iorient % exp_nr_psi;
            
            // Now get the oversampled (rot, tilt, psi) triplets
            // This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
            sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_orientations);
            
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                rot = XX(oversampled_orientations[iover_rot]);
                tilt = YY(oversampled_orientations[iover_rot]);
                psi = ZZ(oversampled_orientations[iover_rot]);
                // Get the Euler matrix
                Euler_angles2matrix(rot, tilt, psi, A);
                
                // Take tilt-series into account
                A = (exp_R_mic * A).inv();
                
#ifdef TIMING
                // Only time one thread, as I also only time one MPI process
                if (thread_id == 0)
                    timer.tic(TIMING_WSUM_PROJ);
#endif
                // Project the reference map (into Fref)
                if (!do_skip_maximization)
                    (mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_INV);
                
#ifdef TIMING
                // Only time one thread, as I also only time one MPI process
                if (thread_id == 0)
                    timer.toc(TIMING_WSUM_PROJ);
#endif
                // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
                // Then outside this loop do the actual backprojection
                Fimg.initZeros(Fref);
                Fweight.initZeros(Fref);
                bool anySignificantWeightTranslations = false;
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        // Which number was this image in the combined array of iseries and part_idpart_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        int group_id = mydata.getGroupId(part_id, exp_iseries);
                        
                        
                        if (!do_skip_maximization)
                        {
                            if (do_map)
                                Minvsigma2 = exp_local_Minvsigma2s[my_image_no];
                            // else Minvsigma2 was initialised to ones
                            
                            // Apply CTF to reference projection
                            if (do_ctf_correction)
                            {
                                Mctf = exp_local_Fctfs[my_image_no];
                                if (refs_are_ctf_corrected)
                                {
                                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
                                    {
                                        DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(Mctf, n);
                                    }
                                }
                                else
                                {
                                    Frefctf = Fref;
                                }
                            }
                            else
                            {
                                // initialise because there are multiple particles and Mctf gets selfMultiplied for scale_correction
                                Mctf.initConstant(1.);
                                Frefctf = Fref;
                            }
                            
                            if (do_scale_correction)
                            {
                                // TODO: implemenent B-factor as well...
                                double myscale = mymodel.scale_correction[group_id];
                                if (myscale > 10000.)
                                {
                                    std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << " my_image_no= " << my_image_no << std::endl;
                                    REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
                                }
                                else if (myscale < 0.001)
                                {
                                    
                                    if (!have_warned_small_scale)
                                    {
                                        std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << myscale <<
                                        "); Use larger groups for more stable scale estimates." << std::endl;
                                        have_warned_small_scale = true;
                                    }
                                    myscale = 0.001;
                                }
                                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
                                {
                                    DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
                                }
                                // For CTF-terms in BP
                                Mctf *= myscale;
                            }
                        } // end if !do_skip_maximization
                        
                        //                        if (currentOrients == 0) {
                        //                            std::cerr << "Projections for orient: " << currentOrients << std::endl;
                        //                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
                        //                            {
                        //                                double Frefreal = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
                        //                                double Frefimag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
                        //                                std::cout << " i: " << n << " Fref: " << Frefreal << " " << Frefimag << std::endl;
                        //                            }
                        //                        }
                        
                        long int ihidden = iorientclass * exp_nr_trans;
                        for (long int itrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
                        {
                            
                            sampling.getTranslations(itrans, adaptive_oversampling, oversampled_translations);
                            
                            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                            {
                                
                                
                                
                                // Only deal with this sampling point if its weight was significant
                                long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                iover_rot * exp_nr_oversampled_trans + iover_trans;
                                
                                
                                double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                
                                // Only sum weights for non-zero weights
                                if (weight >= exp_significant_weight[ipart])
                                {
                                    anySignificantWeightTranslations = true;
#ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (thread_id == 0)
                                        timer.tic(TIMING_WSUM_DIFF2);
#endif
                                    // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                    weight /= exp_sum_weight[ipart];
                                    
                                    if (!do_skip_maximization)
                                    {
                                        // Get the shifted image
                                        long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
                                        itrans * exp_nr_oversampled_trans + iover_trans;
                                        Fimg_shift = exp_local_Fimgs_shifted[ishift];
                                        Fimg_shift_nomask = exp_local_Fimgs_shifted_nomask[ishift];
                                        
                                        // Store weighted sum of squared differences for sigma2_noise estimation
                                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                                        {
                                            int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                                            if (ires > -1)
                                            {
                                                // Use FT of masked image for noise estimation!
                                                double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                                double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                                double wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);
                                                
                                                // group-wise sigma2_noise
                                                DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[group_id], ires) += wdiff2;
                                                // For norm_correction
                                                thr_wsum_norm_correction[ipart] += wdiff2;
                                                
                                                //                                                if (weight > 1e-3) {
                                                //                                                    std::cerr << "Weight: " << weight << " ipart: " << ipart << " iOrient: " << iorient << " itrans: " << itrans << " Fref: " << (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag << " Fimg: " << (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag << " Diff: " << diff_real << " " << diff_imag << " wdiff2: " << wdiff2 << std::endl;
                                                //                                                }
                                            }
                                        }
                                        
                                        // Store the weighted sums of the norm_correction terms
                                        if (do_scale_correction)
                                        {
                                            double sumXA = 0.;
                                            double sumA2 = 0.;
                                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
                                            {
                                                int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
                                                
                                                // Once the reference becomes strongly regularised one does no longer want to store XA and AA!
                                                if (ires > -1 && DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
                                                {
                                                    sumXA = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
                                                    sumXA += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
                                                    DIRECT_A1D_ELEM(thr_wsum_scale_correction_XA[ipart], ires) += weight * sumXA;
                                                    
                                                    // This could be pre-calculated above...
                                                    sumA2 = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
                                                    sumA2 += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
                                                    DIRECT_A1D_ELEM(thr_wsum_scale_correction_AA[ipart], ires) += weight * sumA2;
                                                }
                                            }
                                        }
                                        // Store sum of weights for this group - done above
                                        thr_sumw_group[group_id] += weight;
                                        
                                        // Store weights for this class and orientation - done above
                                        thr_wsum_pdf_class[exp_iclass] += weight;
                                        
                                        if (mymodel.ref_dim ==2)
                                        {
                                            //                                             Also store weighted offset differences for prior_offsets of each class
                                            thr_wsum_prior_offsetx_class[exp_iclass] += weight * XX(exp_old_offset[my_image_no] + oversampled_translations[iover_trans]);
                                            thr_wsum_prior_offsety_class[exp_iclass] += weight * YY(exp_old_offset[my_image_no] + oversampled_translations[iover_trans]);
                                            
                                            // Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
                                            thr_wsum_sigma2_offset += weight * ((mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no] - oversampled_translations[iover_trans]).sum2());
                                            
                                        }
                                        else
                                        {
                                            // Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
                                            thr_wsum_sigma2_offset += weight * ((exp_prior[my_image_no] - exp_old_offset[my_image_no] - oversampled_translations[iover_trans]).sum2());
                                        }
                                        
                                        // Store weight for this direction of this class - done above
                                        if (mymodel.orientational_prior_mode == NOPRIOR)
                                        {
                                            DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
                                        }
                                        else
                                        {
                                            // In the case of orientational priors, get the original number of the direction back
                                            long int mydir = sampling.getDirectionNumberAlsoZeroPrior(idir);
                                            DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
                                        }
                                        
                                        
#ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (thread_id == 0)
                                            timer.toc(TIMING_WSUM_DIFF2);
                                        // Only time one thread, as I also only time one MPI process
                                        if (thread_id == 0)
                                            timer.tic(TIMING_WSUM_SUMSHIFT);
#endif
                                        
                                        // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
                                        // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
                                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                                        {
                                            double myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
                                            // Note that weightxinvsigma2 already contains the CTF!
                                            double weightxinvsigma2 = weight * myctf * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
                                            // now Fimg stores sum of all shifted w*Fimg
                                            (DIRECT_MULTIDIM_ELEM(Fimg, n)).real += (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).real * weightxinvsigma2;
                                            (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).imag * weightxinvsigma2;
                                            // now Fweight stores sum of all w
                                            // Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
                                            DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf;
                                            //                                            if ((currentOrients == 0) && (n == 1)) {
                                            //                                                std::cerr << "iPart: " << ipart << " Ctf: " << myctf << " image no: " << my_image_no << " Weight: " << weight << " MinSigmav2: " << DIRECT_MULTIDIM_ELEM(Minvsigma2, n) << " Fimg: " << (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).imag << " Sum: " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag << " Weight: " << DIRECT_MULTIDIM_ELEM(Fweight, n) << std::endl;
                                            //                                            }
                                        }
                                        
#ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (thread_id == 0)
                                            timer.toc(TIMING_WSUM_SUMSHIFT);
#endif
                                        
                                    } // end if !do_skip_maximization
                                    
                                    // Keep track of max_weight and the corresponding optimal hidden variables
                                    if (weight > thr_max_weight[ipart])
                                    {
                                        // Store optimal image parameters
                                        thr_max_weight[ipart] = weight;
                                        
                                        // Calculate the angles back from the Euler matrix because for tilt series exp_R_mic may have changed them...
                                        //std::cerr << " ORI rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
                                        Euler_matrix2angles(A.inv(), rot, tilt, psi);
                                        //std::cerr << " BACK rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
                                        
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_ROT) = rot;
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_TILT) = tilt;
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PSI) = psi;
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_XOFF) = XX(exp_old_offset[my_image_no]) + XX(oversampled_translations[iover_trans]);
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_YOFF) = YY(exp_old_offset[my_image_no]) + YY(oversampled_translations[iover_trans]);
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_CLASS) = (double)exp_iclass + 1;
                                        DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PMAX) = thr_max_weight[ipart];
                                    }
                                    
                                } // end if weight >= exp_significant_weight
                            } // end loop iover_trans
                        } // end loop itrans
                    }// end loop part_id (i)
                } // end loop ori_part_id
                
                if (!do_skip_maximization)
                {
#ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (thread_id == 0)
                        timer.tic(TIMING_WSUM_BACKPROJ);
#endif
                    
                    MultidimArray<Complex > FimgTest;
                    MultidimArray<double> FweightTest;
                    
                    // Initialising...
                    FweightTest.resize(exp_Fimgs[0]);
                    FimgTest.initZeros(Fref);
                    FweightTest.initZeros(Fref);
                    
                    if (anySignificantWeightTranslations) {
                        bool wrongFimgFound = false;
                        bool wrongFweightFound = false;
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                        {
                            double diff_real = fabs((DIRECT_MULTIDIM_ELEM(Fimg, n)).real - FimgOut[currentOrients * projSize + n].x);
                            double diff_imag = fabs((DIRECT_MULTIDIM_ELEM(Fimg, n)).imag - FimgOut[currentOrients * projSize + n].y);
                            if ((diff_real < 0.000001) && (diff_imag < 0.000001)){
                            } else {
                                wrongFimgFound = true;
                            }
                            
                            //                            double diff_weight = fabs((DIRECT_MULTIDIM_ELEM(Fweight, n)) - FweightOut[currentOrients * projSize + n]);
                            //                            if (diff_weight < 0.000001) {
                            //                            } else {
                            //                                wrongFweightFound = true;
                            //                            }
                            
                        }
                        if (wrongFimgFound) {
                            std::cerr << "Current Orient: " << currentOrients << std::endl;
                            std::cerr << "Wrong - Fimg shift" << std::endl;
                        } else {
                            //                            std::cerr << "Current Orient: " << currentOrients << std::endl;
                            //                            std::cerr << "Right - Fimg shift" << std::endl;
                        }
                        //                        if (wrongFweightFound) {
                        //                            std::cerr << "Wrong - Fweight" << std::endl;
                        //                        } else {
                        //                            std::cerr << "Right - Fweight" << std::endl;
                        //                        }
                        
                        
                        //                        memcpy(FimgTest.data, &FimgOut[currentOrients * projSize], sizeof(cl_double2) * projDim.x * projDim.y);
                        memcpy(FweightTest.data, &FweightOut[currentOrients * projSize], sizeof(cl_double) * projDim.x * projDim.y);
                        
                        //                        if (Fimg.equal(FimgTest)) {
                        //                            std::cerr << "Current Orient: " << currentOrients << " Fimg = FimgTest" << std::endl;
                        //                        } else {
                        //                            std::cerr << "Current Orient: " << currentOrients << " Fimg != FimgTest" << std::endl;
                        //                        }
                        if (Fweight.equal(FweightTest)) {
                            //                            std::cerr << "Current Orient: " << currentOrients << " Fweight = FweightTest" << std::endl;
                        } else {
                            std::cerr << "Current Orient: " << currentOrients << " Fweight != FweightTest" << std::endl;
                        }
                        
                    }
                    
#ifdef CL_PRINT_SAMPLE_VALUES
                    if (anySignificantWeightTranslations) {
                        std::cerr << "Current Orient: " << currentOrients << std::endl;
                        int numSampleValuesPrinted = 0;
                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
                        {
                            if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                                double diff_real = fabs((DIRECT_MULTIDIM_ELEM(Fimg, n)).real - FimgOut[currentOrients * projSize + n].x);
                                double diff_imag = fabs((DIRECT_MULTIDIM_ELEM(Fimg, n)).imag - FimgOut[currentOrients * projSize + n].y);
                                if ((diff_real < 0.000001) && (diff_imag < 0.000001)){
                                    std::cerr << "Fimg Right - CPU: " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag << " CL: " << FimgOut[currentOrients * projSize + n].x << " " << FimgOut[currentOrients * projSize + n].y << std::endl;
                                } else {
                                    std::cerr << "Fimg Wrong - CPU: " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).real << " " << (DIRECT_MULTIDIM_ELEM(Fimg, n)).imag << " CL: " << FimgOut[currentOrients * projSize + n].x << " " << FimgOut[currentOrients * projSize + n].y << std::endl;
                                }
                                
                                double diff_weight = fabs((DIRECT_MULTIDIM_ELEM(Fweight, n)) - FweightOut[currentOrients * projSize + n]);
                                if (diff_weight < 0.000001) {
                                    std::cerr << "Fweight Right - CPU: " << (DIRECT_MULTIDIM_ELEM(Fweight, n)) << " CL: " << FweightOut[currentOrients * projSize + n] << std::endl;
                                } else {
                                    std::cerr << "Fweight Wrong - CPU: " << (DIRECT_MULTIDIM_ELEM(Fweight, n)) << " CL: " << FweightOut[currentOrients * projSize + n] << std::endl;
                                }
                                numSampleValuesPrinted++;
                            }
                            
                        }
                    }
#endif
                    
                    // Perform the actual back-projection.
                    // This is done with the sum of all (in-plane) shifted Fimg's
                    // Perform this inside a mutex
                    global_mutex2.lock();
                    //                    if (anySignificantWeightTranslations && (currentOrients == 0)) {
                    //                        (wsum_model.BPref[exp_iclass]).set2DFourierTransformverbose(Fimg, A, IS_INV, &Fweight);
                    //                    } else {
                    (wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_INV, &Fweight);
                    //                    }
                    global_mutex2.unlock();
                    
#ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (thread_id == 0)
                        timer.toc(TIMING_WSUM_BACKPROJ);
#endif
                } // end if !do_skip_maximization
                
                if (anySignificantWeightTranslations) {
                    currentOrients++;
                }
            }// end if iover_rot
        }// end loop do_proceed
        
    } // end loop ipsi
    
#ifdef CL_PRINT_SPEED
    gettimeofday(&endTV, NULL);
    startMicros = 1000000*startTV.tv_sec + startTV.tv_usec;
    endMicros = 1000000*endTV.tv_sec + endTV.tv_usec;
    
    double elapsedCPUTime = (double)(endMicros - startMicros) / 1000000;
    std::cerr << "SWS CPU Time: " << elapsedCPUTime << "s GPU total time: " << selapsedTime + elapsedTime + felapsedTime << " s, equivalent to: " << elapsedCPUTime / (selapsedTime + elapsedTime + felapsedTime) << " cores" << std::endl;
#endif
#pragma mark CL Verification Comparisons
    
    //Verify that the model weights and values are correct!
    MultidimArray<Complex > modelData = (wsum_model.BPref[exp_iclass]).data;
    MultidimArray<double > modelWeight = (wsum_model.BPref[exp_iclass]).weight;
    bool valuesCorrect = true;
    bool weightsCorrect = true;
    int numSampleValuesPrinted = 0;
#ifdef CL_PRINT_SAMPLE_VALUES
    std::cerr << "Model sample values" << std::endl;
#endif
    
    for (int i = 0; i < model_size; i++) {
        double diffWeight = fabs(modelWeights[i] - DIRECT_MULTIDIM_ELEM(modelWeight, i));
        if ((diffWeight > 0.0000001) || std::isnan(modelWeights[i])) {
            weightsCorrect = false;
        }
        double diffValueReal = fabs(modelValues[i].x - DIRECT_MULTIDIM_ELEM(modelData, i).real);
        double diffValueImag = fabs(modelValues[i].y - DIRECT_MULTIDIM_ELEM(modelData, i).imag);
        if ((diffValueReal > 0.0000001) || std::isnan(modelValues[i].x)) {
            valuesCorrect = false;
        }
        if ((diffValueImag > 0.0000001) || std::isnan(modelValues[i].y)) {
            valuesCorrect = false;
        }
#ifdef CL_PRINT_SAMPLE_VALUES
        if (DIRECT_MULTIDIM_ELEM(modelWeight, i) > 0)
            if (numSampleValuesPrinted < CL_NUMBER_OF_SAMPLES_TO_PRINT) {
                std::cerr << "Model CPU: " << DIRECT_MULTIDIM_ELEM(modelData, i).real << " " << DIRECT_MULTIDIM_ELEM(modelData, i).imag << " weight CPU: " << DIRECT_MULTIDIM_ELEM(modelWeight, i) << " Model CL: " << modelValues[i].x << " " << modelValues[i].y << " weight CL: " << modelWeights[i] << std::endl;
                numSampleValuesPrinted++;
            }
#endif
    }
    if (valuesCorrect) {
        //        std::cerr << "Model Values correct" << std::endl;
    } else {
        std::cerr << "Model Values wrong!" << std::endl;
    }
    if (weightsCorrect) {
        //        std::cerr << "Model Weights correct" << std::endl;
    } else {
        std::cerr << "Model Weights wrong!" << std::endl;
    }
    
    //verify scales - XA and AA
    if (do_scale_correction)
    {
        for (int n = 0; n < exp_nr_particles; n++)
        {
            if (thr_wsum_scale_correction_XA[n].equal(cl_wsum_scale_correction_XA[n])) {
                //                    std::cerr << "Right - xa array " << n << std::endl;
#ifdef CL_PRINT_SAMPLE_VALUES
                std::cerr << "XA Sample values" << std::endl;
                for (int i = 0; i < CL_NUMBER_OF_SAMPLES_TO_PRINT; i++) {
                    std::cerr << "CPU: " << DIRECT_MULTIDIM_ELEM(thr_wsum_scale_correction_XA[n], i) << " CL: " << DIRECT_MULTIDIM_ELEM(cl_wsum_scale_correction_XA[n], i) << std::endl;
                }
#endif
            } else {
                std::cerr << "Wrong - xa array " << n << std::endl;
            }
            if (thr_wsum_scale_correction_AA[n].equal(cl_wsum_scale_correction_AA[n])) {
                //                    std::cerr << "Right - aa array " << n << std::endl;
#ifdef CL_PRINT_SAMPLE_VALUES
                std::cerr << "AA Sample values" << std::endl;
                for (int i = 0; i < CL_NUMBER_OF_SAMPLES_TO_PRINT; i++) {
                    std::cerr << "CPU: " << DIRECT_MULTIDIM_ELEM(thr_wsum_scale_correction_AA[n], i) << " CL: " << DIRECT_MULTIDIM_ELEM(cl_wsum_scale_correction_AA[n], i) << std::endl;
                }
#endif
            } else {
                std::cerr << "Wrong - aa array " << n << std::endl;
            }
            
        }
    }
    
    //verify norm_correction
    for (int n = 0; n < exp_nr_particles; n++)
    {
        double diff = fabs(cl_wsum_norm_correction[n] - thr_wsum_norm_correction[n]);
        //                std::cerr << "Diff: " << diff << std::endl;
        if ((diff > 0.0000000001) || (std::isnan(exp_wsum_norm_correction[n]))) {
            std::cerr << "Wrong - CPU norm: " << thr_wsum_norm_correction[n] << " OpenCL norm: " << exp_wsum_norm_correction[n] << std::endl;
        } else {
            //                    std::cerr << "Right - CPU norm: " << thr_wsum_norm_correction[n] << " OpenCL norm: " << exp_wsum_norm_correction[n] << std::endl;
        }
    }
    
    //Sigma2_noise and sumw_group
    for (int n = 0; n < mymodel.nr_groups; n++)
    {
        if (thr_wsum_sigma2_noise[n].equal(cl_wsum_sigma2_noise[n])) {
            //                std::cerr << "Right - sigma2noise array " << n << std::endl;
#ifdef CL_PRINT_SAMPLE_VALUES
            std::cerr << "Sig 2 noise Sample values" << std::endl;
            for (int i = 0; i < CL_NUMBER_OF_SAMPLES_TO_PRINT; i++) {
                std::cerr << "CPU: " << DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[n], i) << " CL: " << DIRECT_MULTIDIM_ELEM(cl_wsum_sigma2_noise[n], i) << std::endl;
            }
#endif
        } else {
            std::cerr << "Wrong - sigma2noise array " << n << std::endl;
        }
        
        if (fabs(thr_sumw_group[n] - cl_sumw_group[n]) < 0.0000000001) {
            //            std::cerr << "Right - sumw_group " << n << std::endl;
        } else {
            std::cerr << "Wrong - sumw_group " << n << std::endl;
        }
    }
    
    //Check pdf_class, prior_offset_class, pdf_direction, sigma2_offset
    for (int n = 0; n < mymodel.nr_classes; n++)
    {
        if (fabs(thr_wsum_pdf_class[n] - cl_wsum_pdf_class[n]) < 0.0000000001) {
            //            std::cerr << "Right - pdf_class " << n << std::endl;
        } else {
            std::cerr << "Wrong - pdf_class " << n << std::endl;
        }
        
        if (mymodel.ref_dim == 2)
        {
            if (fabs(thr_wsum_prior_offsetx_class[n] - cl_wsum_prior_offsetx_class[n]) < 0.0000000001) {
                //                std::cerr << "Right - wsum_prior_offsetx_class " << n << std::endl;
            } else {
                std::cerr << "Wrong - wsum_prior_offsetx_class " << n << std::endl;
            }
            
            if (fabs(thr_wsum_prior_offsety_class[n] - cl_wsum_prior_offsety_class[n]) < 0.0000000001) {
                //                std::cerr << "Right - wsum_prior_offsety_class " << n << std::endl;
            } else {
                std::cerr << "Wrong - wsum_prior_offsety_class " << n << std::endl;
            }
        }
        if (thr_wsum_pdf_direction[n].equal(cl_wsum_pdf_direction[n])) {
            //                std::cerr << "Right - wsum_pdf_direction array " << n << std::endl;
#ifdef CL_PRINT_SAMPLE_VALUES
            std::cerr << "wsum_pdf_direction Sample values" << std::endl;
            for (int i = 0; i < CL_NUMBER_OF_SAMPLES_TO_PRINT; i++) {
                std::cerr << "CPU: " << DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[n], i) << " CL: " << DIRECT_MULTIDIM_ELEM(cl_wsum_pdf_direction[n], i) << std::endl;
            }
#endif
        } else {
            std::cerr << "Wrong - wsum_pdf_direction array " << n << std::endl;
        }
        
    }
    
    if (fabs(thr_wsum_sigma2_offset - cl_wsum_sigma2_offset) < 0.0000000001) {
        //        std::cerr << "Right - wsum_sigma2_offset " << std::endl;
    } else {
        std::cerr << "Wrong - wsum_sigma2_offset, CPU: " << thr_wsum_sigma2_offset << " CL: " << cl_wsum_sigma2_offset << std::endl;
    }
    
    for (int n = 0; n < exp_nr_particles; n++)
    {
        if (fabs(thr_max_weight[n] - cl_max_weight[n]) < 0.0000000001) {
            //            std::cerr << "Right max_weight " << n << std::endl;
        } else {
            std::cerr << "Wrong max_weight " << n << " CPU: " << thr_max_weight[n] << " CL: " << cl_max_weight[n] << std::endl;
        }
        
        long int my_image_no = exp_starting_image_no[n] + exp_iseries;
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_ROT) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_ROT)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_ROT " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_ROT " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_ROT) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_ROT) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_TILT) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_TILT)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_TILT " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_TILT " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_TILT) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_TILT) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PSI) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PSI)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_PSI " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_PSI " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PSI) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PSI) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_XOFF) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_XOFF)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_XOFF " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_XOFF " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_XOFF) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_XOFF) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_YOFF) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_YOFF)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_YOFF " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_YOFF " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_YOFF) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_YOFF) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_CLASS) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_CLASS)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_CLASS " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_CLASS " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_CLASS) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_CLASS) << std::endl;
        }
        if (fabs(DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PMAX) - DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PMAX)) < 0.0000000001) {
            //            std::cerr << "Right METADATA_PMAX " << n << std::endl;
        } else {
            std::cerr << "Wrong METADATA_PMAX " << n << " CPU: " << DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PMAX) << " CL: " << DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PMAX) << std::endl;
        }
    }
    
#endif
    
#pragma mark Copy data to final place
    
    // Now, inside a global_mutex, update the weighted sums among all threads
    //    global_mutex.lock();
    
    
    if (!do_skip_maximization)
    {
        void *modelData = (void *)(wsum_model.BPref[exp_iclass]).data.data;
        void *modelWeight = (void *)(wsum_model.BPref[exp_iclass]).weight.data;
        memcpy(modelData, modelValues, sizeof(cl_double2) * model_size);
        memcpy(modelWeight, modelWeights, sizeof(cl_double) * model_size);
        
        
        if (do_scale_correction)
        {
            for (int n = 0; n < exp_nr_particles; n++)
            {
                exp_wsum_scale_correction_XA[n] += cl_wsum_scale_correction_XA[n];
                exp_wsum_scale_correction_AA[n] += cl_wsum_scale_correction_AA[n];
            }
        }
        for (int n = 0; n < exp_nr_particles; n++)
        {
            exp_wsum_norm_correction[n] += cl_wsum_norm_correction[n];
        }
        for (int n = 0; n < mymodel.nr_groups; n++)
        {
            wsum_model.sigma2_noise[n] += cl_wsum_sigma2_noise[n];
            wsum_model.sumw_group[n] += cl_sumw_group[n];
        }
        for (int n = 0; n < mymodel.nr_classes; n++)
        {
            wsum_model.pdf_class[n] += cl_wsum_pdf_class[n];
            
            if (mymodel.ref_dim == 2)
            {
                XX(wsum_model.prior_offset_class[n]) += cl_wsum_prior_offsetx_class[n];
                YY(wsum_model.prior_offset_class[n]) += cl_wsum_prior_offsety_class[n];
            }
            wsum_model.pdf_direction[n] += cl_wsum_pdf_direction[n];
        }
        wsum_model.sigma2_offset += cl_wsum_sigma2_offset;
    } // end if !do_skip_maximization
    
    // Check max_weight for each particle and set exp_metadata
    for (int n = 0; n < exp_nr_particles; n++)
    {
        // Equal-to because of the series: the nth images in a series will have the same maximum as the first one
        if (cl_max_weight[n] >= exp_max_weight[n])
        {
            // Set max_weight
            exp_max_weight[n] = cl_max_weight[n];
            
            // Set metadata
            long int my_image_no = exp_starting_image_no[n] + exp_iseries;
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT)  = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_ROT);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT) = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_TILT);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI)  = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PSI);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF) = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_XOFF);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF) = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_YOFF);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS)= DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_CLASS);
            DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX) = DIRECT_A2D_ELEM(cl_metadata, my_image_no, METADATA_PMAX);
        }
    }
    
    //    global_mutex.unlock();
    
    for (int i = 0; i < partCount; i++) {
        clReleaseEvent(clEvents[i]);
    }
    clReleaseMemObject(CL_iter_Proj);
    clReleaseMemObject(cl_modelData);
    clReleaseMemObject(cl_Ainv);
    clReleaseMemObject(cl_eulerAngles);
    clReleaseMemObject(cl_eulerAnglesHost);
    clReleaseMemObject(cl_exp_R_mic);
    clReleaseMemObject(cl_local_Fctfs);
    clReleaseMemObject(cl_FimgShiftAll);
    clReleaseMemObject(cl_FimgShiftAllHost);
    clReleaseMemObject(cl_FimgShiftNoMaskAll);
    clReleaseMemObject(cl_FimgShiftNoMaskAllHost);
    clReleaseMemObject(cl_FrefCTF);
    clReleaseMemObject(cl_Minvsigma2);
    clReleaseMemObject(cl_Mresol_fine);
    clReleaseMemObject(cl_sigma2_noise);
    clReleaseMemObject(cl_dataVsPrior);
    clReleaseMemObject(cl_sumXA);
    clReleaseMemObject(cl_sumAA);
    clReleaseMemObject(cl_FimgOut);
    clReleaseMemObject(cl_FweightOut);
    //    clReleaseMemObject(cl_outputModelWeights);
    clReleaseMemObject(cl_modelValues);
    clReleaseMemObject(cl_modelWeights);
    clReleaseMemObject(cl_Mctfs);
    clReleaseMemObject(cl_weights);
    clReleaseMemObject(cl_weightsHost);
    clReleaseMemObject(cl_sigma2_noise2);
    clReleaseMemObject(cl_sumXA2);
    clReleaseMemObject(cl_sumAA2);
    if (CL_atomicSupport == false) {
        clReleaseMemObject(cl_outputModelValues);
    }
    clReleaseMemObject(cl_oversampled_translations);
    clReleaseMemObject(cl_oversampled_translationsHost);
    clReleaseMemObject(cl_sigma2_offset_return);
    if (mymodel.ref_dim == 2) {
        clReleaseMemObject(cl_prior_offset_return);
    }
    
    free(scales);
    free(myImageArray);
    free(FimgShiftArray);
    free(iorientClassArray);
    free(weights);
    free(groupIdArray);
    free(modelWeights);
    free(modelValues);
    
#ifdef CL_VERIFY_ON_CPU
    free(FimgOut);
    free(FweightOut);    
#endif
    
#ifdef TIMING
    timer.toc(TIMING_WSUM_BACKPROJ);
#endif
    
#ifdef DEBUG_THREAD
    std::cerr << "leaving doOpenCLStoreWeightedSumsAllOrientations" << std::endl;
#endif

}

void MlOptimiser::doOpenCLStoreWeightedSumsAllOrientations()
{
#ifdef DEBUG_THREAD
    std::cerr << "entering doOpenCLStoreWeightedSumsAllOrientations" << std::endl;
#endif
    
    MultidimArray<Complex > Fref;
    
    // Initialising...
    Fref.resize(exp_Fimgs[0]);
    
    long int nr_orients = sampling.NrDirections() * sampling.NrPsiSamplings();
    long int iorientclass_offset = exp_iclass * exp_nr_rot;
    
    //Get number of orients for openCL save data...
    std::vector<bool> translationUsedInAnyRotation;
    translationUsedInAnyRotation.resize(exp_nr_trans * exp_nr_oversampled_trans, false);
    
    int maxOrients = CL_maxMemAlloc / (sizeof(cl_double2) * XSIZE(Fref) * YSIZE(Fref));

    std::vector<int>startIOrients, endIOrients;
    startIOrients.push_back(0);
    
    int currentCalcIter = 1;
    long int totalOrients = 0;
    for (long int iorient = 0; iorient < nr_orients; iorient++)
    {
        
        long int iorientclass = iorientclass_offset + iorient;
        
        // Only proceed if any of the particles had any significant coarsely sampled translation
        if (isSignificantAnyParticleAnyTranslation(iorientclass))
        {
            // Loop over all oversampled orientations (only a single one in the first pass)
            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
            {
                bool anySignificantWeightTranslations = false;
                /// Now that reference projection has been made loop over someParticles!
                for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
                {
                    // loop over all particles inside this ori_particle
                    for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
                    {
                        long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
                        
                        // Which number was this image in the combined array of iseries and part_idpart_id
                        long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
                        int group_id = mydata.getGroupId(part_id, exp_iseries);
                        
                        
                        long int ihidden = iorientclass * exp_nr_trans;
                        for (long int itrans = 0, currentTrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
                        {
                            
                            for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
                            {
                                
                                // Only deal with this sampling point if its weight was significant
                                long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                iover_rot * exp_nr_oversampled_trans + iover_trans;
                                
                                double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
                                if (weight >= exp_significant_weight[ipart]) {
                                    anySignificantWeightTranslations = true;
                                    translationUsedInAnyRotation[currentTrans] = true;
                                }
                                
                                currentTrans++;
                                                                 
                            } // end loop iover_trans
                        } // end loop itrans
                    }// end loop part_id (i)
                } // end loop ori_part_id
                if (anySignificantWeightTranslations) {
                    totalOrients++;
                    if (totalOrients > currentCalcIter * maxOrients) {
                        currentCalcIter++;
                        endIOrients.push_back(iorient);
                        startIOrients.push_back(iorient);
                    }
                }
            }// end if iover_rot
        }// end loop do_proceed
    }
    
    endIOrients.push_back(nr_orients);
    
//    std::cerr << "Total orients: " << totalOrients << " maxOrients: " << maxOrients << std::endl;
    if (totalOrients == 0) {
        return;
    } else {
        for (int i = 0; i < currentCalcIter; i++) {
//            std::cerr << "Iteration " << i << " start ori: " << startIOrients[i] << " end ori: " << endIOrients[i] << std::endl;
            doOpenCLStoreWeightedSumsAllOrientationsCalculations(startIOrients[i], endIOrients[i]);
        }
    }

#ifdef DEBUG_THREAD
    std::cerr << "leaving doOpenCLStoreWeightedSumsAllOrientations" << std::endl;
#endif
    
}

void MlOptimiser::doThreadStoreWeightedSumsAllOrientations(int thread_id)
{
#ifdef DEBUG_THREAD
    std::cerr << "entering doThreadStoreWeightedSumsAllOrientations" << std::endl;
#endif

	std::vector< Matrix1D<double> > oversampled_orientations, oversampled_translations;
	Matrix2D<double> A;
	MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_shift, Fimg_shift_nomask;
	MultidimArray<double> Minvsigma2, Mctf, Fweight;
	double rot, tilt, psi;
	bool have_warned_small_scale = false;

	// Initialising...
	Fref.resize(exp_Fimgs[0]);
	Frefctf.resize(exp_Fimgs[0]);
	Fweight.resize(exp_Fimgs[0]);

	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(exp_Fimgs[0]);
	Mctf.initConstant(1.);

	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(exp_Fimgs[0]);
	Minvsigma2.initConstant(1.);

	// Make local copies of weighted sums (excepts BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_scale_correction_XA, thr_wsum_scale_correction_AA, thr_wsum_pdf_direction;
	std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class, thr_max_weight;
	double thr_wsum_sigma2_offset;
	MultidimArray<double> thr_metadata, zeroArray;

	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(mymodel.nr_groups);
	for (int n = 0; n < mymodel.nr_groups; n++)
	{
		thr_wsum_sigma2_noise[n] = zeroArray;
	}
	// scale-correction terms are a spectrum for each particle
	thr_wsum_scale_correction_XA.resize(exp_nr_particles);
	thr_wsum_scale_correction_AA.resize(exp_nr_particles);
	for (int n = 0; n < exp_nr_particles; n++)
	{
		thr_wsum_scale_correction_XA[n] = zeroArray;
		thr_wsum_scale_correction_AA[n] = zeroArray;
	}
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections(0, true)) for each class
	zeroArray.initZeros(sampling.NrDirections(0, true));
	thr_wsum_pdf_direction.resize(mymodel.nr_classes);
	for (int n = 0; n < mymodel.nr_classes; n++)
	{
		thr_wsum_pdf_direction[n] = zeroArray;
	}
	// wsum_norm_correction is a double for each particle
	thr_wsum_norm_correction.resize(exp_nr_particles, 0.);
	// sumw_group is a double for each group
	thr_sumw_group.resize(mymodel.nr_groups, 0.);
	// wsum_pdf_class is a double for each class
	thr_wsum_pdf_class.resize(mymodel.nr_classes, 0.);
	if (mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.);
	}
	// max_weight is a double for each particle
	thr_max_weight.resize(exp_nr_particles, 0.);
	// wsum_sigma2_offset is just a double
	thr_wsum_sigma2_offset = 0.;
	// metadata is a 2D array of nr_particles x METADATA_LINE_LENGTH
	thr_metadata.initZeros(exp_metadata);

	// exp_iclass loop does not always go from 0 to nr_classes!
	long int iorientclass_offset = exp_iclass * exp_nr_rot;
	size_t first_iorient = 0, last_iorient = 0;
	while (exp_iorient_ThreadTaskDistributor->getTasks(first_iorient, last_iorient))
	{
		for (long int iorient = first_iorient; iorient <= last_iorient; iorient++)
		{

			long int iorientclass = iorientclass_offset + iorient;

			// Only proceed if any of the particles had any significant coarsely sampled translation
			if (isSignificantAnyParticleAnyTranslation(iorientclass))
			{

				long int idir = iorient / exp_nr_psi;
				long int ipsi = iorient % exp_nr_psi;

				// Now get the oversampled (rot, tilt, psi) triplets
				// This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
				sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_orientations);

				// Loop over all oversampled orientations (only a single one in the first pass)
				for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
				{
					rot = XX(oversampled_orientations[iover_rot]);
					tilt = YY(oversampled_orientations[iover_rot]);
					psi = ZZ(oversampled_orientations[iover_rot]);
					// Get the Euler matrix
					Euler_angles2matrix(rot, tilt, psi, A);

					// Take tilt-series into account
					A = (exp_R_mic * A).inv();

#ifdef TIMING
					// Only time one thread, as I also only time one MPI process
					if (thread_id == 0)
						timer.tic(TIMING_WSUM_PROJ);
#endif
					// Project the reference map (into Fref)
					if (!do_skip_maximization)
						(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_INV);

#ifdef TIMING
					// Only time one thread, as I also only time one MPI process
					if (thread_id == 0)
						timer.toc(TIMING_WSUM_PROJ);
#endif
					// Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
					// Then outside this loop do the actual backprojection
					Fimg.initZeros(Fref);
					Fweight.initZeros(Fref);

					/// Now that reference projection has been made loop over someParticles!
					for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
					{
						// loop over all particles inside this ori_particle
						for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
						{
							long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
#ifdef DEBUG_CHECKSIZES
							if (ipart >= exp_starting_image_no.size())
							{
								std::cerr<< "ipart= "<<ipart<<" starting_image_no.size()= "<< exp_starting_image_no.size() <<std::endl;
								REPORT_ERROR("ipart >= starting_image_no.size()");
							}
#endif
							// Which number was this image in the combined array of iseries and part_idpart_id
							long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
							int group_id = mydata.getGroupId(part_id, exp_iseries);

#ifdef DEBUG_CHECKSIZES
							if (group_id >= mymodel.nr_groups)
							{
								std::cerr<< "group_id= "<<group_id<<" ml_model.nr_groups= "<< mymodel.nr_groups <<std::endl;
								REPORT_ERROR("group_id >= ml_model.nr_groups");
							}
#endif

							if (!do_skip_maximization)
							{
								if (do_map)
									Minvsigma2 = exp_local_Minvsigma2s[my_image_no];
								// else Minvsigma2 was initialised to ones

								// Apply CTF to reference projection
								if (do_ctf_correction)
								{
									Mctf = exp_local_Fctfs[my_image_no];
									if (refs_are_ctf_corrected)
									{
										FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
										{
											DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(Mctf, n);
										}
									}
									else
									{
										Frefctf = Fref;
									}
								}
								else
								{
									// initialise because there are multiple particles and Mctf gets selfMultiplied for scale_correction
									Mctf.initConstant(1.);
									Frefctf = Fref;
								}

								if (do_scale_correction)
								{
									// TODO: implemenent B-factor as well...
									double myscale = mymodel.scale_correction[group_id];
									if (myscale > 10000.)
									{
										std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << " my_image_no= " << my_image_no << std::endl;
										REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
									}
									else if (myscale < 0.001)
									{

										if (!have_warned_small_scale)
										{
											std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << myscale <<
													"); Use larger groups for more stable scale estimates." << std::endl;
											have_warned_small_scale = true;
										}
										myscale = 0.001;
									}
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
									{
										DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
									}
									// For CTF-terms in BP
									Mctf *= myscale;
								}
							} // end if !do_skip_maximization

							long int ihidden = iorientclass * exp_nr_trans;
							for (long int itrans = 0; itrans < exp_nr_trans; itrans++, ihidden++)
							{

								sampling.getTranslations(itrans, adaptive_oversampling, oversampled_translations);

								for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
								{

#ifdef DEBUG_CHECKSIZES
									if (iover_trans >= oversampled_translations.size())
									{
										std::cerr<< "iover_trans= "<<iover_trans<<" oversampled_translations.size()= "<< oversampled_translations.size() <<std::endl;
										REPORT_ERROR("iover_trans >= oversampled_translations.size()");
									}
#endif

									// Only deal with this sampling point if its weight was significant
									long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
											iover_rot * exp_nr_oversampled_trans + iover_trans;

#ifdef DEBUG_CHECKSIZES
									if (ihidden_over >= XSIZE(exp_Mweight))
									{
										std::cerr<< "ihidden_over= "<<ihidden_over<<" XSIZE(exp_Mweight)= "<< XSIZE(exp_Mweight) <<std::endl;
										REPORT_ERROR("ihidden_over >= XSIZE(exp_Mweight)");
									}
									if (ipart >= exp_significant_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_significant_weight.size()= "<< exp_significant_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= significant_weight.size()");
									}
									if (ipart >= exp_max_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_max_weight.size()= "<< exp_max_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= exp_max_weight.size()");
									}
									if (ipart >= exp_sum_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_max_weight.size()= "<< exp_sum_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= exp_sum_weight.size()");
									}
#endif
									double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);

									// Only sum weights for non-zero weights
									if (weight >= exp_significant_weight[ipart])
									{
#ifdef TIMING
										// Only time one thread, as I also only time one MPI process
										if (thread_id == 0)
											timer.tic(TIMING_WSUM_DIFF2);
#endif
										// Normalise the weight (do this after the comparison with exp_significant_weight!)
										weight /= exp_sum_weight[ipart];

										if (!do_skip_maximization)
										{
											// Get the shifted image
											long int ishift = my_image_no * exp_nr_oversampled_trans * exp_nr_trans +
													itrans * exp_nr_oversampled_trans + iover_trans;
											Fimg_shift = exp_local_Fimgs_shifted[ishift];
											Fimg_shift_nomask = exp_local_Fimgs_shifted_nomask[ishift];

											// Store weighted sum of squared differences for sigma2_noise estimation
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
											{
												int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
												if (ires > -1)
												{
													// Use FT of masked image for noise estimation!
													double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
													double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
													double wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);

													// group-wise sigma2_noise
													DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[group_id], ires) += wdiff2;
													// For norm_correction
													thr_wsum_norm_correction[ipart] += wdiff2;
												}
											}

											// Store the weighted sums of the norm_correction terms
											if (do_scale_correction)
											{
												double sumXA = 0.;
												double sumA2 = 0.;
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
												{
													int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
#ifdef DEBUG_CHECKSIZES
													if (ires >= XSIZE(thr_wsum_scale_correction_XA[ipart]))
													{
														std::cerr<< "ires= "<<ires<<" XSIZE(thr_wsum_scale_correction_XA[ipart])= "<< XSIZE(thr_wsum_scale_correction_XA[ipart]) <<std::endl;
														REPORT_ERROR("ires >= XSIZE(thr_wsum_scale_correction_XA[ipart])");
													}
													if (ires >= XSIZE(thr_wsum_scale_correction_AA[ipart]))
													{
														std::cerr<< "ires= "<<ires<<" XSIZE(thr_wsum_scale_correction_AA[ipart])= "<< XSIZE(thr_wsum_scale_correction_AA[ipart]) <<std::endl;
														REPORT_ERROR("ires >= XSIZE(thr_wsum_scale_correction_AA[ipart])");
													}
#endif

													// Once the reference becomes strongly regularised one does no longer want to store XA and AA!
													if (ires > -1 && DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
													{
														sumXA = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).real;
														sumXA += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Fimg_shift, n)).imag;
														DIRECT_A1D_ELEM(thr_wsum_scale_correction_XA[ipart], ires) += weight * sumXA;

														// This could be pre-calculated above...
														sumA2 = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
														sumA2 += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
														DIRECT_A1D_ELEM(thr_wsum_scale_correction_AA[ipart], ires) += weight * sumA2;
													}
												}
											}

											// Store sum of weights for this group
											thr_sumw_group[group_id] += weight;

											// Store weights for this class and orientation
											thr_wsum_pdf_class[exp_iclass] += weight;

											if (mymodel.ref_dim ==2)
											{
												// Also store weighted offset differences for prior_offsets of each class
												thr_wsum_prior_offsetx_class[exp_iclass] += weight * XX(exp_old_offset[my_image_no] + oversampled_translations[iover_trans]);
												thr_wsum_prior_offsety_class[exp_iclass] += weight * YY(exp_old_offset[my_image_no] + oversampled_translations[iover_trans]);

												// Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
												thr_wsum_sigma2_offset += weight * ((mymodel.prior_offset_class[exp_iclass] - exp_old_offset[my_image_no] - oversampled_translations[iover_trans]).sum2());

											}
											else
											{
												// Store weighted sum2 of origin offsets (in Angstroms instead of pixels!!!)
												thr_wsum_sigma2_offset += weight * ((exp_prior[my_image_no] - exp_old_offset[my_image_no] - oversampled_translations[iover_trans]).sum2());
											}

#ifdef DEBUG_CHECKSIZES
											if (idir >= XSIZE(thr_wsum_pdf_direction[exp_iclass]))
											{
												std::cerr<< "idir= "<<idir<<" XSIZE(thr_wsum_pdf_direction[exp_iclass])= "<< XSIZE(thr_wsum_pdf_direction[exp_iclass]) <<std::endl;
												REPORT_ERROR("idir >= XSIZE(thr_wsum_pdf_direction[iclass])");
											}
#endif

											// Store weight for this direction of this class
											if (mymodel.orientational_prior_mode == NOPRIOR)
											{
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
											}
											else
											{
												// In the case of orientational priors, get the original number of the direction back
												long int mydir = sampling.getDirectionNumberAlsoZeroPrior(idir);
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
											}

#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (thread_id == 0)
												timer.toc(TIMING_WSUM_DIFF2);
											// Only time one thread, as I also only time one MPI process
											if (thread_id == 0)
												timer.tic(TIMING_WSUM_SUMSHIFT);
#endif

											// Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
											// Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg_shift)
											{
												double myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
												// Note that weightxinvsigma2 already contains the CTF!
												double weightxinvsigma2 = weight * myctf * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
												// now Fimg stores sum of all shifted w*Fimg
												(DIRECT_MULTIDIM_ELEM(Fimg, n)).real += (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).real * weightxinvsigma2;
												(DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += (DIRECT_MULTIDIM_ELEM(Fimg_shift_nomask, n)).imag * weightxinvsigma2;
												// now Fweight stores sum of all w
												// Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
												DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf;
											}

#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (thread_id == 0)
												timer.toc(TIMING_WSUM_SUMSHIFT);
#endif

										} // end if !do_skip_maximization

										// Keep track of max_weight and the corresponding optimal hidden variables
										if (weight > thr_max_weight[ipart])
										{
											// Store optimal image parameters
											thr_max_weight[ipart] = weight;

											// Calculate the angles back from the Euler matrix because for tilt series exp_R_mic may have changed them...
											//std::cerr << " ORI rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
											Euler_matrix2angles(A.inv(), rot, tilt, psi);
											//std::cerr << " BACK rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;

											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_ROT) = rot;
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_TILT) = tilt;
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PSI) = psi;
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_XOFF) = XX(exp_old_offset[my_image_no]) + XX(oversampled_translations[iover_trans]);
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_YOFF) = YY(exp_old_offset[my_image_no]) + YY(oversampled_translations[iover_trans]);
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_CLASS) = (double)exp_iclass + 1;
											DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PMAX) = thr_max_weight[ipart];
										}

									} // end if weight >= exp_significant_weight
								} // end loop iover_trans
							} // end loop itrans
						}// end loop part_id (i)
					} // end loop ori_part_id

					if (!do_skip_maximization)
					{
#ifdef TIMING
						// Only time one thread, as I also only time one MPI process
						if (thread_id == 0)
							timer.tic(TIMING_WSUM_BACKPROJ);
#endif
						// Perform the actual back-projection.
						// This is done with the sum of all (in-plane) shifted Fimg's
						// Perform this inside a mutex
						global_mutex2.lock();
						(wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_INV, &Fweight);
						global_mutex2.unlock();

#ifdef TIMING
						// Only time one thread, as I also only time one MPI process
						if (thread_id == 0)
							timer.toc(TIMING_WSUM_BACKPROJ);
#endif
					} // end if !do_skip_maximization

				}// end if iover_rot
			}// end loop do_proceed

		} // end loop ipsi
	} // end loop idir

	// Now, inside a global_mutex, update the weighted sums among all threads
	global_mutex.lock();

	if (!do_skip_maximization)
	{
		if (do_scale_correction)
		{
			for (int n = 0; n < exp_nr_particles; n++)
			{
				exp_wsum_scale_correction_XA[n] += thr_wsum_scale_correction_XA[n];
				exp_wsum_scale_correction_AA[n] += thr_wsum_scale_correction_AA[n];
			}
		}
		for (int n = 0; n < exp_nr_particles; n++)
		{
			exp_wsum_norm_correction[n] += thr_wsum_norm_correction[n];
		}
		for (int n = 0; n < mymodel.nr_groups; n++)
		{
			wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			wsum_model.sumw_group[n] += thr_sumw_group[n];
		}
		for (int n = 0; n < mymodel.nr_classes; n++)
		{
			wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];

			if (mymodel.ref_dim == 2)
			{
				XX(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
				YY(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
			}
#ifdef CHECKSIZES
			if (XSIZE(wsum_model.pdf_direction[n]) != XSIZE(thr_wsum_pdf_direction[n]))
			{
				std::cerr << " XSIZE(wsum_model.pdf_direction[n])= " << XSIZE(wsum_model.pdf_direction[n]) << " XSIZE(thr_wsum_pdf_direction[n])= " << XSIZE(thr_wsum_pdf_direction[n]) << std::endl;
				REPORT_ERROR("XSIZE(wsum_model.pdf_direction[n]) != XSIZE(thr_wsum_pdf_direction[n])");
			}
#endif
			wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
	} // end if !do_skip_maximization

	// Check max_weight for each particle and set exp_metadata
	for (int n = 0; n < exp_nr_particles; n++)
	{
		// Equal-to because of the series: the nth images in a series will have the same maximum as the first one
		if (thr_max_weight[n] >= exp_max_weight[n])
		{
			// Set max_weight
			exp_max_weight[n] = thr_max_weight[n];

			// Set metadata
			long int my_image_no = exp_starting_image_no[n] + exp_iseries;
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT)  = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_ROT);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT) = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_TILT);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI)  = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PSI);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF) = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_XOFF);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF) = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_YOFF);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS)= DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_CLASS);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX) = DIRECT_A2D_ELEM(thr_metadata, my_image_no, METADATA_PMAX);
		}
	}
	global_mutex.unlock();

	// Wait until all threads have finished
	global_barrier->wait();
#ifdef DEBUG_THREAD
    std::cerr << "leaving doThreadStoreWeightedSumsAllOrientations" << std::endl;
#endif

}

void MlOptimiser::storeWeightedSums()
{

#ifdef TIMING
	timer.tic(TIMING_ESP_WSUM);
#endif

	// Initialise the maximum of all weights to a negative value
	exp_max_weight.resize(exp_nr_particles);
	for (int n = 0; n < exp_nr_particles; n++)
		exp_max_weight[n] = -1.;

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
	{
		// loop over all particles inside this ori_particle
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

			for (exp_iseries = 0; exp_iseries < mydata.getNrImagesInSeries(part_id); exp_iseries++)
			{
				// Re-get all shifted versions of the (current_sized) images, their (current_sized) CTFs and their inverted Sigma2 matrices
				// This may be necessary for when using --strict_highres_exp. Otherwise norm estimation may become unstable!!
                if (do_use_opencl) {
                    doOpenCLPrecalculateShiftedImagesCtfsAndInvSigma2s();
                } else {
                    exp_ipart_ThreadTaskDistributor->reset(); // reset thread distribution tasks
                    global_ThreadManager->run(globalThreadPrecalculateShiftedImagesCtfsAndInvSigma2s);
                }

				int group_id = mydata.getGroupId(part_id, exp_iseries);
				int my_image_no = exp_starting_image_no[ipart] + exp_iseries;
				DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[my_image_no], 0) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], 0));
			}
		}
	}

	for (exp_iseries = 0; exp_iseries < mydata.getNrImagesInSeries((mydata.ori_particles[exp_my_first_ori_particle]).particles_id[0]); exp_iseries++)
	{
		// TODO: check this!!!
		// I think this is just done for the first ipart
		int my_image_no = exp_starting_image_no[0] + exp_iseries;
		// Get micrograph transformation matrix
		exp_R_mic.resize(3,3);
		exp_R_mic(0,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_0);
		exp_R_mic(0,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_1);
		exp_R_mic(0,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_2);
		exp_R_mic(1,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_0);
		exp_R_mic(1,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_1);
		exp_R_mic(1,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_2);
		exp_R_mic(2,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_0);
		exp_R_mic(2,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_1);
		exp_R_mic(2,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_2);

		// For norm_correction of this iseries image:
		exp_wsum_norm_correction.resize(exp_nr_particles);
		for (int n = 0; n < exp_nr_particles; n++)
			exp_wsum_norm_correction[n] = 0.;

		// For scale_correction of this iseries image:
		if (do_scale_correction)
		{
			MultidimArray<double> aux;
			aux.initZeros(mymodel.ori_size/2 + 1);
			exp_wsum_scale_correction_XA.resize(exp_nr_particles);
			exp_wsum_scale_correction_AA.resize(exp_nr_particles);
			for (int n = 0; n < exp_nr_particles; n++)
			{
				exp_wsum_scale_correction_XA[n] = aux;
				exp_wsum_scale_correction_AA[n] = aux;
			}
		}

		// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
		for (exp_iclass = iclass_min; exp_iclass <= iclass_max; exp_iclass++)
		{

			// The loops over all orientations are parallelised using threads
            if (do_use_opencl) {
//                struct timeval curTV;
//                gettimeofday(&curTV, NULL);
//                std::cerr << "Enter sws time: " << curTV.tv_sec << "." << curTV.tv_usec << std::endl;
                doOpenCLStoreWeightedSumsAllOrientations();
//                gettimeofday(&curTV, NULL);
//                std::cerr << "Exit sws time: " << curTV.tv_sec << "." << curTV.tv_usec << std::endl;
            } else {
                exp_iorient_ThreadTaskDistributor->reset(); // reset thread distribution tasks
                global_ThreadManager->run(globalThreadStoreWeightedSumsAllOrientations);
            }

		} // end loop iclass

		// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
		for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
		{
			// loop over all particles inside this ori_particle
			for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
			{
				long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

				// Which number was this image in the combined array of exp_iseries and exp_part_id
				long int my_image_no = exp_starting_image_no[ipart] + exp_iseries;

				// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
				int group_id = mydata.getGroupId(part_id, exp_iseries);
				for (int ires = mymodel.current_size/2 + 1; ires < mymodel.ori_size/2 + 1; ires++)
				{
					DIRECT_A1D_ELEM(wsum_model.sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(exp_power_imgs[my_image_no], ires);
					// Also extend the weighted sum of the norm_correction
					exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(exp_power_imgs[my_image_no], ires);
				}

				// Store norm_correction
				// Multiply by old value because the old norm_correction term was already applied to the image
				if (do_norm_correction)
				{
					double old_norm_correction = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM);
					old_norm_correction /= mymodel.avg_norm_correction;
					// Now set the new norm_correction in the relevant position of exp_metadata
					// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
					// The variance of the total image (on which one normalizes) is twice this value!
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
					wsum_model.avg_norm_correction += old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);

					if (!(iter == 1 && do_firstiter_cc) && DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) > 10.)
					{
						std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) << " for particle " << part_id << " in group " << group_id + 1 << "; Are your groups large enough?" << std::endl;
						std::cout << " mymodel.current_size= " << mymodel.current_size << " mymodel.ori_size= " << mymodel.ori_size << " part_id= " << part_id << std::endl;
						std::cout << " coarse_size= " << coarse_size << std::endl;
						std::cout << " DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM)= " << DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) << std::endl;
						std::cout << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
						std::cout << " exp_wsum_norm_correction[ipart]= " << exp_wsum_norm_correction[ipart] << std::endl;
						std::cout << " old_norm_correction= " << old_norm_correction << std::endl;
						std::cout << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
						std::cout << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
						std::cout << " mymodel.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
						std::cout << " wsum_model.sigma2_noise[group_id]= " << wsum_model.sigma2_noise[group_id] << std::endl;
						std::cout << " exp_power_imgs[my_image_no]= " << exp_power_imgs[my_image_no] << std::endl;
						std::cout << " exp_wsum_scale_correction_XA[ipart]= " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart]= " << exp_wsum_scale_correction_AA[ipart] << std::endl;
						std::cout << " wsum_model.wsum_signal_product_spectra[group_id]= " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id]= " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
						std::cout << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
						std::cout << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
						std::cout << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
						std::cout << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;

					}
					//TMP DEBUGGING
					/*
					if (!(iter == 1 && do_firstiter_cc) && DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) > 10.)
					{
						std::cerr << " mymodel.current_size= " << mymodel.current_size << " mymodel.ori_size= " << mymodel.ori_size << " part_id= " << part_id << std::endl;
						std::cerr << " DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM)= " << DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) << std::endl;
						std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
						std::cerr << " exp_wsum_norm_correction[ipart]= " << exp_wsum_norm_correction[ipart] << std::endl;
						std::cerr << " old_norm_correction= " << old_norm_correction << std::endl;
						std::cerr << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
						std::cerr << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
						std::cerr << " mymodel.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
						std::cerr << " wsum_model.sigma2_noise[group_id]= " << wsum_model.sigma2_noise[group_id] << std::endl;
						std::cerr << " exp_power_imgs[my_image_no]= " << exp_power_imgs[my_image_no] << std::endl;
						std::cerr << " exp_wsum_scale_correction_XA[ipart]= " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart]= " << exp_wsum_scale_correction_AA[ipart] << std::endl;
						std::cerr << " wsum_model.wsum_signal_product_spectra[group_id]= " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id]= " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
						std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
						std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
						std::cerr << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
						std::cerr << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;
						mymodel.write("debug");
						std::cerr << "written debug_model.star" << std::endl;
						REPORT_ERROR("MlOptimiser::storeWeightedSums ERROR: normalization is larger than 10");
					}
					*/

				}

				// Store weighted sums for scale_correction
				if (do_scale_correction)
				{
					// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
					exp_wsum_scale_correction_XA[ipart] /= mymodel.scale_correction[group_id];
					exp_wsum_scale_correction_AA[ipart] /= mymodel.scale_correction[group_id] * mymodel.scale_correction[group_id];

					wsum_model.wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
					wsum_model.wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
				}

			} // end loop part_id (i)
		} // end loop ori_part_id
	} // end loop exp_iseries


#ifdef DEBUG_OVERSAMPLING
	std::cerr << " max_weight= " << max_weight << " nr_sign_sam= "<<nr_significant_samples<<" sign w= "<<exp_significant_weight<<std::endl;
#endif

	// Some analytics...
	// Calculate normalization constant for dLL
	for (long int ori_part_id = exp_my_first_ori_particle, ipart = 0; ori_part_id <= exp_my_last_ori_particle; ori_part_id++)
	{
		// loop over all particles inside this ori_particle
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			double logsigma2 = 0.;
			for (long int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++)
			{
				int group_id = mydata.getGroupId(part_id, iseries);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
				{
					int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
					// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
					// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
					if (ires > 0)
						logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
				}

			}

			if (exp_sum_weight[ipart]==0)
			{
				std::cerr << " part_id= " << part_id << std::endl;
				std::cerr << " ipart= " << ipart << std::endl;
				std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
				std::cerr << " logsigma2= " << logsigma2 << std::endl;
				int group_id = mydata.getGroupId(part_id, 0);
				std::cerr << " group_id= " << group_id << std::endl;
				std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
				std::cerr << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
				std::cerr << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;
				std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
				REPORT_ERROR("ERROR: exp_sum_weight[ipart]==0");
			}

			double dLL;

			if ((iter==1 && do_firstiter_cc) || do_always_cc)
				dLL = -exp_min_diff2[ipart];
			else
				dLL = log(exp_sum_weight[ipart]) - exp_min_diff2[ipart] - logsigma2;

			wsum_model.LL += dLL;
			wsum_model.ave_Pmax += DIRECT_A2D_ELEM(exp_metadata, exp_starting_image_no[ipart], METADATA_PMAX);

			// Also store dLL of each image in the output array
			for (long int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++)
			{
				long int my_image_no = exp_starting_image_no[ipart] + iseries;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL) = dLL;
			}
		} // end loop part_id
	} // end loop ori_part_id
#ifdef TIMING
	timer.toc(TIMING_ESP_WSUM);
#endif

}

/** Monitor the changes in the optimal translations, orientations and class assignments for some particles */
void MlOptimiser::monitorHiddenVariableChanges(long int my_first_ori_particle, long int my_last_ori_particle)
{

	for (long int ori_part_id = my_first_ori_particle, my_image_no = 0; ori_part_id <= my_last_ori_particle; ori_part_id++)
	{

#ifdef DEBUG_CHECKSIZES
		if (ori_part_id >= mydata.ori_particles.size())
		{
			std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
			REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
		}
#endif

		// loop over all particles inside this ori_particle
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++, my_image_no++)
			{
				long int img_id = mydata.getImageId(part_id, iseries);

#ifdef DEBUG_CHECKSIZES
				if (img_id >= mydata.MDimg.numberOfObjects())
				{
					std::cerr<< "img_id= "<<img_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
					REPORT_ERROR("img_id >= mydata.MDimg.numberOfObjects()");
				}
				if (my_image_no >= YSIZE(exp_metadata))
				{
					std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
					REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
				}
#endif

				// Old optimal parameters
				double old_rot, old_tilt, old_psi, old_xoff, old_yoff, old_zoff = 0.;
				int old_iclass;
				mydata.MDimg.getValue(EMDL_ORIENT_ROT,  old_rot, img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_TILT, old_tilt, img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_PSI,  old_psi, img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, old_xoff, img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, old_yoff, img_id);
				mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, old_iclass, img_id);

				// New optimal parameters
				double rot = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT);
				double tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT);
				double psi = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI);
				double xoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF);
				double yoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF);
				double zoff = 0.;
				int iclass = (int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS);

				// Some orientational distance....
				sum_changes_optimal_orientations += sampling.calculateAngularDistance(rot, tilt, psi, old_rot, old_tilt, old_psi);
				sum_changes_optimal_offsets += (xoff-old_xoff)*(xoff-old_xoff) + (yoff-old_yoff)*(yoff-old_yoff) + (zoff-old_zoff)*(zoff-old_zoff);
				if (iclass != old_iclass)
					sum_changes_optimal_classes += 1.;
				sum_changes_count += 1.;
			} // end loop iseries
		} // end loop part_id (i)
	} //end loop ori_part_id


}

void MlOptimiser::updateOverallChangesInHiddenVariables()
{

	// Calculate hidden variable changes
	current_changes_optimal_classes = sum_changes_optimal_classes / sum_changes_count;
	current_changes_optimal_orientations = sum_changes_optimal_orientations / sum_changes_count;
	current_changes_optimal_offsets = sqrt(sum_changes_optimal_offsets / (2. * sum_changes_count));

	// Reset the sums
	sum_changes_optimal_classes = 0.;
	sum_changes_optimal_orientations = 0.;
	sum_changes_optimal_offsets = 0.;
	sum_changes_count = 0.;

	// Update nr_iter_wo_large_hidden_variable_changes if all three assignment types are within 3% of the smallest thus far
	if (1.03 * current_changes_optimal_classes >= smallest_changes_optimal_classes &&
		1.03 * current_changes_optimal_offsets >= smallest_changes_optimal_offsets &&
		1.03 * current_changes_optimal_orientations >= smallest_changes_optimal_orientations)
		nr_iter_wo_large_hidden_variable_changes++;
	else
		nr_iter_wo_large_hidden_variable_changes = 0;

	// Update smallest changes in hidden variables thus far
	if (current_changes_optimal_classes < smallest_changes_optimal_classes)
		smallest_changes_optimal_classes = ROUND(current_changes_optimal_classes);
	if (current_changes_optimal_offsets < smallest_changes_optimal_offsets)
		smallest_changes_optimal_offsets = current_changes_optimal_offsets;
	if (current_changes_optimal_orientations < smallest_changes_optimal_orientations)
		smallest_changes_optimal_orientations = current_changes_optimal_orientations;


}


void MlOptimiser::calculateExpectedAngularErrors(long int my_first_ori_particle, long int my_last_ori_particle)
{

	long int n_trials = 0;
	exp_starting_image_no.clear();
	exp_nr_images = 0;
	for (long int ori_part_id = my_first_ori_particle, my_image_no = 0, ipart = 0; ori_part_id <= my_last_ori_particle; ori_part_id++)
    {
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			exp_starting_image_no.push_back(exp_nr_images);
			exp_nr_images += mydata.getNrImagesInSeries(part_id);
			n_trials++;
		}
    }

	// Set exp_current_image_size to the coarse_size to calculate exepcted angular errors
	if (strict_highres_exp > 0. && !do_acc_currentsize_despite_highres_exp)
	{
		// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
		exp_current_image_size = coarse_size;
	}
	else
	{
		// Use smaller images in the first pass, but larger ones in the second pass
		exp_current_image_size = mymodel.current_size;
	}

	// Separate angular error estimate for each of the classes
	acc_rot = acc_trans = 999.; // later XMIPP_MIN will be taken to find the best class...

	// P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
	// exp(-4.60517) = 0.01
	double pvalue = 4.60517;

	std::cout << " Estimating accuracies in the orientational assignment ... " << std::endl;
	init_progress_bar(n_trials * mymodel.nr_classes);
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{

		// Don't do this for (almost) empty classes
		if (mymodel.pdf_class[iclass] < 0.01)
		{
			mymodel.acc_rot[iclass]   = 999.;
			mymodel.acc_trans[iclass] = 999.;
			continue;
		}

		// Initialise the orientability arrays that will be written out in the model.star file
		// These are for the user's information only: nothing will be actually done with them
#ifdef DEBUG_CHECKSIZES
		if (iclass >= (mymodel.orientability_contrib).size())
		{
			std::cerr<< "iclass= "<<iclass<<" (mymodel.orientability_contrib).size()= "<< (mymodel.orientability_contrib).size() <<std::endl;
			REPORT_ERROR("iclass >= (mymodel.orientability_contrib).size()");
		}
#endif
		(mymodel.orientability_contrib)[iclass].initZeros(mymodel.ori_size/2 + 1);

		double acc_rot_class = 0.;
		double acc_trans_class = 0.;
		// Particles are already in random order, so just move from 0 to n_trials
		for (long int ori_part_id = my_first_ori_particle, my_image_no = 0, ipart = 0; ori_part_id <= my_last_ori_particle; ori_part_id++)
	    {
			for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++, ipart++)
			{
				long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

				// Search 2 times: ang and off
				// Don't estimate rotational accuracies if we're doing do_skip_rotate (for faster movie-frame alignment)
				int imode_start = (do_skip_rotate) ? 1 : 0;
				for (int imode = imode_start; imode < 2; imode++)
				{
					double ang_error = 0.;
					double sh_error = 0.;
					double ang_step;
					double sh_step;
					double my_snr = 0.;

					// Search for ang_error and sh_error where there are at least 3-sigma differences!
					// 13feb12: change for explicit probability at P=0.01
					while (my_snr <= pvalue)
					{
						// Graduallly increase the step size
						if (ang_error < 0.2)
							ang_step = 0.05;
						else if (ang_error < 1.)
							ang_step = 0.1;
						else if (ang_error < 2.)
							ang_step = 0.2;
						else if (ang_error < 5.)
							ang_step = 0.5;
						else if (ang_error < 10.)
							ang_step = 1.0;
						else if (ang_error < 20.)
							ang_step = 2;
						else
							ang_step = 5.0;

						if (sh_error < 0.2)
							sh_step = 0.05;
						else if (sh_error < 1.)
							sh_step = 0.1;
						else if (sh_error < 2.)
							sh_step = 0.2;
						else if (sh_error < 5.)
							sh_step = 0.5;
						else if (sh_error < 10.)
							sh_step = 1.0;
						else
							sh_step = 2.0;

						ang_error += ang_step;
						sh_error += sh_step;

						// Prevent an endless while by putting boundaries on ang_error and sh_error
						if ( (imode == 0 && ang_error > 30.) || (imode == 1 && sh_error > 10.) )
							break;

						init_random_generator(random_seed + part_id);

						// Loop over all images in the series
						// TODO: check this for series!!
						// Initialise the my_snr value (accumulate its sum for all images in the series!!)
						my_snr = 0.;
						for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++)
						{

							int my_image_no = exp_starting_image_no.at(ipart) + iseries;

							Matrix2D<double> R_mic(3,3);
							R_mic(0,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_0);
							R_mic(0,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_1);
							R_mic(0,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_2);
							R_mic(1,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_0);
							R_mic(1,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_1);
							R_mic(1,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_2);
							R_mic(2,0) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_0);
							R_mic(2,1) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_1);
							R_mic(2,2) = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_2);

							int group_id = mydata.getGroupId(part_id, iseries);
#ifdef DEBUG_CHECKSIZES
							if (group_id  >= mymodel.sigma2_noise.size())
							{
								std::cerr<< "group_id = "<<group_id <<" mymodel.sigma2_noise.size()= "<< mymodel.sigma2_noise.size() <<std::endl;
								REPORT_ERROR("group_id  >= mymodel.sigma2_noise.size()");
							}
#endif
							MultidimArray<Complex > F1, F2;
							MultidimArray<double> Fctf;
							Matrix2D<double> A1, A2;


							// TODO: get values through exp_metadata?!
							double rot1 = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT);
							double tilt1 = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT);
							double psi1 = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI);
							double xoff1 = 0.;
							double yoff1 = 0.;

							// Get the FT of the first image
							F1.initZeros(exp_current_image_size, exp_current_image_size/ 2 + 1);

							Euler_angles2matrix(rot1, tilt1, psi1, A1);
							A1 = R_mic * A1.inv();
							(mymodel.PPref[iclass]).get2DFourierTransform(F1, A1, IS_INV);

							// Apply the angular or shift error
							double rot2 = rot1;
							double tilt2 = tilt1;
							double psi2 = psi1;
							Matrix1D<double> shift(2);
							XX(shift) = xoff1;
							YY(shift) = yoff1;
							// Perturb psi or xoff , depending on the mode
							if (imode == 0)
							{
								if (mymodel.ref_dim == 3)
								{
									// Randomly change rot, tilt or psi
									double ran = rnd_unif();
									if (ran < 0.3333)
										rot2 = rot1 + ang_error;
									else if (ran < 0.6667)
										tilt2 = tilt1 + ang_error;
									else
										psi2  = psi1 + ang_error;
								}
								else
								{
									psi2  = psi1 + ang_error;
								}
							}
							else
							{
								// Randomly change xoff or yoff
								double ran = rnd_unif();
								if (ran < 0.5)
									XX(shift) = xoff1 + sh_error;
								else
									YY(shift) = yoff1 + sh_error;
							}
							// Get the FT of the second image
							F2.initZeros(exp_current_image_size, exp_current_image_size / 2 + 1);
							Euler_angles2matrix(rot2, tilt2, psi2, A2);
							A2 = R_mic * A2.inv();
							(mymodel.PPref[iclass]).get2DFourierTransform(F2, A2, IS_INV);
							if (ABS(XX(shift)) > 0. || ABS(YY(shift)) > 0.)
								// shiftImageInFourierTransform takes shifts in pixels!
								shiftImageInFourierTransform(F2, F2, (double) mymodel.ori_size, -shift);

							// Apply CTF to F1 and F2 if necessary
							if (do_ctf_correction)
							{
								CTF ctf;

								ctf.setValues(DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_U),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_V),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_ANGLE),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_VOLTAGE),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_CS),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_Q0),
											  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_BFAC));

								Fctf.resize(F1);
								ctf.getFftwImage(Fctf, mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
#ifdef DEBUG_CHECKSIZES
								if (!Fctf.sameShape(F1) || !Fctf.sameShape(F2))
								{
									std::cerr<<" Fctf: "; Fctf.printShape(std::cerr);
									std::cerr<<" F1:   "; F1.printShape(std::cerr);
									std::cerr<<" F2:   "; F2.printShape(std::cerr);
									REPORT_ERROR("ERROR: Fctf has a different shape from F1 and F2");
								}
#endif
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
								{
									DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
									DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								}
							}

							MultidimArray<int> * myMresol = (YSIZE(F1) == coarse_size) ? &Mresol_coarse : &Mresol_fine;
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
							{
								int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
								if (ires > 0)
								{
									my_snr += norm(DIRECT_MULTIDIM_ELEM(F1, n) - DIRECT_MULTIDIM_ELEM(F2, n)) / (2 * sigma2_fudge * mymodel.sigma2_noise[group_id](ires));
								}
							}

							// Only for the psi-angle and the translations, and only when my_prob < 0.01 calculate a histogram of the contributions at each resolution shell
							if (my_snr > pvalue && imode == 0)
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
								{
									int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
									if (ires > 0)
										mymodel.orientability_contrib[iclass](ires) +=
												norm(DIRECT_MULTIDIM_ELEM(F1, n) - DIRECT_MULTIDIM_ELEM(F2, n)) / ( (2 * sigma2_fudge * mymodel.sigma2_noise[group_id](ires)) );
								}

							}

						} // end for iseries

					} // end while my_snr >= pvalue
					if (imode == 0)
						acc_rot_class += ang_error;
					else if (imode == 1)
						acc_trans_class += sh_error;
				} // end for imode

			}// end for part_id

			progress_bar(n_trials*iclass + ipart);
		} // end for ori_part_id

		mymodel.acc_rot[iclass]   = acc_rot_class / (double)n_trials;
		mymodel.acc_trans[iclass] = acc_trans_class / (double)n_trials;

		// Store normalised spectral contributions to orientability
		if (mymodel.orientability_contrib[iclass].sum() > 0.)
			mymodel.orientability_contrib[iclass]   /= mymodel.orientability_contrib[iclass].sum();

		// Keep the orientational accuracy of the best class for the auto-sampling approach
		acc_rot     = XMIPP_MIN(mymodel.acc_rot[iclass], acc_rot);
		acc_trans   = XMIPP_MIN(mymodel.acc_trans[iclass], acc_trans);


		// Richard's formula with Greg's constant
		//double b_orient = (acc_rot_class*acc_rot_class* particle_diameter*particle_diameter) / 3000.;
		//std::cout << " + expected B-factor from the orientational errors = "
		//		<< b_orient<<std::endl;
		// B=8 PI^2 U^2
		//std::cout << " + expected B-factor from the translational errors = "
		//		<< 8 * PI * PI * mymodel.pixel_size * mymodel.pixel_size * acc_trans_class * acc_trans_class << std::endl;

	} // end loop iclass
	progress_bar(n_trials * mymodel.nr_classes);


	std::cout << " Auto-refine: Estimated accuracy angles= " << acc_rot<< " degrees; offsets= " << acc_trans << " pixels" << std::endl;
	// Warn for inflated resolution estimates
	if (acc_rot > 10.)
	{
		std::cout << " Auto-refine: WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles (yet)!" << std::endl;
		std::cout << " Auto-refine: WARNING: You probably need not worry if the accuracy improves during the next few iterations." << std::endl;
		std::cout << " Auto-refine: WARNING: However, if the problem persists it may lead to spurious FSC curves, so be wary of inflated resolution estimates..." << std::endl;
		std::cout << " Auto-refine: WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
	}

}

void MlOptimiser::updateAngularSampling(bool verb)
{

	if (!do_split_random_halves)
		REPORT_ERROR("MlOptimiser::updateAngularSampling: BUG! updating of angular sampling should only happen for gold-standard (auto-) refinements.");

	if (do_realign_movies)
	{

		// A. Adjust translational sampling to 75% of estimated accuracy
		double new_step = XMIPP_MIN(1.5, 0.75 * acc_trans) * std::pow(2., adaptive_oversampling);

		// Search ranges are three times the estimates std.dev. in the offsets
		double new_range = 3. * sqrt(mymodel.sigma2_offset);

		// Prevent too narrow searches: always at least 3x3 pixels in the coarse search
		if (new_range < 1.5 * new_step)
			new_range = 1.5 * new_step;

		// Also prevent too wide searches: that will lead to memory problems:
		// Just use coarser step size and hope things will settle down later...
		if (new_range > 4. * new_step)
			new_step = new_range / 4.;

		sampling.setTranslations(new_step, new_range);

		if (!do_skip_rotate)
		{
			// B. Find the healpix order that corresponds to at least 50% of the estimated rotational accuracy
			double angle_range = sqrt(mymodel.sigma2_rot) * 3.;
			double new_ang_step, new_ang_step_wo_over;
			int new_hp_order;
			for (new_hp_order = 0; new_hp_order < 8; new_hp_order++)
			{

				new_ang_step = 360. / (6 * ROUND(std::pow(2., new_hp_order + adaptive_oversampling)));
				new_ang_step_wo_over = 2. * new_ang_step;
				// Only consider healpix orders that gives at least more than one (non-oversampled) samplings within the local angular searches
				if (new_ang_step_wo_over > angle_range)
					continue;
				// If sampling is at least twice as fine as the estimated rotational accuracy, then use this sampling
				if (new_ang_step < 0.50 * acc_rot)
					break;
			}

			if (new_hp_order != sampling.healpix_order)
			{
				// Set the new sampling in the sampling-object
				sampling.setOrientations(new_hp_order, new_ang_step * std::pow(2., adaptive_oversampling));
				// Resize the pdf_direction arrays to the correct size and fill with an even distribution
				mymodel.initialisePdfDirection(sampling.NrDirections(0, true));
				// Also reset the nr_directions in wsum_model
				wsum_model.nr_directions = mymodel.nr_directions;
				// Also resize and initialise wsum_model.pdf_direction for each class!
				for (int iclass=0; iclass < mymodel.nr_classes; iclass++)
					wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);
			}
		}
	}
	else
	{

		if (do_skip_rotate)
			REPORT_ERROR("ERROR: --skip_rotate can only be used in classification or in movie-frame refinement ...");

		// Only change the sampling if the resolution has not improved during the last 2 iterations
		// AND the hidden variables have not changed during the last 2 iterations
		double old_rottilt_step = sampling.getAngularSampling(adaptive_oversampling);

		// Only use a finer angular sampling is the angular accuracy is still above 75% of the estimated accuracy
		// If it is already below, nothing will change and eventually nr_iter_wo_resol_gain or nr_iter_wo_large_hidden_variable_changes will go above MAX_NR_ITER_WO_RESOL_GAIN
		if (nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES)
		{
			// Old rottilt step is already below 75% of estimated accuracy: have to stop refinement
			if (old_rottilt_step < 0.75 * acc_rot)
			{
				// don't change angular sampling, as it is already fine enough
				has_fine_enough_angular_sampling = true;

			}
			else
			{
				has_fine_enough_angular_sampling = false;

				// A. Use translational sampling as suggested by acc_trans

				// Prevent very coarse translational samplings: max 1.5
				// Also stay a bit on the safe side with the translational sampling: 75% of estimated accuracy
				double new_step = XMIPP_MIN(1.5, 0.75 * acc_trans) * std::pow(2., adaptive_oversampling);
				// Search ranges are five times the last observed changes in offsets
				double new_range = 5. * current_changes_optimal_offsets;
				// New range can only become 30% bigger than the previous range (to prevent very slow iterations in the beginning)
				new_range = XMIPP_MIN(1.3*sampling.offset_range, new_range);
				// Prevent too narrow searches: always at least 3x3 pixels in the coarse search
				if (new_range < 1.5 * new_step)
					new_range = 1.5 * new_step;
				// Also prevent too wide searches: that will lead to memory problems:
				// If steps size < 1/4th of search range, then decrease search range by 50%
				if (new_range > 4. * new_step)
					new_range /= 2.;
				//If even that was not enough: use coarser step size and hope things will settle down later...
				if (new_range > 4. * new_step)
					new_step = new_range / 4.;
				sampling.setTranslations(new_step, new_range);

				// B. Use twice as fine angular sampling
				int new_hp_order;
				double new_rottilt_step, new_psi_step;
				if (mymodel.ref_dim == 3)
				{
					new_hp_order = sampling.healpix_order + 1;
					new_rottilt_step = new_psi_step = 360. / (6 * ROUND(std::pow(2., new_hp_order + adaptive_oversampling)));
				}
				else if (mymodel.ref_dim == 2)
				{
					new_hp_order = sampling.healpix_order;
					new_psi_step = sampling.getAngularSampling() / 2.;
				}
				else
					REPORT_ERROR("MlOptimiser::autoAdjustAngularSampling BUG: ref_dim should be two or three");

				// Set the new sampling in the sampling-object
				sampling.setOrientations(new_hp_order, new_psi_step * std::pow(2., adaptive_oversampling));

				// Resize the pdf_direction arrays to the correct size and fill with an even distribution
				mymodel.initialisePdfDirection(sampling.NrDirections(0, true));

				// Also reset the nr_directions in wsum_model
				wsum_model.nr_directions = mymodel.nr_directions;

				// Also resize and initialise wsum_model.pdf_direction for each class!
				for (int iclass=0; iclass < mymodel.nr_classes; iclass++)
					wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);

				// Reset iteration counters
				nr_iter_wo_resol_gain = 0;
				nr_iter_wo_large_hidden_variable_changes = 0;

				// Reset smallest changes hidden variables
				smallest_changes_optimal_classes = 9999999;
				smallest_changes_optimal_offsets = 999.;
				smallest_changes_optimal_orientations = 999.;

				// If the angular sampling is smaller than autosampling_hporder_local_searches, then use local searches of +/- 6 times the angular sampling
				if (new_hp_order >= autosampling_hporder_local_searches)
				{
					// Switch ON local angular searches
					mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
					sampling.orientational_prior_mode = PRIOR_ROTTILT_PSI;
					mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 2. * 2. * new_rottilt_step * new_rottilt_step;
					nr_pool = max_nr_pool = 1;
				}

			}

		}
	}

	// Print to screen
	if (verb)
	{
		std::cout << " Auto-refine: Angular step= " << sampling.getAngularSampling(adaptive_oversampling) << " degrees; local searches= ";
		if (sampling.orientational_prior_mode == NOPRIOR)
			std:: cout << "false" << std::endl;
		else
			std:: cout << "true" << std::endl;
		std::cout << " Auto-refine: Offset search range= " << sampling.offset_range << " pixels; offset step= " << sampling.getTranslationalSampling(adaptive_oversampling) << " pixels"<<std::endl;
	}

}

void MlOptimiser::checkConvergence()
{

	if (do_realign_movies)
	{
		// only resolution needs to be stuck
		// Since there does not seem to be any improvement (and sometimes even the opposite)
		// of performing more than one iteration with the movie frames, just perform a single iteration
		//if (nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN)
		//{
		//	has_converged = true;
		//	do_join_random_halves = true;
		//	// movies were already use all data until Nyquist
		//}
	}
	else
	{
		has_converged = false;
		if ( has_fine_enough_angular_sampling && nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES )
		{
			has_converged = true;
			do_join_random_halves = true;
			// In the last iteration, include all data until Nyquist
			do_use_all_data = true;
		}
	}

}

void MlOptimiser::printConvergenceStats()
{

	std::cout << " Auto-refine: Iteration= "<< iter<< std::endl;
	std::cout << " Auto-refine: Resolution= "<< 1./mymodel.current_resolution<< " (no gain for " << nr_iter_wo_resol_gain << " iter) "<< std::endl;
	std::cout << " Auto-refine: Changes in angles= " << current_changes_optimal_orientations << " degrees; and in offsets= " << current_changes_optimal_offsets
			<< " pixels (no gain for " << nr_iter_wo_large_hidden_variable_changes << " iter) "<< std::endl;

	if (has_converged)
	{
		std::cout << " Auto-refine: Refinement has converged, entering last iteration where two halves will be combined..."<<std::endl;
		if (!do_realign_movies)
			std::cout << " Auto-refine: The last iteration will use data to Nyquist frequency, which may take more CPU and RAM."<<std::endl;
	}

}

void MlOptimiser::setMetaDataSubset(int first_ori_particle_id, int last_ori_particle_id)
{

	for (long int ori_part_id = first_ori_particle_id, my_image_no = 0; ori_part_id <= last_ori_particle_id; ori_part_id++)
    {

#ifdef DEBUG_CHECKSIZES
		if (ori_part_id >= mydata.ori_particles.size())
		{
			std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
			REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
		}
#endif

		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];

			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++, my_image_no++)
			{

				long int img_id = mydata.getImageId(part_id, iseries);

#ifdef DEBUG_CHECKSIZES
				if (img_id >= mydata.MDimg.numberOfObjects())
				{
					std::cerr<< "img_id= "<<img_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
					REPORT_ERROR("img_id >= mydata.MDimg.numberOfObjects()");
				}
				if (my_image_no >= YSIZE(exp_metadata))
				{
					std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
					REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
				}
#endif
				mydata.MDimg.setValue(EMDL_ORIENT_ROT,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT), img_id);
				mydata.MDimg.setValue(EMDL_ORIENT_TILT, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT), img_id);
				mydata.MDimg.setValue(EMDL_ORIENT_PSI,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI), img_id);
				mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_X, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF), img_id);
				mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF), img_id);
				mydata.MDimg.setValue(EMDL_PARTICLE_CLASS, (int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS) , img_id);
				mydata.MDimg.setValue(EMDL_PARTICLE_DLL,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL), img_id);
				mydata.MDimg.setValue(EMDL_PARTICLE_PMAX, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX), img_id);
				mydata.MDimg.setValue(EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES,(int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NR_SIGN), img_id);
				mydata.MDimg.setValue(EMDL_IMAGE_NORM_CORRECTION, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM), img_id);

				// For the moment, CTF, prior and transformation matrix info is NOT updated...
				double prior_x = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR);
				double prior_y = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR);
				if (prior_x < 999.)
					mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_X_PRIOR, prior_x, img_id);
				if (prior_y < 999.)
					mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, prior_y, img_id);
			}
		}
	}

}

void MlOptimiser::getMetaAndImageDataSubset(int first_ori_particle_id, int last_ori_particle_id, bool do_also_imagedata)
{

	int nr_images = 0;
	for (long int ori_part_id = first_ori_particle_id; ori_part_id <= last_ori_particle_id; ori_part_id++)
	{

#ifdef DEBUG_CHECKSIZES
		if (ori_part_id >= mydata.ori_particles.size())
		{
			std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
			REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
		}
#endif

		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			nr_images += mydata.getNrImagesInSeries(part_id);
		}
	}

	exp_metadata.initZeros(nr_images, METADATA_LINE_LENGTH);
	if (has_converged && do_use_reconstruct_images)
		exp_imagedata.resize(2*nr_images, mymodel.ori_size, mymodel.ori_size);
	else
		exp_imagedata.resize(nr_images, mymodel.ori_size, mymodel.ori_size);

	for (long int ori_part_id = first_ori_particle_id, my_image_no = 0; ori_part_id <= last_ori_particle_id; ori_part_id++)
    {
		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			for (int iseries = 0; iseries < mydata.getNrImagesInSeries(part_id); iseries++, my_image_no++)
			{
				long int img_id = mydata.getImageId(part_id, iseries);

#ifdef DEBUG_CHECKSIZES
				if (img_id >= mydata.MDimg.numberOfObjects())
				{
					std::cerr<< "img_id= "<<img_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
					REPORT_ERROR("img_id >= mydata.MDimg.numberOfObjects()");
				}
				if (my_image_no >= YSIZE(exp_metadata))
				{
					std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
					REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
				}
				if (my_image_no >= nr_images)
				{
					std::cerr<< "my_image_no= "<<my_image_no<<" nr_images= "<< nr_images <<std::endl;
					REPORT_ERROR("my_image_no >= nr_images");
				}
#endif
				// First read the image from disc
				FileName fn_img, fn_rec_img;
				mydata.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, img_id);
				Image<double> img, rec_img;
				img.read(fn_img);
				if (XSIZE(img()) != XSIZE(exp_imagedata) || YSIZE(img()) != YSIZE(exp_imagedata) )
				{
					std::cerr << " fn_img= " << fn_img << " XSIZE(img())= " << XSIZE(img()) << " YSIZE(img())= " << YSIZE(img()) << std::endl;
					REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect image size");
				}
				if (has_converged && do_use_reconstruct_images)
				{
					mydata.MDimg.getValue(EMDL_IMAGE_RECONSTRUCT_NAME, fn_rec_img, img_id);
					rec_img.read(fn_rec_img);
					if (XSIZE(rec_img()) != XSIZE(exp_imagedata) || YSIZE(rec_img()) != YSIZE(exp_imagedata) )
					{
						std::cerr << " fn_rec_img= " << fn_rec_img << " XSIZE(rec_img())= " << XSIZE(rec_img()) << " YSIZE(rec_img())= " << YSIZE(rec_img()) << std::endl;
						REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect reconstruct_image size");
					}
				}
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
				{
					DIRECT_A3D_ELEM(exp_imagedata, my_image_no, i, j) = DIRECT_A2D_ELEM(img(), i, j);
				}

				if (has_converged && do_use_reconstruct_images)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
					{
						DIRECT_A3D_ELEM(exp_imagedata, nr_images + my_image_no, i, j) = DIRECT_A2D_ELEM(rec_img(), i, j);
					}
				}

				// Now get the metadata
				int iaux;
				mydata.MDimg.getValue(EMDL_ORIENT_ROT,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT), img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_TILT, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT), img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_PSI,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI), img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF), img_id);
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF), img_id);
				mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, iaux, img_id);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS) = (double)iaux;
				mydata.MDimg.getValue(EMDL_PARTICLE_DLL,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL), img_id);
				mydata.MDimg.getValue(EMDL_PARTICLE_PMAX, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX), img_id);
				mydata.MDimg.getValue(EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES, iaux, img_id);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NR_SIGN) = (double)iaux;
				if (!mydata.MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) = 1.;
				if (do_ctf_correction)
				{
					long int mic_id = mydata.getMicrographId(part_id, iseries);
					double kV, DeltafU, DeltafV, azimuthal_angle, Cs, Bfac, Q0;
					if (!mydata.MDimg.getValue(EMDL_CTF_VOLTAGE, kV, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_VOLTAGE, kV, mic_id))
							kV=200;

					if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, DeltafU, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUSU, DeltafU, mic_id))
							DeltafU=0;

					if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, DeltafV, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUSV, DeltafV, mic_id))
							DeltafV=DeltafU;

					if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, mic_id))
							azimuthal_angle=0;

					if (!mydata.MDimg.getValue(EMDL_CTF_CS, Cs, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_CS, Cs, mic_id))
							Cs=0;

					if (!mydata.MDimg.getValue(EMDL_CTF_BFACTOR, Bfac, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_BFACTOR, Bfac, mic_id))
							Bfac=0;

					if (!mydata.MDimg.getValue(EMDL_CTF_Q0, Q0, img_id))
						if (!mydata.MDmic.getValue(EMDL_CTF_Q0, Q0, mic_id))
							Q0=0;

					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_VOLTAGE) = kV;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_U) = DeltafU;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_V) = DeltafV;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_ANGLE) = azimuthal_angle;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_CS) = Cs;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_BFAC) = Bfac;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_Q0) = Q0;

				}

				// beamtilt
				double beamtilt_x = 0., beamtilt_y = 0.;
				if (mydata.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_X))
					mydata.MDimg.getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x, img_id);
				if (mydata.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_Y))
					mydata.MDimg.getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y, img_id);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_X) = beamtilt_x;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_Y) = beamtilt_y;

				// If the priors are NOT set, then set their values to 999.
				if (!mydata.MDimg.getValue(EMDL_ORIENT_ROT_PRIOR,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT_PRIOR), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT_PRIOR) = 999.;
				if (!mydata.MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT_PRIOR), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT_PRIOR) = 999.;
				if (!mydata.MDimg.getValue(EMDL_ORIENT_PSI_PRIOR,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR) = 999.;
				if (!mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR) = 999.;
				if (!mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR), img_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR) = 999.;

				// Pass the transformation matrix (even if it is the Identity matrix...
				Matrix2D<double> R_mic;
				R_mic = mydata.getMicrographTransformationMatrix(part_id, iseries);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_0) = R_mic(0,0);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_1) = R_mic(0,1);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_0_2) = R_mic(0,2);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_0) = R_mic(1,0);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_1) = R_mic(1,1);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_1_2) = R_mic(1,2);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_0) = R_mic(2,0);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_1) = R_mic(2,1);
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_MAT_2_2) = R_mic(2,2);

			}
		}
    }

}

