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

#include "src/ml_model.h"

void MlModel::initialise()
{

	// Auxiliary vector with relevant size in Fourier space
	MultidimArray<double > aux;
    aux.initZeros(ori_size / 2 + 1);

	// Now resize all relevant vectors
    Iref.resize(nr_classes);
    pdf_class.resize(nr_classes, 1./(double)nr_classes);
    pdf_direction.resize(nr_classes);
    group_names.resize(nr_groups, "");
    sigma2_noise.resize(nr_groups, aux);
    nr_particles_group.resize(nr_groups);
    tau2_class.resize(nr_classes, aux);
    fsc_halves_class.resize(nr_classes, aux);
    sigma2_class.resize(nr_classes, aux);
    data_vs_prior_class.resize(nr_classes, aux);
    // TODO handle these two correctly.
    bfactor_correction.resize(nr_groups, 0.);
    scale_correction.resize(nr_groups, 1.);

	acc_rot.resize(nr_classes, 0);
	acc_trans.resize(nr_classes, 0);

	if (ref_dim==2)
	{
		Matrix1D<double> empty(2);
		prior_offset_class.resize(nr_classes, empty);
	}
	// These arrays will be resized when they are filled
	orientability_contrib.resize(nr_classes);

	Projector ref(ori_size, interpolator, padding_factor, r_min_nn);
    PPref.clear();
    // Now fill the entire vector with instances of "ref"
    PPref.resize(nr_classes, ref);

}

// Reading from a file
void MlModel::read(FileName fn_in)
{

	// Clear current model
    clear();

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MlModel::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MDclass, MDgroup, MDlog, MDsigma;

    // Read general stuff
    MDlog.readStar(in, "model_general");

	if (!MDlog.getValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim) ||
		!MDlog.getValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size) ||
		!MDlog.getValue(EMDL_MLMODEL_CURRENT_RESOLUTION, current_resolution) ||
		!MDlog.getValue(EMDL_MLMODEL_CURRENT_SIZE, current_size) ||
		!MDlog.getValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor) ||
		!MDlog.getValue(EMDL_MLMODEL_INTERPOLATOR, interpolator) ||
		!MDlog.getValue(EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, r_min_nn) ||
		!MDlog.getValue(EMDL_MLMODEL_PIXEL_SIZE, pixel_size) ||
		!MDlog.getValue(EMDL_MLMODEL_NR_CLASSES, nr_classes) ||
		!MDlog.getValue(EMDL_MLMODEL_NR_GROUPS, nr_groups) ||
		!MDlog.getValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor) ||
		!MDlog.getValue(EMDL_MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_OFFSET, sigma2_offset) ||
		!MDlog.getValue(EMDL_MLMODEL_PRIOR_MODE, orientational_prior_mode) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_ROT, sigma2_rot) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_TILT, sigma2_tilt) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_PSI, sigma2_psi) ||
		!MDlog.getValue(EMDL_MLMODEL_LL, LL) ||
		!MDlog.getValue(EMDL_MLMODEL_AVE_PMAX, ave_Pmax) )
		REPORT_ERROR("MlModel::readStar: incorrect model_general table");

    // Take inverse again of current resolution:
    current_resolution = 1. / current_resolution;

    sigma2_offset *= sigma2_offset;
	sigma2_rot *= sigma2_rot;
	sigma2_tilt *= sigma2_tilt;
	sigma2_psi *= sigma2_psi;

	// Resize vectors
	initialise();

	// Read classes
	FileName fn_tmp;
	Image<double> img;
	MDclass.readStar(in, "model_classes");
	int iclass = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
	{
		if (!MDclass.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp) ||
			!MDclass.getValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[iclass])||
			!MDclass.getValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]) ||
			!MDclass.getValue(EMDL_MLMODEL_ACCURACY_TRANS, acc_trans[iclass]) 	)
			REPORT_ERROR("MlModel::readStar: incorrect model_classes table");
		if (ref_dim==2)
		{
			if (!MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass])) ||
				!MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass])) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes table: no offset priors for 2D classes");
		}

		// Read in actual reference image
		img.read(fn_tmp);
		Iref[iclass] = img();
		iclass++;
	}

	// Read group stuff
	MDgroup.readStar(in, "model_groups");
	long int igroup;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroup)
	{
        if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_NO, igroup))
                REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
        //Start counting of groups at 1, not at 0....
        if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup-1]) ||
                !MDgroup.getValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_group[igroup-1]) ||
                !MDgroup.getValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup-1]))
                REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
	}

	// Read SSNR, noise reduction, tau2_class spectra for each class
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		MDsigma.readStar(in, "model_class_" + integerToString(iclass + 1));
		int idx;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
		{
			if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
				REPORT_ERROR("MlModel::readStar: incorrect table model_class_"+integerToString(iclass));
			if (!MDsigma.getValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_TAU2_REF, tau2_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_SIGMA2_REF, sigma2_class[iclass](idx)))
				REPORT_ERROR("MlModel::readStar: incorrect table model_class_"+integerToString(iclass));
		}
	}

	// Read sigma models for each group
	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		if (nr_particles_group[igroup] > 0)
		{
			MDsigma.readStar(in, "model_group_" + integerToString(igroup + 1));
			int idx;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
			{
				if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
					REPORT_ERROR("MlModel::readStar: incorrect table model_group_"+integerToString(igroup));
				if (!MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, sigma2_noise[igroup](idx)))
					REPORT_ERROR("MlModel::readStar: incorrect table model_group_"+integerToString(igroup));
			}
		}
		else
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
			{
				DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = 0.;
			}
		}
	}

	// Read pdf_direction models for each class
	if (ref_dim == 3)
	{
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			MDclass.readStar(in, "model_pdf_orient_class_" + integerToString(iclass + 1));
			pdf_direction[iclass].clear();
			double aux;
			std::vector<double> vaux;
			vaux.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
			{
				if (!MDclass.getValue(EMDL_MLMODEL_PDF_ORIENT, aux))
					REPORT_ERROR("MlModel::readStar: incorrect table model_pdf_orient_class"+integerToString(iclass));
				vaux.push_back(aux);
			}
			pdf_direction[iclass].resize(vaux.size());
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(pdf_direction[iclass])
			{
				DIRECT_A1D_ELEM(pdf_direction[iclass], i) = vaux[i];
			}
			nr_directions = vaux.size();
		}
	}
	else
	{
		// For 2D case, just fill pdf_direction with ones.
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			pdf_direction[iclass].clear();
			pdf_direction[iclass].resize(1);
			DIRECT_A1D_ELEM(pdf_direction[iclass], 0) = 1.;
		}
		nr_directions = 1;
	}

	// Close file handler
	in.close();

}

void MlModel::write(FileName fn_out, HealpixSampling &sampling)
{

	MetaDataTable MDclass, MDgroup, MDlog, MDsigma;
    FileName fn_tmp, fn_tmp2;
    double aux;
    std::ofstream  fh;

    // A. Write images
    if (ref_dim == 2)
    {
    	Image<double> img(XSIZE(Iref[0]), YSIZE(Iref[0]), 1, nr_classes);
    	for (int iclass = 0; iclass < nr_classes; iclass++)
    	{
    		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iref[iclass])
			{
    			DIRECT_NZYX_ELEM(img(), iclass, 0, i, j) = DIRECT_A2D_ELEM(Iref[iclass], i, j);
			}
    	}
    	img.write(fn_out + "_classes.mrcs");
    }
    else
    {
    	Image<double> img;
    	// Set correct voxel size in the header
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size);
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size);
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size);
    	for (int iclass = 0; iclass < nr_classes; iclass++)
    	{
    		fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);
    		img() = Iref[iclass];
    		img.write(fn_tmp);


    	}

    	// Also write out bild files with the orientational distribution of each class
		// Also write out angular distributions
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			FileName fn_bild;
			fn_bild.compose(fn_out+"_class",iclass+1,"", 3);
			fn_bild += "_angdist.bild";
			double offset = ori_size * pixel_size / 2.;
			sampling.writeBildFileOrientationalDistribution(pdf_direction[iclass], fn_bild, offset, offset);
		}

	}

    // B. Write STAR file with metadata
    fn_tmp = fn_out + "_model.star";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MlModel::write: Cannot write file: " + fn_tmp);

    // Write the output STAR file
	MDlog.setIsList(true);
	MDlog.addObject();
	MDlog.setName("model_general");
	MDlog.setValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim);
	MDlog.setValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size);
	MDlog.setValue(EMDL_MLMODEL_CURRENT_RESOLUTION, 1./current_resolution);
	MDlog.setValue(EMDL_MLMODEL_CURRENT_SIZE, current_size);
	MDlog.setValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor);
	MDlog.setValue(EMDL_MLMODEL_INTERPOLATOR, interpolator);
	MDlog.setValue(EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, r_min_nn);
	MDlog.setValue(EMDL_MLMODEL_PIXEL_SIZE, pixel_size);
	MDlog.setValue(EMDL_MLMODEL_NR_CLASSES, nr_classes);
	MDlog.setValue(EMDL_MLMODEL_NR_GROUPS, nr_groups);
	MDlog.setValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor);
	MDlog.setValue(EMDL_MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction);
	MDlog.setValue(EMDL_MLMODEL_SIGMA_OFFSET, sqrt(sigma2_offset));
	MDlog.setValue(EMDL_MLMODEL_PRIOR_MODE, orientational_prior_mode);
	MDlog.setValue(EMDL_MLMODEL_SIGMA_ROT, sqrt(sigma2_rot));
	MDlog.setValue(EMDL_MLMODEL_SIGMA_TILT, sqrt(sigma2_tilt));
	MDlog.setValue(EMDL_MLMODEL_SIGMA_PSI, sqrt(sigma2_psi));
	MDlog.setValue(EMDL_MLMODEL_LL, LL);
	MDlog.setValue(EMDL_MLMODEL_AVE_PMAX, ave_Pmax);
	MDlog.write(fh);

	// Write metadata and images for all classes
	MDclass.setName("model_classes");
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		MDclass.addObject();
		if (ref_dim==2)
		{
			fn_tmp = fn_out + "_classes.mrcs";
			fn_tmp.compose(iclass+1, fn_tmp); // fn_tmp = integerToString(iclass) + "@" + fn_tmp;
		}
		else
		{
			fn_tmp.compose(fn_out+"_class",iclass+1,"mrc", 3); // class number from 1 to K!
		}
		MDclass.setValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);
		MDclass.setValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_TRANS, acc_trans[iclass]);

		if (ref_dim==2)
		{
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass]));
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass]));
		}

	}
	MDclass.write(fh);

	// Write radial_average of tau2_class and data_vs_prior_class for each reference
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		MDsigma.clear();
		MDsigma.setName("model_class_"+integerToString(iclass+1));
		for (int ii = 0; ii < XSIZE(tau2_class[iclass]); ii++)
		{
			MDsigma.addObject();
			MDsigma.setValue(EMDL_SPECTRAL_IDX, ii);
			MDsigma.setValue(EMDL_RESOLUTION, getResolution(ii));
			MDsigma.setValue(EMDL_RESOLUTION_ANGSTROM, getResolutionAngstrom(ii));
			MDsigma.setValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_SIGMA2_REF, sigma2_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_TAU2_REF, tau2_class[iclass](ii));
			// Only write orientabilities if they have been determined
			if (XSIZE(orientability_contrib[iclass]) == XSIZE(tau2_class[iclass]))
				MDsigma.setValue(EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION, orientability_contrib[iclass](ii));
		}
		MDsigma.write(fh);
	}

    // Write scale-correction for all groups
    MDgroup.setName("model_groups");
    for (long int igroup = 0; igroup < nr_groups; igroup++)
    {
		MDgroup.addObject();
		//Start counting of groups at 1, not at 0....
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NO, igroup+1);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup]);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_group[igroup]);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup]);
    }
    MDgroup.write(fh);

	// Write sigma models for each group
	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		MDsigma.clear();
		MDsigma.setName("model_group_"+integerToString(igroup+1));
		for (int ii = 0; ii < XSIZE(sigma2_noise[igroup]); ii++)
		{
			MDsigma.addObject();
			// Some points in sigma2_noise arrays are never used...
			aux = sigma2_noise[igroup](ii);
			if (aux > 0.)
			{
				MDsigma.setValue(EMDL_SPECTRAL_IDX, ii);
				MDsigma.setValue(EMDL_RESOLUTION, getResolution(ii));
				MDsigma.setValue(EMDL_MLMODEL_SIGMA2_NOISE, aux);
			}
		}
		MDsigma.write(fh);
	}

	// Write pdf_direction models for each class
	if (ref_dim == 3)
	{
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			MDclass.clear();
			MDclass.setName("model_pdf_orient_class_"+integerToString(iclass+1));
			for (int ii=0; ii < XSIZE(pdf_direction[iclass]); ii++)
			{
				MDclass.addObject();
				MDclass.setValue(EMDL_MLMODEL_PDF_ORIENT, pdf_direction[iclass](ii));
			}
			MDclass.write(fh);
		}
	}

}




void  MlModel::readTauSpectrum(FileName fn_tau, int verb)
{
	MetaDataTable MDtau;
	double val;
	int idx;
	MDtau.read(fn_tau);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDtau)
	{
		MDtau.getValue(EMDL_SPECTRAL_IDX, idx);
		MDtau.getValue(EMDL_MLMODEL_TAU2_REF, val);
		if (idx < XSIZE(tau2_class[0]))
			tau2_class[0](idx) = tau2_fudge_factor * val;
	}
	if (idx < XSIZE(tau2_class[0]) - 1)
	{
		if (verb > 0) std::cerr<< " Warning: provided tau2-spectrum has fewer entries ("<<idx+1<<") than needed ("<<XSIZE(tau2_class[0])<<"). Set rest to zero..."<<std::endl;
	}
	// Use the same spectrum for all classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
		tau2_class[iclass] =  tau2_class[0];

}

// Reading images from disc
void MlModel::readImages(FileName fn_ref, int _ori_size, Experiment &_mydata,
			bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected)
{

	// Set some stuff
	nr_groups = _mydata.groups.size();
	ori_size = _ori_size;
	double avg_norm_correction = 1.;

	// Read references into memory
	Image<double> img;
	FileName fn_tmp;
	if (fn_ref != "None")
	{
		// Read the references into memory
		do_average_unaligned = false;
		// If this is a STAR file, ignore nr_classes and read all references from this file
		if (fn_ref.isStarFile())
		{
			MetaDataTable MDref;
			MDref.read(fn_ref);
			do_generate_seeds = false;
			// ignore nr_classes from the command line, use number of entries in STAR file
			nr_classes = 0;
			Iref.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
			{
				MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);
				img.read(fn_tmp);
				ref_dim = img().getDim();
				if (ori_size != XSIZE(img()) || ori_size != YSIZE(img()))
				{
					std::cerr << " ori_size= " << ori_size << " XSIZE(img())= " << XSIZE(img()) << std::endl;
					REPORT_ERROR("MlOptimiser::read: size of reference images is not the same as the experimental images!");
				}
				Iref.push_back(img());
				nr_classes++;
			}
		}
		// For a single image, read this image as reference and set it in all nr_classes Irefs
		else
		{
			img.read(fn_ref);
			ref_dim = img().getDim();
			if (ori_size != XSIZE(img()) || ori_size != YSIZE(img()))
			{
				std::cerr << " ori_size= " << ori_size << " XSIZE(img())= " << XSIZE(img()) << std::endl;
				REPORT_ERROR("MlOptimiser::read: size of reference image is not the same as the experimental images!");
			}
			Iref.clear();
			for (int iclass = 0; iclass < nr_classes; iclass++)
				Iref.push_back(img());
			if (nr_classes > 1)
				do_generate_seeds = true;
			else
				do_generate_seeds = false;
		}
	}
	else
	{
		// If no -ref is given, assume this is a 2D refinement and calculate average of all unaligned images later on.
		do_average_unaligned = true;
		do_generate_seeds = true;
		refs_are_ctf_corrected = false;
		img().initZeros(ori_size, ori_size);
		ref_dim = 2;
		for (int iclass = 0; iclass < nr_classes; iclass++)
			Iref.push_back(img());
	}

	initialise();

	// Now set the group names from the Experiment groups list
	for (int i=0; i<nr_groups; i++)
		group_names[i] = _mydata.groups[i].name;

}
void MlModel::expandToMovieFrames(Experiment &moviedata, int running_avg_side)
{
	MlModel moviemodel;

	// Start from all equal model (nr_classes, references etc)
	moviemodel = (*this);
	// Reset the average norm correction to one
	moviemodel.avg_norm_correction = 1.;

	// Then delete the current groups: their numbers, names, scale_corrections and sigma2noise spectra
	moviemodel.nr_groups = 0;
	moviemodel.nr_particles_group.clear();
	moviemodel.group_names.clear();
	moviemodel.sigma2_noise.clear();
	moviemodel.scale_correction.clear();
	moviemodel.bfactor_correction.clear();

	// Now go and look in the (already expanded) moviedata to find the unique groups
	for (int i=0; i<moviedata.groups.size(); i++)
	{

		FileName curr_name, movie_name;
		movie_name = moviedata.groups[i].name;

		// Find the corresponding group_name in the current model
		int my_curr_group_nr = -1;
		for (int j=0; j<group_names.size(); j++)
		{
			curr_name = group_names[j].withoutExtension();
			// The moviename should be the current name plus something else...
			if (movie_name.contains(curr_name))
			{
				my_curr_group_nr = j;
				break;
			}
		}
		if (my_curr_group_nr < 0)
			REPORT_ERROR("MlModel::expandToMovieFrames ERROR: cannot find rlnMicrographName or rlnGroupName for movie frame: " + movie_name);

		moviemodel.sigma2_noise.push_back(sigma2_noise[my_curr_group_nr]);
		moviemodel.scale_correction.push_back(scale_correction[my_curr_group_nr]);
		moviemodel.bfactor_correction.push_back(bfactor_correction[my_curr_group_nr]);
		moviemodel.group_names.push_back(movie_name);
		moviemodel.nr_groups++;
	}

	// Also find the number of particles in this group
	moviemodel.nr_particles_group.resize(moviemodel.nr_groups);
	std::vector<int> nr_frames_in_group;
	nr_frames_in_group.resize(moviemodel.nr_groups, -1);

	for (long int ipart = 0; ipart < moviedata.particles.size(); ipart++)
	{
		for (int iimg = 0; iimg < (moviedata.particles[ipart]).images.size(); iimg++)
		{
			long int group_id = ((moviedata.particles[ipart]).images[iimg]).group_id;
			long int img_id = ((moviedata.particles[ipart]).images[iimg]).id;
			// count the number of particles in this group
			moviemodel.nr_particles_group[group_id]++;
			// Get and check the number of frames is constant within one group
			int nframes = -1;
			if (!moviedata.MDimg.getValue(EMDL_PARTICLE_NR_FRAMES, nframes, img_id))
				REPORT_ERROR("MlModel::expandToMovieFrames ERROR: cannot find rlnNrOfFrames in moviedata.");
			if (nr_frames_in_group[group_id] < 0)
				nr_frames_in_group[group_id] = nframes;
			else if (nr_frames_in_group[group_id] != nframes)
				REPORT_ERROR((std::string)"MlModel::expandToMovieFrames ERROR: unequal number of frames in group" + moviemodel.group_names[group_id]);
		}
	}

	// Correct the input sigma2noise spectra by a factor of nframes
	for (int i=0; i < moviemodel.nr_groups; i++)
	{
		moviemodel.sigma2_noise[i] *= (double)nr_frames_in_group[i]/((double)(2 * running_avg_side + 1));
	}

	// Now replace the current model with the expanded moviemodel
	(*this) = moviemodel;

}

void MlModel::initialisePdfDirection(int newsize)
{

	// If the pdf_direction were already filled (size!=0), and newsize=oldsize then leave them as they were
	// If they were still empty, or if the size changes, then initialise them with an even distribution
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		int oldsize = MULTIDIM_SIZE(pdf_direction[iclass]);
		if (oldsize == 0 || oldsize != newsize)
		{
			pdf_direction[iclass].resize(newsize);
			pdf_direction[iclass].initConstant(1./((double) nr_classes * newsize));
		}
	}
	nr_directions = newsize;

}

void MlModel::setFourierTransformMaps(bool update_tau2_spectra, int nr_threads)
{

	for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        if (update_tau2_spectra)
        {
			PPref[iclass].computeFourierTransformMap(Iref[iclass], tau2_class[iclass], current_size, nr_threads);
        }
        else
        {
        	MultidimArray<double> dummy;
        	PPref[iclass].computeFourierTransformMap(Iref[iclass], dummy, current_size, nr_threads);
        }
    }

}

void MlModel::initialiseDataVersusPrior(bool fix_tau)
{

    // Get total number of particles
	double nr_particles = 0.;
	for (int igroup = 0; igroup < nr_particles_group.size(); igroup++)
		nr_particles += (double)nr_particles_group[igroup];

	// Calculate average sigma2_noise over all image groups
	MultidimArray<double> avg_sigma2_noise;
	avg_sigma2_noise.initZeros(sigma2_noise[0]);
	for (int igroup = 0; igroup < nr_particles_group.size(); igroup++)
	{
		avg_sigma2_noise += (double)(nr_particles_group[igroup]) * sigma2_noise[igroup];
	}
	avg_sigma2_noise /= nr_particles;

	// Get the FT of all reference structures
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    // And spectrum is squared, so ori_size*ori_size in the 3D case!
	double normfft = (ref_dim == 3) ? (double)(ori_size * ori_size) : 1.;

    for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		// Initialise output arrays to correct size
		tau2_class[iclass].resize(sigma2_noise[0]);

		// Get the power spectrum of the reference
		MultidimArray<double> spectrum(sigma2_noise[0]);
		getSpectrum(Iref[iclass], spectrum, POWER_SPECTRUM);

		// Factor two because of two-dimensionality of the complex plane
		// (just like sigma2_noise estimates, the power spectra should be divided by 2)
		spectrum *= normfft / 2.;

		// Update the tau2_class spectrum for this reference
		// This is only for writing out in the it000000_model.star file
		if (!fix_tau)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(tau2_class[iclass])
			{
				DIRECT_A1D_ELEM(tau2_class[iclass], i) = tau2_fudge_factor * DIRECT_A1D_ELEM(spectrum, i);
			}
		}

		// Calculate data_vs_prior_class as spectral_nr_observations_per_class/sigma2_noise vs 1/tau2_class
		data_vs_prior_class[iclass].resize(sigma2_noise[0]);
		fsc_halves_class[iclass].initZeros(sigma2_noise[0]);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(tau2_class[iclass])
		{
			double evidence = nr_particles * pdf_class[iclass] / DIRECT_A1D_ELEM(avg_sigma2_noise, i);
			// empirical accounting for ratio of pixels in 3D shells compared to 2D shells
			if (ref_dim == 3 && i > 0)
				evidence /= (2. * (double)i);
			double prior = 1. /  DIRECT_A1D_ELEM(tau2_class[iclass], i);
			double myssnr = evidence / prior;
			DIRECT_A1D_ELEM(data_vs_prior_class[iclass], i ) = myssnr;
			// Also initialise FSC-halves here (...)
			//DIRECT_A1D_ELEM(fsc_halves_class[iclass], i ) = myssnr / (myssnr + 1);
		}
	} // end loop iclass

}

/////////// MlWsumModel
void MlWsumModel::initialise(MlModel &_model, FileName fn_sym)
{
	nr_classes = _model.nr_classes;
    nr_groups = _model.nr_groups;
    nr_directions = _model.nr_directions;
    ref_dim = _model.ref_dim;
    ori_size = _model.ori_size;
    pdf_class = _model.pdf_class;
    if (ref_dim == 2)
    	prior_offset_class = _model.prior_offset_class;
    pdf_direction = _model.pdf_direction;
    sigma2_offset = _model.sigma2_offset;
    sigma2_noise = _model.sigma2_noise;
    sigma2_rot = _model.sigma2_rot;
    sigma2_tilt = _model.sigma2_tilt;
    sigma2_psi = _model.sigma2_psi;
    padding_factor = _model.padding_factor;
    interpolator = _model.interpolator;
    r_min_nn = _model.r_min_nn;

    // Don't need forward projectors in MlWsumModel!
    PPref.clear();
    // Don't need scale_correction and bfactor_correction, keep wsum_signal_product_spectra and wsum_reference_power_spectra instead
    scale_correction.clear();
    bfactor_correction.clear();
    tau2_class.clear();
    data_vs_prior_class.clear();
	acc_rot.clear();
	acc_trans.clear();
    orientability_contrib.clear();


    MultidimArray<double> aux(ori_size / 2 + 1);
    wsum_signal_product_spectra.resize(nr_groups, aux);
    wsum_reference_power_spectra.resize(nr_groups, aux);

    // Resize MlWsumModel-specific vectors
    BackProjector BP(ori_size, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn,
    		         ML_BLOB_ORDER, ML_BLOB_RADIUS, ML_BLOB_ALPHA);
    BPref.clear();
    BPref.resize(nr_classes, BP);
    sumw_group.resize(nr_groups);

}

void MlWsumModel::initZeros()
{

    LL = 0.;
    ave_Pmax = 0.;
    sigma2_offset = 0.;
    avg_norm_correction = 0.;
    sigma2_rot = 0.;
    sigma2_tilt = 0.;
    sigma2_psi = 0.;

    // Set all weighted sums to zero
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	BPref[iclass].initZeros(current_size);
        pdf_class[iclass] = 0.;
        // Assume pdf_direction is already of the right size...
        pdf_direction[iclass].initZeros();
        if (ref_dim == 2)
        	prior_offset_class[iclass].initZeros();
    }

    // Initialise sigma2_noise spectra and sumw_group
    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
        sumw_group[igroup] = 0.;
        sigma2_noise[igroup].initZeros();
        wsum_signal_product_spectra[igroup].initZeros();
        wsum_reference_power_spectra[igroup].initZeros();
    }
}

//#define DEBUG_PACK
#ifdef DEBUG_PACK
#define MAX_PACK_SIZE     100000
#else
// Approximately 1024*1024*1024/8/2 ~ 0.5 Gb
#define MAX_PACK_SIZE 671010000
#endif

void MlWsumModel::pack(MultidimArray<double> &packed)
{
	// for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
	long long int packed_size = 0;
    int spectral_size = (ori_size / 2) + 1;

    packed_size += 7 ;
    // for all group-related stuff
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    // for sumw_group
    packed_size += nr_groups;
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes * 2 * BPref[0].getSize();
    packed_size += nr_classes * BPref[0].getSize();
    packed_size += nr_classes * nr_directions;
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim==2)
    	packed_size += nr_classes*2;

    // Get memory for the packed array
    packed.clear();
    packed.resize(packed_size);

    // Start packing
    long long int idx = 0;

    DIRECT_MULTIDIM_ELEM(packed, idx++) = LL;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = ave_Pmax;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_offset;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = avg_norm_correction;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_rot;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_tilt;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_psi;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
        }
    	sigma2_noise[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n);
        }
        wsum_signal_product_spectra[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n);
        }
        wsum_reference_power_spectra[igroup].clear();

        DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];

    }
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real;
            DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag;
        }
    	BPref[iclass].data.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n);
        }
        BPref[iclass].weight.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
        }
        pdf_direction[iclass].clear();

        DIRECT_MULTIDIM_ELEM(packed, idx++) = pdf_class[iclass];

        if (ref_dim==2)
        {
        	DIRECT_MULTIDIM_ELEM(packed, idx++) = XX(prior_offset_class[iclass]);
        	DIRECT_MULTIDIM_ELEM(packed, idx++) = YY(prior_offset_class[iclass]);
        }
    }
#ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
#endif

    // Just to check whether we went outside our memory...
    if (idx != packed_size)
    {
       	std::cerr << "idx= " << idx << "packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != packed_size");
    }

}
void MlWsumModel::unpack(MultidimArray<double> &packed)
{
    int spectral_size = (ori_size / 2) + 1;

    long long int idx = 0;

    LL = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ave_Pmax = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_offset = DIRECT_MULTIDIM_ELEM(packed, idx++);
    avg_norm_correction = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_rot = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_tilt = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_psi = DIRECT_MULTIDIM_ELEM(packed, idx++);

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	sigma2_noise[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        wsum_signal_product_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        wsum_reference_power_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
    }

    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	BPref[iclass].initialiseDataAndWeight(current_size);
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
    		(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real = DIRECT_MULTIDIM_ELEM(packed, idx++);
    		(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
    		DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    	pdf_direction[iclass].resize(nr_directions);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
        	DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        pdf_class[iclass] = DIRECT_MULTIDIM_ELEM(packed, idx++);

        if (ref_dim==2)
        {
        	XX(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	YY(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    }

    long long int packed_size = MULTIDIM_SIZE(packed);
    packed.clear();

    // Just to check whether we went outside our memory...
    if (idx != packed_size)
    {
       	std::cerr << "idx= " << idx << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }
}


void MlWsumModel::pack(MultidimArray<double> &packed, int &piece, int &nr_pieces, bool do_clear)
{


    // Determine size of the packed array
    int nr_groups = sigma2_noise.size();
    int nr_classes = BPref.size();
    int spectral_size = (ori_size / 2) + 1;
    long long int packed_size = 0;
    long long int idx_start, idx_stop;

	// for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
    packed_size += 7 ;
    // for all group-related stuff
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    // for sumw_group
    packed_size += nr_groups;
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes * 2 * BPref[0].getSize();
    packed_size += nr_classes * BPref[0].getSize();
    packed_size += nr_classes * nr_directions;
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim==2)
    	packed_size += nr_classes*2;

    if (piece < 0 && nr_pieces < 0)
    {
    	// Special case: prevent making multiple pieces if input piece and nr_pieces are both negative
        idx_start = 0;
        idx_stop = packed_size;
    }
    else if (packed_size > MAX_PACK_SIZE)
    {
        idx_start = piece * MAX_PACK_SIZE;
        idx_stop = XMIPP_MIN(idx_start + MAX_PACK_SIZE, packed_size);
        nr_pieces = CEIL((double)packed_size/(double)MAX_PACK_SIZE);
    }
    else
    {
        idx_start = 0;
        idx_stop = packed_size;
        nr_pieces = 1;
    }

    // increment piece so that pack will be called again
    piece++;
#ifdef DEBUG_PACK
    std::cerr << " PACK: idx_start= " << idx_start << " idx_stop= " << idx_stop << " piece= " << piece << " nr_pieces= " << nr_pieces <<" packed_size= "<<packed_size<< std::endl;
    std::cerr << " nr_classes= " << nr_classes << " nr_groups= " << nr_groups << " packed_size= " << packed_size << std::endl;
    std::cerr << " MULTIDIM_SIZE(sigma2_noise[0])= " << MULTIDIM_SIZE(sigma2_noise[0]) << " MULTIDIM_SIZE(wsum_signal_product_spectra[0])= " << MULTIDIM_SIZE(wsum_signal_product_spectra[0]) << " MULTIDIM_SIZE(wsum_reference_power_spectra[0])= " << MULTIDIM_SIZE(wsum_reference_power_spectra[0]) << std::endl;
    std::cerr << " sigma2_noise.size()= " << sigma2_noise.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product_spectra.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product_spectra.size() << std::endl;
    std::cerr << " MULTIDIM_SIZE(pdf_direction[0])= " << MULTIDIM_SIZE(pdf_direction[0]) << " pdf_direction.size()= " << pdf_direction.size()<<std::endl;
#endif

    // Get memory for the packed array
    packed.clear();
    packed.resize(idx_stop - idx_start);

    long long int idx = 0;
    long long int ori_idx = 0;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = LL;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = ave_Pmax;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_offset;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = avg_norm_correction;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_rot;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_tilt;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_psi;
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
            ori_idx++;
        }
    	if (idx == ori_idx && do_clear)
            sigma2_noise[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n);
        	ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            wsum_signal_product_spectra[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            wsum_reference_power_spectra[igroup].clear();

        if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];
        ori_idx++;

    }
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real;
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag;
            ori_idx++;
        }
        // Only clear after the whole array has been packed... i.e. not when we reached the pack_size halfway through
        if (idx == ori_idx && do_clear)
            BPref[iclass].data.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            BPref[iclass].weight.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
        	pdf_direction[iclass].clear();

        if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = pdf_class[iclass];
        ori_idx++;

        if (ref_dim==2)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = XX(prior_offset_class[iclass]);
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = YY(prior_offset_class[iclass]);
            ori_idx++;
        }
    }
#ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
#endif

    // Just to check whether we went outside our memory...
    //std::cerr << " PACK piece= " << piece-1 << " nr_pieces= " << nr_pieces << " ori_idx= " << ori_idx<< " packed_size= " << packed_size << std::endl;
    //std::cerr << " PACK idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop    << std::endl;
    if (idx != idx_stop-idx_start)
    {
       	std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != idx_stop-idx_start");

    }

}

void MlWsumModel::unpack(MultidimArray<double> &packed, int piece, bool do_clear)
{


    int nr_groups = sigma2_noise.size();
    int nr_classes = BPref.size();
    int spectral_size = (ori_size / 2) + 1;
    long long int idx_start;
    long long int idx_stop;
    if (piece < 0)
    {
    	// Special case: prevent making multiple pieces if input piece is negative
        idx_start = 0;
        idx_stop  = MULTIDIM_SIZE(packed);
    }
    else
    {
    	idx_start = piece * MAX_PACK_SIZE;
    	idx_stop  = idx_start + MULTIDIM_SIZE(packed);
    }
    long long int ori_idx = 0;
    long long int idx = 0;
#ifdef DEBUG_PACK
    std::cerr << " UNPACK piece= " << piece << " idx_start= " << idx_start << " idx_stop= " << idx_stop << std::endl;
#endif

    if (ori_idx >= idx_start && ori_idx < idx_stop) LL = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) ave_Pmax = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_offset = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) avg_norm_correction = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_rot = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_tilt = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_psi = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {

    	if (idx == ori_idx)
    		sigma2_noise[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (idx == ori_idx)
    		wsum_signal_product_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (idx == ori_idx)
    		wsum_reference_power_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (ori_idx >= idx_start && ori_idx < idx_stop)
        	sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
        ori_idx++;

    }

    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	if (idx == ori_idx)
    		BPref[iclass].initialiseDataAndWeight(current_size);

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
        	if (ori_idx >= idx_start && ori_idx < idx_stop)
				(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	ori_idx++;

        	if (ori_idx >= idx_start && ori_idx < idx_stop)
            	(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	ori_idx++;
            //DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n) = Complex(re, im);
        }

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

    	if (idx == ori_idx)
    		pdf_direction[iclass].resize(nr_directions);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (ori_idx >= idx_start && ori_idx < idx_stop)
        	pdf_class[iclass] = DIRECT_MULTIDIM_ELEM(packed, idx++);
        ori_idx++;

        if (ref_dim == 2)
        {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				XX(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				YY(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
        }
    }


    long long int packed_size = MULTIDIM_SIZE(packed);
    // Free memory
    if (do_clear)
        packed.clear();

    // Just to check whether we went outside our memory...
    //std::cerr << " UNPACK piece= " << piece << " idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop    << std::endl;
    if (idx != idx_stop-idx_start)
    {
       	std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }


}



