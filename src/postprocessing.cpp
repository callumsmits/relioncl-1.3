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

#include "src/postprocessing.h"

void Postprocessing::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Input rootname, e.g. run1");
	fn_out = parser.getOption("--o", "Output rootname", "postprocess");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms"));

	int mask_section = parser.addSection("Masking options");
	do_auto_mask = parser.checkOption("--auto_mask", "Perform automated masking, based on a density threshold");
	ini_mask_density_threshold = textToFloat(parser.getOption("--inimask_threshold", "Density at which to threshold the map for the initial seed mask", "0.02"));
	extend_ini_mask = textToFloat(parser.getOption("--extend_inimask", "Number of pixels to extend the initial seed mask", "3."));
	width_soft_mask_edge  = textToFloat(parser.getOption("--width_mask_edge", "Width for the raised cosine soft mask edge (in pixels)", "6."));
	fn_mask = parser.getOption("--mask", "Filename of a user-provided mask (1=protein, 0=solvent, all values in range [0,1])", "");

	int sharp_section = parser.addSection("Sharpening options");
	fn_mtf = parser.getOption("--mtf", "User-provided STAR-file with the MTF-curve of the detector", "");
	do_auto_bfac = parser.checkOption("--auto_bfac", "Perform automated B-factor determination (Rosenthal and Henderson, 2003)");
	fit_minres = textToFloat(parser.getOption("--autob_lowres", "Lowest resolution (in A) to include in fitting of the B-factor", "10."));
	fit_maxres = textToFloat(parser.getOption("--autob_highres", "Highest resolution (in A) to include in fitting of the B-factor", "0."));
	adhoc_bfac =  textToFloat(parser.getOption("--adhoc_bfac", "User-provided B-factor (in A^2) for map sharpening, e.g. -400", "0."));

	int filter_section = parser.addSection("Filtering options");
	do_fsc_weighting = !parser.checkOption("--skip_fsc_weighting", "Do not use FSC-weighting (Rosenthal and Henderson, 2003) in the sharpening process");
	// include low-pass filter option in the program? This could be useful for structurally heterogeneous reconstructions (instead of FSC-weighting)
	low_pass_freq = textToFloat(parser.getOption("--low_pass", "Resolution (in Angstroms) at which to low-pass filter the final map (by default at final resolution)", "0."));

	int expert_section = parser.addSection("Expert options");
	randomize_fsc_at = textToFloat(parser.getOption("--randomize_at_fsc", "Randomize phases from the resolution where FSC drops below this value", "0.8"));
	filter_edge_width = textToInteger(parser.getOption("--filter_edge_width", "Width of the raised cosine on the low-pass filter edge (in resolution shells)", "2"));
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void Postprocessing::usage()
{
	parser.writeUsage(std::cerr);
}

void Postprocessing::clear()
{
	fn_in = "";
	fn_out="postprocess";
	angpix = 1.;
	do_auto_mask = false;
	ini_mask_density_threshold = 0.02;
	width_soft_mask_edge = 6.;
	fn_mask = "";
	fn_mtf = "";
	do_auto_bfac = false;
	fit_minres = 10.;
	fit_maxres = 0.;
	adhoc_bfac = 0.;
	do_fsc_weighting = true;
	low_pass_freq = 0.;
	randomize_fsc_at = 0.8;
	filter_edge_width = 2.;
	verb = 1;

}

void Postprocessing::initialise()
{
	// Read in the input maps
	fn_I1 = fn_in + "_half1_class001_unfil.mrc";
	fn_I2 = fn_in + "_half2_class001_unfil.mrc";

	if (verb > 0)
	{
		std::cout <<"== Reading input half-reconstructions: " <<std::endl;
		std::cout.width(35); std::cout << std::left <<"  + half1-map: "; std::cout << fn_I1 << std::endl;
		std::cout.width(35); std::cout << std::left <<"  + half2-map: "; std::cout << fn_I2 << std::endl;
	}

	I1.read(fn_I1);
	I2.read(fn_I2);
	I1().setXmippOrigin();
	I2().setXmippOrigin();

	if (!I1().sameShape(I2()))
	{
		std::cerr << " Size of half1 map: "; I1().printShape(std::cerr); std::cerr << std::endl;
		std::cerr << " Size of half2 map: "; I2().printShape(std::cerr); std::cerr << std::endl;
		REPORT_ERROR("Postprocessing::initialise ERROR: The two half reconstructions are not of the same size!");
	}

	if (do_auto_mask && fn_mask != "")
		REPORT_ERROR("Postprocessing::initialise ERROR: provide either --auto_mask OR --mask, but not both!");

	if (do_auto_bfac && ABS(adhoc_bfac) > 0.)
		REPORT_ERROR("Postprocessing::initialise ERROR: provide either --auto_bfac OR --adhoc_bfac, but not both!");
}

void Postprocessing::getAutoMask()
{

	if (verb > 0)
	{
		std::cout << "== Perform auto-masking ..." << std::endl;
		std::cout.width(35); std::cout << std::left  << "  + density threshold: "; std::cout  << ini_mask_density_threshold << std::endl;
		std::cout.width(35); std::cout << std::left  << "  + extend ini mask: "; std::cout  << extend_ini_mask << " pixels" << std::endl;
		std::cout.width(35); std::cout << std::left  << "  + width soft edge: "; std::cout  << width_soft_mask_edge << " pixels" << std::endl;
	}

	// Store sum of both masks in Im
	I1() += I2();
	I1() /= 2.;
	autoMask(I1(), Im(), ini_mask_density_threshold, extend_ini_mask, width_soft_mask_edge, true); // true sets verbosity

	// Re-read original I1 into memory
	I1.read(fn_I1);
	I1().setXmippOrigin();

}

bool Postprocessing::getMask()
{

	// A. Check whether a user-provided mask is to be used
	if (do_auto_mask)
	{
		getAutoMask();

	}
	else if (fn_mask != "")
	{
		if (verb > 0)
		{
			std::cout << "== Using a user-provided mask ... " <<std::endl;
			std::cout.width(35); std::cout << std::left   << "  + input mask: "; std::cout  << fn_mask <<std::endl;
		}

		// Read the mask in memory
		Im.read(fn_mask);
		Im().setXmippOrigin();

		// Check values are between 0 and 1
		double avg, stddev, minval, maxval;
		Im().computeStats(avg, stddev, minval, maxval);
		if (minval < -1e-6 || maxval - 1. > 1.e-6)
		{
			std::cerr << " minval= " << minval << " maxval= " << maxval << std::endl;
			REPORT_ERROR("Postprocessing::mask ERROR: mask values not in range [0,1]!");
		}

		// Also check the mask is the same size as the input maps
		if (!Im().sameShape(I2()))
		{
			std::cerr << " Size of input mask: "; Im().printShape(std::cerr); std::cerr<< std::endl;
			std::cerr << " Size of input maps: "; I1().printShape(std::cerr); std::cerr<< std::endl;
			REPORT_ERROR("Postprocessing::mask ERROR: mask and input maps do not have the same size!");
		}

	}
	else
	{
		if (verb > 0)
		{
			std::cout << "== Not performing any masking ... " << std::endl;
		}
		return false;
	}

	return true;
}

void Postprocessing::divideByMtf(MultidimArray<Complex > &FT)
{

	if (fn_mtf != "")
	{
		if (verb > 0)
		{
			std::cout << "== Dividing map by the MTF of the detector ..." << std::endl;
			std::cout.width(35); std::cout << std::left <<"  + mtf STAR-file: "; std::cout << fn_mtf << std::endl;
		}

		correctMapForMTF(FT, XSIZE(I1()), fn_mtf);
	}

}


void Postprocessing::sharpenMap()
{

	MultidimArray<Complex > FT;
	FourierTransformer transformer;
	transformer.FourierTransform(I1(), FT, true);

	makeGuinierPlot(FT, guinierin);

	// A. If MTF curve is given, first divide by the MTF
	divideByMtf(FT);
	makeGuinierPlot(FT, guinierinvmtf);

	// B. Then perform B-factor sharpening
	if (do_fsc_weighting)
	{
		if (verb > 0)
		{
			std::cout <<"== Applying sqrt(2*FSC/(FSC+1)) weighting (as in Rosenthal & Henderson, 2003) ..." <<std::endl;
		}
		applyFscWeighting(FT, fsc_true);
	}
	makeGuinierPlot(FT, guinierweighted);

	global_bfactor = 0.;
	if (do_auto_bfac)
	{
		if (verb > 0)
		{
			std::cout <<"== Fitting straight line through Guinier plot to find B-factor ..." <<std::endl;
			std::cout.width(35); std::cout << std::left <<"  + fit from resolution: "; std::cout << fit_minres << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + fit until resolution: "; std::cout << fit_maxres << std::endl;
		}

		fitStraightLine(guinierweighted, global_slope, global_intercept, global_corr_coeff);
		global_bfactor = 4. * global_slope;
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left  <<"  + slope of fit: "; std::cout << global_slope << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + intercept of fit: "; std::cout << global_intercept << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + correlation of fit: "; std::cout << global_corr_coeff << std::endl;
		}
	}
	else if (ABS(adhoc_bfac) > 0.)
	{
		if (verb > 0)
		{
			std::cout <<"== Using a user-provided (ad-hoc) B-factor ..." <<std::endl;
		}
		if (adhoc_bfac > 0.)
			std::cout <<" WARNING: using a positive B-factor. This will effectively dampen your map. Use negative value to sharpen it!" << std::endl;
		global_bfactor = adhoc_bfac;
	}

	// Now apply the B-factor
	if (ABS(global_bfactor) > 0.)
	{
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left  <<"  + apply b-factor of: "; std::cout << global_bfactor << std::endl;
		}
		applyBFactorToMap(FT, XSIZE(I1()), global_bfactor, angpix);
	}

	makeGuinierPlot(FT, guiniersharpen);

	if (verb > 0)
	{
		std::cout << "== Low-pass filtering final map ... " << std::endl;
	}
	double my_filter = (low_pass_freq > 0.) ? low_pass_freq : global_resol;
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left  <<"  + filter frequency: "; std::cout << my_filter << std::endl;
	}
	lowPassFilterMap(FT, XSIZE(I1()), my_filter, angpix, filter_edge_width);

	transformer.inverseFourierTransform(FT, I1());

}

void Postprocessing::applyFscWeighting(MultidimArray<Complex > &FT, MultidimArray<double> my_fsc)
{
	// Find resolution where fsc_true drops below zero for the first time
	// Set all weights to zero beyond that resolution
	int ires_max = 0 ;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(my_fsc)
	{
		if (DIRECT_A1D_ELEM(my_fsc, i) < 1e-10)
			break;
		ires_max = i;
	}
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
    	int ires = ROUND(sqrt((double)kp * kp + ip * ip + jp * jp));
		if (ires <= ires_max)
		{
	        double fsc = DIRECT_A1D_ELEM(my_fsc, ires);
	        if (fsc < 1e-10)
	        	REPORT_ERROR("Postprocessing::applyFscWeighting BUG: fsc <= 0");
	        DIRECT_A3D_ELEM(FT, k, i, j) *= sqrt((2 * fsc) / (1 + fsc));
		}
		else
		{
			DIRECT_A3D_ELEM(FT, k, i, j) = 0.;
		}
	}

}

void Postprocessing::makeGuinierPlot(MultidimArray<Complex > &FT, std::vector<fit_point2D> &guinier)
{

	MultidimArray<int> radial_count(XSIZE(FT));
	MultidimArray<double> lnF(XSIZE(FT));
	fit_point2D      onepoint;

	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
    	int r2 = kp * kp + ip * ip + jp * jp;
    	int ires = ROUND(sqrt((double)r2));
		if (ires < XSIZE(radial_count))
		{

	        lnF(ires) += abs(DIRECT_A3D_ELEM(FT, k, i, j));
	        radial_count(ires)++;
		}
	}

	double xsize = XSIZE(I1());
	guinier.clear();
	FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_count)
	{

		double res = (xsize * angpix)/(double)i; // resolution in Angstrom
		if (res >= angpix * 2.) // Apply B-factor sharpening until Nyquist, then low-pass filter later on (with a soft edge)
        {
            onepoint.x = 1. / (res * res);
            if (DIRECT_A1D_ELEM(lnF, i) > 0.)
            {
                onepoint.y = log ( DIRECT_A1D_ELEM(lnF, i) / DIRECT_A1D_ELEM(radial_count, i) );
                if (res <= fit_minres && res >= fit_maxres)
                {
                    onepoint.w = 1.;
                }
                else
                {
                    onepoint.w = 0.;
                }
            }
            else
            {
                onepoint.y = -99.;
                onepoint.w = 0.;
            }
            //std::cerr << " onepoint.x= " << onepoint.x << " onepoint.y= " << onepoint.y << " onepoint.w= " << onepoint.w << std::endl;
            guinier.push_back(onepoint);
        }
	}

}

void Postprocessing::writeOutput()
{

	if (verb > 0)
	{
		std::cout <<"== Writing out put files ..." <<std::endl;
	}
	// Write all relevant output information
	FileName fn_tmp;
	fn_tmp = fn_out + ".mrc";

	double avg, stddev, minval, maxval;
    I1().computeStats(avg, stddev, minval, maxval);
    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
    I1.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, angpix);
    I1.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, angpix);
    I1.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, angpix);
	I1.write(fn_tmp);
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left   <<"  + Processed map: "; std::cout << fn_tmp<< std::endl;
	}

	// Also write the masked postprocessed map
	if (do_auto_mask || fn_mask != "")
	{
		fn_tmp = fn_out + "_masked.mrc";
		I1() *= Im();
	    I1().computeStats(avg, stddev, minval, maxval);
	    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
	    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
	    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
	    I1.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
		I1.write(fn_tmp);
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left   <<"  + Processed masked map: "; std::cout << fn_tmp<< std::endl;
		}
	}

	// Also write mask
	if (do_auto_mask)
	{
		fn_tmp = fn_out + "_automask.mrc";
		Im().computeStats(avg, stddev, minval, maxval);
		Im.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
		Im.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
		Im.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
		Im.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
		Im.write(fn_tmp);
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left   <<"  + Auto-mask: "; std::cout << fn_tmp<< std::endl;
		}
	}

	// Write an output STAR file with FSC curves, Guinier plots etc
	std::ofstream  fh;
	fn_tmp = fn_out + ".star";
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left <<"  + Metadata file: "; std::cout << fn_tmp<< std::endl;
	}

	fh.open((fn_tmp).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"MlOptimiser::write: Cannot write file: " + fn_tmp);

	// Write the command line as a comment in the header
	fh << "# RELION postprocess" << std::endl;
	fh << "# ";
	parser.writeCommandLine(fh);

	MetaDataTable MDlist, MDfsc, MDguinier;

	MDlist.setIsList(true);
	MDlist.setName("general");
	MDlist.addObject();
	MDlist.setValue(EMDL_POSTPROCESS_FINAL_RESOLUTION, global_resol);
	MDlist.setValue(EMDL_POSTPROCESS_BFACTOR, global_bfactor );
	if (do_auto_bfac)
	{
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_SLOPE, global_slope);
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, global_intercept);
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION, global_corr_coeff);
	}
	MDlist.write(fh);

	MDfsc.setName("fsc");
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
	{
		MDfsc.addObject();
		double res = (i > 0) ? (XSIZE(I1()) * angpix / (double)i) : 999.;
		MDfsc.setValue(EMDL_SPECTRAL_IDX, (int)i);
		MDfsc.setValue(EMDL_RESOLUTION, 1./res);
		MDfsc.setValue(EMDL_RESOLUTION_ANGSTROM, res);
		if (do_mask)
		{
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_TRUE, DIRECT_A1D_ELEM(fsc_true, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_UNMASKED, DIRECT_A1D_ELEM(fsc_unmasked, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_MASKED, DIRECT_A1D_ELEM(fsc_masked, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_RANDOM_MASKED, DIRECT_A1D_ELEM(fsc_random_masked, i) );
		}
		else
		{
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_UNMASKED, DIRECT_A1D_ELEM(fsc_true, i) );
		}
	}
	MDfsc.write(fh);

	// Also write XML file with FSC_true curve for EMDB submission
	writeFscXml(MDfsc);

	MDguinier.setName("guinier");
	for (int i = 0; i < guinierin.size(); i++)
	{
		MDguinier.addObject();
		MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, guinierin[i].x);
		MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_IN, guinierin[i].y);
		if (fn_mtf != "")
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF, guinierinvmtf[i].y);
		if (do_fsc_weighting)
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED, guinierweighted[i].y);
		if (do_auto_bfac || ABS(adhoc_bfac) > 0.)
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED, guiniersharpen[i].y);
		if (do_auto_bfac)
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT, global_intercept);
	}
	MDguinier.write(fh);

	fh.close();

	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left   <<"  + FINAL RESOLUTION: "; std::cout << global_resol<< std::endl;
	}

}

void Postprocessing::writeFscXml(MetaDataTable &MDfsc)
{

    FileName fn_fsc = fn_out + "_fsc.xml";
	std::ofstream  fh;
    fh.open((fn_fsc).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_fsc);

    fh << "<fsc title=\"RELION masked-corrected FSC\" xaxis=\"Resolution (A-1)\" yaxis=\"Correlation Coefficient\">"<<std::endl;

    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDfsc)
    {
    	double xx, yy;
    	MDfsc.getValue(EMDL_RESOLUTION, xx);
    	MDfsc.getValue(EMDL_POSTPROCESS_FSC_TRUE, yy);
    	fh << "  <coordinate>" << std::endl;
    	fh << "    <x>" << xx << "</x>" << std::endl;
    	fh << "    <y>" << yy << "</y>" << std::endl;
    	fh << "  </coordinate>" << std::endl;
    }
    fh << "</fsc>" << std::endl;
    fh.close();
    std::cout.width(35); std::cout << std::left <<"  + EMDB xml-format FSC curve: "; std::cout << fn_fsc << std::endl;
}

void Postprocessing::run()
{

	// Read inout maps and perform some checks
	initialise();

	// Calculate FSC of the unmask maps
	getFSC(I1(), I2(), fsc_unmasked);

	// Check whether we'll do masking
	do_mask = getMask();
	if (do_mask)
	{
		if (verb > 0)
		{
			std::cout <<"== Masking input maps ..." <<std::endl;
		}
		// Mask I1 and I2 and calculated fsc_masked
		I1() *= Im();
		I2() *= Im();
		getFSC(I1(), I2(), fsc_masked);

		// To save memory re-read the same input maps again and randomize phases before masking
		I1.read(fn_I1);
		I2.read(fn_I2);
		I1().setXmippOrigin();
		I2().setXmippOrigin();

		// Check at which resolution shell the FSC drops below randomize_fsc_at
		int randomize_at = -1;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_unmasked)
		{
			if (i > 0 && DIRECT_A1D_ELEM(fsc_unmasked, i) < randomize_fsc_at)
			{
				randomize_at = i;
				break;
			}
		}
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left << "  + randomize phases beyond: "; std::cout << XSIZE(I1())* angpix / randomize_at << " Angstroms" << std::endl;
		}
		if (randomize_at > 0)
		{
			randomizePhasesBeyond(I1(), randomize_at);
			randomizePhasesBeyond(I2(), randomize_at);
			// Mask randomized phases maps and calculated fsc_random_masked
			I1() *= Im();
			I2() *= Im();
			getFSC(I1(), I2(), fsc_random_masked);

		}
		else
			REPORT_ERROR("Postprocessing::run ERROR: FSC curve never drops below randomize_fsc_at");

		// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
		// FSC_true = FSC_t - FSC_n / ( )
		fsc_true.resize(fsc_masked);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
		{
			if (i < randomize_at)
			{
				DIRECT_A1D_ELEM(fsc_true, i) = DIRECT_A1D_ELEM(fsc_masked, i);
			}
			else
			{
				double fsct = DIRECT_A1D_ELEM(fsc_masked, i);
				double fscn = DIRECT_A1D_ELEM(fsc_random_masked, i);
				if (fscn > fsct)
					DIRECT_A1D_ELEM(fsc_true, i) = 0.;
				else
					DIRECT_A1D_ELEM(fsc_true, i) = (fsct - fscn) / (1. - fscn);
			}
		}

		// Now re-read the original maps yet again into memory
		I1.read(fn_I1);
		I2.read(fn_I2);
		I1().setXmippOrigin();
		I2().setXmippOrigin();

	}
	else
	{
		fsc_true = fsc_unmasked;
	}

	global_resol = 999.;
	// See where corrected FSC drops below 0.143
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
	{
		if ( DIRECT_A1D_ELEM(fsc_true, i) < 0.143)
			break;
		global_resol = XSIZE(I1())*angpix/(double)i;
	}

	// Add the two half-maps together for subsequent sharpening
	I1() += I2();

	// Divide by MTF and perform FSC-weighted B-factor sharpening, as in Rosenthal and Henderson, 2003
	// also low-pass filters...
	sharpenMap();

	// Write original and corrected FSC curve, Guinier plot, etc.
	writeOutput();

}
