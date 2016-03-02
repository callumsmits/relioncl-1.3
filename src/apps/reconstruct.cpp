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

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>

class reconstruct_parameters
{
	public:
   	FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_fsc, fn_debug;
	MetaDataTable DF;
	int r_max, r_min_nn, blob_order, padding_factor, ref_dim, interpolator, iter, nr_threads, debug_ori_size, debug_size, ctf_dim;
	double blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres, beamtilt_x, beamtilt_y;
	bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, do_fom_weighting, do_beamtilt;
	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		fn_debug = parser.getOption("--debug", "Rootname for debug reconstruction files", "");
		debug_ori_size =  textToInteger(parser.getOption("--debug_ori_size", "Rootname for debug reconstruction files", "1"));
		debug_size =  textToInteger(parser.getOption("--debug_size", "Rootname for debug reconstruction files", "1"));

		fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
	    fn_out = parser.getOption("--o", "Name for output reconstruction","relion.mrc");
	    fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
       	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
       	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
       	padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));
    	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to use for FFTs", "1"));

	    int ctf_section = parser.addSection("CTF options");
       	do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
    	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
    	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
    	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
    	beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
    	beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
    	do_beamtilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

       	int expert_section = parser.addSection("Expert options");
       	fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");
        if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
       		interpolator = NEAREST_NEIGHBOUR;
       	else
       		interpolator = TRILINEAR;
        blob_radius   = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
        blob_order    = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
        blob_alpha    = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
       	iter = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
       	ref_dim = textToInteger(parser.getOption("--refdim", "Dimension of the reconstruction (2D or 3D)", "3"));
    	angular_error = textToFloat(parser.getOption("--angular_error", "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles", "0."));
    	shift_error = textToFloat(parser.getOption("--shift_error", "Apply random deviations with this standard deviation (in pixels) to each of the 2 translations", "0."));
    	do_fom_weighting = parser.checkOption("--fom_weighting", "Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)");
    	fn_fsc = parser.getOption("--fsc", "FSC-curve for regularized reconstruction", "");

    	// Hidden
       	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

    	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    	// Read MetaData file, which should have the image names and their angles!
    	if (fn_debug == "")
    		DF.read(fn_sel);

     	randomize_random_generator();

     	if (do_beamtilt && ! do_ctf)
     		REPORT_ERROR("ERROR: one can only correct for beamtilt in combination with CTF correction!");

	}

	void reconstruct()
	{

		if (fn_debug != "")
		{
			BackProjector backprojector(debug_ori_size, 3, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha);
			backprojector.initialiseDataAndWeight(debug_size);
			backprojector.data.printShape();
			backprojector.weight.printShape();
			Image<double> It;
			It.read(fn_debug+"_data_real.mrc");
			It().setXmippOrigin();
			It().xinit=0;

			It().printShape();
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.data, k, i, j).real = A3D_ELEM(It(), k, i, j);
			}
			It.read(fn_debug+"_data_imag.mrc");
			It().setXmippOrigin();
			It().xinit=0;
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.data, k, i, j).imag = A3D_ELEM(It(), k, i, j);
			}
			It.read(fn_debug+"_weight.mrc");
			It().setXmippOrigin();
			It().xinit=0;
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.weight, k, i, j) = A3D_ELEM(It(), k, i, j);
			}

			MultidimArray<double> dummy;
			backprojector.reconstruct(It(), iter, false, 1., dummy, dummy, dummy, dummy, 1., false, true, nr_threads, -1);
	    	It.write(fn_out);
	    	std::cerr<<" Done writing map in "<<fn_out<<std::endl;
                exit(1);

		}
		else
		{

		double rot, tilt, psi, fom;
		Matrix2D<double> A3D;
		MultidimArray<Complex > Faux, F2D, Fsub;
		MultidimArray<double> Fweight, Fctf, dummy;
		Image<double> vol, img, sub;
		FourierTransformer transformer;
		Matrix1D< double > trans(2);
		Projector proj;
		int mysize;
//#define DEBUG_WW
#ifdef DEBUG_WW

   		// Get dimension of the images
   		(DF).firstObject();
   		DF.getValue(EMDL_IMAGE_NAME, fn_img);
   		img.read(fn_img);
   		mysize=(int)XSIZE(img());
		BackProjector backprojectort(mysize, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha);
		backprojectort.initZeros(2 * r_max);

		Image<double> Imagn, Iphas, Iw, tvol;
		Imagn.read("FEW_it24_rank2_data_magn.spi");
		Iphas.read("FEW_it24_rank2_data_phas.spi");
		Iw.read("FEW_it24_rank2_weight.spi");
        Iw().setXmippOrigin();
        Iw().xinit=0;

		// Write out invw
		Image<double> oo;
		oo=Iw;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iw())
		{
			if (DIRECT_MULTIDIM_ELEM(Iw(), n) > 1e-2)
				DIRECT_MULTIDIM_ELEM(oo(), n) = 1./ DIRECT_MULTIDIM_ELEM(Iw(), n);
		}
		oo.write("invw.spi");

		Imagn().printShape();
		backprojectort.data.printShape();
		backprojectort.weight.printShape();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imagn())
		{
			double realval = sin(DIRECT_MULTIDIM_ELEM(Iphas(), n)) * DIRECT_MULTIDIM_ELEM(Imagn(), n);
			double imagval = cos(DIRECT_MULTIDIM_ELEM(Iphas(), n)) * DIRECT_MULTIDIM_ELEM(Imagn(), n);
			DIRECT_MULTIDIM_ELEM(backprojectort.data, n) = (Complex)(realval, imagval);
		}
		backprojectort.weight = Iw();
  		std::cerr << "Starting the reconstruction ..." << std::endl;
   		backprojectort.reconstruct(tvol(), iter, false, 1., dummy, dummy, dummy, dummy, 1., false, false, nr_threads, -1);
    	tvol.write(fn_out);
    	std::cerr<<" Done writing TMPPPPPPPPPPPPPPPPP debugging!!!c map in "<<fn_out<<std::endl;
		exit(0);
#endif


   		// Get dimension of the images
		(DF).firstObject();
		DF.getValue(EMDL_IMAGE_NAME, fn_img);
		img.read(fn_img);
		mysize=(int)XSIZE(img());

   		if (DF.containsLabel(EMDL_CTF_MAGNIFICATION) && DF.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		double mag, dstep;
   			DF.getValue(EMDL_CTF_MAGNIFICATION, mag);
   			DF.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
   			angpix = 10000. * dstep / mag;
   			std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
    	}

   		BackProjector backprojector(mysize, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha);

   		if (maxres < 0.)
   			r_max = -1;
   		else
   			r_max = CEIL(mysize * angpix / maxres);

   		backprojector.initZeros(2 * r_max);
   		Projector projector(mysize, interpolator, padding_factor, r_min_nn);

   		if (fn_sub != "")
   		{
   			sub.read(fn_sub);
   			projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
   		}

		// Check for beam-tilt parameters in the input star file
   		if (do_beamtilt && ( DF.containsLabel(EMDL_IMAGE_BEAMTILT_X) || DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y) ) )
   				std::cout << " + Using the beamtilt parameters in the input STAR file" << std::endl;

		std::cerr << "Back-projecting all images ..." << std::endl;
   		int imgno = 0;
		time_config();
   		init_progress_bar(DF.size());
   		FOR_ALL_OBJECTS_IN_METADATA_TABLE(DF)
   		{

   			DF.getValue(EMDL_IMAGE_NAME, fn_img);
   			img.read(fn_img);
   			img().setXmippOrigin();

			// Rotations
			if (ref_dim==2)
   			{
   				rot = tilt = 0.;
   			}
   			else
   			{
   				DF.getValue( EMDL_ORIENT_ROT, rot);
   				DF.getValue( EMDL_ORIENT_TILT, tilt);
   			}
  			psi = 0.;
  			DF.getValue( EMDL_ORIENT_PSI, psi);
  			if (angular_error > 0.)
  			{
  				rot += rnd_gaus(0., angular_error);
  				tilt += rnd_gaus(0., angular_error);
  				psi += rnd_gaus(0., angular_error);
  				//std::cout << rnd_gaus(0., angular_error) << std::endl;
  			}

   			Euler_angles2matrix(rot, tilt, psi, A3D);

   			// Translations (either through phase-shifts or in real space
   			trans.initZeros();
   			DF.getValue( EMDL_ORIENT_ORIGIN_X, XX(trans));
   			DF.getValue( EMDL_ORIENT_ORIGIN_Y, YY(trans));
  			if (shift_error > 0.)
   			{
   				XX(trans) += rnd_gaus(0., shift_error);
   				YY(trans) += rnd_gaus(0., shift_error);
   			}

   			if (do_fom_weighting)
   				DF.getValue( EMDL_PARTICLE_FOM, fom);


   			// Use either selfTranslate OR shiftImageInFourierTransform!!
   			//selfTranslate(img(), trans, WRAP);
   			CenterFFT(img(), true);
   			transformer.FourierTransform(img(), F2D);
   			if (ABS(XX(trans)) > 0. || ABS(YY(trans)) > 0.)
   				shiftImageInFourierTransform(F2D, F2D, XSIZE(img()), trans );

			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			// Apply CTF if necessary
			if (do_ctf)
			{
				CTF ctf;
				ctf.read(DF, DF);
				ctf.getFftwImage(Fctf, mysize, mysize, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);

				if (do_beamtilt || (DF.containsLabel(EMDL_IMAGE_BEAMTILT_X) && DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y) ))
				{
					if (DF.containsLabel(EMDL_IMAGE_BEAMTILT_X))
						DF.getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x);
					if (DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y))
						DF.getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y);
					selfApplyBeamTilt(F2D, beamtilt_x, beamtilt_y, ctf.lambda, ctf.Cs, angpix, mysize);
				}
			}

			// Subtract reference projection
			if (fn_sub != "")
			{
				Fsub.resize(F2D);
				projector.get2DFourierTransform(Fsub, A3D, IS_NOT_INV);

				// Apply CTF if necessary
				if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
					{
						DIRECT_MULTIDIM_ELEM(Fsub, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n) -= DIRECT_MULTIDIM_ELEM(Fsub, n);
				}
				// Back-project difference image
				backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV);
			}
			else
			{
				// "Normal" reconstruction, multiply X by CTF, and W by CTF^2
				if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				// Do the following after squaring the CTFs!
				if (do_fom_weighting)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  *= fom;
						DIRECT_MULTIDIM_ELEM(Fctf, n)  *= fom;
					}
				}

//#define DEBUG_RECONSTRUCT_ONLY
#ifdef DEBUG_RECONSTRUCT_ONLY
					if (fn_img == "/lmb/home/scheres/data/betaGal_rh_withnoise/betaGal_2010_all_p_2x2_unflipped/img00001.win100")
					//if (part_id == my_first_particle_id)
					{
						std::cerr << " fn_img= " << fn_img << std::endl;
						//std::cerr << " myscale= " << myscale << std::endl;
						//std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << " normcorr= " << normcorr << std::endl;
						//std::cerr << " sigma2_fudge= " << sigma2_fudge << " mymodel.tau2_fudge_factor= " << mymodel.tau2_fudge_factor<< std::endl;
						//std::cerr << " A3D= " << A3D << std::endl;
						std::cerr << " A3D= " << A3D << std::endl;
						//std::cerr << " exp_R_mic= " << exp_R_mic << std::endl;
						std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << " xoff= "<< XX(trans)<< " yoff= "<<YY(trans)<<std::endl;
						//std::cerr << "mic_id= "<<mic_id<<" mymodel.sigma2_noise[mic_id]= " << mymodel.sigma2_noise[mic_id] << std::endl;
						Image<double> It;
						It()=Fctf;
						It.write("reconstruct_Fctf.spi");
						It().resize(mysize, mysize);
						MultidimArray<Complex > Faux = F2D;
						FourierTransformer transformer;
						transformer.inverseFourierTransform(Faux, It());
						CenterFFT(It(), false);
						It.write("reconstruct_Mimg.spi");
					}
#endif


				backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
			}



   			if (imgno++%60==0)
   				progress_bar(imgno);
		}
   		progress_bar(DF.size());


   		bool do_map = false;
   		bool do_use_fsc = false;
   		MultidimArray<double> fsc;
   		fsc.resize(mysize/2+1);
   		if (fn_fsc != "")
   		{
   			do_map = true;
   			do_use_fsc =true;
   			MetaDataTable MDfsc;
   			MDfsc.read(fn_fsc);
   			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDfsc)
   			{
   				int idx;
   				double val;
   				MDfsc.getValue(EMDL_SPECTRAL_IDX, idx);
   				MDfsc.getValue(EMDL_MLMODEL_FSC_HALVES_REF, val);
   				fsc(idx) =  val;
   			}
   		}
   		std::cerr << "Starting the reconstruction ..." << std::endl;
   		backprojector.reconstruct(vol(), iter, do_map, 1., dummy, dummy, dummy, fsc, 1., do_use_fsc, true, nr_threads, -1);

   		vol.write(fn_out);
    	std::cerr<<" Done writing map in "<<fn_out<<std::endl;

	}
	}


};


int main(int argc, char *argv[])
{
	reconstruct_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.reconstruct();

    }
    catch (RelionError XE)
    {
        prm.usage();
        std::cout << XE;
        exit(1);
    }
    return 0;
}


