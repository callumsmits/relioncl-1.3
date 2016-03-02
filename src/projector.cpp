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
#include "src/projector.h"
//#define DEBUG

void Projector::initialiseData(int current_size)
{
	// By default r_max is half ori_size
	if (current_size < 0)
		r_max = ori_size / 2;
	else
		r_max = current_size / 2;

	// Never allow r_max beyond Nyquist...
	r_max = XMIPP_MIN(r_max, ori_size / 2);

	// Set pad_size
	pad_size = 2 * (padding_factor * r_max + 1) + 1;

	// Short side of data array
	switch (ref_dim)
	{
	case 2:
	   data.resize(pad_size, pad_size / 2 + 1);
	   break;
	case 3:
	   data.resize(pad_size, pad_size, pad_size / 2 + 1);
	   break;
	default:
	   REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// Set origin in the y.z-center, but on the left side for x.
	data.setXmippOrigin();
	data.xinit=0;

}

void Projector::initZeros(int current_size)
{
	initialiseData(current_size);
	data.initZeros();
}

long int Projector::getSize()
{
	// Short side of data array
	switch (ref_dim)
	{
		case 2:
		   return pad_size * (pad_size / 2 + 1);
		   break;
		case 3:
		   return pad_size * pad_size * (pad_size / 2 + 1);
		   break;
		default:
		   REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");
	}

}

// Fill data array with oversampled Fourier transform, and calculate its power spectrum
void Projector::computeFourierTransformMap(MultidimArray<double> &vol_in, MultidimArray<double> &power_spectrum, int current_size, int nr_threads, bool do_gridding)
{

	MultidimArray<double> Mpad;
	MultidimArray<Complex > Faux;
    FourierTransformer transformer;
    // DEBUGGING: multi-threaded FFTWs are giving me a headache?
	// For a long while: switch them off!
	//transformer.setThreadsNumber(nr_threads);
    double normfft;

	// Size of padded real-space volume
	int padoridim = padding_factor * ori_size;

	// Initialize data array of the oversampled transform
	ref_dim = vol_in.getDim();

	// Make Mpad
	switch (ref_dim)
	{
	case 2:
	   Mpad.initZeros(padoridim, padoridim);
	   normfft = (double)(padding_factor * padding_factor);
	   break;
	case 3:
	   Mpad.initZeros(padoridim, padoridim, padoridim);
	   normfft = (double)(padding_factor * padding_factor * padding_factor * ori_size);
	   break;
	default:
	   REPORT_ERROR("Projector::get2DSlice%%ERROR: Dimension of the data array should be 2 or 3");
	}

	// First do a gridding pre-correction on the real-space map:
	// Divide by the inverse Fourier transform of the interpolator in Fourier-space
	// 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
	// Still to check the 3D case...
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	/////////////////// TODO: removed this for subtomo averaging!!!! Go back and check 3d-> 2d projection!!!!!
	if (do_gridding)
		griddingCorrect(vol_in);

	// Pad translated map with zeros
	vol_in.setXmippOrigin();
	Mpad.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
		A3D_ELEM(Mpad, k, i, j) = A3D_ELEM(vol_in, k, i, j);

	// Translate padded map to put origin of FT in the center
	CenterFFT(Mpad, true);

	// Calculate the oversampled Fourier transform
	transformer.FourierTransform(Mpad, Faux, false);

	// Free memory: Mpad no longer needed
	Mpad.clear();

	// Resize data array to the right size and initialise to zero
	initZeros(current_size);

	// Fill data only for those points with distance to origin less than max_r
	// (other points will be zero because of initZeros() call above
	// Also calculate radial power spectrum
	power_spectrum.initZeros(ori_size / 2 + 1);
	MultidimArray<double> counter(power_spectrum);
	counter.initZeros();

	int max_r2 = r_max * r_max * padding_factor * padding_factor;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) // This will also work for 2D
	{
		int r2 = kp*kp + ip*ip + jp*jp;
		// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
		if (r2 <= max_r2)
		{
			// Set data array
			A3D_ELEM(data, kp, ip, jp) = DIRECT_A3D_ELEM(Faux, k, i, j) * normfft;

			// Calculate power spectrum
			int ires = ROUND( sqrt((double)r2) / padding_factor );
			// Factor two because of two-dimensionality of the complex plane
			DIRECT_A1D_ELEM(power_spectrum, ires) += norm(A3D_ELEM(data, kp, ip, jp)) / 2.;
			DIRECT_A1D_ELEM(counter, ires) += 1.;
		}
	}

	// Calculate radial average of power spectrum
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(power_spectrum)
	{
		if (DIRECT_A1D_ELEM(counter, i) < 1.)
			DIRECT_A1D_ELEM(power_spectrum, i) = 0.;
		else
			DIRECT_A1D_ELEM(power_spectrum, i) /= DIRECT_A1D_ELEM(counter, i);
	}

	transformer.cleanup();

}

void Projector::griddingCorrect(MultidimArray<double> &vol_in)
{
	// Correct real-space map by dividing it by the Fourier transform of the interpolator(s)
	vol_in.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in)
	{
		double r = sqrt((double)(k*k+i*i+j*j));
		// if r==0: do nothing (i.e. divide by 1)
		if (r > 0.)
		{
			double rval = r / (ori_size * padding_factor);
			double sinc = sin(PI * rval) / ( PI * rval);
			//double ftblob = blob_Fourier_val(rval, blob) / blob_Fourier_val(0., blob);
			// Interpolation (goes with "interpolator") to go from arbitrary to fine grid
			if (interpolator==NEAREST_NEIGHBOUR && r_min_nn == 0)
			{
				// NN interpolation is convolution with a rectangular pulse, which FT is a sinc function
            	A3D_ELEM(vol_in, k, i, j) /= sinc;
			}
			else if (interpolator==TRILINEAR || (interpolator==NEAREST_NEIGHBOUR && r_min_nn > 0) )
			{
				// trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
				A3D_ELEM(vol_in, k, i, j) /= sinc * sinc;
			}
			else
				REPORT_ERROR("BUG Projector::griddingCorrect: unrecognised interpolator scheme.");
//#define DEBUG_GRIDDING_CORRECT
#ifdef DEBUG_GRIDDING_CORRECT
			if (k==0 && i==0 && j > 0)
				std::cerr << " j= " << j << " sinc= " << sinc << std::endl;
#endif
		}
	}
}

void Projector::project(MultidimArray<Complex > &f2d, Matrix2D<double> &A, bool inv)
{
	double fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1, y, y2, r2;
	bool is_neg_x;
	Complex d000, d001, d010, d011, d100, d101, d110, d111;
	Complex dx00, dx01, dx10, dx11, dxy0, dxy1;
	Matrix2D<double> Ainv;

    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside r_max should already be zero...
    // f2d.initZeros();

	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f2d) - 1);

    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (double)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;

//#define DEBUG
#ifdef DEBUG
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif

	for (int i=0; i < YSIZE(f2d); i++)
	{
		// Dont search beyond square with side max_r
		if (i <= my_r_max)
		{
			y = i;
		}
		else if (i >= YSIZE(f2d) - my_r_max)
		{
			y = i - YSIZE(f2d);
		}
		else
			continue;

		y2 = y * y;
		for (int x=0; x <= my_r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get logical coordinates in the 3D map
			xp = Ainv(0,0) * x + Ainv(0,1) * y;
			yp = Ainv(1,0) * x + Ainv(1,1) * y;
			zp = Ainv(2,0) * x + Ainv(2,1) * y;

			if (interpolator == TRILINEAR || r2 < min_r2_nn)
			{

				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					zp = -zp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
    			x0 = FLOOR(xp);
				fx = xp - x0;
				x1 = x0 + 1;

				y0 = FLOOR(yp);
				fy = yp - y0;
				y0 -=  STARTINGY(data);
				y1 = y0 + 1;

				z0 = FLOOR(zp);
				fz = zp - z0;
				z0 -= STARTINGZ(data);
				z1 = z0 + 1;

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
				d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
				d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
				d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
				d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
				d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
				d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
				d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

				// Set the interpolated value in the 2D output array
				dx00 = LIN_INTERP(fx, d000, d001);
				dx01 = LIN_INTERP(fx, d100, d101);
				dx10 = LIN_INTERP(fx, d010, d011);
				dx11 = LIN_INTERP(fx, d110, d111);
				dxy0 = LIN_INTERP(fy, dx00, dx10);
				dxy1 = LIN_INTERP(fy, dx01, dx11);
				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fz, dxy0, dxy1);

				// Take complex conjugated for half with negative x
				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));

			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				x0 = ROUND(xp);
				y0 = ROUND(yp);
				z0 = ROUND(zp);
				if (x0 < 0)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(A3D_ELEM(data, -z0, -y0, -x0));
				else
					DIRECT_A2D_ELEM(f2d, i, x) = A3D_ELEM(data, z0, y0, x0);

			} // endif NEAREST_NEIGHBOUR
			else
				REPORT_ERROR("Unrecognized interpolator in Projector::project");

		} // endif x-loop
	} // endif y-loop


#ifdef DEBUG
    std::cerr << "done with project..." << std::endl;
#endif
}


void Projector::projectD(MultidimArray<Complex > &f2d, Matrix2D<double> &A, bool inv)
{
    double fx, fy, fz, xp, yp, zp;
    int x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    Complex d000, d001, d010, d011, d100, d101, d110, d111;
    Complex dx00, dx01, dx10, dx11, dxy0, dxy1;
    Matrix2D<double> Ainv;
    
    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside r_max should already be zero...
    // f2d.initZeros();
    
    // Use the inverse matrix
    if (inv)
        Ainv = A;
    else
        Ainv = A.transpose();
    
    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f2d) - 1);
    
    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (double)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
    //#define DEBUG
#ifdef DEBUG
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif
    
    for (int i=0; i < YSIZE(f2d); i++)
    {
        // Dont search beyond square with side max_r
        if (i <= my_r_max)
        {
            y = i;
        }
        else if (i >= YSIZE(f2d) - my_r_max)
        {
            y = i - YSIZE(f2d);
        }
        else
            continue;
        
        y2 = y * y;
        for (int x=0; x <= my_r_max; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get logical coordinates in the 3D map
            xp = Ainv(0,0) * x + Ainv(0,1) * y;
            yp = Ainv(1,0) * x + Ainv(1,1) * y;
            zp = Ainv(2,0) * x + Ainv(2,1) * y;
            
            
            if (interpolator == TRILINEAR || r2 < min_r2_nn)
            {
                
                // Only asymmetric half is stored
                if (xp < 0)
                {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                    is_neg_x = true;
                }
                else
                {
                    is_neg_x = false;
                }
                
                // Trilinear interpolation (with physical coords)
                // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                x0 = FLOOR(xp);
                fx = xp - x0;
                x1 = x0 + 1;
                
                y0 = FLOOR(yp);
                fy = yp - y0;
                y0 -=  STARTINGY(data);
                y1 = y0 + 1;
                
                z0 = FLOOR(zp);
                fz = zp - z0;
                z0 -= STARTINGZ(data);
                z1 = z0 + 1;
                
                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
                d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
                d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
                d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
                d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
                d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
                d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
                d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);
                
//                if (i == 1) {
//                    std::cerr << "x0:" << x0 << " y0: " << y0 << " z0: " << z0 << " d000: " << d000.real << " " << d000.imag << " d001: " << d001.real << " " << d001.imag << std::endl;
//                }

                // Set the interpolated value in the 2D output array
                dx00 = LIN_INTERP(fx, d000, d001);
                dx01 = LIN_INTERP(fx, d100, d101);
                dx10 = LIN_INTERP(fx, d010, d011);
                dx11 = LIN_INTERP(fx, d110, d111);
                dxy0 = LIN_INTERP(fy, dx00, dx10);
                dxy1 = LIN_INTERP(fy, dx01, dx11);
                DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fz, dxy0, dxy1);

//                if (i == 1) {
//                    std::cerr << "x0:" << x0 << " y0: " << y0 << " z0: " << z0 << " dx00: " << dx00.real << " " << dx00.imag << " dx01: " << dx01.real << " " << dx01.imag << " dx10: " << dx10.real << " " << dx10.imag << " dx11: " << dx11.real << " " << dx11.imag << " dxy0: " << dxy0.real << " " << dxy0.imag << " dxy1: " << dxy1.real << " " << dxy1.imag << std::endl;
//                    std::cerr << "x0: " << x0 << " y0: " << y0 << " z0: " << z0 << " i: " << i << " x " << x << " result: " << LIN_INTERP(fz, dxy0, dxy1).real << " " << LIN_INTERP(fz, dxy0, dxy1).imag << std::endl;
//                }

                
                // Take complex conjugated for half with negative x
                if (is_neg_x)
                    DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
                
            } // endif TRILINEAR
            else if (interpolator == NEAREST_NEIGHBOUR )
            {
                x0 = ROUND(xp);
                y0 = ROUND(yp);
                z0 = ROUND(zp);
                if (x0 < 0)
                    DIRECT_A2D_ELEM(f2d, i, x) = conj(A3D_ELEM(data, -z0, -y0, -x0));
                else
                    DIRECT_A2D_ELEM(f2d, i, x) = A3D_ELEM(data, z0, y0, x0);
                
            } // endif NEAREST_NEIGHBOUR
            else
                REPORT_ERROR("Unrecognized interpolator in Projector::project");
            
        } // endif x-loop
    } // endif y-loop
    
    
#ifdef DEBUG
    std::cerr << "done with project..." << std::endl;
#endif
}

void Projector::rotate2D(MultidimArray<Complex > &f2d, Matrix2D<double> &A, bool inv)
{
	double fx, fy, xp, yp;
	int x0, x1, y0, y1, y, y2, r2;
	bool is_neg_x;
	Complex d00, d01, d10, d11, dx0, dx1;
	Matrix2D<double> Ainv;

    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    // f2d.initZeros();
	// Use the inverse matrix
    if (inv)
    	Ainv = A;
    else
    	Ainv = A.transpose();

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    int my_r_max = XMIPP_MIN(r_max, XSIZE(f2d) - 1);

    // Go from the 2D slice coordinates to the map coordinates
    Ainv *= (double)padding_factor;  // take scaling into account directly
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
#ifdef DEBUG
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
#endif
	for (int i=0; i < YSIZE(f2d); i++)
	{
		// Don't search beyond square with side max_r
		if (i <= my_r_max)
		{
			y = i;
		}
		else if (i >= YSIZE(f2d) - my_r_max)
		{
			y = i - YSIZE(f2d);
		}
		else
			continue;
		y2 = y * y;
		for (int x=0; x <= my_r_max; x++)
		{
	    	// Only include points with radius < max_r (exclude points outside circle in square)
			r2 = x * x + y2;
			if (r2 > max_r2)
				continue;

			// Get logical coordinates in the 3D map
			xp = Ainv(0,0) * x + Ainv(0,1) * y;
			yp = Ainv(1,0) * x + Ainv(1,1) * y;
			if (interpolator == TRILINEAR || r2 < min_r2_nn)
			{
				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

				// Trilinear interpolation (with physical coords)
				// Subtract STARTINGY to accelerate access to data (STARTINGX=0)
				// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
    			x0 = FLOOR(xp);
				fx = xp - x0;
				x1 = x0 + 1;

				y0 = FLOOR(yp);
				fy = yp - y0;
				y0 -=  STARTINGY(data);
				y1 = y0 + 1;

				// Matrix access can be accelerated through pre-calculation of z0*xydim etc.
				d00 = DIRECT_A2D_ELEM(data, y0, x0);
				d01 = DIRECT_A2D_ELEM(data, y0, x1);
				d10 = DIRECT_A2D_ELEM(data, y1, x0);
				d11 = DIRECT_A2D_ELEM(data, y1, x1);

				// Set the interpolated value in the 2D output array
				dx0 = LIN_INTERP(fx, d00, d01);
				dx1 = LIN_INTERP(fx, d10, d11);
				DIRECT_A2D_ELEM(f2d, i, x) = LIN_INTERP(fy, dx0, dx1);
				// Take complex conjugated for half with negative x
				if (is_neg_x)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(DIRECT_A2D_ELEM(f2d, i, x));
			} // endif TRILINEAR
			else if (interpolator == NEAREST_NEIGHBOUR )
			{
				x0 = ROUND(xp);
				y0 = ROUND(yp);
				if (x0 < 0)
					DIRECT_A2D_ELEM(f2d, i, x) = conj(A2D_ELEM(data, -y0, -x0));
				else
					DIRECT_A2D_ELEM(f2d, i, x) = A2D_ELEM(data, y0, x0);
			} // endif NEAREST_NEIGHBOUR
			else
				REPORT_ERROR("Unrecognized interpolator in Projector::project");
		} // endif x-loop
	} // endif y-loop
}
