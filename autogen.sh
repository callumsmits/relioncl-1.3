#Some path variables
setenv RELION_HOME $PWD
setenv PREFIX $RELION_HOME
autoreconf --install
fltk_cxx=`fltk-config --cxxflags`
fltk_ld=`fltk-config --ldflags`
 ./configure --enable-mpi CPPFLAGS="-I/lmb/home/scheres/app/fftw-3.2.2/include -I$PREFIX/include $fltk_cxx" LDFLAGS="-L/lmb/home/scheres/app/fftw-3.2.2/lib -L$PREFIX/lib $fltk_ld" 

