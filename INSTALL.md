Installation
============

Please edit the file src/Makefile to include the appropriate paths and compilers for your system. 

You will need the following packages installed: 

- CFITSIO http://heasarc.gsfc.nasa.gov/fitsio/
- CLAPACK http://www.netlib.org/clapack/
- GSL http://www.gnu.org/software/gsl/
- FFTW http://www.fftw.org

Then edit the include/shear3d.h file so that the data directories and products directories point to you local directories where you store your data 
and where you want to store the data vectors respectively.

Test/Demo 
=========

To test the data vector code please download the following file 

	https://goo.gl/KToz70

place this in your data directory. Warning this is a large file of ~2Gb. This data is based on the CFHTLenS (http://cfhtlens.org) W1 field. 

Then run 

	/bin/smk_transform_data <name> 

<name> refers to the characters that occur before *passcat.cat. The column structure of the catalogue can be determined 
by looking at src/smk_transform_data_v1.c

This should produce a series of outputs in the products/ directory. These are: 

The E and B mode real and imaginary data vectors: 

	products/3dtransform_1iB.<name>.fits   
	products/3dtransform_1iE.<name>.fits   
	products/3dtransform_1rB.<name>.fits   
	products/3dtransform_1rE.<name>.fits

Abscia of the k and l-modes: 

	products/3dtransform_kll.<name>.fits   
	products/3dtransform_ell.<name>.fits   
	
Real and imaginary l-mode signal weights: 

	products/3dtransform_drll.T1.fits     
	products/3dtransform_dill.<name>.fits  

Real and imaginary l-mode noise weights: 

	products/3dtransform_tll.<name>.fits
	products/3dtransform_sll.<name>.fits

Comoving distance values and redshifts for each object: 

	products/3dtransform_rll.<name>.fits   
	products/3dtransform_zll.<name>.fits

Data vector weighting (if needed):

	products/3dtransform_xll.<name>.fits

2-element vector with sigmae and area (can be supplemented with other ancillary information): 

	products/3dtransform_values.<name>.fits

