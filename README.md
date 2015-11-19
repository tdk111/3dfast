3Dfast: 3D Cosmic Shear 
===================================================

3Dfast is open-source software for calculating the 3D harmonic-space 
power spectrum of a spin-2 field. This is developed within the context 
of weak gravitational lensing calculations. 

It was used to generate the 3D cosmic shear power spectra for these papers:

- 3D Cosmic Shear: Cosmology from CFHTLenS http://arxiv.org/abs/1401.6842

Details fo the code are as follows:

- Language: C
- Libraries used: CFITSIO, GSL, CLAPACK
- Other public code used: CAMB
- Architecture: The code is parallelisable and uses in-code POSIX threads

The code has the following functionality, it can

- generate 3D cosmic shear power spectra from data
- create Fisher matrix predictions for cosmic shear surveys
- search cosmological parameter likelihood using an MCMC metropolis-hastings algorithm

Distribution
------------

The current released version of 3dfast is version 1.0.  

Also, feel free to fork the repository:

    https://github.com/tdk111/3dfast

Or clone the repository with either of the following:

    git clone https://github.com/tdk111/3dfast

See INSTALL.md for more information.

The code is licensed under a BSD-style license.  See the file LICENSE for more
details.

Keeping up-to-date with 3dfast
------------------------------

There is a Slack Channel here on which questions can be posted

      https://astro-informatics.slack.com

Please email the address below with questions

      tom.kitching-3dfast@gmail.com


Installation
------------

For installation instructions, please see the file `INSTALL.md` in the main
repository directory.


Getting started
---------------

* Install the code as in `INSTALL.md`.

* Download the test data from here <link> 

Reference documentation
-----------------------

Documentaiton is currently under construction, however theory and development of the 3D cosmic shear 
method can be found in the following papers: 

- Kitching et al., 2015, http://arxiv.org/pdf/1408.7052
- Kitching et al., 2015, http://arxiv.org/pdf/1401.6842.pdf 	
- Merkel & Schaefer (2013) http://arxiv.org/pdf/1306.6466
- Grassi & Schaefer (2013) http://arxiv.org/pdf/1303.1024
- Kitching, Heavens, Miller (2010) http://arxiv.org/pdf/1007.2953.pdf 
- Camera et al. (2010) http://arxiv.org/pdf/1002.4740
- Jimenez et al. (2010) http://arxiv.org/pdf/1003.5918
- Kitching, Taylor, Heavens (2008) http://arxiv.org/pdf/0801.3270
- Kitching et al. (2008) http://arxiv.org/pdf/0801.4565
- Kitching (2007) http://www-astro.physics.ox.ac.uk/~tdk/Thomas%20Kitching's%20Home%20Page/Thomas%20Kitching's%20Home%20Page_files/thesis.pdf
- Kitching et al. (2007) http://arxiv.org/pdf/astro-ph/0610284.pdf
- Heavens et al. (2006) http://arxiv.org/pdf/astro-ph/0606568
- Heavens (2003) http://arxiv.org/pdf/astro-ph/0304151

Repository directory structure
------------------------------

The repository has a number of subdirectories. Below is a guide to their
contents:

* bin/ :      executables (after the compilation procedure is done).
* include/ :  the .h header files for the C parts of 3dfast.
* src/ :      the source code for the purely C++ parts of 3dfast.
* products/ : will contain FITS output from the data vector code, used in the theory code
* data/ :     contains input data in form of ascii file

