#ifndef shear3d_h
#define shear3d_h
#define NUM_THREAD 80       /*number of threads*/
#define MAX_THREAD 100000   /*maximum number of threads allowed by the system*/

#define PLOTTING   0 /*if want to plot in later codes*/
#define PRINTVALS  0 /*1=verbose*/

#define MACHINE_P  1e-30 /*machine precision*/

#define DATADIR    "/home/tdk/code/shear3d/data/"
#define PRODUCTDIR "/home/tdk/code/shear3d/products/"

/*standard variables*/
#define pi 3.1415926
#define light 299792.458 //Kms^-1
#define allsky 41253.  //all sky in square degress

/*safe l ranges*/
#define safelmin 500. 
#define safelmax 5000. 
#define fn_l     128   //specified number of ell modes 128                                                                                               
/*define k range*/
#define fkmin 1e-3
#define fkmax 5.0    //5.
#define fn_k  50    //50
#define fn_ri 10000  //internal U matrix r integral sampling needs to be set to >=2000 

/*flags*/
#define logl_flag 0        //0=linear, 1=log
#define logk_flag 0        //0=linear, 1=log
#define pzs       0        //full p(z) or not

/*fiducial cosmology*/
#define flatcosmo 1     //==1 if flat cosmology
#define cosmo_om  0.3175
#define cosmo_od  0.6825
#define cosmo_ob  0.046 
#define cosmo_w0 -1.00  //w(z)=w0+wa(1-a)
#define cosmo_wa  0.00  
#define cosmo_s8  0.834 //set s8=-1 for CMB normalisation
#define cosmo_ns  0.9624
#define cosmo_h0  0.701 //Reiss et al. 2011 1103.2976
#define cosmo_gm  1.00
#define cosmo_nm  0.00  //eV

/*define some structures*/
/*cosmology*/
struct cosmo {double om; double od; 
  double ob; double w0; 
  double wa; double s8; 
  double ns; double h0; 
  double gm; double norm;
  double nm;};
/*fiducial*/
struct fid {int fit_nl; int fit_tk; int normalise; 
  double kmin; double kmax; int n_k; int n_ki; int k_flag; 
  double rmin; double rmax; int n_r; int n_ri;
  double zmin; double zmax; int n_z;
  int    lmin; int    lmax; int n_l; int l_flag; double dl;
  int n_coarse; 
  double zmedian; double sigmaz0; double area; double n0; double sigmae; 
  int    plot; int    print; int n_linterp;};

/* common global thread variables*/
int *nmodl0;
double **thread_Xsum; 
double ***Plkr;
double **glvalues,**gkvalues;

#endif
