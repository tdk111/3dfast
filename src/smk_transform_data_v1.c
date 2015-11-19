#include <stdio.h>
#include <unistd.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>

/*for these include files also need to change Makefile libraries*/
#include <fitsio2.h>
#include "fitsio.h"

/*3d shear header file*/
#include "shear3d.h"

/*cosmological routines*/
double R0(struct cosmo);
double drdz(struct cosmo, double);
double dist(struct cosmo, double);

/*bessels*/
double sjl(int,double);

/*outputs*/
int write3Dmatrix_fits(char*, double***, int, int, int, int);
int write2Dmatrix_fits(char*, double**, int, int, int, int);

/*transform thread*/
void* stransform_thread(void*);
/* structures to pass to transform threads*/
typedef struct {struct fid fid; int obji; int objf; int nt;} threadindextransform;
pthread_mutex_t transformlock;

/* global variables for the transform loop*/
double ***thread_numkvalues;
int    ***thread_nmodl1,**thread_nmodl0;
double ***thread_sum1rE,***thread_sum1iE;
double ***thread_sum1rB,***thread_sum1iB;
double ***thread_summr,***thread_summi;
double ***thread_sum0r,***thread_sum0i;

int nummodl,rcs;
double *e1,*e2,*e1local,*e2local,***be,*thetax,*thetay,*RA,*DEC,*rg,*eweight,*m,*c2,*c1nb,*c1dp,*c2nb,*c2dp,eweightnorm,mweightnorm;
#define sqrt2overpi 0.79788456

/*data*/
double **ppz,**ppzr,*zpz,****thread_Plkr; 

/******************************************************************/

int main(int argc, char *argv[])
{
  
  /*other variables*/
  int status,nobj;
  int i,j,k,l,t,ostep,ir,ir2,ik,**nmodl1,*lxvals,*lyvals,minmodl,lnum;
  double rvaluec,rvalue1,rvalue2,pnorm;

  int n_z,n_k,n_ki,lmaxk,n_l,n_r,ilx,ily;
  int lmin,lmax;
  double zmin,zmax,zpzdz,kmin,kmax,maxk,rmin,rmax,dz,dk,dr,dl,nobj2;
  double area,xdeg,n0,sigmae;

  time_t t1,t2,t3,t4,t5,t6,t7,t8;

  /*structures*/
  struct cosmo cosmo;
  struct cosmo cosmof;
  struct fid fid;

  /*fits cube output*/
  int yboxno,pixel;
  char *outputname,*idstring;

  /*input arrays*/  
  double *redshift;
  double m1,v1,m2,v2;
  double *zr,rvalue,kvalue,**numkvalues,denom,r0value,**numlvalues;
   
  /*data*/
  int    npz, nvalues=3;
  double pzmax,pznorm;
  char *catname,*scommand;
  double **printvalues;

  int lmod;
  double lx,ly,kr,angle,D1,D2,lnvalue;
  double sum1r,sum1i,sum2r,sum2i;
  double **sum1rE,**sum1iE;
  double **sum1rB,**sum1iB;
  double **summr,**summi;
  double **sum0r,**sum0i;

  double addnoise,pu,pv,theta,en1,en2,sigmar=0.32;
  double *numofz,**numzvalues,**numrvalues,z1,z2;

  FILE *ret;

  /*thread variables*/
  int NUM_THREADs=NUM_THREAD;
  pthread_t threads[NUM_THREADs];
  pthread_attr_t attr;
  void *tstatus;
  int thread_return,nt,max_num_thread,thread_returnd;
  pthread_mutex_init (&transformlock, NULL);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (argc==1) {
    printf(" arguments: ID (field number) \n");
    exit(2);
  }
  
  /************************************************************/

  rcs=0;  //rcs(1) format or cfht(0)

  /*add noise to the ellipticity values?*/
  addnoise=0; //1=yes

  /* create the catalogue name, uses format *passcat.cat*/
  catname=(char *)calloc(500,sizeof(char));
  strcpy(catname,DATADIR);
  strcat(catname,argv[1]);
  strcat(catname,"passcat.cat");
  printf(" %s\n",catname);
 
  /* set an arbitrary ID strings*/
  idstring=(char *)calloc(100,sizeof(char));
  strcpy(idstring,argv[1]);
  printf(" %s\n",idstring);
  printvalues=(double **)calloc(1,sizeof(double*));
  printvalues[0]=(double *)calloc(nvalues,sizeof(double));

  /************************************************************/
  
  /* read catalogue size*/
  nobj=-1; /*if not -1 this is hardwired otherwise uses all objects in the catalogue*/
  if (nobj==-1) {
    scommand=(char *)calloc(500,sizeof(char));
    strcat(scommand,"wc -l ");
    strcat(scommand,catname);
    strcat(scommand," > lines.dat");
    system(scommand);
    FILE *ffile;
    char sread[500];
    ffile=fopen("lines.dat","r");
    fscanf(ffile,"%s",sread);
    fclose(ffile);
    double ndouble;
    ndouble = strtod (sread,NULL);
    nobj=(int)ndouble;
    nobj--; //C convention of starting at zero
    printf(" nobj=%d\n",nobj);
  }

  /**************************************************************************/

  /*fill values from the header file*/

  n_l =fn_l; //specified number of ell modes
  
  /*define k range*/
  kmin=fkmin;
  kmax=fkmax;
  n_k =fn_k;
  
  /*define r number*/
  n_r=n_k; //same for computational speed, but can be different if needed
  
  /*data-hardwired into the p(z)*/
  zmin=0.025;//minimum set in the README for any p(z) 
  zmax=3.500;//zpzdz*npz; //set to a pre-determined maximum
  npz=70;    //hardwired to matach the CFHT p(z) sampling
  zpzdz=(zmax-zmin)/((float)70-1.0); //bin width in z
  n_z=npz;

  if (logk_flag==0) {
    //linear ell spacing:
    dk = ((double)kmax-(double)kmin)/(double)n_k;
  } else if (logk_flag==1) {
    //log spacing:
    dk = (log((double)kmax/(double)kmin))/(double)n_k;
  } else {
    printf("error: log or linear k spacing?\n");
    exit(2);
  }

  /**************************************************************************/

  /*fill in the structures with a fiducial cosmology*/
  cosmo.om=cosmo_om;
  cosmo.od=cosmo_od;
  cosmo.ob=cosmo_ob;
  cosmo.w0=cosmo_w0;
  cosmo.wa=cosmo_wa;
  cosmo.h0=cosmo_h0;
  cosmo.s8=cosmo_s8;
  cosmo.ns=cosmo_ns;
  cosmo.gm=cosmo_gm;
  cosmo.nm=cosmo_nm;
  cosmo.norm=1.;

  cosmof.om=cosmo_om;
  cosmof.od=cosmo_od;
  cosmof.ob=cosmo_ob;
  cosmof.w0=cosmo_w0;
  cosmof.wa=cosmo_wa;
  cosmof.h0=cosmo_h0;
  cosmof.s8=cosmo_s8;
  cosmof.ns=cosmo_ns;
  cosmof.gm=cosmo_gm;
  cosmof.nm=cosmo_nm;
  cosmof.norm=1.;

  /* fill in fiducial structure*/
  fid.plot =PLOTTING;
  fid.print=PRINTVALS;
  fid.kmin=kmin;
  fid.kmax=kmax;
  fid.n_k=n_k;
  fid.k_flag=logk_flag;
  fid.n_r=n_r;
  fid.area   =area;
  fid.n0     =n0;
  fid.zmin=zmin;
  fid.zmax=zmax;
  fid.n_z=n_z;
  fid.lmin=lmin;
  fid.lmax=lmax;
  fid.n_l=n_l;
  fid.dl=dl;
  
  /**************************************************************************/

  /*have to make the coeefficient data*/
    
  /*create the n(z) binned in redshift*/
  numofz=(double *)calloc(n_z,sizeof(double));

  /*read the catalogue*/
  printf(" %d objects in catalogue\n",nobj);
  //assign memory                                                                                                                
  e1=(double *)calloc(nobj,sizeof(double));
  e2=(double *)calloc(nobj,sizeof(double));
  e1local=(double *)calloc(nobj,sizeof(double));
  e2local=(double *)calloc(nobj,sizeof(double));
  thetax=(double *)calloc(nobj,sizeof(double));
  thetay=(double *)calloc(nobj,sizeof(double));
  RA=(double *)calloc(nobj,sizeof(double));
  DEC=(double *)calloc(nobj,sizeof(double));
  eweight=(double *)calloc(nobj,sizeof(double));
  redshift=(double *)calloc(nobj,sizeof(double));
  m=(double *)calloc(nobj,sizeof(double));
  c2=(double *)calloc(nobj,sizeof(double));
  
  //for rcs format
  if (rcs==1) {
    c1dp=(double *)calloc(nobj,sizeof(double));
    c1nb=(double *)calloc(nobj,sizeof(double));
    c2dp=(double *)calloc(nobj,sizeof(double));
    c2nb=(double *)calloc(nobj,sizeof(double));
  }
  
  //the p(z) for each galaxy                                                                                                     
  ppz=(double **)calloc(nobj,sizeof(double*));
  for (ik=0; ik<nobj; ik++) {
    ppz[ik]=(double *)calloc(npz+1,sizeof(double));}
  
  //the p(z) for each galaxy--interpolated onto fine r grid                                                                      
  ppzr=(double **)calloc(nobj,sizeof(double*));
  for (ik=0; ik<nobj; ik++) {
    ppzr[ik]=(double *)calloc(n_r,sizeof(double));}
  
  //abscia of the p(z)                                                                                                           
  zpz=(double *)calloc(npz+1,sizeof(double));
  for (ik=0; ik<npz;ik++) {
    zpz[ik]=zmin+(double)ik*zpzdz;
  }
 
  rg=(double *)calloc(nobj,sizeof(double)); //distance          

  //memory assigned, now read the data

  /**************************************************************************/

  //read arrays                                                                                                                  
  FILE *ufile;
  double minthetax=9999.,minthetay=9999.,maxthetax=-9999.,maxthetay=-9999.,ID,maxlikez,meanRA,minRA,meanDEC,minDEC,phix;
  char fread[500];

  ufile=fopen(catname,"r");
  if (ufile==NULL) {
    printf(" error opening file %s\n",ufile);
    exit(2);
  }
  printf(" file opened\n");
  
  //cuts should already be made at the catalogue level                                                                           
  ir=0;
  eweightnorm=0.;
  mweightnorm=0.;
  meanRA=0.;
  minRA=99999.;
  meanDEC=0.;
  minDEC=99999.;
  for (i=0; i<nobj; i++) { //loop rows                                                                                           
    
    fscanf(ufile,"%s",fread);
    ID = strtod (fread,NULL);
    fscanf(ufile,"%s",fread);
    RA[ir] = strtod (fread,NULL);
    fscanf(ufile,"%s",fread);
    DEC[ir] = strtod (fread,NULL);
    fscanf(ufile,"%s",fread);
    e1[ir] = strtod (fread,NULL);
    fscanf(ufile,"%s",fread);
    e2[ir] = strtod (fread,NULL);
    fscanf(ufile,"%s",fread); //weight                                                                                           
    eweight[ir] = strtod(fread,NULL);
    fscanf(ufile,"%s",fread); //weight                                                                                            
    maxlikez = strtod(fread,NULL);
    fscanf(ufile,"%s",fread); //m (1 and 2)
    m[ir] = strtod(fread,NULL);
    if (rcs==0) {
      fscanf(ufile,"%s",fread); //c2  
      c2[ir] = strtod(fread,NULL);
    } 
    if (rcs==1) {
      fscanf(ufile,"%s",fread); //c2  
      c1dp[ir] = strtod(fread,NULL);
      fscanf(ufile,"%s",fread); //c2  
      c1nb[ir] = strtod(fread,NULL);
      fscanf(ufile,"%s",fread); //c2  
      c2dp[ir] = strtod(fread,NULL);
      fscanf(ufile,"%s",fread); //c2  
      c2nb[ir] = strtod(fread,NULL);
    }
    
    //find peak of the posterior redshifts 
    redshift[ir]=0.;
    pzmax=0.;
    pznorm=0.;
    for (ik=0;ik<npz;ik++) {
      fscanf(ufile,"%s",fread);
      ppz[ir][ik] = strtod (fread,NULL);
      if (ppz[ir][ik] >= pzmax) {
	pzmax=ppz[ir][ik];
	redshift[ir]=zpz[ik];
      }
      pznorm+=ppz[ir][ik];
      }
    if (pznorm>1.1 || pznorm <0.9) {
      //renormalise (assume pdf is in the high z bin)
      for (ik=0;ik<npz;ik++) ppz[ir][ik]=ppz[ir][ik]/pznorm;
    }
    
    if (redshift[ir] <= zmax && redshift[ir] >= zmin && pznorm!=0.) {

      /*n(z)*/
      for (ik=0;ik<npz;ik++) numofz[ik]+=ppz[ir][ik]; //stack the p(z) to get an n(z)-prior for Noise 
      
      /*weights*/
      eweightnorm+=eweight[ir];
      mweightnorm+=(1.0+m[ir])*eweight[ir];
      
      /*angles*/
      meanRA+=RA[ir];
      meanDEC+=DEC[ir];
      if (RA[ir]<=minRA) minRA=RA[ir];
      if (DEC[ir]<=minDEC) minDEC=DEC[ir];
      
      ir++;
    }

  }
  fclose(ufile);
  nobj=ir;
  printf(" read catalogue with %d objects\n",nobj);
  for (ik=0;ik<npz;ik++) {
    nobj2+=numofz[ik];
    printf(" %d %f\n",ik,numofz[ik]);
  }
  printf(" n(z) has %f objects c.f. catalogue n(z) %d\n",nobj2,nobj);    
  
  /*mean of ellipticity and m weights*/
  eweightnorm=eweightnorm/(double)nobj;
  mweightnorm=mweightnorm/(double)nobj; //need to correct the weight norm
  printf(" eweightnorm=%e %e\n",eweightnorm,mweightnorm);

  /*need to convert from spherical to planar coordinates (spherical law of cosines)*/
  meanRA=meanRA/(double)nobj;
  meanDEC=meanDEC/(double)nobj;
  for (i=0;i<nobj;i++) {
    
    /*theta wrt mean center of the field*/
    thetay[i]=DEC[i]-minDEC;
    phix=(90.-DEC[i])*pi/180.;
    thetax[i]=cos(phix)*cos(phix)+sin(phix)*sin(phix)*cos((RA[i]-minRA)*pi/180.);
    thetax[i]=acos(thetax[i])*180./pi;
    
    if (thetax[i] < minthetax) minthetax=thetax[i];
    if (thetay[i] < minthetay) minthetay=thetay[i];
    if (thetax[i] > maxthetax) maxthetax=thetax[i];
    if (thetay[i] > maxthetay) maxthetay=thetay[i];
  }
  printf(" RA range: %f %f, DEC range %f %f\n",minthetax,maxthetax,minthetay,maxthetay);

  //c correction:                               
  if (rcs==0) {
    for (i=0;i<nobj;i++) e2[i]=e2[i]-c2[i];
  } 
  if (rcs==1) {
      for (i=0;i<nobj;i++) {
	e1[i]=e1[i]-c1dp[i]-c1nb[i];
	e2[i]=e2[i]-c2dp[i]-c2nb[i];
      }
  }
 
  /* not compute some simple stats*/
  area   = (maxthetax-minthetax)*(maxthetay-minthetay); //area
  xdeg   = sqrt(area);
  n0     = (double)nobj/(area*(60.*60.));  //galaxies per square arcminute                                                    
  fid.n0 = n0;
  printf(" n0=%f , area=%f sqrt(area)=%f\n",n0,area,xdeg);

  /*define l range*/
  lmin = safelmin;
  printf(" using fundamental mode of %d\n",lmin);
  dl   = 50.; //spacing in lx and ly, can changes this if required

  if (logl_flag==0) {
    //linear ell spacing:                                                                                                  
    lmax = lmin +(int)((double)n_l*dl);
  } else if (logl_flag==1) {
    printf(" log l spacing not yet implemented\n");
    exit(2);
    //log spacing:                                                          
    lmax=(int)((double)lmin*pow(10,(double)n_l*dl));
  } else {
    printf("error: log or linear l spacing?\n");
    exit(2);
  }
  
  fid.kmin=kmin;
  fid.kmax=kmax;
  fid.n_k=n_k;
  fid.k_flag=logk_flag;
  fid.n_r=n_r;
  fid.area   =area;
  fid.n0     =n0;
  fid.zmin=zmin;
  fid.zmax=zmax;
  fid.n_z=n_z;
  fid.lmin=lmin;
  fid.lmax=lmax;
  fid.n_l=n_l;
  fid.dl=dl;

  /*quick test on catalogue*/
  printf(" print first 100 objects for test\n");
  for (i=0;i<100;i++){
    printf(" %f %f %f %f %f\n",e1[i],e2[i],redshift[i],thetax[i],thetay[i]);
  }

  //mean and variance of ellipticities
  m1=m2=v1=v2=0.;
  for (i=0;i<nobj;i++) {
    m1+=e1[i]*eweight[i]/(mweightnorm*(double)nobj);
    m2+=e2[i]*eweight[i]/(mweightnorm*(double)nobj);
  }
  for (i=0;i<nobj;i++) v1+=(m1-e1[i])*(m1-e1[i])*eweight[i]/(mweightnorm*(double)nobj);
  for (i=0;i<nobj;i++) v2+=(m2-e2[i])*(m2-e2[i])*eweight[i]/(mweightnorm*(double)nobj);
  v1=sqrt(v1);
  v2=sqrt(v2);
  printf(" mean e1 = %f std e1 = %f\n",m1,v1);
  printf(" mean e2 = %f std e2 = %f\n",m2,v2);
  sigmae =sqrt(v1*v1+v2*v2); 
  fid.sigmae = sigmae;
  printf(" sum sigmae^2=%f\n",sigmae);

  if (addnoise!=1 && sigmae <=0.1) {
    printf(" WARNING : ellipticity noise is spuriously low\n");
    printf(" suggest changing addnoise=1 and sigmar=0.3\n");
    exit(2);                                                                                                                   
  }

  /*add noise if needed -- maybe needed for shear-only simulations*/
  if (addnoise==1) {
    /*add some shot noise to the ellipticity values*/
    srand48( (long)t1 );
    drand48();
    //using box-muller method for Gaussian random variable
    for (i=0;i<nobj;i++) {
    pu: pu=drand48();
      if (sqrt(-2.*sigmar*sigmar*log(pu)) > 1.0) goto pu;
      pv=drand48();
      en1=sqrt(-2.*sigmar*sigmar*log(pu))*cos(2.*pi*pv);
      en2=sqrt(-2.*sigmar*sigmar*log(pu))*sin(2.*pi*pv);
      e1[i]=e1[i]+en1;
      e2[i]=e2[i]+en2;
    }
    m1=m2=v1=v2=0;
    for (i=0;i<nobj;i++) {
      m1+=e1[i]/(float)nobj;
      m2+=e2[i]/(float)nobj;}
    for (i=0;i<nobj;i++) v1+=(m1-e1[i])*(m1-e1[i])/(float)nobj;
    for (i=0;i<nobj;i++) v2+=(m2-e2[i])*(m2-e2[i])/(float)nobj;
    v1=sqrt(v1);
    v2=sqrt(v2);
    printf(" mean e1 = %f variance e1 = %f\n",m1,v1);
    printf(" mean e2 = %f variance e2 = %f\n",m2,v2);
    fid.sigmae = sigmar; //per component
  }

  /**************************************************************************/
    
  /*calculate fiducial distances to the galaxies*/ //~~~~~not parallelised~~~~~
  ir=0;
  printf(" calculating fiducial r values\n");
  rmax=-99999.;
  rmin= 99999.;
  zmax=-99999.;
  zmin= 99999.;
  numrvalues=(double **)calloc(3,sizeof(double*));
  for (i=0;i<3;i++) numrvalues[i]=(double *)calloc(nobj,sizeof(double)); //for fits write
  for (i=0;i<nobj;i++) { 
    rg[i]=dist(cosmof,redshift[i]);
    if (rg[i]>rmax) rmax=rg[i];
    if (rg[i]<rmin) rmin=rg[i];
    if (redshift[i] > zmax) zmax=redshift[i];
    if (redshift[i] < zmin) zmin=redshift[i];
    numrvalues[0][i]=rg[i];
    numrvalues[1][i]=e1[i];
    numrvalues[2][i]=e2[i];
  }
  printf(" \n done rmin, rmax=%f, %f zmin, zmax=%f, %f n_r=%d\n",rmin,rmax,zmin,zmax,n_r);
  fid.rmin=rmin; //readjust fid values
  fid.rmax=rmax;
  fid.n_r=n_r;
  dr =(rmax-rmin)/((double)n_r); //only linear spacing implemented

  /*compute redshifts that match the r range in fid*/
  zr=(double *)calloc(n_r,sizeof(double));
  zr[0]=fid.rmin*cosmof.h0/3000.;
  for (ir=1;ir<fid.n_r;ir++) {
    rvalue = fid.rmin+(double)ir*dr;
    denom=drdz(cosmof,zr[ir-1]); 
    zr[ir]=0.;
    if ((cosmof.om+cosmof.od)==1.0) {
      if (denom != 0.) zr[ir] = zr[ir-1]+(dr/denom);
    } else if ((cosmof.om+cosmof.od) > 1.0) {
      r0value=R0(cosmof);
      if (denom != 0. && r0value !=0.) zr[ir] = zr[ir-1]+
	((sin(rvalue/r0value)-sin((rvalue-dr)/r0value))/
	 (denom*cos(rvalue/r0value)));
    } else if ((cosmof.om+cosmof.od) < 1.0) {
      r0value=R0(cosmof);
      if (denom != 0. && r0value !=0.) zr[ir] = zr[ir-1]+
	((sinh(rvalue/r0value)-sinh((rvalue-dr)/r0value))/
	 (denom*cosh(rvalue/r0value)));
    }
  } 

  /*interpolate coarse p(z) onto the fine r grid -- z space interpolation*/
  printf(" interpolating p(z) onto the fine r grid\n");
  for (i=0;i<nobj;i++) {
    pnorm=0.;
    for (ir=0;ir<fid.n_r-1;ir++) {

      rvalue = zr[ir];

      j=0;
      for (ir2=0;ir2<npz;ir2++) {
        rvaluec = zpz[ir2];
        if (rvalue>rvaluec && rvalue<=rvaluec+zpzdz) j=ir2;
      }
      if (j<0 || j >=npz) {
        printf(" interpolation error out of bounds\n");
        exit(2);
      }
      rvalue1 = zpz[j];
      rvalue2 = zpz[j+1];

      ppzr[i][ir]=interpolate_1D(ppz[i][j],ppz[i][j+1],rvalue1,rvalue2,rvalue);
      pnorm+=ppzr[i][ir];
    }
    if (pnorm!=0.) { //normalise the interpolated posterior                                                        
      for (ir=0;ir<fid.n_r-1;ir++) ppzr[i][ir]=ppzr[i][ir]/pnorm;
    }
  }
  printf(" done\n");
  
  /*create the n(z) binned in redshift*/
  numzvalues=(double **)calloc(2,sizeof(double*)); //for fits output                                
  numzvalues[0]=(double *)calloc(n_z,sizeof(double));
  numzvalues[1]=(double *)calloc(n_z,sizeof(double));
  
  j=n_z;
  for (i=0;i<n_z;i++) {
    z1=zmin+(double)i*(zmax-zmin)/(double)n_z;
    if (numofz[i]<=0.) {
	printf(" warning: n(z) at z=%f is zero\n",z1);
	exit(2);
    }
    numzvalues[0][i]=numofz[i];
    numzvalues[1][i]=z1;
  }
  fid.zmin=zmin; //re-adjust fid variables 
  fid.zmax=zmin+(double)(j-1)*(zmax-zmin)/(double)n_z;
  zmax=fid.zmax;
  n_z=j;
  fid.n_z=n_z;
  /**************************************************************************/
  
  /* loop over lx, ly space and calculate coefficients*/
  //find number of mod(l) modes 
  nmodl0=(int *)calloc(10000,sizeof(int)); 
  nmodl1=(int **)calloc(1,sizeof(int*)); 
  nmodl1[0]=(int *)calloc(10000,sizeof(int)); 
  for (ilx=0;ilx<10000;ilx++) {
    nmodl0[ilx]=0.;
    nmodl1[0][ilx]=0.;
  }
  k=0;
  for (ilx=0;ilx<fid.n_l/2+1;ilx++) {
    lx=(double)(ilx*dl);
    for (ily=-fid.n_l/2;ily<fid.n_l/2+1;ily++) {
      ly=(double)(ily*dl);
      lmod=(int)(sqrt(lx*lx+ly*ly));   
      nmodl0[k]=lmod;
      k=k+1;
    }
  }  
  for (ily=0;ily<k;ily++) {
    minmodl=99999;
    for (ilx=0;ilx<k;ilx++) {
      if (nmodl0[ilx]<minmodl && nmodl0[ilx]>nmodl1[0][ily-1]) {
	minmodl=nmodl0[ilx];
      }
    }
    nmodl1[0][ily]=minmodl;
    if (nmodl1[0][ily]>99998) {
      nummodl=ily;
      break;
    }
  }

  numkvalues=(double **)calloc(1,sizeof(double*));
  numkvalues[0]=(double *)calloc(fid.n_k,sizeof(double));
  numlvalues=(double **)calloc(1,sizeof(double*));
  numlvalues[0]=(double *)calloc(nummodl,sizeof(double));

  for (ik=0;ik<fid.n_k;ik++) {  
    if (logk_flag==0) {
      //linear k spacing:
      kvalue = kmin+(double)ik*dk;
    } else if (logk_flag==1) {
      //log spacing:
      kvalue = kmin*pow(10,(double)ik*dk);
    } else {
      printf("error: log or linear k spacing?\n");
      exit(2);
    }
    numkvalues[0][ik]=kvalue;
  }
  ret=fopen("lmodes.dat","w");
  printf(" %d unique l-modes\n",nummodl);
  for (ilx=0;ilx<nummodl;ilx++) {
    numlvalues[0][ilx]=(double)nmodl1[0][ilx];
    if (numlvalues[0][ilx]>=safelmin && numlvalues[0][ilx]<=safelmax) { 
      fprintf(ret," %f\n",numlvalues[0][ilx]);
    } 
  }
  fclose(ret);
  printf(" printed lmodes to lmodes.dat\n");

  /********************************************************************/

  /*noise correction matrix*/
  double **xSigma,*xSigmaCount,**ySigma,*ySigmaCount,ddd,**xSigma1,*xSigma2,**ySigma1,*ySigma2;
  xSigma=(double **)calloc(nummodl,sizeof(double*));
  xSigma1=(double **)calloc(1,sizeof(double *));
  xSigma1[0]=(double *)calloc(nummodl,sizeof(double));
  for (ilx=0;ilx<nummodl;ilx++) xSigma[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  xSigmaCount=(double *)calloc(nummodl,sizeof(double));
  
  ySigma=(double **)calloc(nummodl,sizeof(double*));
  ySigma1=(double **)calloc(1,sizeof(double *));
  ySigma1[0]=(double *)calloc(nummodl,sizeof(double));
  for (ilx=0;ilx<nummodl;ilx++) ySigma[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  ySigmaCount=(double *)calloc(nummodl,sizeof(double));

  //initialise
  for (j=0;j<nummodl;j++) {
    xSigmaCount[j]=0.;
    ySigmaCount[j]=0.;
    xSigma1[0][j]=0.;
    ySigma1[0][j]=0.;
    for (i=0;i<fid.n_k;i++) {
      xSigma[j][i]=0.;
      ySigma[j][i]=0.;
    }
  }

  //find the l-multipliers for 1 quandrant
  for (ilx=0;ilx<fid.n_l/2+1;ilx++) {
    lx=(double)(ilx*dl);
    for (ily=-fid.n_l/2;ily<fid.n_l/2+1;ily++) {
      ly=(double)(ily*dl);

      D1=((ly*ly)-(lx*lx))*0.5;
      D2=-lx*ly;
      ddd=1./(D1*D1+D2*D2);
      lmod=(int)(sqrt(lx*lx+ly*ly));

      if (lmod>=safelmin && lmod <=safelmax) {

	for (j=0;j<nummodl;j++) {
	  if ((int)numlvalues[0][j]==lmod) {
	    xSigma[j][0]+=(D1*D1*D1*D1+D2*D2*D1*D1)/((D1*D1+D2*D2)*(D1*D1+D2*D2)); //noise bit 
	    if (ddd*ddd*(D1*D1*D1*D1+D2*D2*D1*D1)!=0.) xSigmaCount[j]+=1.;
	    ySigma[j][0]+=(D2*D2*D2*D2+D2*D2*D1*D1)/((D1*D1+D2*D2)*(D1*D1+D2*D2)); 
	    if (ddd*ddd*(D2*D2*D2*D2+D2*D2*D1*D1)!=0.) ySigmaCount[j]+=1.;
	    xSigma1[0][j]+=D1*D1; //signal bit
	    ySigma1[0][j]+=D2*D2;
	  }
	}

      }

    }
  }
  //avoid infinities, due to numerics 
  for (j=0;j<nummodl;j++) {
    if (xSigmaCount[j]==0.) xSigmaCount[j]=1.;
    if (ySigmaCount[j]==0.) ySigmaCount[j]=1.;
  }

  //average and print
  for (j=0;j<nummodl;j++) {
    xSigma[j][0]=xSigma[j][0]/xSigmaCount[j];
    for (i=1;i<fid.n_k;i++) xSigma[j][i]=xSigma[j][0];
    ySigma[j][0]=ySigma[j][0]/ySigmaCount[j];
    for (i=1;i<fid.n_k;i++) ySigma[j][i]=ySigma[j][0];
    
    xSigma1[0][j]=xSigma1[0][j]/xSigmaCount[j];
    ySigma1[0][j]=ySigma1[0][j]/ySigmaCount[j];
    if (j<100) printf(" %f %f %f %f\n",xSigma1[0][j],ySigma1[0][j],xSigma[0][j],ySigma[0][j]);
  }

  /********************************************************************/
  /********************************************************************/
  /* now read to make the coefficients*/


  //assign memory for the output coefficient arrays
  sum1rE=(double **)calloc(nummodl,sizeof(double*));
  sum1iE=(double **)calloc(nummodl,sizeof(double*));
  for (ilx=0;ilx<nummodl;ilx++) {
    sum1rE[ilx]=(double *)calloc(fid.n_k,sizeof(double));
    sum1iE[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  }
  sum1rB=(double **)calloc(nummodl,sizeof(double*));
  sum1iB=(double **)calloc(nummodl,sizeof(double*));
  for (ilx=0;ilx<nummodl;ilx++) {
    sum1rB[ilx]=(double *)calloc(fid.n_k,sizeof(double));
    sum1iB[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  }
  summr=(double **)calloc(nummodl,sizeof(double*));
  summi=(double **)calloc(nummodl,sizeof(double*));
  for (ilx=0;ilx<nummodl;ilx++) {
    summr[ilx]=(double *)calloc(fid.n_k,sizeof(double));
    summi[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  }
  sum0r=(double **)calloc(nummodl,sizeof(double*));
  sum0i=(double **)calloc(nummodl,sizeof(double*));
  for (ilx=0;ilx<nummodl;ilx++) {
    sum0r[ilx]=(double *)calloc(fid.n_k,sizeof(double));
    sum0i[ilx]=(double *)calloc(fid.n_k,sizeof(double));
  }

  
  if (pzs==1) {
    Plkr=(double ***)calloc(nummodl,sizeof(double**));
    for (ilx=0;ilx<nummodl;ilx++) {
      Plkr[ilx]=(double **)calloc(fid.n_k,sizeof(double*));
      for (ik=0;ik<fid.n_k;ik++) {
	Plkr[ilx][ik]=(double *)calloc(fid.n_r,sizeof(double));
      }
    }
  }
  /********************************************************************/

  /*assign transform thread variables*/
  threadindextransform *threadargtransform;
  threadargtransform = (threadindextransform *)calloc(NUM_THREADs,sizeof(threadindextransform)); 

  /* ell and k values in thread-protected arrays*/
  thread_numkvalues=(double ***)calloc(NUM_THREADs,sizeof(double**));
  for (nt=0;nt<NUM_THREADs;nt++) {
    thread_numkvalues[nt]=(double **)calloc(1,sizeof(double*));
    thread_numkvalues[nt][0]=(double *)calloc(fid.n_k,sizeof(double));

    for (ik=0;ik<fid.n_k;ik++) {  
      thread_numkvalues[nt][0][ik]=numkvalues[0][ik];
    }
  }
  thread_nmodl1=(int ***)calloc(NUM_THREADs,sizeof(int**)); 
  thread_nmodl0=(int **)calloc(NUM_THREADs,sizeof(int*));
  thread_Xsum=(double **)calloc(NUM_THREADs,sizeof(double*));
  for (nt=0;nt<NUM_THREADs;nt++) {
    thread_nmodl1[nt]=(int **)calloc(1,sizeof(int*)); 
    
    thread_nmodl1[nt][0]=(int *)calloc(10000,sizeof(int)); 
    thread_nmodl0[nt]=(int *)calloc(10000,sizeof(int)); 
    thread_Xsum[nt]=(double *)calloc(10000,sizeof(double));

    for (i=0;i<10000;i++) {
      thread_nmodl1[nt][0][i]=nmodl1[0][i];
    }
  }

  thread_sum1rE=(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_sum1iE=(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_sum1rB=(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_sum1iB=(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_summr =(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_summi =(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_sum0r =(double ***)calloc(NUM_THREADs,sizeof(double**));
  thread_sum0i =(double ***)calloc(NUM_THREADs,sizeof(double**));
  for (nt=0;nt<NUM_THREADs;nt++) {
    thread_sum1rE[nt]=(double **)calloc(nummodl,sizeof(double*));
    thread_sum1iE[nt]=(double **)calloc(nummodl,sizeof(double*));
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_sum1rE[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
      thread_sum1iE[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
    }
    thread_sum1rB[nt]=(double **)calloc(nummodl,sizeof(double*));
    thread_sum1iB[nt]=(double **)calloc(nummodl,sizeof(double*));
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_sum1rB[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
      thread_sum1iB[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
    }
    thread_summr[nt]=(double **)calloc(nummodl,sizeof(double*));
    thread_summi[nt]=(double **)calloc(nummodl,sizeof(double*));
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_summr[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
      thread_summi[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
    }
    thread_sum0r[nt]=(double **)calloc(nummodl,sizeof(double*));
    thread_sum0i[nt]=(double **)calloc(nummodl,sizeof(double*));
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_sum0r[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
      thread_sum0i[nt][ilx]=(double *)calloc(fid.n_k,sizeof(double));
    }
  }

  if (pzs==1) {
  //big matrix... careful not to increase n_k or n_r too much                                                      
  thread_Plkr=(double ****)calloc(NUM_THREADs,sizeof(double***));
  for (nt=0;nt<NUM_THREADs;nt++) {
    thread_Plkr[nt]=(double ***)calloc(nummodl,sizeof(double**));
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_Plkr[nt][ilx]=(double **)calloc(fid.n_k,sizeof(double*));
      for (ik=0;ik<fid.n_k;ik++) {
        thread_Plkr[nt][ilx][ik]=(double *)calloc(fid.n_r,sizeof(double));
      }
    }
  }
  }
    
  //Bessel function array 
  be=(double ***)calloc(nobj,sizeof(double**));
  for (i=0;i<NUM_THREADs;i++) {
    be[i]=(double **)calloc(nummodl,sizeof(double*));
    for (ilx=0;ilx<nummodl;ilx++) {
      be[i][ilx]=(double *)calloc(fid.n_k,sizeof(double));
    }
  }

  /**************************************************************************/

  //initialise 
  for (ilx=0;ilx<nummodl;ilx++) nmodl0[ilx]=0.;
  for (nt=0;nt<NUM_THREADs;nt++) {
    for (ilx=0;ilx<nummodl;ilx++) {
      thread_nmodl0[nt][ilx]=0.;
      thread_Xsum[nt][ilx]=0.; //reset
    }
  }
  for (nt=0;nt<NUM_THREADs;nt++) {
    for (ilx=0;ilx<nummodl;ilx++) {
      for (ik=0;ik<fid.n_k;ik++) {
	thread_sum1rE[nt][ilx][ik]=0.;
	thread_sum1iE[nt][ilx][ik]=0.;
	thread_sum1rB[nt][ilx][ik]=0.;
	thread_sum1iB[nt][ilx][ik]=0.;
	thread_summr[nt][ilx][ik]=0.;
        thread_summi[nt][ilx][ik]=0.;
	thread_sum0r[nt][ilx][ik]=0.;
        thread_sum0i[nt][ilx][ik]=0.;

	if (nt==0) {
	  sum1rE[ilx][ik]=0.;
	  sum1iE[ilx][ik]=0.;
	  
	  sum1rB[ilx][ik]=0.;
	  sum1iB[ilx][ik]=0.;

	  summr[ilx][ik]=0.;
          summi[ilx][ik]=0.;

	  sum0r[ilx][ik]=0.;
          sum0i[ilx][ik]=0.;
	}

	if (pzs==1) {
	  for (ir=0;ir<fid.n_k;ir++) {
	    thread_Plkr[nt][ilx][ik][ir]=0.;}
	}
   
      }
    }
  }

  //loop over the transform function in a POSIX-threaded calculation
  j=0;
  ostep=(int)((float)nobj/(float)MAX_THREAD);
  printf(" object steps of %d\n",ostep);
  while (j<nobj) {
    for (nt=0; nt<NUM_THREADs; nt++)
      {

      	if (j-ostep < nobj) {
	 
	threadargtransform[nt].fid  = fid;
	threadargtransform[nt].obji = j;
	if (j+ostep<nobj) {
	  threadargtransform[nt].objf = j+ostep;
	} else {
	  threadargtransform[nt].objf = nobj;
	}
	threadargtransform[nt].nt   = nt;
	printf(" analysing objects %d to %d\n",threadargtransform[nt].obji,threadargtransform[nt].objf);
	  	
	thread_return = pthread_create(&threads[nt], 
				       &attr, 
				       stransform_thread, 
				       (void *)&threadargtransform[nt]);

	if (thread_return)
	  {
	    printf(" error %d from pthread_create %d \n",thread_return,nt);
	    printf(" %s \n",strerror(thread_return));
	    exit(2);
	  }
	j=j+ostep;
     	max_num_thread=nt;
	}
 
      }
    /* wait for all threads to rejoin before going on*/
    for (nt=0; nt<max_num_thread; nt++)
      {
	thread_return = pthread_join(threads[nt], &tstatus);
	
	if (thread_return !=0 && thread_return != ESRCH)
	  {
	    printf(" error %d from pthread_join %d \n",thread_return,nt);
	    printf(" %s \n",strerror(thread_return));
	    exit(2);
	  }
      }
  }
  printf(" done\n");

  /*fill used arrays with thread argument versions*/
  double mMd,mRd,mId,ERd,EId,BRd,BId;
  for (ilx=0;ilx<nummodl;ilx++) {

    for (nt=0;nt<NUM_THREADs;nt++) nmodl0[ilx]+=thread_nmodl0[nt][ilx];

    for (ik=0;ik<fid.n_k;ik++) {

      for (nt=0;nt<NUM_THREADs;nt++) {
	sum1rE[ilx][ik]+=thread_sum1rE[nt][ilx][ik];
	sum1iE[ilx][ik]+=thread_sum1iE[nt][ilx][ik];
	
	sum1rB[ilx][ik]+=thread_sum1rB[nt][ilx][ik];
	sum1iB[ilx][ik]+=thread_sum1iB[nt][ilx][ik];

	summr[ilx][ik]+=thread_summr[nt][ilx][ik];
	summi[ilx][ik]+=thread_summi[nt][ilx][ik];

	sum0r[ilx][ik]+=thread_sum0r[nt][ilx][ik];
        sum0i[ilx][ik]+=thread_sum0i[nt][ilx][ik];

	if (pzs==1) {
	  for (ir=0;ir<fid.n_r;ir++) {
	    Plkr[ilx][ik][ir]+=thread_Plkr[nt][ilx][ik][ir];}
	}
      }
      mMd=(summr[ilx][ik]*summr[ilx][ik])+(summi[ilx][ik]*summi[ilx][ik]);
      mRd=(sum0r[ilx][ik]*summr[ilx][ik])+(sum0i[ilx][ik]*summi[ilx][ik]);
      mId=(sum0i[ilx][ik]*summr[ilx][ik])-(sum0r[ilx][ik]*summi[ilx][ik]);
     if (mMd != 0.) {
	ERd=(mRd/mMd)*sum1rE[ilx][ik]-(mId/mMd)*sum1rB[ilx][ik];
	EId=(mRd/mMd)*sum1iE[ilx][ik]-(mId/mMd)*sum1iB[ilx][ik];
	BRd=(mRd/mMd)*sum1rB[ilx][ik]+(mId/mMd)*sum1rE[ilx][ik];
	BId=(mRd/mMd)*sum1iB[ilx][ik]+(mId/mMd)*sum1iE[ilx][ik];
	if (isnan((mRd*mRd+mId*mId)/(mMd*mMd))!=1) { 
	  xSigma[ilx][ik]=xSigma[ilx][ik]*(mRd*mRd+mId*mId)/(mMd*mMd);
	  ySigma[ilx][ik]=ySigma[ilx][ik]*(mRd*mRd+mId*mId)/(mMd*mMd);
	}
     } else {
       ERd=0.;
       EId=0.;
       BRd=0.;
       BId=0.;
     }
     sum1rE[ilx][ik]=ERd;
     sum1iE[ilx][ik]=EId;
     sum1rE[ilx][ik]=BRd;
     sum1iB[ilx][ik]=BId;
    }

    thread_Xsum[0][ilx]=thread_Xsum[0][ilx]*thread_Xsum[0][ilx]; //|X|^2 per l mode
  }

  for (ilx=0;ilx<nummodl;ilx++) {
    for (ik=0;ik<fid.n_k;ik++) {
      sum1rE[ilx][ik]=sum1rE[ilx][ik]/xSigmaCount[ilx];
      sum1iE[ilx][ik]=sum1iE[ilx][ik]/ySigmaCount[ilx]; 
      sum1rB[ilx][ik]=sum1rB[ilx][ik]/xSigmaCount[ilx];
      sum1iB[ilx][ik]=sum1iB[ilx][ik]/ySigmaCount[ilx];
    }
  }
  printf(" done transforms\n");

  /**************************************************************************/

  //now output the arrays as FITS tables 

  outputname = (char *)calloc(300,sizeof(char));

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_1rE.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,sum1rE,nummodl,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_1iE.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,sum1iE,nummodl,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_1rB.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,sum1rB,nummodl,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_1iB.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,sum1iB,nummodl,yboxno,fid.n_k,1);
  
  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_ell.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,numlvalues,1,yboxno,nummodl,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_kll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable
  write2Dmatrix_fits(outputname,numkvalues,1,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_zll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                       
  write2Dmatrix_fits(outputname,numzvalues,2,yboxno,fid.n_z,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_rll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                         
  write2Dmatrix_fits(outputname,numrvalues,3,yboxno,nobj,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_xll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                         
  write2Dmatrix_fits(outputname,thread_Xsum,1,yboxno,nummodl,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_sll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                                                            
  write2Dmatrix_fits(outputname,xSigma,nummodl,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_tll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                                                             
  write2Dmatrix_fits(outputname,ySigma,nummodl,yboxno,fid.n_k,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_drll.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                                                             
  write2Dmatrix_fits(outputname,xSigma1,1,yboxno,nummodl,1);

  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_dill.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable                                                                                                    
  write2Dmatrix_fits(outputname,ySigma1,1,yboxno,nummodl,1);

  if (pzs==1) {
    strcpy(outputname,PRODUCTDIR);
    strcat(outputname,"3dtransform_plkr.");
    strcat(outputname,idstring);
    strcat(outputname,".fits");
    remove(outputname);
    yboxno=1; //dummy variable                                                                             
    write3Dmatrix_fits(outputname,Plkr,nummodl,yboxno,fid.n_k,fid.n_k);
  }

  //used to store some doubles needed
  strcpy(outputname,PRODUCTDIR);
  strcat(outputname,"3dtransform_values.");
  strcat(outputname,idstring);
  strcat(outputname,".fits");
  remove(outputname);
  yboxno=1; //dummy variable 
  printvalues[0][0]=sigmae;
  printvalues[0][1]=area;
  write2Dmatrix_fits(outputname,printvalues,1,yboxno,nvalues,1);

  /**************************************************************************/

  return 0;
}

/*****************************************************/
//the thread function where the coefficients are calculated 
void *stransform_thread(void *threadargs) {

  int ilx,ily,tn,k,ik,ir,lmod,objn,tobji,tobjf,lnum;
  double D1,D2,ddd,sum1r,sum1i,sum2r,sum2i,summr,summi,sum0r,sum0i,angle,lx,ly;
  double gammaR,gammaI,kvalue,kr,averageDsig1,averageDsig2,nDsig;
  double mammaR,mammaI,mamm0R,mamm0I; 

  struct fid fidp;
  struct cosmo cosmop;
  
  threadindextransform *this_thread;
  this_thread = (threadindextransform *)threadargs;

  fidp   = this_thread->fid;
  tobji  = this_thread->obji;
  tobjf  = this_thread->objf;
  tn     = this_thread->nt;

  for (objn=tobji;objn<tobjf;objn++) {
    for (ilx=0;ilx<nummodl;ilx++) {
      lmod=thread_nmodl1[tn][0][ilx];
      
      if (lmod >=safelmin && lmod <=safelmax) {
	for (ik=0;ik<fidp.n_k;ik++) {
	  //calculate kvalues and bessel functions 
	  kvalue=thread_numkvalues[tn][0][ik];
	  kr = kvalue*rg[objn];
	  be[tn][ilx][ik] = sjl(lmod,kr);
	}
      }
    }

    if (pzs==1) {
      //printf("create the P_\ell(k,r) matrix\n");/*(Kitching,Heavens,Miller equation 5)*/                           
      for (ilx=0;ilx<nummodl;ilx++) {
	lmod=thread_nmodl1[tn][0][ilx];      
	if (lmod >=safelmin && lmod <=safelmax) {
	  for (ik=0;ik<fidp.n_k;ik++) {
	    for (ir=0;ir<fidp.n_r;ir++) {
	      thread_Plkr[tn][ilx][ik][ir]+=be[tn][ilx][ik]*ppzr[objn][ir];
	    }
	  }
	}
      }
    }

    /*calculating transform*/
    averageDsig1=averageDsig2=0.;
    for (ilx=0;ilx<fidp.n_l/2+1;ilx++) {
      
      lx=(double)(ilx*fidp.dl);
      
      for (ily=-fidp.n_l/2;ily<fidp.n_l/2+1;ily++) {
	
	ly=(double)(ily*fidp.dl);
	
	lmod=(int)(sqrt(lx*lx+ly*ly));

	if (lmod>=safelmin && lmod <=safelmax) {
	  
	  D1=((ly*ly)-(lx*lx))*0.5;
	  D2=-lx*ly;
	  
	  averageDsig1+=D1*D1/(D1*D1+D2*D2);
          averageDsig2+=D2*D2/(D1*D1+D2*D2);
	  nDsig+=1.;

	  ddd=1./(D1*D1+D2*D2);

	  for (k=0;k<nummodl;k++) {
	    if (lmod==thread_nmodl1[tn][0][k]) lnum=k;
	  }
	  if (objn==0 && tn==0) { //only need to count the modes once so do at obj=0 for first thread  
	    thread_nmodl0[tn][lnum]++;
	    thread_Xsum[tn][lnum]+=2.*sqrt(D1*D1+D2*D2); //sum of mod(X_l)'s
	  }
	  
	  for (ik=0;ik<fidp.n_k;ik++) {
	    
	    sum1r = 0.0;
	    sum1i = 0.0;
	    sum2r = 0.0;
	    sum2i = 0.0;
	    summr = 0.0;
	    summi = 0.0;
	    sum0r = 0.0;
	    sum0i = 0.0;

	    angle = (lx*thetax[objn])+(ly*thetay[objn]);
	    angle = angle*pi/180.; //radians
 
	    sum1r = sum1r + e1[objn]*be[tn][lnum][ik]*cos(angle)*sqrt2overpi*eweight[objn]/eweightnorm;
	    sum1i = sum1i + e1[objn]*be[tn][lnum][ik]*sin(angle)*sqrt2overpi*eweight[objn]/eweightnorm;
	    sum2r = sum2r + e2[objn]*be[tn][lnum][ik]*cos(angle)*sqrt2overpi*eweight[objn]/eweightnorm;
	    sum2i = sum2i + e2[objn]*be[tn][lnum][ik]*sin(angle)*sqrt2overpi*eweight[objn]/eweightnorm;

	    summr = summr + (1.+m[objn])*be[tn][lnum][ik]*cos(angle)*sqrt2overpi*eweight[objn];
	    summi = summi + (1.+m[objn])*be[tn][lnum][ik]*sin(angle)*sqrt2overpi*eweight[objn];
	    sum0r = sum0r + be[tn][lnum][ik]*cos(angle)*sqrt2overpi*eweight[objn];
            sum0i = sum0i + be[tn][lnum][ik]*sin(angle)*sqrt2overpi*eweight[objn];

	    //real and imag parts of the sum
	    gammaR=sum1r+sum2i; //real
	    gammaI=sum2r-sum1i; //imag
	    
	    if (isnan(gammaR)==1) {
	      gammaR=0.;
	    }
	    if(isnan(gammaI)==1) {
	      gammaI=0.;
	    }
	    
	    //rotate out E and B modes 
	    thread_sum1rE[tn][lnum][ik]+=D1*ddd*(D1*gammaR+D2*gammaI); 
	    thread_sum1iE[tn][lnum][ik]+=D2*ddd*(D1*gammaR+D2*gammaI); 
	    
	    thread_sum1rB[tn][lnum][ik]+=D1*ddd*(D1*gammaI-D2*gammaR); 
	    thread_sum1iB[tn][lnum][ik]+=D2*ddd*(D1*gammaI-D2*gammaR); 
	    
	    //shear measurement bias m-correction terms 
	    mammaR=summr+summi; 
	    mammaI=summr-summi;
	    mamm0R=sum0r+sum0i;
            mamm0I=sum0r-sum0i;

	    thread_summr[tn][lnum][ik]+=D1*ddd*(D1*mammaR+D2*mammaI);
            thread_summi[tn][lnum][ik]+=D2*ddd*(D1*mammaR+D2*mammaI);
            thread_sum0r[tn][lnum][ik]+=D1*ddd*(D1*mamm0R+D2*mamm0I);
            thread_sum0i[tn][lnum][ik]+=D2*ddd*(D1*mamm0R+D2*mamm0I);
	    
	  }
	}
      }
    }

    averageDsig1=averageDsig1/nDsig;
    averageDsig2=averageDsig2/nDsig;
  }

  pthread_exit((void*) 0);
}
/*****************************************************/
