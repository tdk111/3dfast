#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"

int write2Dmatrix_fits(char *outputname, double **f, int numl, int yboxnum, int pwidth, int pheight)
{
  fitsfile *afptr;
  int status = 0;
  int atype, anaxis, check = 1;
  long anaxes[3], fpixel[3]={1,1,1};
  int i, n, bitpix, row, psfnum, psfnumx, psfnumy;
  float *psf, *psf3d;
  double sum;
  int psfsize, psfsize3d, dwidth, dfactor, maxboxno;
  int ix, iy, x, y, j, ij, npixels;

  /*  write out matrix data as a fits file */

  if (pwidth<=0 || pheight <=0 || numl<=0 || yboxnum<=0)
    {
      fprintf (stderr," error in subimage dimensions or number \n");
      return 1;
    }

  bitpix = -32;
  anaxis = 3;
  anaxes[0] = pwidth;
  anaxes[1] = pheight;
  maxboxno  = numl*yboxnum;
  anaxes[2] = maxboxno;

  psfsize = anaxes[0]*anaxes[1];
  psf = (float*)calloc(psfsize,sizeof(float));
  if (psf==NULL)
    {
      fprintf (stderr," error allocating memory for psf in writepsf_fits \n");
      fprintf (stderr," memory requested = %d floats \n",psfsize);
      exit (2);
    }


  psfsize3d = maxboxno*psfsize;
  psf3d = (float*)calloc(psfsize3d,sizeof(float));
  if (psf3d==NULL)
    {
      fprintf (stderr," error allocating memory for 3D psf in writepsf_fits \n");
      fprintf (stderr," memory requested = %d floats \n",psfsize3d);
      exit (2);
    }

  fits_create_file(&afptr, outputname, &status); 
  fits_create_img(afptr, bitpix, anaxis, anaxes, &status);

  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }

  for (psfnum=0; psfnum<maxboxno; psfnum++)
    {
      
      npixels= (int)(pwidth*pheight);

      /* write matrix into 3D array */
      for (i=0; i<npixels; i++)
	{
	  j = psfnum*npixels + i;
	  //printf("%d\n",j);
	  psf3d[j] = (float)f[psfnum][i];
	}
    }

  /* write all matrix data into image array */
  if (fits_write_pix(afptr, TFLOAT, fpixel, psfsize3d, psf3d, &status) )
    {
      printf(" error reading pixel data \n");
      exit (2);
    }

  /* close the fits file */

  fits_close_file(afptr, &status);

  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }


  return (0);

}
