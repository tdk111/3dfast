#include <stdio.h>
#include <unistd.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>

#include "shear3d.h"

/*****************************************************/
double interpolate_1D(double y1,double y2,double x1,double x2, double xvalue) {

  //1D interpolation

  double grad, inter, output;
  
  grad=(y2-y1)/(x2-x1);
  inter=y1-(grad*x1);
  output=(grad*xvalue)+inter;

  return output;
}
/*****************************************************/
double hubble(struct cosmo cosmo, double z) {

  //H(z)

  double a,output;

  a=1./(1+z);

  output=cosmo.h0*100.*sqrt(cosmo.om*pow(a,-3.)+
			    cosmo.od*pow(a,-3.*(1+cosmo.w0+cosmo.wa))*
			    exp(-3*cosmo.wa*(1.-a))+
			    (1.-cosmo.om-cosmo.od)*pow(a,-2.));
    
  return output;

}
/*****************************************************/
double dist(struct cosmo cosmo, double z) {

  //D(z) : by linear interpolated integration

  int i,iint=100;
  double zmin=0.,dz,output=0.,zvalue,h;
  
  dz=(z-zmin)/((double)iint);

  for (i=0;i<iint;i++) {
    zvalue=zmin+(double)i*dz;
    h=hubble(cosmo,zvalue);
    if (h!=0.) output+=dz/h;
  }
  output=output*light;
  
  return output;
}

/*****************************************************/
//==============================================================================
// Recursive definition of determinate using expansion by minors.
//
// Notes: 1) arguments:
//             a (double **) pointer to a pointer of an arbitrary square matrix
//             n (int) dimension of the square matrix
//
//        2) Determinant is a recursive function, calling itself repeatedly
//           each time with a sub-matrix of the original till a terminal
//           2X2 matrix is achieved and a simple determinat can be computed.
//           As the recursion works backwards, cumulative determinants are
//           found till untimately, the final determinate is returned to the
//           initial function caller.
//
//        3) m is a matrix (4X4 in example)  and m13 is a minor of it.
//           A minor of m is a 3X3 in which a row and column of values
//           had been excluded.   Another minor of the submartix is also
//           possible etc.
//             m  a b c d   m13 . . . .
//                e f g h       e f . h     row 1 column 3 is elminated
//                i j k l       i j . l     creating a 3 X 3 sub martix
//                m n o p       m n . p
//
//        4) the following function finds the determinant of a matrix
//           by recursively minor-ing a row and column, each time reducing
//           the sub-matrix by one row/column.  When a 2X2 matrix is
//           obtained, the determinat is a simple calculation and the
//           process of unstacking previous recursive calls begins.
//
//                m n
//                o p  determinant = m*p - n*o
//
//        5) this function uses dynamic memory allocation on each call to
//           build a m X m matrix  this requires **  and * pointer variables
//           First memory allocation is ** and gets space for a list of other
//           pointers filled in by the second call to malloc.
//
//        6) C++ implements two dimensional arrays as an array of arrays
//           thus two dynamic malloc's are needed and have corresponsing
//           free() calles.
//
//        7) the final determinant value is the sum of sub determinants
//
//==============================================================================


double determ(double **a,int n)
{
    int i,j,j1,j2 ;                    // general loop and matrix subscripts
    double det = 0 ;                   // init determinant
    double **m;                         // pointer to pointers to implement 2d
                                        // square array

    if (n < 1)    {   }                // error condition, should never get here

    else if (n == 1) {                 // should not get here
        det = a[0][0] ;
        }

    else if (n == 2)  {                // basic 2X2 sub-matrix determinate
                                       // definition. When n==2, this ends the
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1] ;// the recursion series
        }


                                       // recursion continues, solve next sub-matrix
    else {                             // solve the next minor by building a
                                       // sub matrix
        det = 0 ;                      // initialize determinant of sub-matrix

                                           // for each column in sub-matrix
        for (j1 = 0 ; j1 < n ; j1++) {
                                           // get space for the pointer list
	  m = (double **) calloc((n-1),sizeof(double *)) ;

            for (i = 0 ; i < n-1 ; i++)
	      m[i] = (double *) calloc((n-1),sizeof(double)) ;
                       //     i[0][1][2][3]  first malloc
                       //  m -> +  +  +  +   space for 4 pointers
                       //       |  |  |  |          j  second malloc
                       //       |  |  |  +-> _ _ _ [0] pointers to
                       //       |  |  +----> _ _ _ [1] and memory for
                       //       |  +-------> _ a _ [2] 4 doubles
                       //       +----------> _ _ _ [3]
                       //
                       //                   a[1][2]
                      // build sub-matrix with minor elements excluded
            for (i = 1 ; i < n ; i++) {
                j2 = 0 ;               // start at first sum-matrix column position
                                       // loop to copy source matrix less one column
                for (j = 0 ; j < n ; j++) {
                    if (j == j1) continue ; // don't copy the minor column element

                    m[i-1][j2] = a[i][j] ;  // copy source element into new sub-matrix
                                            // i-1 because new sub-matrix is one row
                                            // (and column) smaller with excluded minors
                    j2++ ;                  // move to next sub-matrix column position
                    }
                }

            det += pow(-1.0,1.0 + j1 + 1.0) * a[0][j1] * determ(m,n-1) ;
                                            // sum x raised to y power
                                            // recursively get determinant of next
                                            // sub-matrix which is now one
                                            // row & column smaller

            for (i = 0 ; i < n-1 ; i++) free(m[i]) ;// free the storage allocated to
                                            // to this minor's set of pointers
            free(m) ;                       // free the storage for the original
                                            // pointer to pointer
        }
    }
    return(det) ;
}
/*****************************************************/
double ierfc(double y)
// inverse of the error function erfc from http://www.mathematik.uni-bielefeld.de/~sillke/ALGORITHMS/special-functions/ierfc.c
// Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp)
{
  const double IERFC_LIM = 27;
  double z = (y > 1)? 2-y: y;
  if (z < 1e-300) return (y > 1)? -IERFC_LIM : IERFC_LIM ;
  double w = 0.916461398268964 - log(z);
  double u = sqrt(w);
  double s = (log(u) + 0.488826640273108) / w;
  double t = 1 / (u + 0.231729200323405);
  double x = u * (1 - s * (s * 0.124610454613712 + 0.5)) -
    ((((-0.0728846765585675 * t + 0.269999308670029) * t +
       0.150689047360223) * t + 0.116065025341614) * t +
     0.499999303439796) * t;
  t = 3.97886080735226 / (x + 3.97886080735226);
  u = t - 0.5;
  s = (((((((((0.00112648096188977922 * u +
	       1.05739299623423047e-4) * u - 0.00351287146129100025) * u -
	     7.71708358954120939e-4) * u + 0.00685649426074558612) * u +
	   0.00339721910367775861) * u - 0.011274916933250487) * u -
	 0.0118598117047771104) * u + 0.0142961988697898018) * u +
       0.0346494207789099922) * u + 0.00220995927012179067;
  s = ((((((((((((s * u - 0.0743424357241784861) * u -
		 0.105872177941595488) * u + 0.0147297938331485121) * u +
	       0.316847638520135944) * u + 0.713657635868730364) * u +
	     1.05375024970847138) * u + 1.21448730779995237) * u +
	   1.16374581931560831) * u + 0.956464974744799006) * u +
	 0.686265948274097816) * u + 0.434397492331430115) * u +
       0.244044510593190935) * t -
    z * exp(x * x - 0.120782237635245222);
  x += s * (x * s + 1);
  return (y > 1)? -x: x;
}
/*****************************************************/
double drdz(struct cosmo cosmo,double z) {

  /* Returns dr/dz (in Mpc/h) at redshift z, for a universe with                                                                                                   
     matter, vacuum and radiation density parameters Omega_m,Omega_v,Omega_r*/

  double a,output;

  a =1.0/(1.0+z);

  /* note: h0 independent*/

  /* From Peacock p.89/Huterer and Turner 0012510 for w*/
  output = (3000.0)/(R0(cosmo)*
                     sqrt((1.0-cosmo.om-cosmo.od)/a/a+
                          cosmo.od*(pow(a,(-3.0*(1.0+cosmo.w0+cosmo.wa)))*
                                    exp(-3.0*(cosmo.wa*(1.0-a))))+
                          cosmo.om/a/a/a));
  //3000. from speed of light/100 Kms-1                                                                                                                            

  output = output/cosmo.h0;

  return output;
}
/*****************************************************/
double R0(struct cosmo cosmo) {

  /* Function return curvature scale R0 if k<>0. Units are h^{-1} Mpc                                                                                              
     If k=0, 1 is returned*/

  double k,ot,output;

  if (flatcosmo==1) {
    output = 1.0;
  } else {
    ot=cosmo.om+cosmo.od;
    output = (3000.0/cosmo.h0)/sqrt(sqrt((ot-1.0)*(ot-1.0)));
  }

  return output;
}
/*****************************************************/
