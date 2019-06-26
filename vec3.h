#include <math.h>

#include "math_const.h"

/* /\* 3-vector operations *\/ */
/* static inline double dot( const double *x1, const double *x2); */
/* static inline void cross( const double *x1, const double *x2, double *xout); */
/* static inline double mag( const double *x); */
/* static inline void vec_sum( const double *x1, const double *x2, double *xout); */
/* static inline void vec_dif( const double *x1, const double *x2, double *xout); */
/* static inline double dist( const double *x1, const double *x2); */
/* static inline void hat_vec( const double *x, double *xout); */
/* static inline void dif_hat( const double *x1, const double *x2, double *xout); */

/* /\* Matrix vector operations *\/ */
/* static inline double det3( const double *m); */
/* static inline void solve3( const double *m, const double *b, double *x); */
/* /\* void sym_ev_decomp( const double *s, double *lam, double *vec); *\/ */

/* /\* Triangle operations *\/ */
/* static inline void n_vec( const double *x1, const double *x2, const double *x3, */
/* 		   double *xout); */
/* static inline void n_hat( const double *x1, const double *x2, const double *x3, */
/* 		   double *xout); */
/* static inline void diheadral( const double *x1, const double *x2, */
/* 		       const double *x3, const double *x4); */
/* static inline double angle( const double *x1, const double *x2, const double *x3, */
/* 		     const double *x4, const double *x5); */

static inline double dot( const double *x1, const double *x2)
{
  int i;
  double sum=0;
  for( i=0; i<3; i++){
    sum += x1[i]*x2[i];
  }
  return(sum);
}
static inline double dot_d( const double *x1, const double *x2,
		     double *dx1, double *dx2){
  int i;
  double sum=0;
  for( i=0; i<3; i++){
    sum += x1[i]*x2[i];
    dx1[i] = x2[i];
    dx2[i] = x1[i];
  }
  return( sum);
}

static inline void cross( const double *x1, const double *x2, double *xout)
{
  int i;
  for(i=0; i<3; i++){
    xout[i] = x1[(i+1)%3]*x2[(i+2)%3]-x1[(i+2)%3]*x2[(i+1)%3];
  }
}
static inline void cross_d( const double *x1, const double *x2, double *xout,
		     double *dxout_dx1, double *dxout_dx2)
{
  int i;
  for( i=0; i<3; i++){
    xout[i] = x1[(i+1)%3]*x2[(i+2)%3]-x1[(i+2)%3]*x2[(i+1)%3];
    dxout_dx1[3*i+i] = 0.;
    dxout_dx1[3*i+(i+1)%3] = x2[(i+2)%3];
    dxout_dx1[3*i+(i+2)%3] = -x2[(i+1)%3];
    dxout_dx2[3*i+i] = 0.;
    dxout_dx2[3*i+(i+1)%3] = -x1[(i+2)%3];
    dxout_dx2[3*i+(i+2)%3] = x1[(i+1)%3];
  }
}

static inline double mag( const double *x)
{
  int i;
  double sum;
  sum = 0.;
  for( i=0; i<3; i++){
    sum += x[i]*x[i];
  }
  return(sqrt( sum));
}
static inline double mag_d( const double *x, double *dx)
{
  int i;
  double sum;
  sum = 0.;
  for( i=0; i<3; i++){
    sum += x[i]*x[i];
  }
  sum = sqrt(sum);
  if( sum == 0. ){
    for( i=0; i<3; i++){
      dx[i] = 1.;
    }
  }else{
    for( i=0; i<3; i++){
      dx[i] = x[i]/sum;
    }
  }
  return( sum);
}

static inline void vec_sum( const double *x1, const double *x2, double *xout)
{
  int i;
  for( i=0; i<3; i++){
    xout[i] = x1[i] + x2[i];
  }
}
static inline void vec_sum_d( const double *x1, const double *x2, double *xout,
		double *dxout_dx1, double *dxout_dx2)
{
  int i, j;
  for( i=0; i<3; i++){
    xout[i] = x1[i] + x2[i];
    for( j=0; j<3; j++){
      if( j==i ){
      	dxout_dx1[3*i+j] = 1.;
      	dxout_dx2[3*i+j] = 1.;
      }else{
      	dxout_dx1[3*i+j] = 0;
      	dxout_dx2[3*i+j] = 0;
      }
    }
  }
}

static inline void vec_dif( const double *x1, const double *x2, double *xout)
{
  int i;
  for( i=0; i<3; i++){
    xout[i] = x1[i] - x2[i];
  }
}

static inline void vec_dif_d( const double *x1, const double *x2, double *xout,
		double *dxout_dx1, double *dxout_dx2)
{
  int i, j;
  for( i=0; i<3; i++){
    xout[i] = x1[i] - x2[i];
    for( j=0; j<3; j++){
      if( j==i ){
      	dxout_dx1[3*i+j] = 1.;
      	dxout_dx2[3*i+j] = -1.;
      }else{
      	dxout_dx1[3*i+j] = 0;
	      dxout_dx2[3*i+j] = 0;
      }
    }
  }
}

static inline double dist( const double *x1, const double *x2)
{
  double t[3];
  vec_dif( x1, x2, t);
  return( mag( t));
}
static inline double dist_d( const double *x1, const double *x2,
		      double *dx1, double *dx2)
{
  int i, j;
  double b, db_dt[3], t[3], dt_dx1[9], dt_dx2[9];
  vec_dif_d( x1, x2, t, dt_dx1, dt_dx2);
  b = mag_d( t, db_dt);
  for( i=0; i<3; i++){
    dx1[i] = 0.;
    dx2[i] = 0.;
    for( j=0; j<3; j++){
      dx1[i] += db_dt[j]*dt_dx1[3*j+i];
      dx2[i] += db_dt[j]*dt_dx2[3*j+i];
    }
  }
  return( b);
}
 
static inline double dist2( const double *x1, const double *x2)
{
  double t[3];
  vec_dif( x1, x2, t);
  return( dot(t, t));
}

static inline void hat_vec( const double *x, double *xout)
{
  int i;
  double tmag;
  tmag = mag( x);
  if( tmag == 0. ){
    xout[0] = 1.;
    xout[1] = 0.;
    xout[2] = 0.;
  }else{
    for( i=0; i<3; i++){
      xout[i] = x[i]/tmag;
    }
  }
}
static inline void hat_vec_d( const double *x, double *xout, double *dxout_dx)
{
  int i, j;
  double tmag, dtmag_dx[3];
  tmag = mag_d( x, dtmag_dx);
  for( i=0; i<3; i++){
    xout[i] = x[i]/tmag;
    for( j=0; j<3; j++){
        dxout_dx[3*i+j] = -x[i]/(tmag*tmag)*dtmag_dx[j];
      if( i==j){
      	dxout_dx[3*i+j] += 1/(tmag);
      }
    }
  }
}

static inline void dif_hat( const double *x1, const double *x2, double *xout)
{
  double v[3];
  vec_dif( x1, x2, v);
  hat_vec( v, xout);
}
static inline void dif_hat_d( const double *x1, const double *x2, double *xout,
		       double *dxout_dx1, double *dxout_dx2)
{
  int i, j, k;
  double v[3], dv_dx1[9], dv_dx2[9], dxout_dv[9];;
  vec_dif_d( x1, x2, v, dv_dx1, dv_dx2);
  hat_vec_d( v, xout, dxout_dv);
  for( i=0; i<3; i++){
    for( j=0; j<3; j++){
      dxout_dx1[3*i+j] = 0.;
      dxout_dx2[3*i+j] = 0.;
      for( k=0; k<3; k++){
      	dxout_dx1[3*i+j] += dxout_dv[3*i+k]*dv_dx1[3*k+j];
      	dxout_dx2[3*i+j] += dxout_dv[3*i+k]*dv_dx2[3*k+j];
      }
    }
  }
}

static inline double det3( const double *m)
{
  int i, ip, ipp;
  double sum;
  sum = 0;
  for( i=0; i<3; i++){
    ip = (i+1)%3;
    ipp = (i+2)%3;
    sum += m[3*1+ip]*(m[3*0+i]*m[3*2+ipp] - m[3*2+i]*m[3*0+ipp]);
  }
  return( sum);
}

static inline void solve3( const double *m, const double *b, double *x)
{
  int i, j, k;
  double a, mt[9];
  
  a = det3( m);
  if( a != 0 ){
    for( k=0; k<3; k++){
      for( i=0; i<3; i++){
	for( j=0; j<3; j++){
	  mt[3*i+j] = ((k==j)? b[i] : m[3*i+j]) ;
	}
      }
      x[k] = det3(mt)/a;
    }
  }else{
    for( k=0; k<3; k++){
      x[k] = 0.;
    }
  }
}


static inline void n_vec( const double *x1, const double *x2, const double *x3,
		   double *nout)
{
  int i;
  double t1[3], t2[3], t3[3];
  cross( x1, x2, t1);
  cross( x2, x3, t2);
  cross( x3, x1, t3);
  for( i=0; i<3; i++){
    nout[i] = t1[i] + t2[i] + t3[i];
  }
}
static inline void n_vec_d( const double *x1, const double *x2, 
		     const double *x3, double *nout,
		     double *dnout_dx1, double *dnout_dx2, 
		     double *dnout_dx3)
{
  int i, j;
  double t1[3], t2[3], t3[3];
  double dt1_dx1[9], dt2_dx2[9], dt3_dx3[9];
  double dt3_dx1[9], dt1_dx2[9], dt2_dx3[9];
  cross_d( x1, x2, t1, dt1_dx1, dt1_dx2);
  cross_d( x2, x3, t2, dt2_dx2, dt2_dx3);
  cross_d( x3, x1, t3, dt3_dx3, dt3_dx1);
  for( i=0; i<3; i++){
    nout[i] = t1[i] + t2[i] + t3[i];
    for( j=0; j<3; j++){
      dnout_dx1[3*i+j] = dt1_dx1[3*i+j]+dt3_dx1[3*i+j];
      dnout_dx2[3*i+j] = dt1_dx2[3*i+j]+dt2_dx2[3*i+j];
      dnout_dx3[3*i+j] = dt2_dx3[3*i+j]+dt3_dx3[3*i+j];
    }
  }
}

static inline void n_hat( const double *x1, const double *x2, 
		   const double *x3, double *nout)
{
  double t[3];
  n_vec( x1, x2, x3, t);
  hat_vec( t, nout);
}
static inline void n_hat_d( const double *x1, const double *x2, 
		     const double *x3, double *nout,
		     double *dnout_dx1, double *dnout_dx2, 
		     double *dnout_dx3)
{
  int i, j, k;
  double t[3], dt_dx1[9], dt_dx2[9], dt_dx3[9], dnout_dt[9];
  n_vec_d( x1, x2, x3, t, dt_dx1, dt_dx2, dt_dx3);
  hat_vec_d( t, nout, dnout_dt);
  for( i=0; i<3; i++){
    for( j=0; j<3; j++){
      dnout_dx1[3*i+j] = 0.;
      dnout_dx2[3*i+j] = 0.;
      dnout_dx3[3*i+j] = 0.;
      for( k=0; k<3; k++){
	dnout_dx1[3*i+j] += dnout_dt[3*i+k]*dt_dx1[3*k+j];
	dnout_dx2[3*i+j] += dnout_dt[3*i+k]*dt_dx2[3*k+j];
	dnout_dx3[3*i+j] += dnout_dt[3*i+k]*dt_dx3[3*k+j];
      }
    }
  }
}

static inline void dihedral( const double *x1, const double *x2,
		      const double *x3, const double *x4,
		      double *co, double *si)
{
  double n1[3], n2[3], t[3], siv[3];
  n_hat( x1, x2, x3, n1);
  n_hat( x2, x1, x4, n2);
  dif_hat( x2, x1, t);
  cross( n1, n2, siv);
  *co = dot( n1, n2);
  if( dot( t, siv) >= 0. ){
    *si = mag( siv);
  }else{
    *si = -mag( siv);
  }
}
  
static inline void dihedral_d( const double *x1, const double *x2,
			const double *x3, const double *x4,
			double *co, double *si,
			double *dco_dx1, double *dco_dx2,
			double *dco_dx3, double *dco_dx4,
			double *dsi_dx1, double *dsi_dx2,
			double *dsi_dx3, double *dsi_dx4)
{
  int i, j, k, sign;
  double n1[3], dn1_dx1[9], dn1_dx2[9], dn1_dx3[9];
  double n2[3], dn2_dx2[9], dn2_dx1[9], dn2_dx4[9];
  double t[3], siv[3], dsiv_dn1[9], dsiv_dn2[9];
  double dsiv_dx1[9], dsiv_dx2[9], dsiv_dx3[9], dsiv_dx4[9];
  double dsi_dsiv[3];
  double dco_dn1[3], dco_dn2[3];

  n_hat_d( x1, x2, x3, n1, dn1_dx1, dn1_dx2, dn1_dx3);
  n_hat_d( x2, x1, x4, n2, dn2_dx2, dn2_dx1, dn2_dx4);
  dif_hat( x2, x1, t);
  cross_d( n1, n2, siv, dsiv_dn1, dsiv_dn2);
  for( i=0; i<3; i++){
    for( j=0; j<3; j++){
      dsiv_dx1[3*i+j] = 0.;
      dsiv_dx2[3*i+j] = 0.;
      dsiv_dx3[3*i+j] = 0.;
      dsiv_dx4[3*i+j] = 0.;
      for( k=0; k<3; k++){
      	dsiv_dx1[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx1[3*k+j] +
	        dsiv_dn2[3*i+k]*dn2_dx1[3*k+j];
	      dsiv_dx2[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx2[3*k+j] +
	        dsiv_dn2[3*i+k]*dn2_dx2[3*k+j];
      	dsiv_dx3[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx3[3*k+j];
	      dsiv_dx4[3*i+j] += dsiv_dn2[3*i+k]*dn2_dx4[3*k+j];
      }
    }
  }
  *co = dot_d( n1, n2, dco_dn1, dco_dn2);
  for( i=0; i<3; i++){
    dco_dx1[i] = 0.;
    dco_dx2[i] = 0.;
    dco_dx3[i] = 0.;
    dco_dx4[i] = 0.;
    for( j=0; j<3; j++){
      dco_dx1[i] += dco_dn1[j]*dn1_dx1[3*j+i] +
      	dco_dn2[j]*dn2_dx1[3*j+i];
      dco_dx2[i] += dco_dn1[j]*dn1_dx2[3*j+i] +
      	dco_dn2[j]*dn2_dx2[3*j+i];
      dco_dx3[i] += dco_dn1[j]*dn1_dx3[3*j+i];
      dco_dx4[i] += dco_dn2[j]*dn2_dx4[3*j+i];
    }
  }
  if( dot( t, siv) >= 0. ){
    sign = 1;
  }else{
    sign = -1;
  }
  *si = ((double)sign)*mag_d( siv, dsi_dsiv);
  for( i=0; i<3; i++){
    dsi_dx1[i] = 0.;
    dsi_dx2[i] = 0.;
    dsi_dx3[i] = 0.;
    dsi_dx4[i] = 0.;
    for( j=0; j<3; j++){
      dsi_dx1[i] += dsi_dsiv[j]*dsiv_dx1[3*j+i];
      dsi_dx2[i] += dsi_dsiv[j]*dsiv_dx2[3*j+i];
      dsi_dx3[i] += dsi_dsiv[j]*dsiv_dx3[3*j+i];
      dsi_dx4[i] += dsi_dsiv[j]*dsiv_dx4[3*j+i];
    }
    dsi_dx1[i] *= ((double) sign);
    dsi_dx2[i] *= ((double) sign);
    dsi_dx3[i] *= ((double) sign);
    dsi_dx4[i] *= ((double) sign);
  }
}

/* static inline void sin_vec( const double *x1, const double *x2, */
/* 		     const double *x3, const double *x4, */
/* 		     double *siv) */
/* { */
/*   double n1[3], n2[3]; */
/*   n_hat( x1, x2, x3, n1); */
/*   n_hat( x2, x1, x4, n2); */
/*   cross( n1, n2, siv); */
/* } */
  
/* static inline void sin_vec_d( const double *x1, const double *x2, */
/* 		       const double *x3, const double *x4, */
/* 		       double *siv, */
/* 		       double *dsiv_dx1, double *dsiv_dx2, */
/* 		       double *dsiv_dx3, double *dsiv_dx4) */
/* { */
/*   int i, j, k, sign; */
/*   double n1[3], dn1_dx1[9], dn1_dx2[9], dn1_dx3[9]; */
/*   double n2[3], dn2_dx2[9], dn2_dx1[9], dn2_dx4[9]; */
/*   double dsiv_dn1[9], dsiv_dn2[9]; */

/*   n_hat_d( x1, x2, x3, n1, dn1_dx1, dn1_dx2, dn1_dx3); */
/*   n_hat_d( x2, x1, x4, n2, dn2_dx2, dn2_dx1, dn2_dx4); */
/*   cross_d( n1, n2, siv, dsiv_dn1, dsiv_dn2); */
/*   for( i=0; i<3; i++){ */
/*     for( j=0; j<3; j++){ */
/*       dsiv_dx1[3*i+j] = 0.; */
/*       dsiv_dx2[3*i+j] = 0.; */
/*       dsiv_dx3[3*i+j] = 0.; */
/*       dsiv_dx4[3*i+j] = 0.; */
/*       for( k=0; k<3; k++){ */
/* 	dsiv_dx1[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx1[3*k+j] + */
/* 	  dsiv_dn2[3*i+k]*dn2_dx1[3*k+j]; */
/* 	dsiv_dx2[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx2[3*k+j] + */
/* 	  dsiv_dn2[3*i+k]*dn2_dx2[3*k+j]; */
/* 	dsiv_dx3[3*i+j] += dsiv_dn1[3*i+k]*dn1_dx3[3*k+j]; */
/* 	dsiv_dx4[3*i+j] += dsiv_dn2[3*i+k]*dn2_dx4[3*k+j]; */
/*       } */
/*     } */
/*   } */
/* } */

/*!
 * The `angle` function calculates the angle bewteen two edges i_l and
 * i_r on loop j. It does so by calculating the cosine of the angle
 * between the vectors defined by the edges. It gets the quadrant by
 * comparing with the average of the normals of the triangles that
 * form the edges.
 */
static inline double angle( const double *x1, const double *x2, const double *x3,
		     const double *x4, const double *x5)
{
  double v1[3], v2[3], n1[3], n2[3], t[3], s[3], a;

  dif_hat( x1, x2, v1);
  dif_hat( x3, x2, v2);
  n_hat( x1, x2, x4, n1);
  n_hat( x2, x3, x5, n2);
  vec_sum( n1, n2, t);
  cross( v1, v2, s);
  if( dot( t, s) >= 0.){
    a = acos( dot( v1, v2) );
  }else{
    a = TWOPI- acos( dot( v1, v2) );
  }
  return( a);
}

/*!
 * The `intersect` function: This function test whether a line from va
 * to vb crosses the triangle formed between v0,v1, and v2. To do so,
 * it solves a matrix equation for where the intesection between the
 * line and plane occurs (if they extend to infinity they either cross
 * or are parallel) and test whether this point lies within the
 * triangle. This to see if they interesect, I check if the
 * determinant is equal to zero. If it is not, I solve the system, and
 * check to point.
 */
static inline int intersect( const double *va, const double *vb,
		      const double *v0, const double *v1, const double *v2)
{
  int i;
  double x[3], b[3], m[9], a;
  double epsilon = 0.05;
  for( i=0; i<3; i++){
    m[3*i+0] = va[i]-vb[i];
    m[3*i+1] = v1[i]-v0[i];
    m[3*i+2] = v2[i]-v0[i];
    b[i] = va[i]-v0[i];
  }
  a = det3( m);
  if( a != 0. ){
    solve3( m, b, x);
    if( 0.+epsilon < x[0] && x[0] < 1.-epsilon &&
	0.+epsilon < x[1] && x[1] < 1.-epsilon &&
	0.+epsilon < x[2] && x[2] < 1.-epsilon &&
	x[1] + x[2] < 1.+epsilon){
      return( 1);
    }
  }
  return( 0);
}

  
