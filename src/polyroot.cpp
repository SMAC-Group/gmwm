/* Modified source from complex.c in R's directory
 * 
 * Formerly src/appl/cpoly.c:
 *
 *  Copyright (C) 1997-1998 Ross Ihaka
 *  Copyright (C) 1999-2001 R Core Team
 *
 *	cpoly finds the zeros of a complex polynomial.
 *
 *	On Entry
 *
 *	opr, opi      -	 double precision vectors of real and
 *			 imaginary parts of the coefficients in
 *			 order of decreasing powers.
 *
 *	degree	      -	 int degree of polynomial.
 *
 *
 *	On Return
 *
 *	zeror, zeroi  -	 output double precision vectors of
 *			 real and imaginary parts of the zeros.
 *
 *	fail	      -	 output int parameter,	true  only if
 *			 leading coefficient is zero or if cpoly
 *			 has found fewer than degree zeros.
 *
 *	The program has been written to reduce the chance of overflow
 *	occurring. If it does occur, there is still a possibility that
 *	the zerofinder will work provided the overflowed quantity is
 *	replaced by a large number.
 *
 *	This is a C translation of the following.
 *
 *	TOMS Algorithm 419
 *	Jenkins and Traub.
 *	Comm. ACM 15 (1972) 97-99.
 *
 *	Ross Ihaka
 *	February 1997
 */

#include <RcppArmadillo.h>
#include <cfloat> // Needed for contants DBL_EPSILON, DBL_MAX, M_SQRT2
#include <cmath> // Used for std::isfinite

#include "polyroot.h"


double myfmod_cpp(double x1, double x2)
{
  return x1 - floor(x1 / x2) * x2;
}


double R_pow_cpp(double x, double y) /* = x ^ y */
{
  if(x == 1. || y == 0.0)
    return(1.);
  if(x == 0.0) {
    if(y > 0.0) return(0.0);
    /* y < 0 */return(std::numeric_limits<double>::infinity());
  }
  if (std::isfinite(x) && std::isfinite(y))
    return(pow(x,y));
  if (std::isnan(x) || std::isnan(y)) {
    return(std::numeric_limits<double>::quiet_NaN());
  }
  if(!std::isfinite(x)) {
    if(x > 0)		/* Inf ^ y */
    return((y < 0.0)? 0.0 : std::numeric_limits<double>::infinity() );
    else {			/* (-Inf) ^ y */
    if(std::isfinite(y) && y == floor(y)) /* (-Inf) ^ n */
    return((y < 0.0) ? 0.0 : (myfmod_cpp(y,2.) ? x  : -x));
    }
  }
  if(!std::isfinite(y)) {
    if(x >= 0) {
      if(y > 0)		/* y == +Inf */
    return((x >= 1)? std::numeric_limits<double>::infinity() : 0.0);
      else		/* y == -Inf */
    return((x < 1) ? std::numeric_limits<double>::infinity() : 0.0);
    }
  }
  return(std::numeric_limits<double>::quiet_NaN());		// all other cases: (-Inf)^{+-Inf, non-int}; (neg)^{+-Inf} 
}

double R_pow_di_cpp(double x, int n)
{
  double pow = 1.0;
  
  if (std::isnan(x)) return x;
  if (n != 0) {
    if (!std::isfinite(x)) return R_pow_cpp(x, (double)n);
    if (n < 0) { n = -n; x = 1/x; }
    for(;;) {
      if(n & 01) pow *= x;
      if(n >>= 1) x *= x; else break;
    }
  }
  return pow;
}

//' @title Root Finding C++
//' @description Used to interface with Armadillo
//' @param z A \code{cx_vec} (complex vector) that has 1 in the beginning (e.g. c(1,3i,-3i))
//' @keywords internal
// [[Rcpp::export]]
arma::cx_vec do_polyroot_arma(const arma::cx_vec& z)
{
  
  std::vector<double> zr, zi, rr, ri;
  arma::cx_vec r;
  
  std::complex<double> obs;
  
  bool fail;
  int degree, i, n;
  
  // Find out the size of n
  n = z.size();
  
  // Obtain degree
  degree = 0;
  for(i = 0; i < n; i++) {
    if(real(z(i)) != 0.0 || imag(z(i)) != 0.0) degree = i;
  }
  
  n = degree + 1; // omit trailing zeroes
  if(degree >= 1) {
    rr = std::vector<double>(n);
    ri = rr;
    zr = rr;
    zi = rr;
    
    for(i = 0 ; i < n ; i++) {
      if(!std::isfinite(z[i].real()) || !std::isfinite(z[i].imag())){
        throw std::invalid_argument( "Invalid polynomial coefficient" );
      }
      zr[degree-i] = z[i].real();
      zi[degree-i] = z[i].imag();
    }
    polyroot_cpp(zr, zi, degree, rr, ri, fail);
    if(fail){ throw std::runtime_error("Root finding code failed!"); }
    
    r = arma::cx_vec(degree); 
    for(i = 0 ; i < degree ; i++) {
      obs.real(rr[i]);
      obs.imag(ri[i]);
      r(i) = obs;
    }
  }
  else {
    r = arma::zeros<arma::cx_vec>(0);
  }
  return r;
}


//' @title Root Finding C++
//' @description Vroom Vroom
//' @param z A \code{vec<complex<double>} (complex vector) that has 1 in the beginning (e.g. c(1,3i,-3i))
//' @keywords internal
// [[Rcpp::export]]
std::vector< std::complex<double> > do_polyroot_cpp(const std::vector< std::complex<double> >& z)
{
  
  std::vector<double> zr, zi, rr, ri;
  std::vector<std::complex<double> > r;
  
  bool fail;
  int degree, i, n;

  // Find out the size of n
  n = z.size();
  
  // Obtain degree
  degree = 0;
  for(i = 0; i < n; i++) {
    if(z[i].real() != 0.0 || z[i].imag() != 0.0) degree = i;
  }
  
  n = degree + 1; // omit trailing zeroes
  if(degree >= 1) {
    rr = std::vector<double>(n);
    ri = rr;
    zr = rr;
    zi = rr;
    
    for(i = 0 ; i < n ; i++) {
      if(!std::isfinite(z[i].real()) || !std::isfinite(z[i].imag())){
        throw std::invalid_argument( "Invalid polynomial coefficient" );
      }
      zr[degree-i] = z[i].real();
      zi[degree-i] = z[i].imag();
    }
    polyroot_cpp(zr, zi, degree, rr, ri, fail);
    if(fail){ throw std::runtime_error("Root finding code failed!"); }

    r = std::vector<std::complex<double> >(degree); 
    for(i = 0 ; i < degree ; i++) {
      r[i].real(rr[i]);
      r[i].imag(ri[i]);
    }
  }
  else {
    r = std::vector<std::complex<double> >(0);
  }
  return r;
}

  

// Global Variables (too many!)

static int nn;
static std::vector<double> pr, pi, hr, hi, qpr, qpi, qhr, qhi, shr, shi;

static double sr, si;
static double tr, ti;
static double pvr, pvi;

static const double eta =  DBL_EPSILON;
static const double are = /* eta = */DBL_EPSILON;
static const double mre = 2.0 * M_SQRT2 * /* eta, i.e. */DBL_EPSILON;
static const double infin = DBL_MAX;

void polyroot_cpp(const std::vector<double> &opr, const std::vector<double> &opi, int &degree,
                         std::vector<double> &zeror, std::vector<double> &zeroi, bool &fail)
{
  static const double smalno = DBL_MIN;
  static const double base = (double)FLT_RADIX;
  static int d_n, i, i1, i2;
  static double zi, zr, xx, yy;
  static double bnd, xxx;
  bool conv;
  int d1;
  static const double cosr =/* cos 94 */ -0.06975647374412529990;
  static const double sinr =/* sin 94 */  0.99756405025982424767;
  xx = M_SQRT1_2;/* 1/sqrt(2) = 0.707.... */
  
  yy = -xx;
  fail = false;
  
  nn = degree;
  d1 = nn - 1;
  
  // Algorithm fails if the leading coefficient is zero. 
  if (opr[0] == 0.0 && opi[0] == 0.0) {
    fail = true;
    return;
  }
  
  // remove the zeros at the origin if any.
  while (opr[nn] == 0.0 && opi[nn] == 0.0) {
    d_n = d1-nn+1;
    zeror[d_n] = 0.0;
    zeroi[d_n] = 0.0;
    nn--;
  }
  nn++;
  /*-- Now, global var.  nn := #{coefficients} = (relevant degree)+1 */
  
  if (nn == 1) return;
  
  // Use a single allocation as these as small 
  pr = std::vector<double>(10*nn); 
  pi = pr; 
  hr = pr; 
  hi = pr;
  qpr = pr; 
  qpi = pr; 
  qhr = pr; 
  qhi = pr;
  shr = pr; 
  shi = pr; 
  
  // make a copy of the coefficients and shr[] = | p[] | 
  for (i = 0; i < nn; i++) {
    pr[i] = opr[i];
    pi[i] = opi[i];
    shr[i] = hypot(pr[i], pi[i]);
  }
  
  // scale the polynomial with factor 'bnd'
  bnd = cpoly_scale_cpp(nn, shr, eta, infin, smalno, base);
  if (bnd != 1.0) {
    for (i=0; i < nn; i++) {
      pr[i] *= bnd;
      pi[i] *= bnd;
    }
  }
  
  // Start the algorithm for one zero 
  
  while (nn > 2) {
    
    // Calculate bnd, a lower bound on the modulus of the zeros. 
    
    for (i=0 ; i < nn ; i++){
      shr[i] = hypot(pr[i], pi[i]);
    }
    bnd = cpoly_cauchy_cpp(nn, shr, shi);
    
    // Outer loop to control 2 major passes
    // with different sequences of shifts
    
    for (i1 = 1; i1 <= 2; i1++) {
      
      // First stage calculation, no shift
      
      noshft_cpp(5);
      
      //	Inner loop to select a shift
      for (i2 = 1; i2 <= 9; i2++) {
        
        // Shift is chosen with modulus bnd 
        // and amplitude rotated by 94 degrees 
        // from the previous shift 
        
        xxx= cosr * xx - sinr * yy;
        yy = sinr * xx + cosr * yy;
        xx = xxx;
        sr = bnd * xx;
        si = bnd * yy;
        
        // Second stage calculation, fixed shift
        
        conv = fxshft_cpp(i2 * 10, zr, zi);
        if (conv)
          goto L10;
      }
    }
    
    // The zerofinder has failed on two major passes
    // return empty handed
    
    fail = true;
    return;
    
    // The second stage jumps directly to the third stage iteration.
    // If successful, the zero is stored and the polynomial deflated.

    L10:
      d_n = d1+2 - nn;
    zeror[d_n] = zr;
    zeroi[d_n] = zi;
    --nn;
    for (i=0; i < nn ; i++) {
      pr[i] = qpr[i];
      pi[i] = qpi[i];
    }
  }// End while loop
    
  // Calculate the final zero and return
  cdivid_cpp(-pr[1], -pi[1], pr[0], pi[0], zeror[d1], zeroi[d1]);
  
  return;
}


/*  Computes the derivative polynomial as the initial
*  polynomial and computes l1 no-shift h polynomials.	*/

void noshft_cpp(int l1)
{
  int i, j, jj, n = nn - 1, nm1 = n - 1;
  
  double t1, t2, xni;
  
  for (i=0; i < n; i++) {
    xni = (double)(nn - i - 1);
    hr[i] = xni * pr[i] / n;
    hi[i] = xni * pi[i] / n;
  }
  
  for (jj = 1; jj <= l1; jj++) {
    
    if (hypot(hr[n-1], hi[n-1]) <=
        eta * 10.0 * hypot(pr[n-1], pi[n-1])) {
      /*	If the constant term is essentially zero, */
      /*	shift h coefficients. */
      
      for (i = 1; i <= nm1; i++) {
        j = nn - i;
        hr[j-1] = hr[j-2];
        hi[j-1] = hi[j-2];
      }
      hr[0] = 0.0;
      hi[0] = 0.0;
    }
    else {
      cdivid_cpp(-pr[nn-1], -pi[nn-1], hr[n-1], hi[n-1], tr, ti);
      for (i = 1; i <= nm1; i++) {
        j = nn - i;
        t1 = hr[j-2];
        t2 = hi[j-2];
        hr[j-1] = tr * t1 - ti * t2 + pr[j-1];
        hi[j-1] = tr * t2 + ti * t1 + pi[j-1];
      }
      hr[0] = pr[0];
      hi[0] = pi[0];
    }
  }
}


/*  Computes l2 fixed-shift h polynomials and tests for convergence.
*  initiates a variable-shift iteration and returns with the
*  approximate zero if successful.
*/

/*  l2	  - limit of fixed shift steps
 *  zr,zi - approximate zero if convergence (result true)
 *
 * Return value indicates convergence of stage 3 iteration
 *
 * Uses global (sr,si), nn, pr[], pi[], .. (all args of polyev_cpp() !)
 */
bool fxshft_cpp(const int l2, double &zr, double &zi)
{

  
  bool pasd, bol, test;
  static double svsi, svsr;
  static int i, j, n;
  static double oti, otr;
  
  n = nn - 1;
  
  /* evaluate p at s. */
  
  polyev_cpp(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi);
  
  test = true;
  pasd = false;
  
  /* calculate first t = -p(s)/h(s). */
  
  calct_cpp(bol);
  
  /* main loop for one second stage step. */
  
  for (j=1; j<=l2; j++) {
    
    otr = tr;
    oti = ti;
    
    /* compute next h polynomial and new t. */
    
    nexth_cpp(bol);
    calct_cpp(bol);
    zr = sr + tr;
    zi = si + ti;
    
    /* test for convergence unless stage 3 has */
    /* failed once or this is the last h polynomial. */
    
    if (!bol && test && j != l2) {
      if (hypot(tr - otr, ti - oti) >= hypot(zr, zi) * 0.5) {
        pasd = false;
      }
      else if (! pasd) {
        pasd = true;
      }
      else {
        
        /* the weak convergence test has been */
        /* passed twice, start the third stage */
        /* iteration, after saving the current */
        /* h polynomial and shift. */
        
        for (i = 0; i < n; i++) {
          shr[i] = hr[i];
          shi[i] = hi[i];
        }
        svsr = sr;
        svsi = si;
        if (vrshft_cpp(10, zr, zi)) {
          return true;
        }
        
        /* the iteration failed to converge. */
        /* turn off testing and restore */
        /* h, s, pv and t. */
        
        test = false;
        for (i=1 ; i<=n ; i++) {
          hr[i-1] = shr[i-1];
          hi[i-1] = shi[i-1];
        }
        sr = svsr;
        si = svsi;
        polyev_cpp(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi);
        calct_cpp(bol);
      }
    }
  }
  
  /* attempt an iteration with final h polynomial */
  /* from second stage. */
  
  return(vrshft_cpp(10, zr, zi));
}


/* carries out the third stage iteration.
*/
bool vrshft_cpp(const int l3, double &zr, double &zi)
{
  /*  l3	    - limit of steps in stage 3.
  *  zr,zi   - on entry contains the initial iterate;
  *	      if the iteration converges it contains
  *	      the final iterate on exit.
  * Returns true if iteration converges
  *
  * Assign and uses  GLOBAL sr, si
  */
  bool bol, b;
  static int i, j;
  static double r1, r2, mp, ms, tp, relstp;
  static double omp;
  
  b = false;
  sr = zr;
  si = zi;
  
  /* main loop for stage three */
  
  for (i = 1; i <= l3; i++) {
    
    /* evaluate p at s and test for convergence. */
    polyev_cpp(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi);
    
    mp = hypot(pvr, pvi);
    ms = hypot(sr, si);
    if (mp <=  20.0 * errev_cpp(nn, qpr, qpi, ms, mp, /*are=*/eta, mre)) {
      goto L_conv;
    }
    
    /* polynomial value is smaller in value than */
    /* a bound on the error in evaluating p, */
    /* terminate the iteration. */
    
    if (i != 1) {
      
      if (!b && mp >= omp && relstp < .05) {
        
        /* iteration has stalled. probably a */
        /* cluster of zeros. do 5 fixed shift */
        /* steps into the cluster to force */
        /* one zero to dominate. */
        
        tp = relstp;
        b = true;
        if (relstp < eta)
          tp = eta;
        r1 = sqrt(tp);
        r2 = sr * (r1 + 1.) - si * r1;
        si = sr * r1 + si * (r1 + 1.);
        sr = r2;
        polyev_cpp(nn, sr, si, pr, pi, qpr, qpi, pvr, pvi);
        for (j = 1; j <= 5; ++j) {
          calct_cpp(bol);
          nexth_cpp(bol);
        }
        omp = infin;
        goto L10;
      }
      else {
        
        /* exit if polynomial value */
        /* increases significantly. */
        
        if (mp * .1 > omp)
          return false;
      }
    }
    omp = mp;
    
    /* calculate next iterate. */
    
    L10:
      calct_cpp(bol);
    nexth_cpp(bol);
    calct_cpp(bol);
    if (!bol) {
      relstp = hypot(tr, ti) / hypot(sr, si);
      sr += tr;
      si += ti;
    }
  }
  return false;
  
  L_conv:
    zr = sr;
  zi = si;
  return true;
}

void calct_cpp(bool &bol)
{
  /* computes	 t = -p(s)/h(s).
  * bol   - logical, set true if h(s) is essentially zero.	*/
  
  int n = nn - 1;
  double hvi, hvr;
  
  /* evaluate h(s). */
  polyev_cpp(n, sr, si, hr, hi,
         qhr, qhi, hvr, hvi);
  
  bol = hypot(hvr, hvi) <= are * 10.0 * hypot(hr[n-1], hi[n-1]);
  if (!bol) {
    cdivid_cpp(-pvr, -pvi, hvr, hvi, tr, ti);
  }
  else {
    tr = 0.0;
    ti = 0.0;
  }
}

void nexth_cpp(bool bol)
{
  /* calculates the next shifted h polynomial.
  * bol :	if true  h(s) is essentially zero
  */
  int j, n = nn - 1;
  double t1, t2;
  
  if (!bol) {
    for (j=1; j < n; j++) {
      t1 = qhr[j - 1];
      t2 = qhi[j - 1];
      hr[j] = tr * t1 - ti * t2 + qpr[j];
      hi[j] = tr * t2 + ti * t1 + qpi[j];
    }
    hr[0] = qpr[0];
    hi[0] = qpi[0];
  }
  else {
    /* if h(s) is zero replace h with qh. */
    
    for (j=1; j < n; j++) {
      hr[j] = qhr[j-1];
      hi[j] = qhi[j-1];
    }
    hr[0] = 0.0;
    hi[0] = 0.0;
  }
}

/*--------------------- Independent Complex Polynomial Utilities ----------*/


void polyev_cpp(int n,
              double s_r, double s_i,
              std::vector<double> &p_r, std::vector<double> &p_i,
              std::vector<double> &q_r, std::vector<double> &q_i,
              double &v_r, double &v_i)
  {
    /* evaluates a polynomial  p  at  s	 by the horner recurrence
    * placing the partial sums in q and the computed value in v_.
    */
    int i;
    double t;
    
    q_r[0] = p_r[0];
    q_i[0] = p_i[0];
    v_r = q_r[0];
    v_i = q_i[0];
    for (i = 1; i < n; i++) {
      t = v_r * s_r - v_i * s_i + p_r[i];
      q_i[i] = v_i = v_r * s_i + v_i * s_r + p_i[i];
      q_r[i] = v_r = t;
    }
  }


  double errev_cpp(const int n, 
               const std::vector<double> &qr, std::vector<double> &qi,
               const double ms, const double mp, const double a_re, const double m_re)  {
    /*	bounds the error in evaluating the polynomial by the horner
    *	recurrence.
    *
    *	qr,qi	 - the partial sum vectors
    *	ms	 - modulus of the point
    *	mp	 - modulus of polynomial value
    * a_re,m_re - error bounds on complex addition and multiplication
    */
    double e;
    int i;
    
    e = hypot(qr[0], qi[0]) * m_re / (a_re + m_re);
    for (i=0; i < n; i++)
      e = e*ms + hypot(qr[i], qi[i]);
    
    return e * (a_re + m_re) - mp * m_re;
  }

double cpoly_cauchy_cpp(const int n, std::vector<double> &pot, std::vector<double> &q)
  {
    /* Computes a lower bound on the moduli of the zeros of a polynomial
    * pot[1:nn] is the modulus of the coefficients.
    */
    double f, x, delf, dx, xm;
    int i, n1 = n - 1;
    
    pot[n1] = -pot[n1];
    
    /* compute upper estimate of bound. */
    
    x = exp((log(-pot[n1]) - log(pot[0])) / (double) n1);
    
    /* if newton step at the origin is better, use it. */
    
    if (pot[n1-1] != 0.0) {
      xm = -pot[n1] / pot[n1-1];
      if (xm < x)
        x = xm;
    }
    
    /* chop the interval (0,x) unitl f le 0. */
    
    for(;;) {
      xm = x * 0.1;
      f = pot[0];
      for (i = 1; i < n; i++)
        f = f * xm + pot[i];
      if (f <= 0.0) {
        break;
      }
      x = xm;
    }
    
    dx = x;
    
    /* do Newton iteration until x converges to two decimal places. */
    
    while (fabs(dx / x) > 0.005) {
      q[0] = pot[0];
      for(i = 1; i < n; i++)
        q[i] = q[i-1] * x + pot[i];
      f = q[n1];
      delf = q[0];
      for(i = 1; i < n1; i++)
        delf = delf * x + q[i];
      dx = f / delf;
      x -= dx;
    }
    return x;
  }

double cpoly_scale_cpp(const int n, std::vector<double> &pot,
                          const double eps, const double BIG,
                          const double small, const double base)
  {
    /* Returns a scale factor to multiply the coefficients of the polynomial.
    * The scaling is done to avoid overflow and to avoid
    *	undetected underflow interfering with the convergence criterion.
    * The factor is a power of the base.
    * pot[1:n] : modulus of coefficients of p
    * eps,BIG,
    * small,base - constants describing the floating point arithmetic.
    */
    
    int i, ell;
    double x, high, sc, lo, min_, max_;
    
    /* find largest and smallest moduli of coefficients. */
    high = sqrt(BIG);
    lo = small / eps;
    max_ = 0.0;
    min_ = BIG;
    for (i = 0; i < n; i++) {
      x = pot[i];
      if (x > max_) max_ = x;
      if (x != 0.0 && x < min_)
        min_ = x;
    }
    
    /* scale only if there are very large or very small components. */
    
    if (min_ < lo || max_ > high) {
      x = lo / min_;
      if (x <= 1.0)
        sc = 1.0 / (sqrt(max_) * sqrt(min_));
      else {
        sc = x;
        if (BIG / sc > max_)
          sc = 1.0;
      }
      ell = (int) (log(sc) / log(base) + 0.5);
      return R_pow_di(base, ell);
    }
    else return 1.0;
  }


void cdivid_cpp(const double ar, const double ai, 
                   const double br, const double bi, 
                   double &cr, double &ci) {
      /* complex division c = a/b, i.e., (cr +i*ci) = (ar +i*ai) / (br +i*bi),
     avoiding overflow. */
    
    double d, r;
    
    if (br == 0.0 && bi == 0.0) {
      /* division by zero, c = infinity. */
      cr = ci = std::numeric_limits<double>::infinity();
    }  
    else if (fabs(br) >= fabs(bi)) {
      r = bi / br;
      d = br + r * bi;
      cr = (ar + ai * r) / d;
      ci = (ai - ar * r) / d;
    }  
    else {
      r = br / bi;
      d = bi + r * br;
      cr = (ar * r + ai) / d;
      ci = (ai * r - ar) / d;
    }
}
       