#include <Rcpp.h>

typedef std::complex<long double> cplx;

cplx fromRcplx(Rcomplex zr) {
  return cplx((long double)(zr.r), (long double)(zr.i));
}

Rcomplex toRcplx(cplx z) {
  Rcomplex zr;
  zr.r = (double)(z.real());
  zr.i = (double)(z.imag());
  return zr;
}

// [[Rcpp::export]]
Rcomplex Carlson_RF_(Rcomplex xr, Rcomplex yr, Rcomplex zr, double err) {
  cplx x = fromRcplx(xr);
  cplx y = fromRcplx(yr);
  cplx z = fromRcplx(zr);
  long double dx = 2.0 * (long double)err;
  long double dy = 2.0 * (long double)err;
  long double dz = 2.0 * (long double)err;
  cplx A;
  while(dx > err || dy > err || dz > err) {
    cplx srx = std::sqrt(x);
    cplx sry = std::sqrt(y);
    cplx srz = std::sqrt(z);
    cplx lambda = srx * sry + sry * srz + srz * srx;
    x = (x + lambda) / (long double)4.0;
    y = (y + lambda) / (long double)4.0;
    z = (z + lambda) / (long double)4.0;
    A = (x + y + z) / (long double)3.0;
    dx = std::abs((A - x)/A);
    dy = std::abs((A - y)/A);
    dz = std::abs((A - z)/A);
  }
  long double E2 = dx * dy + dy * dz + dz * dx;
  long double E3 = dy * dx * dz;
  long double h = (
    1.0 - E2/10.0 + E3/14.0 + E2*E2/24.0 - 3.0*E2*E3/44.0 - 5.0*E2*E2*E2/208.0
      + 3.0*E3*E3/104.0 + E2*E2*E3/16.0
  );
  cplx out = h / std::sqrt(A);
  return toRcplx(out);

}
