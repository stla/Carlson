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
  long double dy = dx;
  long double dz = dx;
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
  long double h =
    1.0 - E2/10.0 + E3/14.0 + E2*E2/24.0 - 3.0*E2*E3/44.0 - 5.0*E2*E2*E2/208.0
      + 3.0*E3*E3/104.0 + E2*E2*E3/16.0;
  cplx out = h / std::sqrt(A);
  return toRcplx(out);
}

// [[Rcpp::export]]
Rcomplex Carlson_RD_(Rcomplex xr, Rcomplex yr, Rcomplex zr, double err) {
  cplx x = fromRcplx(xr);
  cplx y = fromRcplx(yr);
  cplx z = fromRcplx(zr);
  long double dx = 2.0 * (long double)err;
  long double dy = dx;
  long double dz = dx;
  cplx s(0.0, 0.0);
  long double fac = 1.0;
  cplx A;
  while(dx > err || dy > err || dz > err) {
    cplx srx = std::sqrt(x);
    cplx sry = std::sqrt(y);
    cplx srz = std::sqrt(z);
    cplx lambda = srx * sry + sry * srz + srz * srx;
    s += fac / (srz * (z + lambda));
    fac /= 4.0;
    x = (x + lambda) / (long double)4.0;
    y = (y + lambda) / (long double)4.0;
    z = (z + lambda) / (long double)4.0;
    A = (x + y + ((long double)3.0)*z) / (long double)5.0;
    dx = std::abs((A - x)/A);
    dy = std::abs((A - y)/A);
    dz = std::abs((A - z)/A);
  }
  long double E2 = dx*dy + 3.0*dy*dz + 3.0*dz*dz + 3.0*dx*dz;
  long double E3 = dz*dz*dz + 3.0*dx*dz*dz + 3.0*dx*dy*dz + 3.0*dy*dz*dz;
  long double E4 = dy*dz*dz*dz + dx*dz*dz*dz + 3.0*dx*dy*dz*dz;
  long double E5 = dx*dy*dz*dz*dz;
  long double h =
    fac * (1.0 - 3.0*E2/14.0 + E3/6.0 + 9.0*E2*E2/88.0 - 3.0*E4/22.0
             - 9.0*E2*E3/52.0 + 3.0*E5/26.0 - E2*E2*E2/16.0 + 3.0*E3*E3/40.0
             + 3.0*E2*E4/20.0 + 45.0*E2*E2*E3/272.0 - 9.0*(E3*E4 + E2*E5)/68.0);
  cplx out = ((long double)3.0)*s + h / (A * std::sqrt(A));
  return toRcplx(out);
}

// [[Rcpp::export]]
Rcomplex Carlson_RJ_(
    Rcomplex xr, Rcomplex yr, Rcomplex zr, Rcomplex pr, double err
) {
  cplx x = fromRcplx(xr);
  cplx y = fromRcplx(yr);
  cplx z = fromRcplx(zr);
  cplx p = fromRcplx(pr);
  cplx A0 = (x + y + z + p + p) / (long double)5.0;
  cplx A = A0;
  cplx delta = (p - x) * (p - y) * (p - z);
  long double M =
    std::max(
      std::max(
        std::max(
          std::abs(A - x),
          std::abs(A - y)
        ),
        std::abs(A-z)
      ),
      std::abs(A-p)
    );
  long double Q = std::pow(4.0/err, 1.0/6.0) * M;
  std::vector<cplx> d(0);
  std::vector<cplx> e(0);
  long double f = 1.0;
  long double fac = 1.0;
  while(std::abs(A) <= Q) {
    cplx srx = std::sqrt(x);
    cplx sry = std::sqrt(y);
    cplx srz = std::sqrt(z);
    cplx srp = std::sqrt(p);
    cplx dnew = (srp + srx) * (srp + sry) * (srp + srz);
    d.push_back(f * dnew);
    e.push_back(fac * delta / (dnew * dnew));
    f *= 4.0;
    fac /= 64.0;
    cplx lambda = srx * sry + sry * srz + srz * srx;
    x = (x + lambda) / (long double)4.0;
    y = (y + lambda) / (long double)4.0;
    z = (z + lambda) / (long double)4.0;
    p = (p + lambda) / (long double)4.0;
    A = (A + lambda) / (long double)4.0;
    Q /= 4.0;
  }
  cplx fA = f * A;
  cplx X = (A0-x) / fA;
  cplx Y = (A0-y) / fA;
  cplx Z = (A0-z) / fA;
  cplx P = -(X + Y + Z) / (long double)2.0;
  cplx P2 = P * P;
  cplx E2 = X * Y + Y * Z + Z * X - ((long double)3.0)*P2;
  cplx E3 = X * Y * Z + (((long double)2.0)*E2 + ((long double)4.0)*P2) * P;
  cplx E4 = (((long double)2.0)*X*Y*Z + (E2 + ((long double)3.0)*P2) * P) * P;
  cplx E5 = X * Y * Z * P2;
  cplx h =
    ((long double)1.0 - (long double)3.0*E2/(long double)14.0
       + E3/(long double)6.0 + (long double)9.0*E2*E2/(long double)88.0
       - (long double)3.0*E4/(long double)22.0
       - (long double)9.0*E2*E3/(long double)52.0
       + (long double)3.0*E5/(long double)26.0) / f;
  cplx s(0.0, 0.0);
  cplx one(1.0, 0.0);
  int n = e.size();
  for(int i = 0; i < n; i++) {
    cplx srei = std::sqrt(e[i]);
    cplx a = srei == (long double)0.0 ? one : std::atan(srei)/srei;
    s += a / d[i];
  }
  cplx out = h / (A * std::sqrt(A)) + (long double)6.0 * s;
  return toRcplx(out);
}
