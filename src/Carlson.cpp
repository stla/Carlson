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


Rcomplex ellE(Rcomplex phi, Rcomplex m, double err) {
  Rcomplex out;
  if(Rcpp::ComplexVector::is_na(phi) || Rcpp::ComplexVector::is_na(m)) {
    out.r = NA_REAL;
    out.i = NA_REAL;
  } else if(phi.r == 0.0 && phi.i == 0.0) {
    out.r = 0.0;
    out.i = 0.0;
  } else if(phi.r >= -M_PI_2 && phi.r <= M_PI_2) {
    if(m.r == 0.0 && m.i == 0.0) {
      out = phi;
    } else if(m.r == 1.0 && m.i == 0.0) {
      out = toRcplx(std::sin(fromRcplx(phi)));
    } else {
      cplx sine = std::sin(fromRcplx(phi));
      if(std::isinf(sine.real()) || std::isinf(sine.imag())) {
        Rcpp::stop("`sin(phi)` is not finite.");
      }
      Rcomplex sine2 = toRcplx(sine * sine);
      Rcomplex one {1.0, 0.0};
      Rcomplex cosine2 = one - sine2;
      Rcomplex oneminusmsine2 = one - m*sine2;
      Rcomplex one_3 {1.0 / 3.0, 0.0};
      out = toRcplx(sine) * (Carlson_RF_(cosine2, oneminusmsine2, one, err) -
        one_3 * m * sine2 * Carlson_RD_(cosine2, oneminusmsine2, one, err));
    }
  } else {
    double k = phi.r > M_PI_2 ?
      ceil(phi.r/M_PI - 0.5) : -floor(0.5 - phi.r/M_PI);
    Rcomplex kpi {M_PI * k, 0.0};
    phi = phi - kpi;
    Rcomplex ktimes2 {2.0 * k, 0.0};
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ktimes2 * ellE(PI_2, m, err) + ellE(phi, m, err);
  }
  return out;
}

//[[Rcpp::export]]
Rcpp::ComplexVector ellEcpp(
    Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err
) {
  int n = phi_.size();
  Rcpp::ComplexVector out(n);
  for(int i = 0; i < n; i++) {
    out(i) = ellE(phi_(i), m_(i), err);
  }
  return out;
}


/*
Rcpp::ComplexVector ellEcpp(
    Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err
) {
  int n = phi_.size();
  Rcpp::ComplexVector out(n);
  Rcpp::LogicalVector phi_NA = Rcpp::is_na(phi_);
  Rcpp::LogicalVector m_NA   = Rcpp::is_na(m_);
  for(int j = 0; j < n; j++) {
    Rcomplex outj;
    Rcomplex phi = phi_(j);
    Rcomplex m   = m_(j);
    if(phi_NA(j) || m_NA(j)) {
      outj = Rcpp::ComplexVector::get_na();
    } else if(phi.r == 0.0 && phi.i == 0.0) {
      outj.r = 0.0;
      outj.i = 0.0;
    } else if(phi.r >= -M_PI_2 && phi.r <= M_PI_2) {
      if(m.r == 0.0 && m.i == 0.0) {
        outj = phi;
      } else if(m.r == 1.0 && m.i == 0.0) {
        outj = toRcplx(std::sin(fromRcplx(phi)));
      } else {
        cplx sine = std::sin(fromRcplx(phi));
        if(std::isinf(sine.real()) || std::isinf(sine.imag())) {
          Rcpp::stop("`sin(phi)` is not finite.");
        }
        Rcomplex sine2 = toRcplx(sine * sine);
        Rcomplex one;
        one.r = 1.0;
        one.i = 0.0;
        Rcomplex cosine2 = one - sine2;
        Rcomplex oneminusmsine2 = one - m*sine2;
        Rcomplex one_3;
        one_3.r = 1.0 / 3.0;
        one_3.i = 0.0;
        outj = toRcplx(sine) * (Carlson_RF_(cosine2, oneminusmsine2, one, err) -
          one_3 * m * sine2 * Carlson_RD_(cosine2, oneminusmsine2, one, err));
      }
    } else {
      double k = phi.r > M_PI_2 ?
        ceil(phi.r/M_PI - 0.5) : -floor(0.5 - phi.r/M_PI);
      Rcomplex kpi;
      kpi.r = M_PI * k;
      kpi.i = 0.0;
      phi = phi - kpi;
      Rcomplex ktimes2;
      ktimes2.r = 2.0 * k;
      ktimes2.i = 0.0;
      Rcomplex PI_2;
      PI_2.r = M_PI_2;
      PI_2.i = 0.0;
      Rcpp::ComplexVector PI_2vec = Rcpp::ComplexVector::create(PI_2);
      Rcpp::ComplexVector mvec    = Rcpp::ComplexVector::create(m);
      Rcpp::ComplexVector phivec  = Rcpp::ComplexVector::create(phi);
      Rcpp::ComplexVector E1 = ellEcpp(PI_2vec, mvec, err);
      Rcpp::ComplexVector E2 = ellEcpp(phivec, mvec, err);
      outj = ktimes2 * E1(0) + E2(0);
    }
    out(j) = outj;
  }
  return out;
}
*/

Rcomplex ellF(Rcomplex phi, Rcomplex m, double err) {
  Rcomplex out;
  if(Rcpp::ComplexVector::is_na(phi) || Rcpp::ComplexVector::is_na(m)) {
    out.r = NA_REAL;
    out.i = NA_REAL;
  } else if(
      (phi.r == 0.0 && phi.i == 0.0) || std::isinf(m.r) || std::isinf(m.i)
    ) {
    out.r = 0.0;
    out.i = 0.0;
  } else if(
      phi.r == 0.0 && std::isinf(phi.i) && m.i == 0.0 && m.r > 0 && m.r < 1
  ) {
    Rcomplex PI_2 {M_PI_2, 0.0};
    Rcomplex minv {1.0 / m.r, 0.0};
    Rcomplex msqrt {std::sqrt(m.r), 0.0};
    out = ellF(PI_2, m, err) - ellF(PI_2, minv, err) / msqrt;
    if(phi.i < 0) {
      out.r = -out.r;
      out.i = -out.i;
    }
  } else if((phi.r == M_PI_2 || phi.r == -M_PI_2) && m.r == 1.0 && m.i == 0.0) {
    out.r = NAN;
    out.i = NAN;
  } else if(phi.r >= -M_PI_2 && phi.r <= M_PI_2) {
    if(m.r == 1.0 && m.i == 0.0 && std::fabs(phi.r) < M_PI_2) {
      out = toRcplx(std::asinh(std::tan(fromRcplx(phi))));
    } else if(m.r == 0.0 && m.i == 0.0) {
      out = phi;
    } else {
      cplx sine = std::sin(fromRcplx(phi));
      if(std::isinf(sine.real()) || std::isinf(sine.imag())) {
        Rcpp::stop("`sin(phi)` is not finite.");
      }
      Rcomplex sine2 = toRcplx(sine * sine);
      Rcomplex one {1.0, 0.0};
      Rcomplex cosine2 = one - sine2;
      Rcomplex oneminusmsine2 = one - m*sine2;
      out = toRcplx(sine) * Carlson_RF_(cosine2, oneminusmsine2, one, err);
    }
  } else {
    double k = phi.r > M_PI_2 ?
      ceil(phi.r/M_PI - 0.5) : -floor(0.5 - phi.r/M_PI);
    Rcomplex kpi {M_PI * k, 0.0};
    phi = phi - kpi;
    Rcomplex ktimes2 {2.0 * k, 0.0};
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ktimes2 * ellF(PI_2, m, err) + ellF(phi, m, err);
  }
  return out;
}

//[[Rcpp::export]]
Rcpp::ComplexVector ellFcpp(
    Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err
) {
  int n = phi_.size();
  Rcpp::ComplexVector out(n);
  for(int i = 0; i < n; i++) {
    out(i) = ellF(phi_(i), m_(i), err);
  }
  return out;
}


Rcomplex ellZ(Rcomplex phi, Rcomplex m, double err) {
  Rcomplex out;
  if(Rcpp::ComplexVector::is_na(phi) || Rcpp::ComplexVector::is_na(m)) {
    out.r = NA_REAL;
    out.i = NA_REAL;
  } else if(std::isinf(m.r) && m.i == 0.0) {
    out.r = NAN;
    out.i = NAN;
  } else if(m.r == 1.0 && m.i == 0.0) {
    if(std::fabs(phi.r) <= M_PI_2) {
      out = toRcplx(std::sin(fromRcplx(phi)));
    } else {
      double k = phi.r > M_PI_2 ?
      ceil(phi.r/M_PI - 0.5) : -floor(0.5 - phi.r/M_PI);
      Rcomplex kpi {M_PI * k, 0.0};
      phi = phi - kpi;
      out = toRcplx(std::sin(fromRcplx(phi)));
    }
  } else {
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ellE(phi, m, err) -
      ellE(PI_2, m, err) / ellF(PI_2, m, err) * ellF(phi, m, err);
  }
  return out;
}

//[[Rcpp::export]]
Rcpp::ComplexVector ellZcpp(
    Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err
) {
  int n = phi_.size();
  Rcpp::ComplexVector out(n);
  for(int i = 0; i < n; i++) {
    out(i) = ellZ(phi_(i), m_(i), err);
  }
  return out;
}


Rcomplex ellPI(Rcomplex phi, Rcomplex n, Rcomplex m, double err) {
  Rcomplex out;
  bool complete = phi.r == M_PI_2 && phi.i == 0.0;
  if(
      Rcpp::ComplexVector::is_na(phi) || Rcpp::ComplexVector::is_na(n)
      || Rcpp::ComplexVector::is_na(m)
  ) {
    out.r = NA_REAL;
    out.i = NA_REAL;
  } else if(
      (phi.r == 0.0 && phi.i == 0.0) ||
        (std::isinf(n.r) && n.i == 0.0) ||
        (std::isinf(m.r) && m.i == 0.0)
  ) {
    out.r = 0.0;
    out.i = 0.0;
  } else if(complete && m.r == 1.0 && m.i == 0.0 && n.r != 1.0 && n.i == 0.0) {
    out.r = n.r > 1.0 ? -std::numeric_limits<double>::infinity() :
                        std::numeric_limits<double>::infinity();
    out.i = 0.0;
  } else if(complete && n.r == 1.0 && n.i == 0.0) {
    out.r = NAN;
    out.i = NAN;
  } else if(complete && m.r == 0.0 && m.i == 0.0) {
    cplx PI_2(M_PI_2, 0.0);
    cplx one(1.0, 0.0);
    out = toRcplx(PI_2 / std::sqrt(one - fromRcplx(n)));
  } else if(complete && n == m) {
    Rcomplex one {1.0, 0.0};
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ellE(PI_2, m, err) / (one - m);
  } else if(complete && n.r == 0.0 && n.i == 0.0) {
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ellF(PI_2, m, err);
  } else if(phi.r >= -M_PI_2 && phi.r <= M_PI_2) {
    cplx sine = std::sin(fromRcplx(phi));
    if(std::isinf(sine.real()) || std::isinf(sine.imag())) {
      Rcpp::stop("`sin(phi)` is not finite.");
    }
    Rcomplex sine2 = toRcplx(sine * sine);
    Rcomplex one {1.0, 0.0};
    Rcomplex cosine2 = one - sine2;
    Rcomplex oneminusmsine2 = one - m*sine2;
    Rcomplex oneminusnsine2 = one - n*sine2;
    Rcomplex one_3 {1.0 / 3.0, 0.0};
    out = toRcplx(sine) * (Carlson_RF_(cosine2, oneminusmsine2, one, err) +
        one_3 * n * sine2 *
        Carlson_RJ_(cosine2, oneminusmsine2, one, oneminusnsine2, err));
  } else {
    double k = phi.r > M_PI_2 ?
      ceil(phi.r/M_PI - 0.5) : -floor(0.5 - phi.r/M_PI);
    Rcomplex kpi {M_PI * k, 0.0};
    phi = phi - kpi;
    Rcomplex ktimes2 {2.0 * k, 0.0};
    Rcomplex PI_2 {M_PI_2, 0.0};
    out = ktimes2 * ellPI(PI_2, n, m, err) + ellPI(phi, n, m, err);
  }
  return out;
}

//[[Rcpp::export]]
Rcpp::ComplexVector ellPIcpp(
    Rcpp::ComplexVector phi_, Rcpp::ComplexVector n_,
    Rcpp::ComplexVector m_, double err
) {
  int n = phi_.size();
  Rcpp::ComplexVector out(n);
  for(int i = 0; i < n; i++) {
    out(i) = ellPI(phi_(i), n_(i), m_(i), err);
  }
  return out;
}
