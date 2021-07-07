// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_ffunctions.hpp"
#include "gm2_dilog.hpp"
#include "gm2_log.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

#include <boost/numeric/odeint.hpp>

namespace gm2calc {

namespace {
   const double eps = 10.0*std::numeric_limits<double>::epsilon();

   /// returns number squared
   template <typename T> T sqr(T x) noexcept { return x*x; }

   /// returns number cubed
   template <typename T> T pow3(T x) noexcept { return x*x*x; }

   /// returns number to the power 4
   template <typename T> T pow4(T x) noexcept { return sqr(sqr(x)); }

   bool is_zero(double a, double prec)
   {
      return std::fabs(a) < prec;
   }

   bool is_equal(double a, double b, double prec)
   {
      const double max = std::max(std::abs(a), std::abs(b));
      return is_zero(a - b, prec*(1.0 + max));
   }

   void sort(double& x, double& y, double& z) noexcept
   {
      if (x > y) { std::swap(x, y); }
      if (y > z) { std::swap(y, z); }
      if (x > y) { std::swap(x, y); }
   }

   template <typename T>
   double integrate(T fun, double start, double stop, double dt) {
      using boost::numeric::odeint::integrate;
      using state_type = std::array<double, 1>;

      const auto f = [fun](const state_type& /* unused */, state_type& dxdt,
                           double x) { dxdt[0] = fun(x); };

      state_type x0 = {start};
      integrate(f, x0, start, stop, dt);

      return x0[0];
   }

   /// lambda^2(u,v)
   double lambda_2(double u, double v) noexcept
   {
      return sqr(1 - u - v) - 4*u*v;
   }

   /// u < 1 && v < 1, lambda^2(u,v) > 0; note: phi_pos(u,v) = phi_pos(v,u)
   double phi_pos(double u, double v) noexcept
   {
      const double eps = 1.0e-7;

      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const double Pi = 3.141592653589793;
      const auto lambda = std::sqrt(lambda_2(u,v));

      if (is_equal(u, v, eps)) {
         return (-(sqr(std::log(u)))
                 + 2*sqr(std::log((1 - lambda)/2.))
                 - 4*dilog((1 - lambda)/2.)
                 + sqr(Pi)/3.)/lambda;
      }

      return (-(std::log(u)*std::log(v))
              + 2*std::log((1 - lambda + u - v)/2.)*std::log((1 - lambda - u + v)/2.)
              - 2*dilog((1 - lambda + u - v)/2.)
              - 2*dilog((1 - lambda - u + v)/2.)
              + sqr(Pi)/3.)/lambda;
   }

   /// lambda^2(u,v) < 0, u = 1
   double phi_neg_1v(double v, double lambda) noexcept
   {
      return 2*(+ clausen_2(2*std::acos((2 - v)/2))
                + 2*clausen_2(2*std::acos(0.5*std::sqrt(v))))/lambda;
   }

   /// lambda^2(u,v) < 0; note: phi_neg(u,v) = phi_neg(v,u)
   double phi_neg(double u, double v) noexcept
   {
      const double eps = 1.0e-7;

      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const auto lambda = std::sqrt(-lambda_2(u,v));

      if (is_equal(u, 1.0, eps)) {
         return phi_neg_1v(v, lambda);
      }

      if (is_equal(v, 1.0, eps)) {
         return phi_neg_1v(u, lambda);
      }

      if (is_equal(u, v, eps)) {
         return 2*(2*clausen_2(2*std::acos(1/(2.*std::sqrt(u))))
                   + clausen_2(2*std::acos((-1 + 2*u)/(2.*std::abs(u)))))/lambda;
      }

      const auto sqrtu = std::sqrt(u);
      const auto sqrtv = std::sqrt(v);

      return 2*(+ clausen_2(2*std::acos(0.5*(1 + u - v)/sqrtu))
                + clausen_2(2*std::acos(0.5*(1 - u + v)/sqrtv))
                + clausen_2(2*std::acos(0.5*(-1 + u + v)/(sqrtu*sqrtv))))/lambda;
   }

   /**
    * Phi(u,v) with u = x/z, v = y/z.
    *
    * The following identities hold:
    * Phi(u,v) = Phi(v,u) = Phi(1/u,v/u)/u = Phi(1/v,u/v)/v
    */
   double phi_uv(double u, double v) noexcept
   {
      const auto lambda = lambda_2(u,v);

      if (is_zero(lambda, 1e-11)) {
         // phi_uv is always multiplied by lambda.  So, in order to
         // avoid nans if lambda == 0, we simply return 0
         return 0.0;
      }

      if (lambda > 0.) {
         if (u <= 1 && v <= 1) {
            return phi_pos(u,v);
         }
         if (u >= 1 && v/u <= 1) {
            return phi_pos(1./u,v/u)/u;
         }
         // v >= 1 && u/v <= 1
         return phi_pos(1./v,u/v)/v;
      }

      return phi_neg(u,v);
   }

} // anonymous namespace

double F1C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 4.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.6 + d*(0.4 + d*(-2.0/7.0
         + d*(3.0/14.0 + d*(-1.0/6.0
         + 2.0/15.0*d)))));
   }

   return 2.0/pow4(d)*(2.0 + x*(3.0 + 6.0*std::log(x) + x*(-6.0 + x)));
}

double F2C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + d*(-0.75 + d*(0.6 + d*(-0.5 + d*(3.0/7.0
         + d*(-0.375 + 1.0/3.0*d)))));
   }

   return 3.0/(2.0*pow3(1.0 - x))*(-3.0 - 2.0*std::log(x) + x*(4.0 - x));
}

double F3C(double x) noexcept {
   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0
         + d*(1059.0/1175.0
         + d*(-4313.0/3525.0
         + d*(70701.0/57575.0
         + d*(-265541.0/230300.0
         + d*(+48919.0/46060.0
         - 80755.0/82908.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return 4.0/(141.0*pow4(d)) * (
      + (1.0 - x) * (151.0 * x2 - 335.0 * x + 592.0)
      + 6.0 * (21.0 * pow3(x) - 108.0 * x2 - 93.0 * x + 50.0) * lx
      - 54.0 * x * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 108.0 * x * (x2 - 2.0 * x + 12.0) * dilog(1.0 - x)
      );
}

double F4C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0
         + d*(-45.0/122.0
         + d*(941.0/6100.0
         + d*(-17.0/305.0
         + d*(+282.0/74725.0
         + d*(+177.0/6832.0 - 47021.0/1076040.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return -9.0/(122.0 * pow3(1.0 - x)) * (
      + 8.0 * (x2 - 3.0 * x + 2.0)
      + (11.0 * x2 - 40.0 * x + 5.0) * lx
      - 2.0 * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 4.0 * (x2 - 2.0 * x + 9.0) * dilog(1.0 - x)
      );
}

double F1N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 2.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.4 + d*(0.2 + d*(-4.0/35.0
         + d*(1.0/14.0 + d*(-1.0/21.0 + 1.0/30.0*d)))));
   }

   return 2.0/pow4(d)*(1.0 + x*(-6.0 + x*(+3.0 - 6.0 * std::log(x) + 2.0 * x)));
}

double F2N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 3.0;
   }

   if (is_equal(x, 1.0, 0.04)) {
      const double d = x - 1.0;

      return 1. + d*(-0.5 + d*(0.3 + d*(-0.2
         + d*(1.0/7.0 + d*(-3.0/28.0 + 1.0/12.0*d)))));
   }

   return 3.0/pow3(1.0 - x) * (1.0 + x*(2.0 * std::log(x) - x));
}

double F3N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 8.0/105.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(76/875.0 + d*(-431/2625.0 + d*(5858/42875.0
         + d*(-3561/34300.0 + d*(23/294.0 - 4381/73500.0*d)))));
   }

   const double x2 = sqr(x);

   return 4.0/(105.0 * pow4(d)) * (
      + (1.0 - x) * (- 97.0 * x2 - 529.0 * x + 2.0)
      + 6.0 * x2 * (13.0 * x + 81.0) * std::log(x)
      + 108.0 * x * (7.0 * x + 4.0) * dilog(1.0 - x)
      );
}

double F4N(double x) noexcept {
   const double PI2 = 9.8696044010893586; // Pi^2

   if (is_zero(x, eps)) {
      return -3.0/4.0*(-9.0 + PI2);
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + sqr(d)*(-111.0/800.0 + d*(59.0/400.0 + d*(-129.0/980.0
         + d*(177.0/1568.0 - 775.0/8064.0*d))));
   }

   return -2.25/pow3(1.0 - x) * (
      + (x + 3.0) * (x * std::log(x) + x - 1.0)
      + (6.0 * x + 2.0) * dilog(1.0 - x)
      );
}

namespace {

/// Fb(1,1)
double Fb11(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      + 1.0/12.0 + (-0.05 + y1/30.)*y1
      + x1*(-0.05 + (1.0/30.0 - y1/42.0)*y1
      + x1*(1.0/30.0 + (-1.0/42.0 + y1/56.0)*y1));
}

/// Fb(x,1), x != 1, x != 0
double Fb1(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (2.0 + x*(3.0 + 6.0*lx + x*(-6.0 + x)))/(6.0*x14)
      + y1*(3.0 + x*(10.0 + 12.0*lx + x*(-18.0 + x*(6.0 - x))))/(12.0*x15)
      + sqr(y1)*(12.0 + x*(65.0 + 60.0*lx + x*(-120.0 + x*(60.0 + x*(-20.0 + 3.0*x)))))/(60.*x16);
}

/// Fb(x,x), x != 1, x != 0
double Fbx(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-5.0 - 2.0*lx + x*(4.0 - 4.0*lx + x))/(2.0*x14)
      - d*(-1.0 + x*(-9.0 - 6.0*lx + x*(9.0 - 6.0*lx + x)))/(2.0*x15*x)
      - sqr(d)*(-1.0 + x*(12.0 + x*(36.0 + 36.0*lx + x*(-44.0 + 24.0*lx - 3.0*x))))/(6.*x16*sqr(x));
}

} // anonymous namespace

double Fb(double x, double y) noexcept {
   if ((is_zero(x, eps) && is_zero(y, eps)) || is_zero(x, eps) ||
       is_zero(y, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.01) && is_equal(y, 1.0, 0.01)) {
      return Fb11(x, y);
   }

   if (is_equal(x, 1.0, 0.01)) {
      return Fb1(y, x);
   }

   if (is_equal(y, 1.0, 0.01)) {
      return Fb1(x, y);
   }

   if (is_equal(x, y, 0.01)) {
      return Fbx(x, y);
   }

   return - (G4(x) - G4(y)) / (x - y);
}

namespace {

/// Fa(1,1)
double Fa11(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      0.25 + (-0.2 + y1/6.)*y1
      + x1*(-0.2 + (1.0/6.0 - y1/7.0)*y1)
      + sqr(x1)*(1.0/6.0 + (-1.0/7.0 + y1/8.0)*y1);
}

/// Fa(x,1), x != 1, x != 0
double Fa1(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-11.0 - 6.0*lx + x*(18.0 + x*(-9.0 + 2.0*x)))/(6.0*x14)
      + y1*(-25.0 - 12.0*lx + x*(48.0 + x*(-36.0 + x*(16.0 - 3.0*x))))/(12.0*x15)
      + sqr(y1)*(-137.0 - 60.0*lx + x*(300.0 + x*(-300.0 + x*(200.0 + x*(-75.0 + 12.0*x)))))/(60.0*x16);
}

/// Fa(x,x), x != 1, x != 0
double Fax(double x, double y) noexcept {
   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;
   const double x2 = sqr(x);
   const double x3 = x2*x;

   return (2.0 + x*(3.0 + 6.0*lx + x*(-6.0 + x)))/(2.0*x14*x)
      - d*(-1.0 + x*(8.0 + x*(12.0*lx + x*(-8.0 + x))))/(2.0*x15*x2)
      - sqr(d)*(-2.0 + x*(15.0 + x*(-60.0 + x*(20.0 - 60.0*lx + x*(30.0 - 3.0*x)))))/(6.0*x16*x3);
}

} // anonymous namespace

double Fa(double x, double y) noexcept {
   if ((is_zero(x, eps) && is_zero(y, eps)) || is_zero(x, eps) ||
       is_zero(y, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.001) && is_equal(y, 1.0, 0.001)) {
      return Fa11(x,y);
   }

   if (is_equal(x, 1.0, 0.001)) {
      return Fa1(y,x);
   }

   if (is_equal(y, 1.0, 0.001)) {
      return Fa1(x,y);
   }

   if (is_equal(x, y, 0.001)) {
      return Fax(x,y);
   }

   return - (G3(x) - G3(y)) / (x - y);
}

double G3(double x) noexcept {
   if (is_equal(x, 1.0, 0.01)) {
      const double d = x - 1.0;
      return 1.0/3.0 + d*(-0.25 + d*(0.2 + (-1.0/6.0 + d/7.0)*d));
   }

   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x - 3.0) + 2.0*std::log(x));
}

double G4(double x) noexcept {
   if (is_equal(x, 1.0, 0.01)) {
      const double d = x - 1.0;
      return 1.0/6.0 + d*(-1.0/12.0 + d*(0.05 + (-1.0/30.0 + d/42.0)*d));
   }

   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x + 1.0) - 2.0*x*std::log(x));
}

namespace {

/// I2abc(a,a,a), squared arguments, a != 0
double I2aaa(double a, double b, double c) noexcept {
   const double ba = b - a;
   const double ca = c - a;
   const double a2 = sqr(a);
   const double a3 = a2*a;

   return 0.5/a + (-ba - ca)/(6.0*a2) + (sqr(ba) + ba*ca + sqr(ca))/(12.0*a3);
}

/// I2abc(a,a,c), squared arguments, a != c
double I2aac(double a, double b, double c) noexcept {
   const double ba = b - a;
   const double ac = a - c;
   const double a2 = sqr(a);
   const double a3 = a2*a;
   const double c2 = sqr(c);
   const double c3 = c2*c;
   const double ac2 = sqr(ac);
   const double ac3 = ac2*ac;
   const double ac4 = ac2*ac2;
   const double lac = std::log(a/c);

   return (ac - c*lac)/ac2
      + ba*(-a2 + c2 + 2*a*c*lac)/(2.0*a*ac3)
      + sqr(ba)*((2*a3 + 3*a2*c - 6*a*c2 + c3 - 6*a2*c*lac)/(6.*a2*ac4));
}

/// I2abc(a,a,0), squared arguments, a != 0
double I2aa0(double a, double b) noexcept {
   const double a2 = sqr(a);
   const double a3 = a2*a;
   const double ba = b - a;
   const double ba2 = sqr(ba);

   return 1.0/a - ba/(2.0*a2) + ba2/(3.0*a3);
}

/// I2abc(0,b,c), squared arguments, b != c
double I20bc(double b, double c) noexcept {
   return std::log(b/c)/(b - c);
}

} // anonymous namespace

double Iabc(double a, double b, double c) noexcept {
   if ((is_zero(a, eps) && is_zero(b, eps) && is_zero(c, eps)) ||
       (is_zero(a, eps) && is_zero(b, eps)) ||
       (is_zero(a, eps) && is_zero(c, eps)) ||
       (is_zero(b, eps) && is_zero(c, eps))) {
      return 0.0;
   }

   const double a2 = sqr(a);
   const double b2 = sqr(b);
   const double c2 = sqr(c);
   const double eps_eq = 0.001;

   if (is_equal(a2, b2, eps_eq) && is_equal(a2, c2, eps_eq)) {
      return I2aaa(a2, b2, c2);
   }

   if (is_equal(a2, b2, eps_eq)) {
      if (is_zero(c, eps)) {
         return I2aa0(a2, b2);
      }
      return I2aac(a2, b2, c2);
   }

   if (is_equal(b2, c2, eps_eq)) {
      if (is_zero(a, eps)) {
         return I2aa0(b2, c2);
      }
      return I2aac(b2, c2, a2);
   }

   if (is_equal(a2, c2, eps_eq)) {
      if (is_zero(b, eps)) {
         return I2aa0(a2, c2);
      }
      return I2aac(a2, c2, b2);
   }

   if (is_zero(a, eps)) {
      return I20bc(b2, c2);
   }

   if (is_zero(b, eps)) {
      return I20bc(c2, a2);
   }

   if (is_zero(c, eps)) {
      return I20bc(a2, b2);
   }

   return (+ a2 * b2 * std::log(a2/b2)
           + b2 * c2 * std::log(b2/c2)
           + c2 * a2 * std::log(c2/a2))
           / ((a2 - b2) * (b2 - c2) * (a2 - c2));
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 */
double f_PS(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_PS: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   } else if (z < 0.25) {
      const double y = std::sqrt(1. - 4. * z);
      return 2.0*z/y*(dilog(1.0 - 0.5*(1.0 - y)/z) - dilog(1.0 - 0.5*(1.0 + y)/z));
   } else if (z == 0.25) {
      return 1.3862943611198906; // Log[4]
   }

   // z > 0.25
   const std::complex<double> y = std::sqrt(std::complex<double>(1.0 - 4.0*z, 0.0));
   return std::real(2.0*z/y*(dilog(1.0 - 0.5*(1.0 - y)/z) - dilog(1.0 - 0.5*(1.0 + y)/z)));
}

/**
 * Calculates \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
 */
double f_S(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_S: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   }

   return (2.0*z - 1.0)*f_PS(z) - 2.0*z*(2.0 + std::log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_sferm: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   }

   return 0.5*z*(2.0 + std::log(z) - f_PS(z));
}

double F1(double w) noexcept {
   if (w < 0.0) {
      ERROR("F1: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.0) {
      return 0.0;
   } else if (w == 0.25) {
      return -0.5;
   }

   const auto y = std::sqrt(std::complex<double>(1 - 4*w));
   const auto lm1my = std::log(-1.0 - y);
   const auto l1my = std::log(1.0 - y);
   const auto lw = std::log(w);
   const auto l16 = 2.7725887222397812; // 4 Log[2]

   const auto res =
      - 2.0*y - y*lw + w*l16*l1my + 2.0*w*lw*l1my
      - 2.0*w*l1my*(l1my + std::log(1.0 + y))
      + (l1my - (1.0 - 2.0*w)*lm1my)*std::log((1.0 - y)*(1.0 + y)/(4.0*w))
      + (1.0 - 2.0*w)*(dilog((1.0 + y)/(-1.0 + y)) - dilog((-1.0 + y)/(1.0 + y)))
      ;

   return std::real(w/y*res);
}

double F1t(double w) noexcept {
   return 0.5*f_PS(w);
}

double F2(double w) noexcept {
   if (w < 0.0) {
      ERROR("F2: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return -0.38629436111989062; // 1 - Log[4]
   }

   const auto y = std::sqrt(std::complex<double>(1 - 4*w));
   const auto lm1my = std::log(-1.0 - y);
   const auto l1my = std::log(1.0 - y);
   const auto lw = std::log(w);
   const auto l2 = 0.69314718055994531; // Log[2]
   const auto l3 = 1.0397207708399180; // 3/2 Log[2]

   auto res =
      + l1my*(l1my - lm1my - l3 - lw + 0.5*std::log(1.0 + y))
      + lm1my*(l2 + lw) + (1.0 + 0.5*lw)*y/w
      - (0.5*l1my - lm1my)*std::log(2.0/(1.0 + y))
      + dilog((1.0 + y)/(-1.0 + y)) - dilog((-1.0 + y)/(1.0 + y))
      ;

   return std::real(w/y*res);
}

double F3(double w) noexcept {
   if (w < 0.0) {
      ERROR("F2: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return 19.0/4.0;
   }

   const auto y = std::sqrt(std::complex<double>(1 - 4*w));
   const auto lm1my = std::log(-1.0 - y);
   const auto l1my = std::log(1.0 - y);
   const auto l1py = std::log(1.0 + y);
   const auto lm1py = std::log(-1.0 + y);
   const auto l2o1py = std::log(2.0/(1.0 + y));
   const auto lw = std::log(w);
   const auto l2 = 0.69314718055994531; // Log[2]
   const auto l8 = 3*l2; // 3 Log[2]
   const auto l12 = 12*l2; // 12 Log[2]

   const auto res =
      + l1my*l1my*(-45.0 - 57.0*y + 36.0*w*(4.0 + y))
      + lw*(6.0*(15.0 + 1.0/w)*y + 3.0*(-19.0 + 12.0*w)*(-1.0 + y)*lm1py)
      + 6.0*(
         + 30.0*y + (2.0*y)/w + 19.0*l2*lm1py - 12.0*w*l2*lm1py
         - 19.0*y*l2*lm1py + 12.0*w*y*l2*lm1py
         + (17.0 - 30.0*w)*(dilog((-1.0 + y)/(1.0 + y)) - dilog((1.0 + y)/(-1.0 + y)))
         )
      + lm1my*(
         - 45.0*l2 + (-19.0 + 12.0*w)*y*l8 + 12.0*w*l12
         + (-45.0 - 57.0*y + 36.0*w*(4.0 + y))*lw
         + 3.0*(-15.0 - 19.0*y + 12.0*w*(4.0 + y))*l2o1py
         )
      + l1my*(
         + (lm1my + lw)*(45.0 + 57.0*y - 36.0*w*(4.0 + y))
         + l1py*(63.0 - 57.0*y + 18.0*w*(1.0 + 2.0*y))
         + 39.0*l2 - 198.0*w*l2 + 19.0*y*l8 - 12.0*w*y*l8
         - 3.0*(-19.0 + 12.0*w)*(-1.0 + y)*lm1py
         + (51.0 + 57.0*y - 18.0*w*(5.0 + 2.0*y))*l2o1py
         )
      + l1py*(19.0 - 12.0*w)*(-1.0 + y)*(3.0*lw + l8 + 3.0*lm1py + 3.0*l2o1py)
      ;

   return std::real(w/(12.0*y)*res);
}

double G(double wa, double wb, double x) noexcept {
   return std::log((wa*x+wb*(1-x))/(x*(1-x)))/(x*(1-x)-wa*x-wb*(1-x));
}

double Gn(double wa, double wb, int n) noexcept {
   const auto fun = [wa, wb, n](double x) {
      return std::pow(x,n)*G(wa, wb, x);
   };

   double res = std::numeric_limits<double>::quiet_NaN();

   try {
      integrate(fun, 0.0 + eps, 1.0 - eps, eps);
   } catch (...) {}

   return res;
}

/**
 * \f$\Phi(x,y,z)\f$ function.  The arguments x, y and z are
 * interpreted as squared masses.
 *
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\Phi(x,y,z)\f$
 */
double Phi(double x, double y, double z) noexcept
{
   sort(x, y, z);
   const auto u = x/z, v = y/z;
   return phi_uv(u,v)*z*lambda_2(u, v)/2;
}

} // namespace gm2calc
