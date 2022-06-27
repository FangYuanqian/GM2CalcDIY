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
#include <cmath>
#include <complex>
#include <limits>
#include <tuple>

namespace gm2calc {

namespace {
   constexpr double eps = 10.0*std::numeric_limits<double>::epsilon();
   const double qdrt_eps = std::pow(eps, 0.25);

   /// returns number squared
   template <typename T> T sqr(T x) noexcept { return x*x; }

   /// returns number cubed
   template <typename T> T pow3(T x) noexcept { return x*x*x; }

   /// returns number to the power 4
   template <typename T> T pow4(T x) noexcept { return sqr(sqr(x)); }

   bool is_zero(double a, double prec) noexcept
   {
      return std::fabs(a) < prec;
   }

   bool is_equal(double a, double b, double prec) noexcept
   {
      const double max = std::max(std::abs(a), std::abs(b));
      return is_zero(a - b, prec*(1.0 + max));
   }

   void sort(double& x, double& y) noexcept
   {
      if (x > y) { std::swap(x, y); }
   }

   void sort(double& x, double& y, double& z) noexcept
   {
      if (x > y) { std::swap(x, y); }
      if (y > z) { std::swap(y, z); }
      if (x > y) { std::swap(x, y); }
   }

   /// lambda^2(u,v)
   double lambda_2(double u, double v) noexcept
   {
      return sqr(1 - u - v) - 4*u*v;
   }

   /// expansion of (1 - lambda + u - v)/2 for u ~ v ~ 0 up to including O(u^3 v^3)
   double l00(double u, double v) noexcept
   {
      return v*(1 + u*(1 + u*(1 + u)) + v*(u*(1 + u*(3 + 6*u)) + u*(1 + u*(6 + 20*u))*v));
   }

   /// expansion of (1 - lambda + u - v)/2 for u ~ 0 and v < 1 up to including O(u^3 v^3)
   double l0v(double u, double v) noexcept
   {
      const double a = 1 - v;
      const double a2 = a*a;
      const double a3 = a2*a;
      return u*(0.5*(1 + (1 + v)/a) + u*(v + u*v*(1 + v)/a2)/a3);
   }

   /// expansion of (1 - lambda - u + v)/2 for u ~ 0 and v < 1 up to including O(u^3 v^3)
   double lv0(double u, double v) noexcept
   {
      const double a = 1 - v;
      const double a2 = a*a;
      const double a3 = a2*a;
      return v + u*(0.5*(-1 + (1 + v)/a) + u*(v + u*v*(1 + v)/a2)/a3);
   }

   /// returns tuple (0.5*(1 - lambda + u - v), 0.5*(1 - lambda - u + v))
   std::tuple<double,double> luv(double lambda, double u, double v) noexcept
   {
      if (v < qdrt_eps) {
         return std::make_tuple(l00(u, v), l00(v, u));
      } else if (u < qdrt_eps) {
         return std::make_tuple(l0v(u, v), lv0(u, v));
      }
      return std::make_tuple(0.5*(1 - lambda + u - v),
                             0.5*(1 - lambda - u + v));
   }

   /// u < 1 && v < 1, lambda^2(u,v) > 0; note: phi_pos(u,v) = phi_pos(v,u)
   double phi_pos(double u, double v) noexcept
   {
      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const double pi23 = 3.2898681336964529; // Pi^2/3
      const auto lambda = std::sqrt(lambda_2(u,v));

      if (is_equal(u, v, eps)) {
         const double x = u < qdrt_eps ? u*(1 + u*(1 + u*(2 + 5*u))) : 0.5*(1 - lambda);

         return (- sqr(std::log(u)) + 2*sqr(std::log(x))
                 - 4*dilog(x) + pi23)/lambda;
      }

      double x = 0, y = 0;
      std::tie(x, y) = luv(lambda, u, v);

      return (- std::log(u)*std::log(v) + 2*std::log(x)*std::log(y)
              - 2*dilog(x) - 2*dilog(y) + pi23)/lambda;
   }

   /// clausen_2(2*acos(x))
   double cl2acos(double x) noexcept
   {
      return clausen_2(2*std::acos(x));
   }

   /// lambda^2(u,v) < 0, u = 1
   double phi_neg_1v(double v) noexcept
   {
      return 2*(cl2acos(1 - 0.5*v) + 2*cl2acos(0.5*std::sqrt(v)));
   }

   /// lambda^2(u,v) < 0; note: phi_neg(u,v) = phi_neg(v,u)
   double phi_neg(double u, double v) noexcept
   {
      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         // -I/9 (Pi^2 - 36 PolyLog[2, (1 - I Sqrt[3])/2])/Sqrt[3]
         return 2.343907238689459;
      }

      const auto lambda = std::sqrt(-lambda_2(u,v));

      if (is_equal(u, v, eps)) {
         return 4*clausen_2(2*std::asin(std::sqrt(0.25/u)))/lambda;
      }

      if (is_equal(u, 1.0, eps)) {
         return phi_neg_1v(v)/lambda;
      }

      if (is_equal(v, 1.0, eps)) {
         return phi_neg_1v(u)/lambda;
      }

      const auto sqrtu = std::sqrt(u);
      const auto sqrtv = std::sqrt(v);

      return 2*(+ cl2acos(0.5*(1 + u - v)/sqrtu)
                + cl2acos(0.5*(1 - u + v)/sqrtv)
                + cl2acos(0.5*(-1 + u + v)/(sqrtu*sqrtv)))/lambda;
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

      if (is_zero(lambda, eps)) {
         // phi_uv is always multiplied by lambda.  So, in order to
         // avoid nans if lambda == 0, we simply return 0
         return 0.0;
      }

      if (lambda > 0.) {
         if (u <= 1 && v <= 1) {
            return phi_pos(u,v);
         }
         const auto vou = v/u;
         if (u >= 1 && vou <= 1) {
            const auto oou = 1/u;
            return phi_pos(oou,vou)*oou;
         }
         // v >= 1 && u/v <= 1
         const auto oov = 1/v;
         return phi_pos(oov,1/vou)*oov;
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

   return 4.0/105.0/pow4(d) * (
      + (1.0 - x) * (-97.0 * x2 - 529.0 * x + 2.0)
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
   if (x < 0 || y < 0) {
      ERROR("Fb: x and y must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   if (is_zero(y, eps)) {
      return 0;
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
   if (x < 0 || y < 0) {
      ERROR("Fa: x and y must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   if (is_zero(y, eps)) {
      return 0;
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

/// Ixy(0,y), squared arguments, y != 0
double I0y(double y) noexcept {
   if (is_equal(y, 1, eps)) {
      const double d = y - 1;
      return 1 + d*(-0.5 + 1./3*d);
   }

   return std::log(y)/(y - 1);
}

/// I(x,y), squared arguments, x == 1, y != 0
double I1y(double x, double y) noexcept {
   const double dy = y - 1;
   const double dy2 = sqr(dy);
   const double dx = (x - 1)/dy2;
   const double y2 = sqr(y);
   const double yly = y*std::log(y);

   return (1 - y + yly)/dy2
      + dx*(0.5 - 0.5*y2 + yly)/dy
      + sqr(dx)*(1./3 + 0.5*y + yly + y2*(1./6*y - 1));
}

/// I(x,y), squared arguments, x == y, x != 0, y != 0
double Ixx(double x, double y) noexcept {
   const double eps_eq = 0.001;

   if (is_equal(y, 1, eps_eq)) {
      const double dx = x - 1;
      const double dy = y - 1;
      const double dy2 = sqr(dy);

      return 0.5 + dx*(-1./6 + 1./12*dy - 1./20*dy2)
         + sqr(dx)*(1./12 - 1./20*dy + 1./30*dy2)
         - 1./6*dy + 1./12*dy2;
   }

   const double y2 = sqr(y);
   const double dy = y - 1;
   const double dy2 = sqr(dy);
   const double dxy = (x - y)/dy2;
   const double ly = std::log(y);

   return (dy - ly)/dy2
      + dxy*(0.5 - 0.5*y2 + y*ly)/(dy*y)
      + sqr(dxy)*(1./6 - y + y2*(0.5 + 1./3*y - ly))/y2;
}

/// I(x,y), x < y, x and y are squared arguments
double Ixy(double x, double y) noexcept {
   const double eps_eq = 0.001;

   if (is_zero(y, eps)) {
      return 0;
   }

   if (is_zero(x, eps)) {
      return I0y(y);
   }

   if (is_equal(x/y, 1, eps_eq)) {
      return Ixx(x, y);
   }

   if (is_equal(x, 1, eps_eq)) {
      return I1y(x, y);
   }

   if (is_equal(y, 1, eps_eq)) {
      return I1y(y, x);
   }

   const double lx = std::log(x);
   const double ly = std::log(y);

   return (x*(y - 1)*lx - y*(x - 1)*ly)/((x - 1)*(x - y)*(y - 1));
}

/// I(x,y,z), x, y and z are squared arguments
double Ixyz(double x, double y, double z) noexcept {
   sort(x, y, z);

   if (is_zero(z, eps)) {
      return 0;
   }

   return Ixy(x/z, y/z)/z;
}

} // anonymous namespace

double Iabc(double a, double b, double c) noexcept {
   return Ixyz(sqr(a), sqr(b), sqr(c));
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 * @author Alexander Voigt
 */
double f_PS(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_PS: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   } else if (z < 0.25) {
      const double y = std::sqrt(1 - 4*z); // 0 < y < 1
      const double c = -9.8696044010893586; // -Pi^2
      const double q = (1 + y)/(1 - y);
      const double lq = std::log(q);
      return z/y*(4*dilog(1 + q) - lq*(2*std::log(z) - lq) + c);
   } else if (z == 0.25) {
      return 1.3862943611198906; // Log[4]
   }

   // z > 0.25
   const double y = std::sqrt(-1 + 4*z);
   const double theta = std::atan2(y, 2*z - 1);
   return 4*z/y*clausen_2(theta);
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

   return (2*z - 1)*f_PS(z) - 2*z*(2 + std::log(z));
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

   return 0.5*z*(2 + std::log(z) - f_PS(z));
}

/**
 * \f$\mathcal{F}_1(\omega)\f$, Eq (25) arxiv:1502.04199
 */
double F1(double w) noexcept {
   if (w < 0.0) {
      ERROR("F1: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.0) {
      return 0.0;
   } else if (w == 0.25) {
      return -0.5;
   }

   return (w - 0.5)*f_PS(w) - w*(2 + std::log(w));
}

/**
 * \f$\tilde{\mathcal{F}}_1(\omega)\f$, Eq (26) arxiv:1502.04199
 */
double F1t(double w) noexcept {
   return 0.5*f_PS(w);
}

/**
 * \f$\mathcal{F}_2(\omega)\f$, Eq (27) arxiv:1502.04199
 */
double F2(double w) noexcept {
   if (w < 0.0) {
      ERROR("F2: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return -0.38629436111989062; // 1 - Log[4]
   }

   return 1 + 0.5*(std::log(w) - f_PS(w));
}

/**
 * \f$\mathcal{F}_3(\omega)\f$, Eq (28) arxiv:1502.04199
 * @author Alexander Voigt
 */
double F3(double w) noexcept {
   if (w < 0.0) {
      ERROR("F3: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return 19.0/4.0;
   }

   return (0.5 + 7.5*w)*(2 + std::log(w)) + (4.25 - 7.5*w)*f_PS(w);
}

/**
 * Källén lambda function \f$\lambda^2(x,y,z) = x^2 + y^2 + z^2 - 2xy - 2yz - 2xz\f$.
 * The arguments u and v are interpreted as squared masses.
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\lambda^2(x,y,z)\f$
 */
double lambda_2(double x, double y, double z) noexcept
{
   return z*z*lambda_2(x/z, y/z);
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
