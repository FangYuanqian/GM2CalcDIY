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

#include "GeneralTHDM/gm2_2loop_helpers.hpp"
#include "gm2_ffunctions.hpp"

/**
 * \file gm2_2loop_B.cpp
 *
 * Contains functions necessary to calculate the bosonic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace general_thdm {

namespace {

const double pi = 3.1415926535897932;

double sqr(double x) noexcept { return x*x; }

/// Eq.(102), arxiv:1607.06292
double YF1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   // @todo(alex) avoid re-calculation of common sub-expressions
   return
      - 72*cw2*(-1 + cw2)*(u + 2*w)/u
      - 36*cw2*(-1 + cw2)*(u + 2*w)/u*std::log(w)
      + 9*(-8*cw4 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*std::log(u)
      - 9*(3 - 10*cw2 + 8*cw4)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi(w,w,1)
      + 9*(8*cw4 + 3*u - 2*cw2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u*u)*Phi(w,w,w)
      ;
}

/// Eq.(99), arxiv:1607.06292
double fb(double u, double w, double al, double cw2) noexcept
{
   return al*pi/(cw2*(-1.0 + cw2))*(u + 2*w);
}

/// Eq.(100), arxiv:1607.06292
double Fm0(double u, double w, double al, double cw2, double mm2, double mz2) noexcept
{
   const auto al2 = al*al;
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;

   // @todo(alex) cancel out al*pi
   return al2/(576*sqr(pi)*cw4*sw4) * mm2/mz2
      * 1.0/(al*pi) * YF1(u,w,cw2);
}

} // anonymous namespace

/**
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607:06292
 */
double amu2L_B_EWadd(double eta, double zetal)
{
   return 2.3e-11 * eta * zetal;
}

/**
 * Calculates 2-loop bosonic non-Yukawa contributions.
 *
 * Eq (71), arxiv:1607:06292
 */
double amu2L_B_nonYuk()
{
   // @todo(alex) implementation missing
   return 0.0;
}

/**
 * Calculates 2-loop bosonic Yukawa contributions.
 *
 * Eq (52), arxiv:1607:06292
 */
double amu2L_B_Yuk()
{
   // @todo(alex) pass parameters to function
   const auto tb = 1.0;
   const auto zetal = 0.0;
   const auto lambda5 = 0.0;
   const auto eta = 0.0;
   const auto al = 0.0;
   const auto mhSM2 = 0.0;
   const auto mHp2 = 0.0;
   const auto mw2 = 0.0;
   const auto mz2 = 0.0;
   const auto mm2 = 0.0;

   const auto cw2 = mw2/mz2;
   const auto xhSM = mhSM2/mz2;
   const auto xHp = mHp2/mz2;

   const auto a000 = fb(xhSM, xHp, al, cw2)*Fm0(xhSM, xHp, al, cw2, mm2, mz2);
   const auto a0z0 = 0.0; // @todo(alex) implementation missing
   const auto a500 = 0.0; // @todo(alex) implementation missing
   const auto a5z0 = 0.0; // @todo(alex) implementation missing
   const auto a001 = 0.0; // @todo(alex) implementation missing
   const auto a0z1 = 0.0; // @todo(alex) implementation missing
   const auto a501 = 0.0; // @todo(alex) implementation missing
   const auto a5z1 = 0.0; // @todo(alex) implementation missing

   const double res =
      + a000
      + a0z0*(tb - 1.0/tb)*zetal
      + a500*lambda5
      + a5z0*(tb - 1.0/tb)*lambda5*zetal
      + (
         + a001*(tb - 1.0/tb)
         + a0z1*zetal
         + a501*(tb - 1.0/tb)*lambda5
         + a5z1*lambda5*zetal
         )*eta;

   return res;
}

} // namespace general_thdm

} // namespace gm2calc
