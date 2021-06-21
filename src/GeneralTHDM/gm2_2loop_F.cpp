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
#include "gm2_dilog.hpp"
#include "gm2_ffunctions.hpp"

/**
 * \file gm2_2loop_F.cpp
 *
 * Contains functions necessary to calculate the fermionic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace {

struct THDM_parameters {
   double alpha{}; ///< alpha_em
   double mm{};    ///< muon mass
   double mw{};    ///< W boson mass
   double mz{};    ///< Z boson mass
};

struct F_neut_pars {
   double mf2;   /// squared mass of fermion f
   double qf;    ///< electromagnetic charge of fermion f
   double ql;    ///< electromagnetic charge of fermion l
   double t3f;   ///< SU(2)_L charge of ferimon f
   double t3l;   ///< SU(2)_L charge of ferimon l
   double nc;    ///< number of colors of fermion f
};

const double pi = 3.1415926535897932;
const double q_u =  2.0/3.0;  ///< electric charge of up-type quarks
const double q_d = -1.0/3.0;  ///< electric charge of down-type quarks
const double q_l = -1.0;      ///< electric charge of down-type leptons
const double t3_u = +1.0;     ///< SU(2)_L charge of up-type quark @todo(alex) check convention
const double t3_d = -1.0;     ///< SU(2)_L charge of down-type quark @todo(alex) check convention
const double t3_l = -1.0;     ///< SU(2)_L charge of charged lepton @todo(alex) check convention

double sqr(double x) noexcept { return x*x; }

/// Eq (56), arxiv:1607.06292, S = h or H
double FS(double ms2, double mf2)
{
   return -2.0 + std::log(ms2/mf2)
      - (ms2 - 2*mf2)/ms2 * Phi(ms2, mf2, mf2)/(ms2 - 4*mf2);
}

/// Eq (57), arxiv:1607.06292, S = A
double FA(double ms2, double mf2)
{
   return Phi(ms2, mf2, mf2)/(ms2 - 4*mf2);
}

/// Eq (54), arxiv:1607.06292, S = h or H
template <typename F>
double fSgamma(double ms2, const F_neut_pars& pars, const THDM_parameters& thdm, F fS) noexcept
{
   const double al2 = sqr(thdm.alpha);
   const double mm2 = sqr(thdm.mm);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double sw2 = 1.0 - mw2/mz2;
   const double mf2 = pars.mf2;
   const double qf2 = sqr(pars.qf);
   const double nc = pars.nc;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * qf2*nc * mf2/ms2 * fS(ms2, mf2);
}

/// Eq (55), arxiv:1607.06292, S = h or H
template <typename F>
double fSZ(double ms2, const F_neut_pars& pars, const THDM_parameters& thdm, F fS) noexcept
{
   const double al2 = sqr(thdm.alpha);
   const double mm2 = sqr(thdm.mm);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double mf2 = pars.mf2;
   const double qf = pars.qf;
   const double ql = pars.ql;
   const double nc = pars.nc;
   const double gvf = 0.5*pars.t3f - qf*sw2;
   const double gvl = 0.5*pars.t3l - ql*sw2;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * (-nc*qf*gvl*gvf)/(sw2*cw2)
      * mf2/(ms2 - mz2) * (fS(ms2, mf2) - fS(mz2, mf2));
}

/// Eq (53), arxiv:1607.06292, S = h or H
template <typename F>
double ffS(double ms2, const F_neut_pars& pars, const THDM_parameters& thdm, F fS) noexcept
{
   return fSgamma(ms2, pars, thdm, fS) + fSZ(ms2, pars, thdm, fS);
}

/// Eq (60), arxiv:1607.06292
double FlHp(double ms2, double mf2) noexcept
{
   const double xl = mf2/ms2;

   return xl + xl*(xl - 1.0)*(dilog(1.0 - 1.0/xl) - sqr(pi)/6)
      + (xl - 0.5)*std::log(xl);
}

/// Eq (61), arxiv:1607.06292
double FdHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double sqrt_xu = std::sqrt(xu);
   const double sqrt_xd = std::sqrt(xd);
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;
   const double s = 0.25*(qu + qd);
   const double c = sqr(xu - xd) - qu*xu + qd*xd;
   const double cbar = (xu - qu)*xu - (xd + qd)*xd;

   return -(xu - xd) + (cbar/y - c*(xu - xd)/y) * Phi(sqrt_xd, sqrt_xu, 1.0)
      + c*(dilog(1.0 - xd/xu) - 0.5*std::log(xu)*std::log(xd/xu)*Phi(sqrt_xd, sqrt_xu, 1.0))
      + (s + xd)*std::log(xd) + (s - xu)*std::log(xu);
}

/// Eq (62), arxiv:1607.06292
double FuHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double sqrt_xu = std::sqrt(xu);
   const double sqrt_xd = std::sqrt(xd);
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;

   return FdHp(ms2, md2, mu2, 2.0 + qd, 2.0 + qu)
      - 4.0/3.0*(xu - xd - 1.0)/y*Phi(sqrt_xd, sqrt_xu, 1.0)
      - 1.0/3.0*(sqr(std::log(xd)) - sqr(std::log(xu)));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, double mf2, const THDM_parameters& thdm) noexcept
{
   const double al2 = sqr(thdm.alpha);
   const double mm2 = sqr(thdm.mm);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double nc = 1.0;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*mf2/(ms2 - mw2)
      * (FlHp(ms2, mf2) - FlHp(mw2, mf2));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = u
double fuHp(double ms2, double md2, double mu2, double qd, double qu, const THDM_parameters& thdm) noexcept
{
   const double al2 = sqr(thdm.alpha);
   const double mm2 = sqr(thdm.mm);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double nc = 3.0;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*mu2/(ms2 - mw2)
      * (FuHp(ms2, md2, mu2, qd, qu) - FuHp(mw2, md2, mu2, qd, qu));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = d
double fdHp(double ms2, double md2, double mu2, double qd, double qu, const THDM_parameters& thdm) noexcept
{
   const double al2 = sqr(thdm.alpha);
   const double mm2 = sqr(thdm.mm);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double nc = 3.0;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*md2/(ms2 - mw2)
      * (FdHp(ms2, md2, mu2, qd, qu) - FdHp(mw2, md2, mu2, qd, qu));
}

} // anonymous namespace

/**
 * \fn amu2L_F
 *
 * Calculates 2-loop fermionic contributions.
 *
 * Eq (63), arxiv:1607:06292
 */
double amu2L_F()
{
   THDM_parameters thdm;
   double mu2{}, md2{}, ml2{};
   double mh2{}, mH2{}, mA2{}, mHp2{}, mhSM2{};
   double yuh{}, ydh{}, ylh{};
   double yuH{}, ydH{}, ylH{};
   double yuA{}, ydA{}, ylA{};

   const F_neut_pars pars_u{mu2, q_u, q_l, t3_u, t3_l, 3.0};
   const F_neut_pars pars_d{md2, q_d, q_l, t3_d, t3_l, 3.0};
   const F_neut_pars pars_l{ml2, q_l, q_l, t3_l, t3_l, 1.0};

   const auto lFS = [] (double ms2, double mf2) { return FS(ms2, mf2); };
   const auto lFA = [] (double ms2, double mf2) { return FA(ms2, mf2); };

   double res = 0.0;

   // h
   res += ffS(mh2, pars_u, thdm, lFS)*yuh*ylh;
   res += ffS(mh2, pars_d, thdm, lFS)*ydh*ylh;
   res += ffS(mh2, pars_l, thdm, lFS)*ylh*ylh;

   // H
   res += ffS(mH2, pars_u, thdm, lFS)*yuH*ylH;
   res += ffS(mH2, pars_d, thdm, lFS)*ydH*ylH;
   res += ffS(mH2, pars_l, thdm, lFS)*ylH*ylH;

   // A
   res += ffS(mA2, pars_u, thdm, lFA)*yuA*ylA;
   res += ffS(mA2, pars_d, thdm, lFA)*ydA*ylA;
   res += ffS(mA2, pars_l, thdm, lFA)*ylA*ylA;

   // H^\pm
   res += fuHp(mHp2, md2, mu2, q_d, q_u, thdm)*yuA*ylA;
   res += fdHp(mHp2, md2, mu2, q_d, q_u, thdm)*ydA*ylA;
   res += flHp(mHp2, ml2, thdm)*ylA*ylA;

   // subtract hSM
   res -= ffS(mhSM2, pars_u, thdm, lFS);
   res -= ffS(mhSM2, pars_d, thdm, lFS);
   res -= ffS(mhSM2, pars_l, thdm, lFS);

   return res;
}

} // namespace gm2calc
