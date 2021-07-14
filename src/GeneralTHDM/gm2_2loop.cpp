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

#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/GeneralTHDM.hpp"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"

namespace gm2calc {

/**
 * Calculates 2-loop bosonic contribution to a_mu in the THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop_bosonic(const GeneralTHDM& model)
{
   general_thdm::THDM_B_parameters pars_b;
   pars_b.alpha_em = model.get_alpha_em();
   pars_b.mm = model.get_MFe(1);
   pars_b.mw = model.get_MVWm();
   pars_b.mz = model.get_MVZ();
   pars_b.mhSM = model.get_sm().get_mh();
   pars_b.mA = model.get_MAh(1);
   pars_b.mHp = model.get_MHm(1);
   pars_b.mh = model.get_Mhh();
   pars_b.tb = model.get_tan_beta();
   pars_b.zetal = model.get_zeta_bar_l();
   pars_b.eta = model.get_eta();
   pars_b.lambda5 = model.get_LambdaFive();

   return amu2L_B(pars_b);
}

/**
 * Calculates fermionic 2-loop contribution to a_mu in the THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop_fermionic(const GeneralTHDM& model)
{
   const double sqrt2 = 1.4142135623730950; // Sqrt[2]
   const double zeta_u = model.get_zeta_bar_u();
   const double zeta_d = model.get_zeta_bar_d();
   const double zeta_l = model.get_zeta_bar_l();
   const double alpha_h = model.get_alpha_h();
   const double beta = model.get_beta();
   const double sba = std::sin(beta - alpha_h);
   const double cba = std::cos(beta - alpha_h);
   const Eigen::Matrix<double,3,1> id = Eigen::Matrix<double,3,1>::Ones();

   // Eq.(18), arxiv:1607.06292
   Eigen::Matrix<double,3,4> yuS;
   yuS.col(0) = id*(sba + cba*zeta_u); // S = h
   yuS.col(1) = id*(cba - sba*zeta_u); // S = H
   yuS.col(2) = id*zeta_u;             // S = A
   yuS.col(3) = sqrt2*yuS.col(2);      // S = H^+

   Eigen::Matrix<double,3,4> ydS;
   ydS.col(0) = id*(sba + cba*zeta_d); // S = h
   ydS.col(1) = id*(cba - sba*zeta_d); // S = H
   ydS.col(2) = id*(-zeta_d);          // S = A
   ydS.col(3) = sqrt2*ydS.col(2);      // S = H^+

   Eigen::Matrix<double,3,4> ylS;
   ylS.col(0) = id*(sba + cba*zeta_l); // S = h
   ylS.col(1) = id*(cba - sba*zeta_l); // S = H
   ylS.col(2) = id*(-zeta_l);          // S = A
   ylS.col(3) = sqrt2*ylS.col(2);      // S = H^+

   general_thdm::THDM_F_parameters pars_f;
   pars_f.alpha_em = model.get_alpha_em();
   pars_f.mm = model.get_MFe(1);
   pars_f.mw = model.get_MVWm();
   pars_f.mz = model.get_MVZ();
   pars_f.mhSM = model.get_sm().get_mh();
   pars_f.mA = model.get_MAh(1);
   pars_f.mHp = model.get_MHm(1);
   pars_f.mh = model.get_Mhh();
   pars_f.ml = model.get_MFe();
   pars_f.mu = model.get_MFu();
   pars_f.md = model.get_MFd();
   pars_f.yuS = yuS;
   pars_f.ydS = ydS;
   pars_f.ylS = ylS;

   return amu2L_F(pars_f);
}

/**
 * Calculates full 2-loop contribution to a_mu in the general THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop(const GeneralTHDM& model)
{
   return calculate_amu_2loop_bosonic(model)
      + calculate_amu_2loop_fermionic(model);
}

} // namespace gm2calc
