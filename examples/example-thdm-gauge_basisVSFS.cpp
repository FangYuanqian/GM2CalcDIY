#include "/home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_1loop.h"
#include "/home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_2loop.h"
#include "/home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_uncertainty.hpp"
#include "/home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_error.hpp"
#include "/home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/THDM.hpp"

#include <cstdio>
#include </home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_2loop.hpp>
#include </home/fyq/DIYbyfyq/GM2CalcDIY/include/gm2calc/gm2_1loop.hpp>

int main()
{
   // define THDM parameters in the gauge basis
   gm2calc ::thdm:: Gauge_basis basis;
   basis.yukawa_type = gm2calc ::thdm:: Yukawa_type :: type_2;
   basis.lambda << 4, 0.2, 0.5, 0.8, -2, 0, 0;    // lambda_ {1,...,7}
   basis.tan_beta = 10;                                   // tan(beta)
   basis.m122 = 1000;                                   // m_ {12}^2 in GeV^2
   basis.zeta_u = 0;                                     // zeta_u
   basis.zeta_d = 0;                                     // zeta_d
   basis.zeta_l = 0;                                     // zeta_l
   basis.Delta_u << 0, 0, 0, 0, 0, 0, 0, 0, 0;           // Delta_u
   basis.Delta_d << 0, 0, 0, 0, 0, 0, 0, 0, 0;           // Delta_d
   basis.Delta_l << 0, 0, 0, 0, 0, 0, 0, 0, 0;           // Delta_l
   basis.Pi_u << 0, 0, 0, 0, 0, 0, 0, 0, 0;              // Pi_u
   basis.Pi_d << 0, 0, 0, 0, 0, 0, 0, 0, 0;              // Pi_d
   basis.Pi_l << 0, 0, 0, 0, 0, 0, 0, 0, 0;              // Pi_l

// define SM parameters
   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/128.94579);  // electromagnetic coupling
   sm.set_mu(2, 173.34);               // top quark mass 
   sm.set_mu(1, 1.28);                 // charm quark mass
   sm.set_md(2, 4.18);                 // bottom quark mass
   sm.set_ml(2, 1.77684);              // tau lepton mass

// define options to customize the calculation
   gm2calc::thdm::Config config;
   config.force_output = false;
//the “running masses” scheme is chosen
   config.running_couplings = false;    // use running couplings
//the 2HDM model should be created within a try block
   try {
      // setup the THDM
      const gm2calc::THDM model(basis, sm, config);
      // calculate a_mu up to (including) the 2-loop level
      const double amu = gm2calc::calculate_amu_1loop(model)
                       + gm2calc::calculate_amu_2loop(model);
                       
      // calculate the uncertainty of the 2-loop a_mu
      const double delta_amu =
         gm2calc::calculate_uncertainty_amu_2loop(model);

      std::printf("amu = %.5e +- %.5e\n", amu, delta_amu);
   } catch (const gm2calc::Error& e) {
      std::printf("%s\n", e.what());
   }

   return 0;
}
