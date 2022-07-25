/* Calculate the GOS by executing the other parts of GOScalc and combine them to obtain the result
    Copyright (C) 2019-2021 The GOScalc devolopers

    This file is part of GOScalc.

    GOScalc is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include"calc_gos.hpp"

#include<assert.h>

#include<armadillo>
#include<wignerSymbols.h>

#include"goscalc_constants.hpp"


std::vector<double> calc_energy_losses(const Goscalc_shared_data& gsd){
  const unsigned no_energy_losses = gsd.cd.energy_free_steps+1;
  std::vector<double> energy_losses(no_energy_losses);
  for (unsigned int i_e = 0; i_e < no_energy_losses; i_e++){
    const double energy_free = gsd.cd.energy_free_start*pow(gsd.cd.energy_free_increase,i_e);
    energy_losses[i_e] = abs(gsd.cd.bound_energy) + energy_free;
  }

  return energy_losses;
}



Value_pairs calculate_gos(unsigned int index_eloss, Hankel_trafo& trafo_plan, const Goscalc_shared_data& gsd){
  const unsigned int l = gsd.cd.l_bound;
  Value_pairs result;

  for (unsigned int lfree=0; lfree<gsd.free_waves[index_eloss].size(); lfree++){
    Value_pairs integrand = gsd.prepd_bound_wave*gsd.free_waves[index_eloss][lfree];

    //see leapman1980 eq. 4 and eq. 7
    int lower_valid_lambda = std::abs(static_cast<int>(lfree-l));
    int upper_valid_lambda = static_cast<int>(lfree+l);

    for (int lambda = lower_valid_lambda; lambda<=upper_valid_lambda; lambda++){
      double wigner_symb = WignerSymbols::wigner3j(lfree,lambda,l,0,0,0); //values seem to be in [-1,1]
      if (std::abs(wigner_symb) <= 1.E-10){continue;}

      const Value_pairs rad_int = trafo_plan.perform_hankel_trafo(integrand.column_as_vector(1), static_cast<unsigned int>(lambda));

      if ((lambda == lower_valid_lambda) && (lfree == 0)){
        result = (2*lfree+1)*(2*lambda+1)*pow(wigner_symb,2)*rad_int*rad_int;
      } else {
        result += (2*lfree+1)*(2*lambda+1)*pow(wigner_symb,2)*rad_int*rad_int;
      }
    }
  }

  const arma::colvec inv_q2 = arma::pow(result.column_as_colvec(0),-2);
  const double energy_loss = gsd.results_energy_losses[index_eloss];
  const arma::colvec prefactor = cst::m_e2dhbarh2*energy_loss*inv_q2; //unitless(!)
  return prefactor*result;
}

