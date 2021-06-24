#ifndef calc_gos_HPP
#define calc_gos_HPP

#include"value_pairs.hpp"
#include"goscalc_shared_data.hpp"
#include"hankel_trafo.hpp"

std::vector<double> calc_energy_losses(const Goscalc_shared_data& gsd);

Value_pairs calculate_gos(unsigned int index_eloss, Hankel_trafo& trafo_plan, const Goscalc_shared_data& gsd);


#endif
