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

#ifndef calc_gos_HPP
#define calc_gos_HPP

#include"value_pairs.hpp"
#include"goscalc_shared_data.hpp"
#include"hankel_trafo.hpp"

std::vector<double> calc_energy_losses(const Goscalc_shared_data& gsd);

Value_pairs calculate_gos(unsigned int index_eloss, Hankel_trafo& trafo_plan, const Goscalc_shared_data& gsd);


#endif
