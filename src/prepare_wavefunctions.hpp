#ifndef prepare_wavefunctions_HPP
#define prepare_wavefunctions_HPP

#include<cmath>


#include"goscalc_shared_data.hpp"
#include"value_pairs.hpp"

/**
 *@file prepare_wavefunctions.hpp
 *@brief Functions to prepare the results of wavegen (header)
 */

///Convert values of the bound wave function of the atomic electron to appropriate units
/*!
 * _wavegen_ outputs the wave function resulting from the DFT calculation (which is assumed to be stored in \p tsd)
 * as \f$u = R/r\f$, where \f$r\f$ is the distance from the nucleus and \f$R\f$ is the radial wave function.
 * This function (_prepare_bound_wave_) obtains \f$R\f$ by multiplying the results of _wavegen_ with \f$r\f$.
 * It then converts \f$R\f$ from atomic units to Angstrom.
 */
void prepare_bound_wave(Goscalc_shared_data& gsd);

///Calculate wave function of the (free) atomic electron after ionization
/*!
 * This function first calculates the wave function of the free atomic electron with the help of the
 * contwave-function from \cite Frigge2011. The results of this function are then converted from
 * atomic units and Hartree to Angstrom and electron Volt. The parameter \p l_free is the azimuthal
 * quantum number of the free atomic electron and the potential resulting from the DFT calculation,
 * the energy of the free atomic electron and the grid parameter \p amesh used in wavegen are
 * required form the \p tsd storage container.
 */
void prepare_free_waves(Goscalc_shared_data& gsd);


#endif
