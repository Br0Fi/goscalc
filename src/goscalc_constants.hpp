/* define useful constants for other parts of GOScalc to use
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

#ifndef goscalc_constants_HPP
#define goscalc_constants_HPP


#include<unordered_map>
#include<string>
#include<math.h>
#include<complex>

#include<boost/math/constants/constants.hpp>

/**
 *@file transpot_constants.hpp
 *@brief Defines physical constants used throughout transpot
 */


///Physical constants and combinations thereof
/*!References used: \cite Mohr2015\n*/
namespace cst {
	constexpr double a_0 = 0.52917721067;										///<Bohr radius in Angstrom
	constexpr double el = 1.6021766208; 										///<\f$\cdot 10^{-19}\rightarrow\f$ electron charge in As
	constexpr double m_e = 9.10938356; 											///<\f$\cdot 10^{-31}\rightarrow\f$ electron rest mass in kg
	constexpr double hbar = 1.054571800; 										///<\f$\cdot 10^{-34}\rightarrow\f$ reduced Planck constant in J
	constexpr double e_ryd = 13.605693009;                                      ///<Rydberg energy in eV  (unused)

	//Derived constants (\AA = Angstrom):
	constexpr double m_e2dhbarh2 = 2.*el*m_e/pow(hbar,2)*0.01;					///<\f$2m_e/\hbar^2\f$ in \f$1/(\T{\AA}^2\T{eV})\f$

	//Conversion constants:
	constexpr double hart_to_eV = 27.21138602;									///<Conversion factor Hartree to eV
	constexpr double au_to_ang = a_0;											///<Conversion factor atomic units to \f$\T{\AA}\f$

	///Map to convert atomic number to abbreviation
	const std::unordered_map<unsigned int, std::string> elemental_abbrv{
		{1,"H"}, {2,"He"}, {3,"Li"}, {4,"Be"}, {5,"B"}, {6,"C"},
		{7,"N"}, {8,"O"}, {9,"F"}, {10,"Ne"}, {11,"Na"}, {12,"Mg"},
		{13,"Al"}, {14,"Si"}, {15,"P"}, {16,"S"}, {17,"Cl"}, {18,"Ar"},
		{19,"K"}, {20,"Ca"}, {21,"Sc"}, {22,"Ti"}, {23,"V"}, {24,"Cr"},
		{25,"Mn"}, {26,"Fe"}, {27,"Co"}, {28,"Ni"}, {29,"Cu"}, {30,"Zn"},
		{31,"Ga"}, {32,"Ge"}, {33,"As"}, {34,"Se"}, {35,"Br"}, {36,"Kr"},
		{37,"Rb"}, {38,"Sr"}, {39,"Y"}, {40,"Zr"}, {41,"Nb"}, {42,"Mo"},
		{43,"Tc"}, {44,"Ru"}, {45,"Rh"}, {46,"Pd"}, {47,"Ag"}, {48,"Cd"},
		{49,"In"}, {50,"Sn"}, {51,"Sb"}, {52,"Te"}, {53,"I"}, {54,"Xe"},
		{55,"Cs"}, {56,"Ba"}, {57,"La"}, {58,"Ce"}, {59,"Pr"}, {60,"Nd"},
		{61,"Pm"}, {62,"Sm"}, {63,"Eu"}, {64,"Gd"}, {65,"Tb"}, {66,"Dy"},
		{67,"Ho"}, {68,"Er"}, {69,"Tm"}, {70,"Yb"}, {71,"Lu"}, {72,"Hf"},
		{73,"Ta"}, {74,"W"}, {75,"Re"}, {76,"Os"}, {77,"Ir"}, {78,"Pt"},
		{79,"Au"}, {80,"Hg"}, {81,"Tl"}, {82,"Pb"}, {83,"Bi"}, {84,"Po"},
		{85,"At"}, {86,"Rn"}, {87,"Fr"}, {88,"Ra"}, {89,"Ac"}, {90,"Th"},
		{91,"Pa"}, {92,"U"}, {93,"Np"}, {94,"Pu"}, {95,"Am"}, {96,"Cm"},
		{97,"Bk"}, {98,"Cf"}, {99,"Es"}, {100,"Fm"}, {101,"Md"}, {102,"No"},
		{103,"Lr"}
	};
	///Map to convert azimuthal quantum number to abbreviation
	const std::unordered_map<unsigned int, std::string> subshell_abbrv{
		{0,"s"}, {1,"p"}, {2,"d"}, {3,"e"}, {5,"f"}, {6,"g"}, {7,"h"}
	};
}



#endif
