/* Class to switch the working path
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

#ifndef PATHSWITCHER_HPP
#define PATHSWITCHER_HPP


#include<boost/filesystem.hpp> 


/**
 *@file path_switcher.hpp
 *@brief Small class that changes the current path
 */

///Class to switch the working path
/*!Class that switches the working path (using boost::filesystem) and then switches it back upon
   its destruction.*/
class Path_switcher{
  public:
	///Constructor that switches the working directory to \p new_path
	explicit Path_switcher(boost::filesystem::path new_path)
	: old_path(boost::filesystem::current_path()){
		boost::filesystem::current_path(new_path);
	}


	///Destructor that switches the working director back to the old path
	~Path_switcher(){
		boost::filesystem::current_path(old_path);	
	}


	///Delete copy-assignement operator (move equivalent implicitly deleted)
	Path_switcher & operator = (const Path_switcher&) = delete;
	///Delete copy operator (move equivalent implicitly deleted)
	Path_switcher(const Path_switcher&) = delete;		

	


  private:
	const boost::filesystem::path old_path;
	

};





#endif
