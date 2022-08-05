cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project("GOScalc" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_VERBOSE_MAKEFILE ON)

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_SOURCE_DIR}/build
 )
include(version.cmake)


## Sources
set(SOURCES main.cpp
			goscalc_shared_data.hpp
			path_switcher.hpp
			logger.hpp
			prepare_wavefunctions.hpp	prepare_wavefunctions.cpp
			goscalc_io.hpp				goscalc_io.cpp
			value_pairs.hpp				value_pairs.cc
			hankel_trafo.hpp			hankel_trafo.cc
			calc_gos.hpp				calc_gos.cpp
      conwave.c)
set(SOURCES_ABS)
foreach (source ${SOURCES})
	list(APPEND SOURCE_ABS ${PROJECT_SOURCE_DIR}/src/${source})
endforeach()

# set path to wignerSymbols library.
SET(CMAKE_EXE_LINKER_FLAGS
          "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib")

## Targets
add_executable(goscalc ${SOURCE_ABS})
set_target_properties(goscalc PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(goscalc fftw3 wignerSymbols armadillo boost_system boost_filesystem)
target_compile_definitions(goscalc PRIVATE VERSION="${VERSION}")


##Compiler flags
set(RELEASE_FLAGS -Ofast -march=native)
set(DEBUG_FLAGS -Wall -Wextra -Wconversion -Wpedantic -Ofast -march=native)
set(DEBUG_GDB_FLAGS -Wall -Wextra -Wconversion -Wpedantic -g -O0) # -fsanitize=address -Wfloat-equal) #gdb; file goscalc; run;
#set(DEBUG_GDB_FLAGS -Wall -Wextra -Wconversion -Wpedantic -g -O0 -fsanitize=address -Wfloat-equal) #gdb; file goscalc; run;
set(PROFILE_FLAGS -Ofast -march=native -pg -no-pie) #after execution: gprof ./goscalc > profile.txt
#set_target_properties(goscalc PROPERTIES LINK_FLAGS "-pg -no-pie") #necessary for profile build
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
    #target_compile_options( goscalc PRIVATE ${RELEASE_FLAGS})
    target_compile_options( goscalc PRIVATE ${DEBUG_GDB_FLAGS})
    #target_link_options(goscalc PRIVATE ${DEBUG_GDB_FLAGS})
endif()
