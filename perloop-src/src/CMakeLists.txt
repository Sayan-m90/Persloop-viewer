#project( Persloop )

cmake_minimum_required( VERSION 3.5 FATAL_ERROR )
#cmake_policy( SET CMP0017 NEW )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" )
link_directories(${CMAKE_SOURCE_DIR}/src/lib)
# Boost
find_package( Boost 1.5 COMPONENTS filesystem program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

include_directories( "." ".." ) 

add_executable( persloop
	persloop.cpp
	AnnotationMatrix.cpp
	ShortestGraph.h
	SimplexNode.h
	SimplicialComplex.cpp
	SimplicialComplex.h
	UnionFindDeletion.cpp
	Legal.h
	ParseCommand.h
	trygraph.h
	) 

target_link_libraries( persloop ${Boost_LIBRARIES})
