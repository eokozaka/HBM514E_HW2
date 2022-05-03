cmake_minimum_required(VERSION 3.10)

# set the project name
project(HBM512HW2 VERSION 0.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(USE_MYSORT "Use crappy sorting" ON)

configure_file(mainConfig.h.in mainConfig.h)

if(USE_MYSORT)
	add_subdirectory(sort)
	list(APPEND EXTRA_LIBS SortMethods)
	list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR/sort}")
endif()

#add the executable
add_executable(gaussElimination main.c++)

target_link_libraries(gaussElimination PUBLIC ${EXTRA_LIBS})

target_include_directories(gaussElimination PUBLIC
							"${PROJECT_BINARY_DIR}"
							${EXTRA_INCLUDES}
							)	
