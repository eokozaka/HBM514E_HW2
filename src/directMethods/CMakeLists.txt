find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
add_library(directMethods gaussElimination.c++ lu.c++ )
target_link_libraries(directMethods PUBLIC generalTools mpi)
target_include_directories(directMethods
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_C_INCLUDE_PATH}
	)
