find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
add_library(iterativeMethods jacobi.c++ )
target_link_libraries(iterativeMethods PUBLIC generalTools mpi)
target_include_directories(iterativeMethods
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_C_INCLUDE_PATH}
	)
