add_library(directMethods gaussElimination.c++ lu.c++ )
target_link_libraries(directMethods PUBLIC generalTools)
target_include_directories(directMethods
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
	)
