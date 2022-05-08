add_library(directMethods directMethods.c++)
target_include_directories(directMethods
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
	)
