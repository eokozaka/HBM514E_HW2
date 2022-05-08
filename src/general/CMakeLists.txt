add_library(generalTools general.c++)
target_include_directories(generalTools
	INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
	)
