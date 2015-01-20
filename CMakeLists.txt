cmake_minimum_required (VERSION 2.6)


project (blur_hdr_cube)


set( PROJECT_SOURCES
	cubicblur/cubicblur.cpp
	cubicblur/cubicblur.h
	cubicblur/blurprocess.cpp
	cubicblur/blurprocess.h
	cubicblur/hdrloader.cpp
	cubicblur/hdrloader.h
	cubicblur/lodepng.cpp
	cubicblur/lodepng.h
	cubicblur/targetver.h
)

include_directories("${PROJECT_SOURCE_DIR}/submodules/tclap/include")

add_executable(blur_hdr_cube ${PROJECT_SOURCES})