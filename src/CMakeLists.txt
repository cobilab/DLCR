cmake_minimum_required(VERSION 2.8.4)

SET(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lpthread -lm -O3 ")
add_executable (DLCR dlcr.c mem.c msg.c common.c buffer.c cache.c levels.c pmodels.c context.c defs.h)

