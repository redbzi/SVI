cmake_minimum_required(VERSION 3.6)
project(AMDP)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp interpolation.cpp black-scholes.cpp reverse-bs.cpp extrapolation.cpp svimodel.cpp surfaceVol.cpp plot.cpp)
SET(LIBRARIES libboost_system.dylib libboost_iostreams.dylib libboost_filesystem.dylib)

add_definitions( "-Wall -O3" )

INCLUDE_DIRECTORIES(/usr/include/)
LINK_DIRECTORIES(/usr/local/lib/)

add_executable(AMDP ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(AMDP ${LIBRARIES})