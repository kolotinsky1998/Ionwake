cmake_minimum_required(VERSION 3.14)
project(Optimization)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Optimization
        src/TArray.cpp
        src/TArray.h
        src/converter.cpp
        src/converter.hpp
        src/main.cpp
        src/scheme.cpp
        src/scheme.hpp
        src/ionwake/IonWake.h
        src/ionwake/IonWake.cpp
        )
