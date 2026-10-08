cmake_minimum_required(VERSION 3.20)
project(stormline LANGUAGES CXX VERSION 0.1.0)

add_executable(stormline
    src/main.cpp
)

target_compile_features(stormline PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(stormline PRIVATE /W4 /permissive-)
else()
    target_compile_options(stormline PRIVATE -Wall -Wextra -Wpedantic)
endif()
