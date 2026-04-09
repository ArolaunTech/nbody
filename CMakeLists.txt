cmake_minimum_required(VERSION 3.25)

project(nbody LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE nbody_src "src/*.h" "src/*.cpp")

add_compile_definitions(LIB="/Users/aaditya/cspice/include/SpiceUsr.h")

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-Wall temp)

if(temp)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

check_cxx_compiler_flag(-O3 temp2)

if(temp2)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

add_executable(nbody ${nbody_src})
target_link_libraries(nbody /Users/aaditya/cspice/lib/cspice.a)

target_include_directories(nbody PUBLIC "glad/include")