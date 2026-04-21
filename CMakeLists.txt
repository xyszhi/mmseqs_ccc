cmake_minimum_required(VERSION 3.10...4.2)
project(mmseqs_ccc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mmseqs_ccc main.cpp)

find_package(Threads REQUIRED)
target_link_libraries(mmseqs_ccc Threads::Threads)

add_executable(extract_cluster extract_cluster.cpp)
