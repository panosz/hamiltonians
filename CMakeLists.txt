cmake_minimum_required(VERSION 3.12)
project(hamiltonians)

set(CMAKE_CXX_STANDARD 17)

find_package( Boost REQUIRED)

add_executable(hamiltonians main.cpp line.cpp line.hpp State2.cpp State2.hpp Hamiltonian.cpp Hamiltonian.hpp dynamic_system.hpp observer.cpp observer.hpp)

target_link_libraries(hamiltonians Boost::boost)