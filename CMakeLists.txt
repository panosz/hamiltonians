cmake_minimum_required(VERSION 3.12)
project(hamiltonians)

set(CMAKE_CXX_STANDARD 17)

find_package( Boost REQUIRED)

add_executable(hamiltonians main.cpp line.cpp line.hpp State.cpp State.hpp Hamiltonian.cpp Hamiltonian.hpp dynamic_system.hpp observer.cpp observer.hpp Integration.cpp Integration.hpp)

target_link_libraries(hamiltonians Boost::boost)