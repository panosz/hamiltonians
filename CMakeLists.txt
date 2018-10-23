cmake_minimum_required(VERSION 3.12)
project(hamiltonians)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

find_package(Armadillo REQUIRED)

find_package(myUtilities REQUIRED)

add_library(hamiltonians line.cpp line.hpp State.cpp State.hpp Hamiltonian.cpp Hamiltonian.hpp dynamic_system.hpp observer.cpp observer.hpp Integration.cpp Integration.hpp PoincareSurface.hpp dynamic_system.cpp action_angle.cpp action_angle.hpp)


target_link_libraries(hamiltonians Boost::boost myUtilities::myUtilities armadillo)

target_compile_options(hamiltonians
        PRIVATE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
        -Weffc++ -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wmisleading-indentation
        -Wduplicated-cond
        #-Wduplicated-branches
        -Wlogical-op
        -Wnull-dereference -Wuseless-cast -Wdouble-promotion
        -Wformat=2

        $<$<CONFIG:Release>:
        -O2
        -march=native
        >
        )

add_executable(duffing_example main.cpp)
target_link_libraries(duffing_example hamiltonians myUtilities::myUtilities Boost::boost)

add_executable(pentulum_numeric pendulum_example.cpp)