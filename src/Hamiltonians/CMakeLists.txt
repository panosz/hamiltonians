find_package(Boost REQUIRED)

message(STATUS "boost include dirs = " ${Boost_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)

find_package(myUtilities REQUIRED)

add_library(
        Hamiltonians src/line.cpp include/line.hpp src/State.cpp include/State.hpp
        src/Hamiltonian.cpp include/Hamiltonian.hpp include/dynamic_system.hpp src/observer.cpp
        include/observer.hpp src/Integration.cpp include/Integration.hpp include/PoincareSurface.hpp
        src/dynamic_system.cpp src/action_angle.cpp include/action_angle.hpp include/periodic_q_surface.hpp)


set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};$ENV{HOME}")


target_link_libraries(Hamiltonians Boost::boost myUtilities::myUtilities armadillo)

target_compile_options(Hamiltonians
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


target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


target_include_directories(
        ${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

install (
        TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
