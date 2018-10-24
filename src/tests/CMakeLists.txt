

set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};$ENV{HOME}")


add_library(gmock ${GOOGLETEST_DIR}/src/gtest-all.cc ${GOOGLEMOCK_DIR}/src/gmock-all.cc)

target_include_directories(gmock PUBLIC ${GOOGLETEST_DIR} ${GOOGLEMOCK_DIR}
  ${GOOGLETEST_DIR}/include ${GOOGLEMOCK_DIR}/include)




find_package(myUtilities REQUIRED)


add_executable(periodic_q_surfaceTest periodic_q_surfaceTest.cpp)

target_link_libraries(periodic_q_surfaceTest  PUBLIC gmock  ${PROJECT_NAME} myUtilities::myUtilities )
