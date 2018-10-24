find_package(Boost REQUIRED)
find_package(myUtilities REQUIRED)

add_executable(duffing_example main.cpp)
target_link_libraries(duffing_example PUBLIC ${PROJECT_NAME}  myUtilities::myUtilities Boost::boost)
set_target_properties(duffing_example PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_executable(pentulum_numeric pendulum_example.cpp)
target_link_libraries(pentulum_numeric PUBLIC ${PROJECT_NAME}  myUtilities::myUtilities Boost::boost)
set_target_properties(pentulum_numeric PROPERTIES EXCLUDE_FROM_ALL TRUE)
