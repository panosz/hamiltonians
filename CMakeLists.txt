cmake_minimum_required(VERSION 3.12)
project(Hamiltonians VERSION 0.10)

set(CMAKE_CXX_STANDARD 17)

# Offer the user the choice of overriding the googlemock directory
set(GOOGLEMOCK_DIR $ENV{GMOCK_DIR} CACHE PATH "googlemock directory")
set(GOOGLETEST_DIR $ENV{GTEST_DIR} CACHE PATH "googletest directory")

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
#set(INSTALL_INCLUDE_DIR include CACHE PATH
#  "Installation directory for header files")





# The interesting stuff goes here
# ===============================

# add the project subdirectories
#


include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})


add_subdirectory(${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME})

## Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}Targets
        FILE
        ${PROJECT_NAME}Targets.cmake
        NAMESPACE
        ${PROJECT_NAME}::
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)


configure_package_config_file (${PROJECT_NAME}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR} )


#Install the config and configversion files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )


#############################################
## Exporting from the build tree


export(EXPORT ${PROJECT_NAME}Targets FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake NAMESPACE ${PROJECT_NAME}::)


#Register package in user's package registry
export(PACKAGE ${PROJECT_NAME})

#add tests
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tests)

#add examples
add_subdirectory(${PROJECT_SOURCE_DIR}/src/examples)

