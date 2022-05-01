cmake_minimum_required(VERSION 3.4) 

project("EnsembleMC")

SET(scalfmm_DIR "${PROJECT_SOURCE_DIR}/include/scalfmm/lib/cmake/scalfmm/")

find_package(scalfmm CONFIG REQUIRED)

if(scalfmm_FOUND)
  message(STATUS "ScalFMM Found")
  add_executable(${PROJECT_NAME} ./src/${PROJECT_NAME}.cpp )
  target_link_libraries(${PROJECT_NAME} scalfmm::scalfmm)
else()
  message(FATAL_ERROR "ScalFMM NOT FOUND")
endif()

target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include/")

# install include files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")
install(FILES ${HEADER_FILES} DESTINATION include)
