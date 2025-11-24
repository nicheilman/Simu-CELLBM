#The version of cmake required
cmake_minimum_required(VERSION 3.0)

FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(LBM ${SRC})


target_include_directories(LBM PUBLIC  "${PROJECT_SOURCE_DIR}/include/"
                                        "${PROJECT_SOURCE_DIR}/include/LBM"
                                        "${PROJECT_SOURCE_DIR}/include/time_integration"
                                        "${PROJECT_SOURCE_DIR}/include/mesh/cell_types"
                                        "${PROJECT_SOURCE_DIR}/include/math_modules"
)

target_link_libraries(LBM PUBLIC 
    math_modules    
)

