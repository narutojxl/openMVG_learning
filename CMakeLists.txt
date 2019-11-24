

project(test1)
cmake_minimum_required(VERSION 2.8)


#set( CMAKE_BUILD_TYPE "Debug" ) vs code 运行时不能识别该语句
set( CMAKE_CXX_FLAGS "-std=c++11" )


set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/lib)

# 添加cmake模块
#list( APPEND CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake_modules )


# OpenCV
#set(OpenCV_DIR /home/jxl/opencv-3.2.0/build_dir)  # 这两种方式都可以
set(CMAKE_PREFIX_PATH "/home/jxl/opencv-3.2.0/installed")

find_package(OpenCV 3.2 REQUIRED )
#find_package(OpenCV REQUIRED ) # use this if in OpenCV2 
include_directories(${OpenCV_INCLUDE_DIRS})



#set(Eigen_DIR  "/home/jxl/Eigen_3.3.3/build")
find_package(Eigen3  REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})


# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


# Ceres  我们的ceres是安装在系统路径下，在build openMVG时，即使选择自带的ceres，最后还是find系统路径下的ceres
#find_package(Ceres REQUIRED )
#include_directories( ${CERES_INCLUDE_DIRS} )



# OpenMVG
set(OpenMVG_DIR  "/home/jxl/openMVG/installed/share/openMVG/cmake")
find_package(OpenMVG REQUIRED)
include_directories(${OPENMVG_INCLUDE_DIRS})





add_executable(test1  src/main.cpp)
target_link_libraries(test1  
        ${OpenCV_LIBS}
        ${Boost_LIBRARIES}      
        #${CERES_LIBRARIES}
        ${OPENMVG_LIBRARIES}
        ) 



add_executable(test2  src/main2.cpp)
target_link_libraries(test2  
        ${OpenCV_LIBS}
        ${Boost_LIBRARIES}      
        #${CERES_LIBRARIES}
        ${OPENMVG_LIBRARIES}
        ) 



