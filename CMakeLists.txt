#--------------------------------------------------------------
#
# Example of file to be used in CMake to build an external 
# project based on Chrono::Engine
#
#--------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


#--------------------------------------------------------------
#         Modify the project name if you want: 

PROJECT(CoronaElectrostaticSeparator)

#--------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

include(cmake/SBELUtils.cmake)
enable_openmp_support()
enable_cuda_support()

FIND_PACKAGE(ChronoEngine)
FIND_PACKAGE(GLM)
FIND_PACKAGE(GLFW)
FIND_PACKAGE(ZLIB)

SET(CHRONO_PARALLEL_INC  ""    CACHE   PATH    "Where is ChronoEngine Source located")
SET(CHRONO_PARALLEL_LIB_PATH  ""   CACHE    PATH    "Where is ChronoParallel lib located")
SET(BLAZE_DIR ""  CACHE PATH   "Where is Blaze located?")

FIND_LIBRARY(CHRONO_LIB_PARALLEL NAMES ChronoEngine_Parallel PATHS ${CHRONO_PARALLEL_LIB_PATH} REQUIRED)
FIND_LIBRARY(CHRONO_LIB_UTILS NAMES ChronoEngine_Utils PATHS ${CHRONO_PARALLEL_LIB_PATH} REQUIRED)
FIND_LIBRARY(CHRONO_LIB_OPENGL NAMES ChronoEngine_Opengl PATHS ${CHRONO_PARALLEL_LIB_PATH} REQUIRED)
#--------------------------------------------------------------
SET(LIBRARIES ${CHRONOENGINE_LIBRARIES} ${CHRONO_LIB_PARALLEL} ${CHRONO_LIB_UTILS}  ${CHRONO_LIB_OPENGL} ${ZLIB_LIBRARIES} )

include_directories(
    ${CHRONOENGINE_INCLUDES}
    ${CHRONO_PARALLEL_INC}
    ${GLFW_INCLUDE_DIR}
    ${BLAZE_DIR}
    "${CHRONO_PARALLEL_LIB_PATH}/../" 
    "${CUDA_TOOLKIT_ROOT_DIR}/include"
    "${CUDA_SDK_ROOT_DIR}/common/inc"
    ${ZLIB_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    )

#--------------------------------------------------------------
#         Set c++ sources for building the exe, as usual in CMake
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

ADD_EXECUTABLE(conveyor
                source/conveyor_Ida.cpp
				source/SimulatorCES.h				
				#source/UserInterfaceEventReceiver.cpp
				#source/UserInterfaceEventReceiver.h
				source/ElectricParticleProperty.h
				source/ElectricForcesCES.h
				source/ProcessFlow.h
				source/ParserEmitter.h
				source/ParserElectricForcesCES.h
				source/ParserChRandomShapeCreator.h
				source/ParserChDistribution.h )

#--------------------------------------------------------------

#--------------------------------------------------------------

#--------------------------------------------------------------

TARGET_LINK_LIBRARIES(conveyor ${LIBRARIES} )
