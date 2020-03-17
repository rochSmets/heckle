cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME HECKLE)
project(${PROJECT_NAME})

set(SOURCES_INC
    inc/trans.h
    inc/setup.h
    inc/maxwell.h
    inc/sources.h
    inc/start.h
    inc/check.h
    inc/ions.h
    inc/scrunch.h
    inc/smooth.h
    inc/defines.h
    inc/memory.h
    inc/fill.h
    inc/smalloc.h
    inc/read.h
    inc/ohm.h
    inc/init.h
    inc/structures.h
    inc/subdomains.h
    inc/misc.h
    inc/wall.h
    inc/push.h
    inc/iamdead.h
    src/Collisions/collisionIonIon.h
    src/Collisions/collision.h
    src/Memory/HeckleMem.h
    src/Particle/particle.h
    src/InitModel/initModelUniform.h
    src/InitModel/initModelHarris.h
    src/InitModel/initModelRecoAsymangle.h
    src/InitModel/initModelDoubleHarris.h
    src/InitModel/initModel.h
    src/InitModel/initModelNBeams.h
    src/InitModel/initModelShearB.h
    src/InitModel/initModelNLasers.h
    src/HeckleMPI/ghostg2p.h
    src/HeckleMPI/particlecomm.h
    src/HeckleMPI/ghostg2mom.h
    src/HeckleMPI/ghostg4mom.h
    src/HeckleMPI/ghosts.h
    src/HeckleMPI/ghostg2fields.h
    src/HeckleMPI/ghostg2pwinske.h
    src/HeckleMPI/ghostg4current.h
    src/HeckleMPI/ghostg4winske.h
    src/HeckleIO/HeckleIORestarts.h
    src/HeckleIO/HeckleIOTime.h
    src/HeckleIO/HeckleIOSpecies.h
    src/HeckleIO/HeckleIOFields.h
    src/HeckleBC/hecklebc.h
    src/HeckleBC/perfcond.h
    src/HeckleBC/periodic.h
    src/HeckleBC/particlebc.h
    src/HeckleBC/bc_constants.h
    src/Species/species.h
    src/CloseModel/closeModelIsotherm.h
    src/CloseModel/closeModelPolytrop.h
    src/CloseModel/closeModelFullPressureCommon.h
    src/CloseModel/closeModelFullPressureImplicit.h
    src/CloseModel/closeModelFullPressureSubCycling.h
    src/CloseModel/closeModel.h
    )

set(SOURCES_C
    src/memory.c
    src/Collisions/collisionIonIon.c
    src/Collisions/collision.c
    src/iamdead.c
    src/ions.c
    src/Memory/HeckleMem.c
    src/check.c
    src/read.c
    src/Particle/particle.c
    src/maxwell.c
    src/ohm.c
    src/InitModel/initModelHarris.c
    src/InitModel/initModelRecoAsymangle.c
    src/InitModel/initModelDoubleHarris.c
    src/InitModel/initModelUniform.c
    src/InitModel/initModel.c
    src/InitModel/initModelNBeams.c
    src/InitModel/initModelShearB.c
    src/InitModel/initModelNLasers.c
    src/smooth.c
    src/start.c
    src/push.c
    src/wall.c
    src/smalloc.c
    src/heckle.c
    src/misc.c
    src/HeckleMPI/ghostg2p.c
    src/HeckleMPI/ghosts.c
    src/HeckleMPI/ghostg2mom.c
    src/HeckleMPI/ghostg4mom.c
    src/HeckleMPI/ghostg2fields.c
    src/HeckleMPI/particlecomm.c
    src/HeckleMPI/ghostg2pwinske.c
    src/HeckleMPI/ghostg4current.c
    src/HeckleMPI/ghostg4winske.c
    src/init.c
    src/subdomains.c
    src/sources.c
    src/fill.c
    src/HeckleIO/HeckleIOFields.c
    src/HeckleIO/HeckleIOSpecies.c
    src/HeckleIO/HeckleIORestarts.c
    src/HeckleIO/HeckleIOTime.c
    src/HeckleBC/hecklebc.c
    src/HeckleBC/particlebc.c
    src/HeckleBC/periodic.c
    src/HeckleBC/perfcond.c
    src/Species/species.c
    src/setup.c
    src/CloseModel/closeModelIsotherm.c
    src/CloseModel/closeModelPolytrop.c
    src/CloseModel/closeModelFullPressureCommon.c
    src/CloseModel/closeModelFullPressureImplicit.c
    src/CloseModel/closeModelFullPressureSubCycling.c
    src/CloseModel/closeModel.c
    )

add_executable(
    ${PROJECT_NAME}
    ${SOURCES_INC}
    ${SOURCES_C}
    )

find_package(MPI REQUIRED)

set(HDF5_USE_STATIC_LIBRARIES OFF)
set(HDF5_PREFER_PARALLEL ON)

#set(HDF5_ROOT "/usr/lib64/openmpi/bin/")

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 11)

find_package(HDF5 COMPONENTS C REQUIRED)

if(${HDF5_IS_PARALLEL})
    #message("hdf5 is parallel")
else()
    message("hoops : hdf5 is not parallel !")
endif()

include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_C_LIBRARIES} ${MPI_C_LIBRARIES} m)

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/CloseModel>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/Collisions>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/HeckleBC>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/HeckleIO>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/HeckleMPI>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/InitModel>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/Memory>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/Particle>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/plasma>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/Species>
  )

#add_definitions( -D NS=2 -D FULLP)
add_definitions( -D NS=2 )
