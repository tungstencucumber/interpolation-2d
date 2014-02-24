cmake_minimum_required (VERSION 2.8)

set (PROJECT gcm3d-interpolation)
project (${PROJECT})

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMakeModules)

if (CMAKE_COMPILER_IS_GNUCXX)
	message (STATUS "Adding -Wno-deprecated to supress warnings from VTK.")
	add_definitions (-Wno-deprecated)
	message (STATUS "Adding external CXXFLAGS and LDFLAGS.")
	add_definitions (${CXXFLAGS})
	add_definitions (${LDFLAGS})
else ()
	message (WARNING "Compiler untested!")
endif ()

set ( SOURCES
	./Mesh.cpp
	./Method.cpp
	./General.cpp
	./Node.cpp
	./Triangle.cpp
	./VTKSnapshotWriter.cpp
	./SnapshotWriter.cpp
	./Utils.cpp
	./main.cpp )

add_executable (${PROJECT} ${SOURCES})
set_target_properties(${PROJECT} PROPERTIES COMPILE_FLAGS -O3)

add_executable (${PROJECT}-profile EXCLUDE_FROM_ALL ${SOURCES})
set_target_properties(${PROJECT}-profile PROPERTIES COMPILE_FLAGS -pg)


find_package (VTK REQUIRED)
if (VTK_FOUND)
	message(STATUS "VTK found")
	include_directories(${VTK_INCLUDE_DIRS})
	set(LINK_LIBRARIES vtkCommon vtkFiltering vtkIO)
endif (VTK_FOUND)


set(LINK_LIBRARIES ${LINK_LIBRARIES} m)
set(LINK_LIBRARIES ${LINK_LIBRARIES} gmsh)

target_link_libraries (${PROJECT} ${LINK_LIBRARIES})
target_link_libraries (${PROJECT}-profile ${LINK_LIBRARIES})
