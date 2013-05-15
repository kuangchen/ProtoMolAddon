project( ProtoMol )
cmake_minimum_required( VERSION 2.6 )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
SET(CMAKE_CXX_FLAGS "-std=c++0x -O3")
add_definitions (-fpermissive -Wall -g -O3 )

find_package(Lua51)
if (LUA51_FOUND)
  include_directories( ${LUA_INCLUDE_DIRS} )
  LIST( APPEND LIBS ${LUA_LIBRARIES} )
endif (LUA51_FOUND)
                   
find_package(GSL)  
if (GSL_FOUND)
  include_directories( ${GSL_INCLUDE_DIRS} )
  LIST( APPEND LIBS ${GSL_LIBRARIES} )
endif (GSL_FOUND)

find_package(ZLIB)  
if (ZLIB_FOUND)
  include_directories( ${ZLIB_INCLUDE_DIRS} )
  LIST( APPEND LIBS ${ZLIB_LIBRARIES} )
endif (ZLIB_FOUND)

set (HDF5_FIND_COMPONENTS CXX)
find_package(HDF5 REQUIRED ${HDF5_FIND_COMPONENTS})  
if (HDF5_FOUND)
  include_directories( ${HDF5_INCLUDE_DIRS} )
  LIST( APPEND LIBS ${HDF5_LIBRARIES} )
endif (HDF5_FOUND)



option( BUILD_LAPACK "Build with Lapack" Off )
if( BUILD_LAPACK )
	set( BUILD_LAPACK_TYPE "" CACHE STRING "Lapack type to compile against. Choices are: lapack simtk mkl" )

	if( ${BUILD_LAPACK_TYPE} MATCHES "lapack" )
		find_library( LIB_BLAS "blas" HINTS "${BLAS_LIBPATH}" )
		find_library( LIB_LAPACK "lapack" HINTS "${LAPACK_LIBPATH}" )

		if( NOT ${LIB_BLAS} MATCHES "LIB_BLAS-NOTFOUND" AND NOT ${LIB_LAPACK} MATCHES "LIB_LAPACK-NOTFOUND" )
			set( LIBS ${LIBS} ${LIB_BLAS} ${LIB_LAPACK} )
			add_definitions( "-DHAVE_LAPACK" )
		else()
			if( ${LIB_BLAS} MATCHES "LIB_BLAS-NOTFOUND" )
				set( BLAS_INCPATH "" CACHE PATH "Blas include directory" )
				set( BLAS_LIBPATH "" CACHE PATH "Blas library directory" )
			endif()

			if( ${LIB_LAPACK} MATCHES "LIB_LAPACK-NOTFOUND" )
				set( LAPACK_INCPATH "" CACHE PATH "Lapack include directory" )
				set( LAPACK_LIBPATH "" CACHE PATH "Lapack library directory" )
			endif()
		endif()
	endif()

	if( ${BUILD_LAPACK_TYPE} MATCHES "simtk" )
		find_library( LIB_LAPACK_SIMTK "SimTKlapack" HINTS "${SIMTK_LAPACK_LIBPATH}" )
		if( NOT ${LIB_LAPACK_SIMTK} MATCHES "LIB_LAPACK_SIMTK-NOTFOUND" )
			list( APPEND LIBS ${LIB_LAPACK_SIMTK} )
			add_definitions( "-DHAVE_SIMTK_LAPACK" )
		else()
			set( SIMTK_LAPACK_INCPATH "" CACHE PATH "SimTK Lapack include directory" )
			set( SIMTK_LAPACK_LIBPATH "" CACHE PATH "SimTK Lapack library directory" )
		endif()
	endif()

	if( ${BUILD_LAPACK_TYPE} MATCHES "mkl" )
		add_definitions( "-DHAVE_MKL_LAPACK" )
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl=parallel")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl=parallel")
		set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -mkl=parallel")
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -mkl=parallel")
	endif()
endif()

option( BUILD_OPENMM "Build with OpenMM" Off )
if( BUILD_OPENMM )
	add_definitions( "-DHAVE_OPENMM" )

	set( OPENMM_DIR "$ENV{OPENMM_HOME}" CACHE PATH "OpenMM Install Directory")
	if( OPENMM_DIR )
		link_directories( "${OPENMM_DIR}/lib/" )
		include_directories( "${OPENMM_DIR}/include/" )
	endif()

	find_library( OPENMM_LIB "OpenMM" HINT "${OPENMM_DIR}/lib" )
	set( LIBS ${LIBS} ${OPENMM_LIB} )

	option( BUILD_OPENMM_LTMD "Build with LTMD OpenMM" On )
	if( BUILD_OPENMM_LTMD )
		add_definitions( "-DHAVE_OPENMM_LTMD" )
		
		set( LTMDPLUGIN_DIR "$ENV{LTMDOPENMM_HOME}" CACHE PATH "LTMD OpenMM Install Directory")
		if( LTMDPLUGIN_DIR )
			link_directories( "${LTMDPLUGIN_DIR}/lib/" )
			include_directories( "${LTMDPLUGIN_DIR}/include/" )
		endif()
	endif()
	
	find_library( LTMDPLUGIN_LIB "OpenMMLTMD" HINT "${LTMDPLUGIN_DIR}/lib" )
	set( LIBS ${LIBS} ${LTMDPLUGIN_LIB} )
endif()

option( BUILD_GUI "Build with GUI Support" Off )
if( BUILD_GUI )
	find_package( Threads REQUIRED )
	add_definitions( "-DHAVE_GUI" )
	set( LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT} )
endif()

option( BUILD_GROMACS "Build with GROMACS Support" Off )
if( BUILD_GROMACS )
	set( GROMACS_DIR "$ENV{GROMACS_HOME}" CACHE PATH "GROMACS Install Directory" )
	if( GROMACS_DIR )
		link_directories( "${GROMACS_DIR}/lib" )
		include_directories( "${GROMACS_DIR}/include" )
	endif( GROMACS_DIR )

	find_library( GROMACS_MD_LIB "md" HINT "${GROMACS_DIR}/lib" )
	find_library( GROMACS_GMX_LIB "gmx" HINT "${GROMACS_DIR}/lib" )
	set( LIBS ${LIBS} ${GROMACS_MD_LIB} ${GROMACS_GMX_LIB} )

	add_definitions( "-DHAVE_GROMACS" )
endif( BUILD_GROMACS )

# OpenMP
option( BUILD_OPENMP "Build with openmp support" Off )
if( BUILD_OPENMP )
    include( FindOpenMP )
    if( OPENMP_FOUND )
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif( OPENMP_FOUND )
endif( BUILD_OPENMP )

# Parallel
option( BUILD_PARALLEL "Build with MPI Support" On )
if( BUILD_PARALLEL )
	include( FindMPI )
	if( MPI_FOUND )
		add_definitions( "-DHAVE_MPI" )
		include_directories( ${MPI_INCLUDE_PATH} )
		set( LIBS ${LIBS} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
		set( CMAKE_EXE_LINKER_FLAGS ${MPI_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS} )
	endif( MPI_FOUND )
endif( BUILD_PARALLEL )

# Code Coverage
option( CODE_COVERAGE "Build library with code coverage instrumentation" Off )
if( CODE_COVERAGE )
	set( COVERAGE_COMMAND "gcov" )

	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage" )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage" )
	set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage" )

	set( CMAKE_CXX_STANDARD_LIBRARIES "-lgcov" )
endif( CODE_COVERAGE )

# Warnings
option( BUILD_WARNINGS "Build with all warnings turned on" Off )
if( BUILD_WARNINGS )
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wstrict-aliasing")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wstrict-aliasing")
endif( BUILD_WARNINGS )

# Compilation
include_directories( src )

# Library Headers
file( GLOB HEADERS "src/protomol/*.h" )
source_group( "include" FILES ${HEADERS} )
file( GLOB BASE_HEADERS "src/protomol/base/*.h" )
source_group( "include\\Base" FILES ${BASE_HEADERS} )
file( GLOB CONFIG_HEADERS "src/protomol/config/*.h" )
source_group( "include\\Config" FILES ${CONFIG_HEADERS} )
file( GLOB FACTORY_HEADERS "src/protomol/factory/*.h" )
source_group( "include\\Factory" FILES ${FACTORY_HEADERS} )
file( GLOB MODIFIER_HEADERS "src/protomol/modifier/*.h" )
source_group( "include\\Modifier" FILES ${MODIFIER_HEADERS} )
file( GLOB MODULE_HEADERS "src/protomol/module/*.h" )
source_group( "include\\Module" FILES ${MODULE_HEADERS} )
file( GLOB OS_HEADERS "src/protomol/os/*.h" )
source_group( "include\\OS" FILES ${OS_HEADERS} )
file( GLOB OUTPUT_HEADERS "src/protomol/output/*.h" )
source_group( "include\\Output" FILES ${OUTPUT_HEADERS} )
file( GLOB PARALLEL_HEADERS "src/protomol/parallel/*.h" )
source_group( "include\\Parallel" FILES ${PARALLEL_HEADERS} )
file( GLOB SWITCH_HEADERS "src/protomol/switch/*.h" )
source_group( "include\\Switch" FILES ${SWITCH_HEADERS} )
file( GLOB TOPOLOGY_HEADERS "src/protomol/topology/*.h" )
source_group( "include\\Topology" FILES ${TOPOLOGY_HEADERS} )
file( GLOB TYPE_HEADERS "src/protomol/type/*.h" )
source_group( "include\\Type" FILES ${TYPE_HEADERS} )

set( LIBRARY_BASE_HEADERS ${HEADERS} ${BASE_HEADERS} ${CONFIG_HEADERS}
	${FACTORY_HEADERS} ${MODIFIER_HEADERS} ${MODULE_HEADERS} ${OS_HEADERS}
	${OUTPUT_HEADERS} ${PARALLEL_HEADERS} ${SWITCH_HEADERS} ${TOPOLOGY_HEADERS}
	${TYPE_HEADERS} )

if( NOT MSVC )
	file( GLOB DEBUG_HEADERS "src/protomol/debug/*.h" )
	source_group( "include\\Debug" FILES ${DEBUG_HEADERS} )

	set( LIBRARY_BASE_HEADERS ${LIBRARY_BASE_HEADERS} ${DEBUG_HEADERS} )
endif( NOT MSVC )

# Library Force Headers
file( GLOB FORCE_HEADERS "src/protomol/force/*.h" )
source_group( "include\\Force" FILES ${FORCE_HEADERS} )
file( GLOB BONDED_FORCE_HEADERS "src/protomol/force/bonded/*.h" )
source_group( "include\\Force\\Bonded" FILES ${BONDED_FORCE_HEADERS} )
file( GLOB BORN_FORCE_HEADERS "src/protomol/force/born/*.h" )
source_group( "include\\Force\\Born" FILES ${BORN_FORCE_HEADERS} )
file( GLOB COULOMB_FORCE_HEADERS "src/protomol/force/coulomb/*.h" )
source_group( "include\\Force\\Coulomb" FILES ${COULOMB_FORCE_HEADERS} )
file( GLOB EXTENDED_FORCE_HEADERS "src/protomol/force/extended/*.h" )
source_group( "include\\Force\\Extended" FILES ${EXTENDED_FORCE_HEADERS} )
file( GLOB GB_FORCE_HEADERS "src/protomol/force/GB/*.h" )
source_group( "include\\Force\\GB" FILES ${GB_FORCE_HEADERS} )
file( GLOB HESSIAN_FORCE_HEADERS "src/protomol/force/hessian/*.h" )
source_group( "include\\Force\\Hessian" FILES ${HESSIAN_FORCE_HEADERS} )
file( GLOB NONBONDED_FORCE_HEADERS "src/protomol/force/nonbonded/*.h" )
source_group( "include\\Force\\NonBonded" FILES ${NONBONDED_FORCE_HEADERS} )
file( GLOB SYSTEM_FORCE_HEADERS "src/protomol/force/system/*.h" )
source_group( "include\\Force\\System" FILES ${SYSTEM_FORCE_HEADERS} )
file( GLOB TABLE_FORCE_HEADERS "src/protomol/force/table/*.h" )
source_group( "include\\Force\\Table" FILES ${TABLE_FORCE_HEADERS} )

set( LIBRARY_FOURCE_HEADERS ${FORCE_HEADERS} ${BONDED_FORCE_HEADERS}
	${BORN_FORCE_HEADERS} ${COULOMB_FORCE_HEADERS} ${EXTENDED_FORCE_HEADERS}
	${GB_FORCE_HEADERS} ${HESSIAN_FORCE_HEADERS} ${NONBONDED_FORCE_HEADERS}
	${SYSTEM_FORCE_HEADERS} ${TABLE_FORCE_HEADERS} )

# Library Integrator Headers
file( GLOB INTEGRATOR_HEADERS "src/protomol/integrator/*.h" )
source_group( "include\\Integrator" FILES ${INTEGRATOR_HEADERS} )
file( GLOB BASE_INTEGRATOR_HEADERS "src/protomol/integrator/base/*.h" )
source_group( "include\\Integrator\\Base" FILES ${BASE_INTEGRATOR_HEADERS} )
file( GLOB HESSIAN_INTEGRATOR_HEADERS "src/protomol/integrator/hessian/*.h" )
source_group( "include\\Integrator\\Hessian" FILES ${HESSIAN_INTEGRATOR_HEADERS} )
file( GLOB LEAPFROG_INTEGRATOR_HEADERS "src/protomol/integrator/leapfrog/*.h" )
source_group( "include\\Integrator\\Leapfrog" FILES ${LEAPFROG_INTEGRATOR_HEADERS} )
file( GLOB NORMAL_INTEGRATOR_HEADERS "src/protomol/integrator/normal/*.h" )
source_group( "include\\Integrator\\Normal" FILES ${NORMAL_INTEGRATOR_HEADERS} )

if( BUILD_OPENMM )
    set( OPENMM_INTEGRATOR_HEADERS "src/protomol/integrator/openMM/OpenMMIntegrator.h" )
    if( BUILD_OPENMM_LTMD )
         set( OPENMM_INTEGRATOR_HEADERS ${OPENMM_INTEGRATOR_HEADERS} src/protomol/integrator/openMM/NormalModeOpenMM.h )
    endif( BUILD_OPENMM_LTMD )
    
    source_group( "include\\Integrator\\OpenMM" FILES ${OPENMM_INTEGRATOR_HEADERS} )
endif( BUILD_OPENMM )

set( LIBRARY_INTEGRATOR_HEADERS ${INTEGRATOR_HEADERS} ${BASE_INTEGRATOR_HEADERS}
	${HESSIAN_INTEGRATOR_HEADERS} ${LEAPFROG_INTEGRATOR_HEADERS}
	${NORMAL_INTEGRATOR_HEADERS} )

if( BUILD_OPENMM )
	set( LIBRARY_INTEGRATOR_HEADERS ${LIBRARY_INTEGRATOR_HEADERS} ${OPENMM_INTEGRATOR_HEADERS} )
endif( BUILD_OPENMM )

# Library IO Headers
file( GLOB IO_HEADERS "src/protomol/io/*.h" )
source_group( "include\\IO" FILES ${IO_HEADERS} )
file( GLOB GROMACS_IO_HEADERS "src/protomol/io/gromacs/*.h" )
source_group( "include\\IO\\Gromacs" FILES ${GROMACS_IO_HEADERS} )

set( LIBRARY_IO_HEADERS ${IO_HEADERS} ${GROMACS_IO_HEADERS} )

# Combined Headers
set( LIBRARY_HEADERS ${LIBRARY_BASE_HEADERS} ${LIBRARY_FOURCE_HEADERS} ${LIBRARY_INTEGRATOR_HEADERS} ${LIBRARY_IO_HEADERS} )

# Library sources
file( GLOB SOURCES "src/protomol/*.cpp" )
source_group( "src" FILES ${SOURCES} )
file( GLOB BASE_SOURCES "src/protomol/base/*.cpp" )
source_group( "src\\Base" FILES ${BASE_SOURCES} )
file( GLOB CONFIG_SOURCES "src/protomol/config/*.cpp" )
source_group( "src\\Config" FILES ${CONFIG_SOURCES} )
file( GLOB FACTORY_SOURCES "src/protomol/factory/*.cpp" )
source_group( "src\\Factory" FILES ${FACTORY_SOURCES} )
file( GLOB MODIFIER_SOURCES "src/protomol/modifier/*.cpp" )
source_group( "src\\Modifier" FILES ${MODIFIER_SOURCES} )
file( GLOB MODULE_SOURCES "src/protomol/module/*.cpp" )
source_group( "src\\Module" FILES ${MODULE_SOURCES} )
file( GLOB OS_SOURCES "src/protomol/os/*.cpp" )
source_group( "src\\OS" FILES ${OS_SOURCES} )
file( GLOB OUTPUT_SOURCES "src/protomol/output/*.cpp" )
source_group( "src\\Output" FILES ${OUTPUT_SOURCES} )
file( GLOB PARALLEL_SOURCES "src/protomol/parallel/*.cpp" )
source_group( "src\\Parallel" FILES ${PARALLEL_SOURCES} )
file( GLOB SWITCH_SOURCES "src/protomol/switch/*.cpp" )
source_group( "src\\Switch" FILES ${SWITCH_SOURCES} )
file( GLOB TOPOLOGY_SOURCES "src/protomol/topology/*.cpp" )
source_group( "src\\Topology" FILES ${TOPOLOGY_SOURCES} )
file( GLOB TYPE_SOURCES "src/protomol/type/*.cpp" )
source_group( "src\\Type" FILES ${TYPE_SOURCES} )

# Remove Application Sources
list( REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/protomol/main.cpp" )
list( REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/protomol/modules.cpp" )

set( LIBRARY_BASE_SOURCES ${SOURCES} ${BASE_SOURCES} ${CONFIG_SOURCES}
	${FACTORY_SOURCES} ${MODIFIER_SOURCES} ${MODULE_SOURCES} ${OS_SOURCES}
	${OUTPUT_SOURCES} ${PARALLEL_SOURCES} ${SWITCH_SOURCES} ${TOPOLOGY_SOURCES}
	${TYPE_SOURCES} )

if( NOT MSVC )
	file( GLOB DEBUG_SOURCES "src/protomol/debug/*.cpp" )
	source_group( "src\\Debug" FILES ${DEBUG_SOURCES} )

	set( LIBRARY_BASE_SOURCES ${LIBRARY_BASE_SOURCES} ${DEBUG_SOURCES} )
endif( NOT MSVC )

# Library Force sources
file( GLOB FORCE_SOURCES "src/protomol/force/*.cpp" )
source_group( "src\\Force" FILES ${FORCE_SOURCES} )
file( GLOB BONDED_FORCE_SOURCES "src/protomol/force/bonded/*.cpp" )
source_group( "src\\Force\\Bonded" FILES ${BONDED_FORCE_SOURCES} )
file( GLOB BORN_FORCE_SOURCES "src/protomol/force/born/*.cpp" )
source_group( "src\\Force\\Born" FILES ${BORN_FORCE_SOURCES} )
file( GLOB COULOMB_FORCE_SOURCES "src/protomol/force/coulomb/*.cpp" )
source_group( "src\\Force\\Coulomb" FILES ${COULOMB_FORCE_SOURCES} )
file( GLOB EXTENDED_FORCE_SOURCES "src/protomol/force/extended/*.cpp" )
source_group( "src\\Force\\Extended" FILES ${EXTENDED_FORCE_SOURCES} )
file( GLOB GB_FORCE_SOURCES "src/protomol/force/GB/*.cpp" )
source_group( "src\\Force\\GB" FILES ${GB_FORCE_SOURCES} )
file( GLOB HESSIAN_FORCE_SOURCES "src/protomol/force/hessian/*.cpp" )
source_group( "src\\Force\\Hessian" FILES ${HESSIAN_FORCE_SOURCES} )
file( GLOB NONBONDED_FORCE_SOURCES "src/protomol/force/nonbonded/*.cpp" )
source_group( "src\\Force\\NonBonded" FILES ${NONBONDED_FORCE_SOURCES} )
file( GLOB SYSTEM_FORCE_SOURCES "src/protomol/force/system/*.cpp" )
source_group( "src\\Force\\System" FILES ${SYSTEM_FORCE_SOURCES} )
file( GLOB TABLE_FORCE_SOURCES "src/protomol/force/table/*.cpp" )
source_group( "src\\Force\\Table" FILES ${TABLE_FORCE_SOURCES} )

set( LIBRARY_FOURCE_SOURCES ${FORCE_SOURCES} ${BONDED_FORCE_SOURCES}
	${BORN_FORCE_SOURCES} ${COULOMB_FORCE_SOURCES} ${EXTENDED_FORCE_SOURCES}
	${GB_FORCE_SOURCES} ${HESSIAN_FORCE_SOURCES} ${NONBONDED_FORCE_SOURCES}
	${SYSTEM_FORCE_SOURCES} ${TABLE_FORCE_SOURCES} )

# Library Integrator sources
file( GLOB INTEGRATOR_SOURCES "src/protomol/integrator/*.cpp" )
source_group( "src\\Integrator" FILES ${INTEGRATOR_SOURCES} )
file( GLOB BASE_INTEGRATOR_SOURCES "src/protomol/integrator/base/*.cpp" )
source_group( "src\\Integrator\\Base" FILES ${BASE_INTEGRATOR_SOURCES} )
file( GLOB HESSIAN_INTEGRATOR_SOURCES "src/protomol/integrator/hessian/*.cpp" )
source_group( "src\\Integrator\\Hessian" FILES ${HESSIAN_INTEGRATOR_SOURCES} )
file( GLOB LEAPFROG_INTEGRATOR_SOURCES "src/protomol/integrator/leapfrog/*.cpp" )
source_group( "src\\Integrator\\Leapfrog" FILES ${LEAPFROG_INTEGRATOR_SOURCES} )
file( GLOB NORMAL_INTEGRATOR_SOURCES "src/protomol/integrator/normal/*.cpp" )
source_group( "src\\Integrator\\Normal" FILES ${NORMAL_INTEGRATOR_SOURCES} )

if( BUILD_OPENMM )
    set( OPENMM_INTEGRATOR_SOURCES "src/protomol/integrator/openMM/OpenMMIntegrator.cpp" )
    if( BUILD_OPENMM_LTMD )
         set( OPENMM_INTEGRATOR_SOURCES ${OPENMM_INTEGRATOR_SOURCES} src/protomol/integrator/openMM/NormalModeOpenMM.cpp )
    endif( BUILD_OPENMM_LTMD )
    
    source_group( "src\\Integrator\\OpenMM" FILES ${OPENMM_INTEGRATOR_SOURCES} )
endif( BUILD_OPENMM )

set( LIBRARY_INTEGRATOR_SOURCES ${INTEGRATOR_SOURCES} ${BASE_INTEGRATOR_SOURCES}
	${HESSIAN_INTEGRATOR_SOURCES} ${LEAPFROG_INTEGRATOR_SOURCES}
	${NORMAL_INTEGRATOR_SOURCES} )

if( BUILD_OPENMM )
	set( LIBRARY_INTEGRATOR_SOURCES ${LIBRARY_INTEGRATOR_SOURCES} ${OPENMM_INTEGRATOR_SOURCES} )
endif( BUILD_OPENMM )

# Library IO sources
file( GLOB IO_SOURCES "src/protomol/io/*.cpp" )
source_group( "src\\IO" FILES ${IO_SOURCES} )
file( GLOB GROMACS_IO_SOURCES "src/protomol/io/gromacs/*.cpp" )
source_group( "src\\IO\\Gromacs" FILES ${GROMACS_IO_SOURCES} )

set( LIBRARY_IO_SOURCES ${IO_SOURCES} ${GROMACS_IO_SOURCES} )

# My add ons
file( GLOB ADDON_SOURCES "src/protomol/addon/*.cpp" )
source_group( "src\\Addons" FILES ${ADDON_SOURCES} )


# Combined sources
set( LIBRARY_SOURCES ${LIBRARY_BASE_SOURCES} ${LIBRARY_FOURCE_SOURCES} ${LIBRARY_INTEGRATOR_SOURCES} ${LIBRARY_IO_SOURCES} ${ADDON_SOURCES} )

# Library Compilation
add_library( Library STATIC ${LIBRARY_SOURCES} ${LIBRARY_HEADERS} )
target_link_libraries( Library ${LIBS} )

set_target_properties( Library PROPERTIES OUTPUT_NAME ProtoMol )

# Application Compilation
add_executable( Executable src/protomol/main.cpp src/protomol/modules.cpp )
target_link_libraries( Executable Library )

set_target_properties( Executable PROPERTIES OUTPUT_NAME ProtoMol )

install( TARGETS Executable RUNTIME DESTINATION bin )
