cmake_minimum_required(VERSION 3.11) 

project(lsm LANGUAGES Fortran)

MESSAGE( STATUS "test?: ${TEST_FLAGS}")

# 编译选项
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  # 累积编译选项
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast -qopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast")

  if (CMAKE_BUILD_TYPE MATCHES Debug) #
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  endif ()

elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  if (CMAKE_BUILD_TYPE MATCHES Debug) # 测试
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()

endif ()

# 打印编译器
MESSAGE( STATUS "compiler is: ${CMAKE_Fortran_COMPILER_ID}")
MESSAGE( STATUS "compiler flag: ${CMAKE_Fortran_FLAGS}")
MESSAGE( STATUS "compiler type: ${CMAKE_BUILD_TYPE}")

# module位置
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
MESSAGE( STATUS "module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")

# 检索库的位置
if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF}) # 鲁棒性更好
endif ()
if (DEFINED ENV{NETCDF_ROOT}) # 依靠环境变量进行检索
  set(NETCDF_INCLUDEDIR "$ENV{NETCDF_ROOT}/include")
  set(NETCDF_LIBDIR "$ENV{NETCDF_ROOT}/lib")
  MESSAGE( STATUS "netcdf lib: ${NETCDF_LIBDIR}")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

# 获取父目录, 为了嵌套在其他的仓中进行编译
get_directory_property(parent_dir PARENT_DIRECTORY)
MESSAGE( STATUS "parent directory: ${parent_dir}")

# 编译库
add_subdirectory(lib/task)
add_subdirectory(lib/flogger)

set(sources
    src/module_cfg.f90
    src/main.f90
)

add_executable(lsm.exe ${sources})
target_link_directories(lsm.exe PUBLIC ${LAPACK_LIBDIR})
target_link_libraries(lsm.exe task flogger)
