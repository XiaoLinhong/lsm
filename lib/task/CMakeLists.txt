set(sources
  src/module_tool.f90
  src/module_constant.f90
  src/module_structure.f90
  src/module_ncio.f90
)

add_library(task STATIC ${sources})
target_include_directories(task PUBLIC ${NETCDF_INCLUDEDIR})
target_link_directories(task PUBLIC ${NETCDF_LIBDIR})
target_link_libraries(task netcdff)
