file(GLOB_RECURSE C_HEADERS "*.h")
foreach(_headerFile ${C_HEADERS})
  get_filename_component(_dir ${_headerFile} PATH)
  list(APPEND C_INCLUDE_DIRS ${_dir})
endforeach()

file(GLOB_RECURSE CXX_HEADERS "*.hpp")
foreach(_headerFile ${CXX_HEADERS})
  get_filename_component(_dir ${_headerFile} PATH)
  list(APPEND C_INCLUDE_DIRS ${_dir})
endforeach()

if(C_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES C_INCLUDE_DIRS)
endif()
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${C_INCLUDE_DIRS}")

set(MOLSIM_SOURCES "molsim.cpp")
add_library(molsim ${MOLSIM_SOURCES})
target_link_libraries(molsim ${EXTRA_LIBRARIES})

set_target_properties(molsim PROPERTIES LINK_FLAGS "${C_LINK_FLAGS}")
