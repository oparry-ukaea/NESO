# Add all subdirectories (must have a CMakeLists)
file(GLOB tool_dirs "*")
foreach(dir ${tool_dirs})
  if(EXISTS ${dir}/CMakeLists.txt)
    add_subdirectory(${dir})
  endif(EXISTS ${dir}/CMakeLists.txt)
endforeach(dir ${tool_dirs})
