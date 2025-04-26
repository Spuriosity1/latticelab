main_sources = files(
  'UnitCellSpecifier.cpp',
  'preset_cellspecs.cpp',
  'rationalmath.cpp'
  )


# the main static lib
lattice_indexing_lib = static_library('lattice_indexing_lib', 
  main_sources, 
  dependencies: main_deps,
  include_directories: main_include,
  install: true
  )

subdir('tools')
