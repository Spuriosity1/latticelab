struct_source = [
'UnitCellSpecifier.cpp',
'preset_cellspecs.cpp',
'rationalmath.cpp'
]

project_sources += struct_source

struct_lib = static_library('struct', struct_source, 
  dependencies: [global_deps],
  )

project_lib += struct_lib
