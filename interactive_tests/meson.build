all_itest_deps = []
all_itest_dep_libs = []

itest_deps = [main_deps]

modtest = executable('modtest', 'modtest.cpp', 
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )

unitcellspectest = executable('unitcellspectest', 'unitcell_spec_test.cpp', 
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )

jsonllitest = executable('jsonllitest', 'json_lli_test.cpp', 
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )

map_bench = executable('bench_map_type', 'bench_map_type.cpp',
  include_directories: g_include
  )


bench_cellspec = executable('bench_cellspec', 'bench_cellspec.cpp',
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )
