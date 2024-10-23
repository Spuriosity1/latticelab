all_itest_deps = []
all_itest_dep_libs = []

itest_deps = [global_deps]

modtest = executable('modtest', 'modtest.cpp', 
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs]
  )

unitcellspectest = executable('unitcellspectest', 'unitcell_spec_test.cpp', 
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs]
  )


jsonllitest = executable('jsonllitest', 'json_lli_test.cpp', 
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs]
  )
