structtest_sources = [
  'structtest.cc'
]

structtest_deps = []
structtest_deps_libs = [struct_lib]


structtest = executable('structtest', [structtest_sources, test_main],
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs]
  )

test('IndexTest', structtest)

all_test_deps += structtest_deps
all_test_dep_libs += structtest_deps_libs
