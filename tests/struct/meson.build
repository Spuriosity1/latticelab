# The general dependencies for everyhting in this folder
structtest_deps = []
structtest_deps_libs = [struct_lib]


cellspectest = executable('cellspectest', ['cellspectest.cc', test_main],
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs],
  cpp_args : ['-DDEBUG']
  )

modulotest = executable('modulotest', ['modtest.cc', test_main],
  include_directories: project_inc,
  dependencies: [structtest_deps, test_deps],
  link_with: [structtest_deps_libs, test_dep_libs]
  )

test('modulotest', modulotest)
test('cellspectest', cellspectest)

all_test_deps += structtest_deps
all_test_dep_libs += structtest_deps_libs
