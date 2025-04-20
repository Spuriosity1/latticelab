all_test_deps = []
all_test_dep_libs = []
project_test_sources = []

gtest_dep = gtest.get_variable('gtest_dep')
thread_dep = dependency('threads')
test_main = files('tests_main.cpp')


test_deps = [
    gtest_dep, 
    thread_dep
]

test_dep_libs = [   
]





cellspectest = executable('cellspectest', ['cellspectest.cpp', test_main],
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs],
  cpp_args : ['-DDEBUG']
  )


modulotest = executable('modulotest', ['modtest.cpp', test_main],
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )

rationaltest = executable('rationaltest', ['rationaltest.cpp', test_main],
  include_directories: g_include,
  dependencies: [main_deps, test_deps],
  link_with: [lattice_indexing_lib, test_dep_libs]
  )

test('modulotest', modulotest)
test('cellspectest', cellspectest)
test('rationaltest', rationaltest)

all_test_deps += main_deps
all_test_dep_libs += lattice_indexing_lib


# This executable contains all the tests
project_test_sources += test_main
all_test_deps += test_deps
all_test_dep_libs += test_dep_libs

all_tests = executable('all_tests', 
    project_test_sources, 
    include_directories : g_include, 
    dependencies : all_test_deps, 
    link_with : all_test_dep_libs)
