project('lattice_lib', ['cpp', 'c'],
  version: '1.0.1',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=release', 'cpp_std=c++2a' ]
)

project_description = 'Indexing library for lattices'

cmake = import('cmake')


# (global) variables tracking sources and libraries
library_path = []
project_sources = []
project_header_files = []
project_test_sources = []
project_benchmark_sources = []

project_inc = ['include']
project_lib = []



# build dependencies
snf_project = subproject('smith_normal_form', default_options: [
  'enable-tests=false'
  ])
snf_dep = snf_project.get_variable('snf_dep')


# Configure the CMake project
# HighFive_subproj = cmake.subproject('HighFive')

# Fetch the dependency object
# highfive_dep = HighFive_subproj.dependency('HighFive')


# import the #define directives
cpp_proj_arguments = snf_project.get_variable('cpp_proj_arguments')
add_project_arguments(cpp_proj_arguments, language : 'cpp')


main_deps = [
  #  dependency('Armadillo', required: true),
  #  dependency('argparse', required: true),
  #  dependency('HDF5', required: true),
  dependency('nlohmann_json', required: true),
  snf_dep
  ]



# builds libraries
main_include = include_directories(project_inc)
g_include = main_include

# collect the sources
subdir('src')
subdir('include')


# build tests
if get_option('enable-tests')
  gtest = subproject('gtest')
  subdir('tests')
  subdir('interactive_tests')
endif


# make project available
pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : lattice_indexing_lib,
                 version : '1.0',
                 name : 'liblatindex',
                 description : 'A library for indexing r-chains')


latlib_dep = declare_dependency(
  link_with : lattice_indexing_lib,
  include_directories : g_include
)
