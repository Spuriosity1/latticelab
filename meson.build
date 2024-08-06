project('qsi_defect_pbc', 'cpp',
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=debugoptimized', 'cpp_std=c++2a' ]
)

# (global) variables tracking sources and libraries
library_path = []
project_sources = []
project_header_files = []
project_test_sources = []
project_benchmark_sources = []

project_inc = []
project_lib = []

# add further include dirs here

# build dependencies
snf_project = subproject('smith_normal_form', default_options: [
  'build-armadillo=true'
  ])
snf_dep = snf_project.get_variable('snf_dep')

# import the #define directives
cpp_proj_arguments = snf_project.get_variable('cpp_proj_arguments')
add_project_arguments(cpp_proj_arguments, language : 'cpp')


global_deps = [
  dependency('Armadillo', required: true),
  snf_dep
  ]

# builds libraries
subdir('src')


# build tests
if get_option('enable-tests')
  gtest = subproject('gtest')
  subdir('tests')
  subdir('interactive_tests')
endif

# builds the main executable
main_bin = executable('generate_spinon_basis', gen_spinon_basis_source,
  include_directories: project_inc, dependencies: global_deps, link_with: project_lib)

#
## This adds the clang format file to the build directory
#configure_file(input : '.clang-format',
#               output : '.clang-format',
#	       copy: true)
#run_target('format',
#  command : ['clang-format','-i','-style=file',project_sources,project_test_sources,project_benchmark_sources,project_header_files])
#
## This regex excludes any sources from the third_party, tests, benchmarks 
## and gtest related files.
#regex = '^((?!(third_party|tests|benchmarks|gtest)).)*$'
## This adds clang tidy support
#configure_file(input : '.clang-tidy',
#               output : '.clang-tidy',
#               copy : true)
#run_target('tidy',
#  command : ['run-clang-tidy.py','-fix', '-j' , '8', 'files', regex ,'-format', '-p='+ meson.build_root()])
