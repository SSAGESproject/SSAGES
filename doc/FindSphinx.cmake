# CMake module to locate Sphinx executable necessary to build User Manual
#
# use as find_package (Sphinx) in CMakeLists.txt
#
# Defines the following variables
#
# SPHINX_EXECUTABLE
# SPHINX_FOUND
# SPHINX_VERSION
#

find_program (SPHINX_EXECUTABLE
              NAMES sphinx-build sphinx-build2
              DOC "Path to sphinx-build executable")

# Handle REQUIRED and QUIET
# This will also set SPHINX_FOUND
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Sphinx
                                   DEFAULT_MSG
                                   SPHINX_EXECUTABLE)

mark_as_advanced (SPHINX_EXECUTABLE)
