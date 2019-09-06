include_directories(${PYTHON_INCLUDE_DIR})

add_library(python-biggles MODULE _biggles.cpp)
target_link_libraries(python-biggles biggles biggles_tools ${Boost_PYTHON_LIBRARY})
set_target_properties(python-biggles PROPERTIES
    OUTPUT_NAME _biggles
    PREFIX ""
)

# Find our where the site packages build destination is
execute_process(COMMAND
    ${PYTHON_EXECUTABLE} -c
    "from distutils.sysconfig import get_python_lib; print(get_python_lib(prefix=''))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# This target installs the python modules proper.
install(TARGETS python-biggles LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES})
install(DIRECTORY biggles DESTINATION ${PYTHON_SITE_PACKAGES} FILES_MATCHING PATTERN *.py)

# Install the biggles wrapper script
configure_file(biggles-tools.in.py tickles @ONLY)
configure_file(biggles-track.in.py biggles @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/biggles DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/tickles DESTINATION bin)

# Set the PYTHONPATH for the biggles and _biggles modules in the build tree
set(BIGGLES_INTERNAL_PYTHONPATH "${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_CURRENT_SOURCE_DIR}")

# Create the internal biggles tool which uses modules in the build tree
configure_file(biggles-internal.in.py biggles-internal.py @ONLY)
set(BIGGLES_INTERNAL_PY "${CMAKE_CURRENT_BINARY_DIR}/biggles-internal.py"
    CACHE INTERNAL "Location of biggles command-line tool which is hard-coded to use build tree paths.")

# Create the internal tickles tool which uses modules in the build tree
configure_file(tickles-internal.in.py tickles-internal.py @ONLY)
set(TICKLES_INTERNAL_PY "${CMAKE_CURRENT_BINARY_DIR}/tickles-internal.py"
    CACHE INTERNAL "Location of tickles command-line tool which is hard-coded to use build tree paths.")

# Install the msmm feature conversion wrapper script
configure_file(msmm-features-to-biggles.in.py msmm-features-to-biggles @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/msmm-features-to-biggles DESTINATION bin)

configure_file(biggles/info.py chunk-info @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/chunk-info DESTINATION bin)

# Create a 'python' target which depends on any custom modules, tools, etc
add_custom_target(python DEPENDS python-biggles ${BIGGLES_INTERNAL_PY})

# Python tests
add_subdirectory(test)
