if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
    # Already in cache, be silent
    set(GMP_FIND_QUIETLY TRUE)
endif()

find_path(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h PATH_SUFFIXES x64)
find_library(GMP_LIBRARY NAMES gmp gmpxx REQUIRED PATH_SUFFIXES x64)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

# NOTE: this has been adapted from CMake's FindPNG.cmake.
if(GMP_FOUND AND NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" IMPORTED_LOCATION "${GMP_LIBRARY}")
endif()
