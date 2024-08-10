# Build

Replace the existing TRIANGLE compilation in CMakeLists.txt with the following:

```
# TRIANGLE (with external test)
# set(TRIANGLE_SOURCE_DIR "${EXT_PROJECTS_DIR}/Triangle/")
# include_directories(${TRIANGLE_SOURCE_DIR})
# message("Included Triangle in: ${TRIANGLE_SOURCE_DIR}")
# set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} -DUSETRILIBRARY")
# add_library(triangle SHARED "${TRIANGLE_SOURCE_DIR}/triangle.c" "${CMAKE_CURRENT_LIST_DIR}/external/Triangle/external_test.cpp")
# target_compile_definitions(triangle PRIVATE -DTRILIBRARY -DEXTERNAL_TEST)# -DNO_TIMER)
# target_link_libraries(shell triangle)
# target_link_libraries(libshell triangle)
```

# Usage

In simulation, add the following snippet:

```
if (parser.parse<bool>("-resFunc", false)) {
    ::setTriangleResolutionFunction([](Real x, Real y) {
        Real originX = 0.625;
        Real originY = -0.54;
        Real squaredDistanceToOrigin = (x - originX) * (x - originX)
                                        + (y - originY) * (y - originY);
        Real sigma = 0.05;
        return 0.04 * (1.0 - 0.95 * exp(-squaredDistanceToOrigin / 2.0 / sigma / sigma));
    });
}
```
