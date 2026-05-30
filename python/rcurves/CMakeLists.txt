# Ref: https://github.com/google-deepmind/mujoco/blob/main/python/mujoco/CMakeLists.txt

project(rcurves_python)
cmake_minimum_required(VERSION 3.15)

# Support new IN_LIST if() operator.
set(CMAKE_POLICY_DEFAULT_CMP0057 NEW)
# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

enable_language(C)
enable_language(CXX)

if (MSVC AND MSVC_VERSION GREATER_EQUAL 1927)
    set(CMAKE_CXX_STANDARD 20) # For forceinline lambdas.
else ()
    set(CMAKE_CXX_STANDARD 17)
endif ()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

if (MSVC)
    add_compile_options(/Gy /Gw /Oi)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fdata-sections -ffunction-sections)
endif ()

find_package(Python3 COMPONENTS Interpreter Development)
find_package(Eigen3 REQUIRED)

option(RCURVES_PYTHON_USE_SYSTEM_PYBIND11 "Use installed pybind11 version." OFF)
set(PYBIND11_VERSION 2.12.0)
if (RCURVES_PYTHON_USE_SYSTEM_PYBIND11)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)
else ()
    include(FetchContent)
    FetchContent_Declare(pybind11 URL https://github.com/pybind/pybind11/archive/refs/tags/v${PYBIND11_VERSION}.tar.gz)
    FetchContent_MakeAvailable(pybind11)
endif ()

# ==================== RCURVES PYTHON BINDINGS ==================================
macro(rcurves_pybind11_module name)
    pybind11_add_module(${name} ${ARGN})
    #target_compile_options(${name} PRIVATE ${AVX_COMPILE_OPTIONS})
    if (NOT MSVC)
        target_compile_options(
                ${name}
                PRIVATE -Werror
                #-Wall
                -Wimplicit-fallthrough
                #-Wunused
                -Wno-int-in-bool-context
                -Wno-sign-compare
        )
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            # Only allow fallthrough annotation via __attribute__.
            target_compile_options(${name} PRIVATE -Wimplicit-fallthrough=5)

            # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=109561
            if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 11)
                target_compile_options(${name} PRIVATE -Wno-maybe-uninitialized)
            endif ()
        endif ()
    endif ()
    set_target_properties(${name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties(${name} PROPERTIES INSTALL_RPATH $ORIGIN BUILD_WITH_INSTALL_RPATH TRUE)
endmacro()

rcurves_pybind11_module(rcurves_python rcurves_pybind11.cpp)
target_link_libraries(rcurves_python PRIVATE geometry-central polyscope Eigen3::Eigen)
target_compile_options(rcurves_python PUBLIC -DGLM_ENABLE_EXPERIMENTAL)
target_compile_options(rcurves_python PUBLIC -fvisibility=hidden)
target_include_directories(rcurves_python PUBLIC ../../release/include
        ../../deps/libgmultigrid/include
        ../../deps/polyscope/deps/imgui/imgui
        ../../deps/polyscope/deps/imgui/ImGuizmo
        ../../deps/polyscope/deps/imgui/implot
        ../../deps/polyscope/deps/glm
        EIGEN3_INCLUDE_DIRS
)
target_link_directories(rcurves_python PUBLIC ../../release/lib)

set(LIBRARIES_FOR_WHEEL
        "$<TARGET_FILE:rcurves_python>"
)

set(RCURVES_PYTHON_MAKE_WHEEL ON)
if (RCURVES_PYTHON_MAKE_WHEEL)
    add_custom_target(
            wheel ALL
            COMMAND "${CMAKE_COMMAND}" -E rm -rf "${CMAKE_CURRENT_BINARY_DIR}/dist"
            COMMAND "${CMAKE_COMMAND}" -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/.."
            "${CMAKE_CURRENT_BINARY_DIR}/dist"
            COMMAND "${CMAKE_COMMAND}" -E copy "${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE"
            "${CMAKE_CURRENT_BINARY_DIR}/dist/LICENSE"
            COMMAND "${CMAKE_COMMAND}" -E copy ${LIBRARIES_FOR_WHEEL}
            "${CMAKE_CURRENT_BINARY_DIR}/dist/rcurves"
            COMMAND "${Python3_EXECUTABLE}" -m pip wheel --wheel-dir "${CMAKE_BINARY_DIR}" --no-deps -vvv
            "${CMAKE_CURRENT_BINARY_DIR}/dist"
    )

    add_dependencies(
            wheel
            rcurves_python
    )
endif ()