cmake_minimum_required(VERSION 3.16)
project(my_wave_solver CXX)

# Tell CMake to use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Path to MFEM
set(MFEM_DIR "${CMAKE_SOURCE_DIR}/../mfem")

# Include MFEM headers
include_directories(${MFEM_DIR})

# Add the executable
add_executable(my_wave_solver my_wave_solver.cpp)

# Link to the MFEM static library
target_link_libraries(my_wave_solver
  ${MFEM_DIR}/libmfem.a
  dl
  pthread
)