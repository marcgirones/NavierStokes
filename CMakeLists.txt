cmake_minimum_required(VERSION 3.16)
project(NavierStokes LANGUAGES CXX)

# --- Estándar y tipo de build -------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo de build" FORCE)
endif()

# --- FFTW3 (sudo apt install libfftw3-dev) ------------------------------------
# Primero intenta pkg-config; si no está instalado, busca la cabecera y la
# librería directamente.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
endif()

if(TARGET PkgConfig::FFTW3)
  set(FFTW3_TARGET PkgConfig::FFTW3)
else()
  find_path(FFTW3_INCLUDE_DIR fftw3.h)
  find_library(FFTW3_LIBRARY NAMES fftw3)
  if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
    message(FATAL_ERROR
      "FFTW3 not found. Use: sudo apt install libfftw3-dev")
  endif()
  add_library(fftw3_dep INTERFACE)
  target_include_directories(fftw3_dep INTERFACE ${FFTW3_INCLUDE_DIR})
  target_link_libraries(fftw3_dep INTERFACE ${FFTW3_LIBRARY})
  set(FFTW3_TARGET fftw3_dep)
endif()

# --- Librería del solver (Mesh, términos FVM, Solver) -------------------------
add_library(nscore
  src/Mesh.cpp
  src/NavierStokesTerms.cpp
  src/Solver.cpp
)
target_include_directories(nscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(nscore PUBLIC ${FFTW3_TARGET})
target_compile_options(nscore PRIVATE -Wall -Wextra)

# --- Ejecutable ---------------------------------------------------------------
# Se llama "NS" porque es el valor por defecto de --exe en
# scripts/convergence_time_check.py
#add_executable(NS main/NS.cpp)
add_executable(case_A main/case_A.cpp)
target_link_libraries(case_A PRIVATE nscore)

add_executable(case_B main/case_A.cpp)
target_link_libraries(case_B PRIVATE nscore)

add_executable(case_C main/case_C.cpp)
target_link_libraries(case_C PRIVATE nscore)

foreach(target nscore case_A case_B case_C)
	target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()

# Copia el ejecutable junto a la raíz del proyecto para poder usar ./NS
add_custom_command(TARGET case_A POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:case_A> ${CMAKE_SOURCE_DIR}/case_A
	COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:case_B> ${CMAKE_SOURCE_DIR}/case_B
	COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:case_C> ${CMAKE_SOURCE_DIR}/case_C
	COMMENT "Copiando case_A a la raíz del proyecto")
