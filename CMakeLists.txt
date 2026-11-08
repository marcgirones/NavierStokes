cmake_minimum_required(VERSION 3.16)
project(NavierStokes LANGUAGES CXX)

file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/results")

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

add_executable(case_B main/case_B.cpp)
target_link_libraries(case_B PRIVATE nscore)

add_executable(case_C main/case_C.cpp)
target_link_libraries(case_C PRIVATE nscore)

foreach(target nscore case_A case_B case_C)
	target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()

# Copia cada ejecutable junto a la raíz del proyecto (./case_A, ./case_B, ...)
foreach(target case_A case_B case_C)
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:${target}> ${CMAKE_SOURCE_DIR}/${target}
    COMMENT "Copiando ${target} a la raíz del proyecto")
endforeach()

# --- Targets para lanzar los scripts de Python --------------------------------
# Uso:
#   cmake --build build --target caseA     # valoracio_c_d.py
#   cmake --build build --target caseB     # case_B + divergence_plot.py
#   cmake --build build --target caseC     # convergence_time_check.py
#   cmake --build build --target caseAll   # los tres, en orden A -> B -> C
#
# ("all" no se puede usar como nombre: CMake lo reserva para "compilar todo".)
#
# Cada target compila antes el ejecutable C++ que necesita. Los scripts se
# ejecutan con el mismo directorio de trabajo para el que fueron escritos
# (A desde scripts/, B y C desde la raíz del proyecto).
#
# Argumentos extra opcionales, por ejemplo:
#   cmake -S . -B build -DCASEA_ARGS="--sizes;10;50;100" -DCASEC_ARGS="--sizes;16;32;--time;2"
find_package(Python3 COMPONENTS Interpreter)

if(Python3_Interpreter_FOUND)
  set(CASEA_ARGS "" CACHE STRING "Argumentos extra para valoracio_c_d.py (separados por ;)")
  set(CASEB_ARGS "" CACHE STRING "Argumentos del solver case_B (N tFinal nu, separados por ;)")
  set(CASEC_ARGS "" CACHE STRING "Argumentos extra para convergence_time_check.py (separados por ;)")

  set(SCRIPTS_DIR ${CMAKE_SOURCE_DIR}/scripts)

  # caseA: valoracio_c_d.py (sus rutas por defecto son relativas a scripts/)
  add_custom_target(caseA
    COMMAND ${Python3_EXECUTABLE} valoracio_c_d.py --exe $<TARGET_FILE:case_A> ${CASEA_ARGS}
    WORKING_DIRECTORY ${SCRIPTS_DIR}
    DEPENDS case_A
    USES_TERMINAL
    COMMENT "caseA: valoracio_c_d.py")

  # caseB: case_B escribe los CSV en results/case_B_results (relativo al cwd)
  # y divergence_plot.py los lee de ahí.
  add_custom_target(caseB
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/results/case_B_results
    COMMAND $<TARGET_FILE:case_B> ${CASEB_ARGS}
    COMMAND ${Python3_EXECUTABLE} ${SCRIPTS_DIR}/divergence_plot.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS case_B
    USES_TERMINAL
    COMMENT "caseB: case_B + divergence_plot.py")

  # caseC: convergence_time_check.py (rutas por defecto relativas a la raíz)
  add_custom_target(caseC
    COMMAND ${Python3_EXECUTABLE} ${SCRIPTS_DIR}/convergence_time_check.py --exe $<TARGET_FILE:case_C> ${CASEC_ARGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS case_C
    USES_TERMINAL
    COMMENT "caseC: convergence_time_check.py")

  # caseAll: los tres, secuencialmente (aunque se use make -j)
  add_custom_target(caseAll
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target caseA
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target caseB
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target caseC
    USES_TERMINAL
    COMMENT "caseAll: A -> B -> C")
else()
  message(WARNING "Python3 no encontrado: no se crean los targets caseA/caseB/caseC/caseAll")
endif()
