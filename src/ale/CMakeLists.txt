
set(alerax_SOURCES ale.cpp
  AleArguments.cpp
  AleEvaluator.cpp
  AleOptimizer.cpp
  AleState.cpp
  Highways.cpp
  OptimizationClasses.cpp
  TrimFamilies.cpp
  )

add_program(alerax "${alerax_SOURCES}")
target_link_libraries(alerax PRIVATE mpi-scheduler)
install(TARGETS alerax DESTINATION bin)

# WGD regression test (not installed): _hasWGD all-false must reproduce the
# original DL likelihood. See tests/wgd_regression.cpp.
add_program(wgd_regression "${CMAKE_SOURCE_DIR}/tests/wgd_regression.cpp")

