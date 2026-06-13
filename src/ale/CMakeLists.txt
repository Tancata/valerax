
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

# LORe resolution-branch marginal readout test (not installed): the U-aware
# backtrace posterior under the fitted global r. See tests/lore_marginal.cpp.
add_program(lore_marginal "${CMAKE_SOURCE_DIR}/tests/lore_marginal.cpp")

# LORe AORe-vs-LORe model-selection sanity check (the old r in {0,1} profile
# proxy, kept as a coarse sanity check). See tests/lore_modelselection_sanity.cpp.
add_program(lore_modelselection_sanity
            "${CMAKE_SOURCE_DIR}/tests/lore_modelselection_sanity.cpp")

