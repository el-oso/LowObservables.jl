# Entry point for the test/ project environment.
# For CI: julia --project=test test/runtests.jl
# During development (faster, uses main project):
#   julia --project=. -e 'using LowObservables, ReTestItems; runtests(LowObservables)'

using LowObservables
using ReTestItems

runtests(LowObservables; testitem_timeout = 120, nworkers = 0)
