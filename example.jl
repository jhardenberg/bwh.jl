include("./bwh_pattern.jl")
using NCDatasets
#set the parameters
#run a simulation with two identical species except for dispersal method. The seed dispersing species has
#a characteristic dispersal distance of 3 non-dimensional units
P=params(p=1.5, ls=3, nt=5000)
bwh(P)