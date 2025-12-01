include("./bwh_pattern.jl")
using NCDatasets
#set the parameters
#ϕ will be interpreted as either ϕ_b or ϕ_w depending on the flag settings in bwh_pattern.jl
P=params(p=1.5, ϕ=0.01, nt=5000)
bwh(p)