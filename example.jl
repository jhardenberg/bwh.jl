include("./bwh_pattern.jl")
using NCDatasets
#set the parameters
#ϕ will be interpreted as either ϕ_b or ϕ_w depending on the flag settings in bwh_pattern.jl
P=params(p_low=1.5, p_high=1.8, p_mode="linear", t_drop=200, tau_rate=0.003, tau_drought=1000, nt=2000, dir_out="/work/users/sfilippini/code/")
bwh(P)