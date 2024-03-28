module Simulations

using Random
using StatsBase, Distributions

export simulate,
    pop_size, chromosome_length,
    DiscreteGensDiploid,
    evolve!,
    getIBDs, joinIBDs, getIBSs, get_paths,
    isfullycoalesced, noncoalesced_length


include("populations.jl")
include("chromosomes.jl")
include("models.jl")

end # module Simulations
