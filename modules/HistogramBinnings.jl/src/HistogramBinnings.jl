module HistogramBinnings




import StatsAPI: fit
import StatsBase: midpoints, Histogram, sturges




export LogEdgeVector, LinEdgeVector, midpoints

abstract type AbstractEdgeVector{T} <:  AbstractVector{T} end


"""
    LogEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Real

Construct a logarithmically spaced edge vector of type `T`. The edges are spaced
evenly in log space. `lo` and `hi` are the lower and upper limits of the histogram.
`nbins` is the number of bins.

# Examples

```julia
using StatsBase, HistogramBinnings

h = Histogram(LogEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
append!(h, vs)
```
"""
struct LogEdgeVector{T} <: AbstractEdgeVector{T}
    edges::Vector{T}
end

"""
    LinEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Real

Construct a linearly spaced edge vector of type `T`. The edges are spaced
evenly in linear space. `lo` and `hi` are the lower and upper limits of the histogram.
`nbins` is the number of bins.

# Examples

```julia
using StatsBase, HistogramBinnings

h = Histogram(LinEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
append!(h, vs)
```
"""
struct LinEdgeVector{T} <: AbstractEdgeVector{T}
    edges::Vector{T}
end

Base.show(io::IO, r::AbstractEdgeVector) = print(io, "$(typeof(r)): $(r.edges)")




function LogEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Integer
    @assert (lo > 0) && (hi > 0) && (nbins > 0) && (hi > lo)
    lo = floor(Int, lo)
    hi = ceil(Int, hi)
    edges = unique(floor.(Int,
            10 .^ (range(log10(lo), stop = log10(hi+1), length = nbins+1))))
    edges[1] = lo
    edges[end] = hi + 1
    LogEdgeVector{T}(edges)
end

function LogEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Real
    @assert (lo > 0) && (hi > 0) && (nbins > 0) && (hi > lo)
    edges = 10 .^ (range(log10(lo), stop = log10(hi), length = nbins+1))
    edges[1] = lo
    edges[end] = hi + eps(T)
    LogEdgeVector{T}(edges)
end

function LogEdgeVector(; lo = 1, hi = 10, nbins::Integer)
    lo, hi = promote(lo, hi)
    LogEdgeVector{typeof(lo)}(lo = lo, hi = hi, nbins = nbins)
end



function LinEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Integer
    @assert (nbins > 0) && (hi > lo)
    lo = floor(Int, lo)
    hi = ceil(Int, hi)
    edges = unique(floor.(Int, range(lo, hi + 1, length = nbins + 1)))
    LinEdgeVector{T}(edges)
end

function LinEdgeVector{T}(; lo = 1, hi = 10, nbins::Integer) where T <: Real
    @assert (nbins > 0) && (hi > lo)
    edges = range(lo, hi, length = nbins + 1)
    LinEdgeVector{T}(edges)
end

function LinEdgeVector(; lo = 1, hi = 10, nbins::Integer)
    lo, hi = promote(lo, hi)
    LinEdgeVector{typeof(lo)}(lo = lo, hi = hi, nbins = nbins)
end



function Base.getindex(r::AbstractEdgeVector{T}, i::Int) where T
    r.edges[i]
end

function Base.getindex(r::AbstractEdgeVector{T}, i::UnitRange{Int}) where T
    r.edges[i]
end

Base.size(r::AbstractEdgeVector{T}) where T = (length(r.edges),)

Base.eltype(r::AbstractEdgeVector{T}) where T = T




"""
    midpoints(r::AbstractEdgeVector{T}) where T

Return the midpoints of the bins in `r`.
"""
function midpoints(r::LogEdgeVector{T}) where T <: Real
    sqrt.(r.edges[1:end-1] .* (r.edges[2:end]))
end

function midpoints(r::LogEdgeVector{T}) where T <: Integer
    sqrt.(r.edges[1:end-1] .* (r.edges[2:end] .- 1))
end

function midpoints(r::LinEdgeVector{T}) where T <: Real
    (r.edges[1:end-1] .+ r.edges[2:end]) ./ 2
end

function midpoints(r::LinEdgeVector{T}) where T <: Integer
    (r.edges[1:end-1] .+ r.edges[2:end] .- 1) ./ 2
end

"""
    fit(::Type{Histogram}, ::Type{T}, vs::AbstractVector;
        nbins = sturges(length(vs)),
        lo = minimum(vs),
        hi = maximum(vs)
    ) where T <: AbstractEdgeVector

Fit a histogram with edges `T` to the data `vs`.

`T` is either `LogEdgeVector` or `LinEdgeVector`. The edges are `Int` if 
the `vs` are `Int` otherwise they are `Float64`. `lo` and `hi` are the lower and upper
limits of the histogram. `nbins` is the number of bins. 

# Examples

```julia
using StatsBase, HistogramBinnings
using Distributions

vs = floor.(Int, rand(Pareto(), 10000000))

h = fit(Histogram, LogEdgeVector, vs)

# or 

h = Histogram(LogEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
append!(h1, vs)
```
"""
function fit(::Type{Histogram}, ::Type{T}, vs::AbstractVector; 
        nbins = sturges(length(vs)), 
        lo = minimum(vs), 
        hi = maximum(vs)
    )                             where T <: AbstractEdgeVector
    e = T(lo = lo, hi = hi, nbins = nbins)
    fit(Histogram, vs, e)
end




end # module HistogramBinnings
