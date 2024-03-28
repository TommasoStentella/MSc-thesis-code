# HistogramBinnings.jl

A package to generate log-binned histograms.


```julia
using StatsBase, HistogramBinnings


data = [1, 2, 10, 1000, 1000]

h1 = Histogram(LogEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
append!(h1, data)

# or

h2 = fit(Histogram, LogEdgeVector, data, nbins=100, lo = 1, hi = 1_000_000)
```

Similarly one may use `LinEdgeVector` for linear binning.


For plotting purposes one may define
```julia
using LinearAlgebra

function xy(h::Histogram{T, 1, E}; mode = :density) where {T, E}
    hn = normalize(h; mode)
    return midpoints(h.edges[1]), hn.weights
end
```

and then plot

```julia
plot(xy(h)...)
```