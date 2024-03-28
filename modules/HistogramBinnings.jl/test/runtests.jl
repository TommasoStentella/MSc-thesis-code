using Test
using LinearAlgebra
using StatsBase
using HistogramBinnings


@testset "Edges" begin
    e = LogEdgeVector([1,3,10])
    @test size(e)[1] == 3
    @test eltype(e) == Int
    @test e[1] == 1
    @test midpoints(e) == [sqrt(2), sqrt(27)]
    

    e = LogEdgeVector([1.0,3,10])
    @test size(e)[1] == 3
    @test eltype(e) == Float64
    @test e[1] == 1
    @test midpoints(e) == [sqrt(3), sqrt(30)]
    

    e = LinEdgeVector([1,3,10])
    @test size(e)[1] == 3
    @test eltype(e) == Int
    @test e[1] == 1
    @test midpoints(e) == [1.5, 6]


    e = LinEdgeVector([1.0,3,10])
    @test size(e)[1] == 3
    @test eltype(e) == Float64
    @test e[1] == 1
    @test midpoints(e) == [2.0, 6.5]


    e = LinEdgeVector([1, 10, 100])
    for i in e
        @test i isa Int
    end
    @test length(collect(e)) == 3
end



@testset "Edge generation lin" begin
    e = LinEdgeVector{Int}(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 11

    e = LinEdgeVector{Float64}(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 101

    e = LinEdgeVector(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 11
    @test eltype(e) == Int

    e = LinEdgeVector(lo = 1, hi = 10.0, nbins = 100)
    @test size(e)[1] == 101
    @test eltype(e) == Float64

end



@testset "Edge generation log" begin
    e = LogEdgeVector{Int}(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 11

    e = LogEdgeVector{Float64}(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 101

    e = LogEdgeVector(lo = 1, hi = 10, nbins = 100)
    @test size(e)[1] == 11
    @test eltype(e) == Int

    e = LogEdgeVector(lo = 1, hi = 10.0, nbins = 100)
    @test size(e)[1] == 101
    @test eltype(e) == Float64

end


@testset "Histogram" begin
    data = [1,2,10,100,1000]
    
    h = Histogram(LogEdgeVector{Int}(lo = 1, hi = 1_000_000, nbins = 60))
    append!(h, data)
    @test sum(h.weights) == length(data)

    h = Histogram(LogEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
    append!(h, data)
    @test sum(h.weights) == length(data)
    push!(h, 1_000_000)
    @test sum(h.weights) == length(data) + 1

    push!(h, 1, 6)
    @test sum(h.weights) == length(data) + 1 + 6
    @test h.weights[1] == 7

    append!(h, [1,2,3], [10, 10, 10])
    @test sum(h.weights) == length(data) + 1 + 6 + 30
    @test h.weights[1] == 17
    @test h.weights[2] == 11


    h = Histogram(LinEdgeVector{Int}(lo = 1, hi = 1_000_000, nbins = 60))
    append!(h, data)
    @test sum(h.weights) == length(data)

    h = Histogram(LinEdgeVector(lo = 1, hi = 1_000_000, nbins = 60))
    append!(h, data)
    @test sum(h.weights) == length(data)
    push!(h, 1_000_000)
    @test sum(h.weights) == length(data) + 1

    @test StatsBase.binvolume(h, 2) == 16667
end


@testset "fit" begin
    data = [1,2,10,100,1000]

    # StatsBase way to do things
    h = fit(Histogram, data)
    @test sum(h.weights) == length(data)

    h = fit(Histogram, data, LogEdgeVector{Int}([1,2,10]))
    @test h.weights == [1,1]
    @test sum(h.weights) == sum(<(10), data)



    # HistogramBinnings way to do things
    h = fit(Histogram, LogEdgeVector{Int}, data)
    @test sum(h.weights) == length(data)

    h = fit(Histogram, LogEdgeVector, data)
    @test sum(h.weights) == length(data)

    h = fit(Histogram, LogEdgeVector, data, nbins=100, hi = 1_000_000)
    @test sum(h.weights) == length(data)
end

@testset "xy" begin
    function xy(h::Histogram{T, 1, E}; mode = :density) where {T, E}
        hn = normalize(h; mode)
        return midpoints(h.edges[1]), hn.weights
    end

    data = [1,2,10,2,1]
    h = fit(Histogram, LogEdgeVector{Int}, data, nbins = 10)

    x, y = xy(h)

    @test x[1] == 1
    @test x[2] == 2
end