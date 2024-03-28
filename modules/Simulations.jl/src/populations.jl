

struct Chromosome{Ti <: Integer, Tl <: Integer, Tt}
    id::Ti
    length::Tl
    breakpoints::Vector{Tl}
    pid1::Ti
    pid2::Ti
    time::Tt
end


struct Individual{Ti <: Integer}
    ids::Vector{Ti}
end

Base.show(io::IO, indv::Individual) = 
    print(io, "Individual(", indv.ids, ")")


mutable struct ChromosomePool{Ti <: Integer, Tl <: Integer, Tt}
    chromosomes::Vector{Chromosome{Ti,Tl,Tt}}
end

Base.show(io::IO, pool::ChromosomePool) = 
    print(io, "ChromosomePool(", length(pool.chromosomes), " chromosomes, ", chromosome_length(pool), " bp)")

# ChromosomePool() = ChromosomePool{Int32, Int32, Float64}(Chromosome{Int32, Int32, Float64}[], nothing)


function Base.show(io::IO, c::Chromosome) 
    print(io, "Chromosome(", c.id, ", ", c.length, " bp, ")
    if 0 < length(c.breakpoints) <= 10
        print(io, "breakpoints: ", c.breakpoints, ", ")
    else
        print(io, length(c.breakpoints), " breakpoints, ")
    end
    print(io, "pids: ", c.pid1, ", ", c.pid2, ", birth at: ", c.time, ")")
end


mutable struct Demography{Tn <: Function}
    population_size::Tn
end

struct RecombinationMap{Tl <: Integer}
    λ::Vector{Poisson{Float32}}
    ran::Vector{StepRange{Tl,Tl}}
    low_h::Vector{Tl}
    lim_h::Vector{Tl}
    low_c::Vector{Tl}
    lim_c::Vector{Tl}
end

mutable struct Population{Ti <: Integer, Tl <: Integer, Tt <: Real, Tn <: Function}
    ploidy::Ti
    alive::Vector{Individual{Ti}}
    time::Tt
    demography::Demography{Tn}
    rmap::RecombinationMap{Tl}
end


Base.show(io::IO, p::Population) = 
    print(io, "Population(", p.ploidy, "-ploid, ", pop_size(p), " individuals)")


function pop_size(pop::Population)
    length(pop.alive)
end


function offspring_size(pop::Population)
    Int(round(pop.demography.population_size(pop.time)))
end

function chromosome_length(pool::ChromosomePool)
    length(pool.chromosomes) == 0 ? 0 : pool.chromosomes[1].length
end

get_chromosomes_ids(pop::Population, indv::Integer) = 
    pop.alive[indv].ids

function map_break(breakpoint; low=[], lim=[])
    index = searchsortedfirst(lim, breakpoint)
    return low[index] + breakpoint - 1
end


# function Population(pool; L::Int = 1000, N::Int = 1000)
#     Ti = Int32
#     Tl = Int64
#     Tt = Float64
#     ploidy = 2
#     chromosomes = map(1:ploidy * N) do i
#         Chromosome{Ti, Tl, Tt}(newID(pool), L, Int32[], 0, 0, 0.0)
#     end
#     alive = map(1:N) do i
#         Individual{Ti}(collect((i-1)*ploidy+1:i*ploidy))
#     end

#     # default is WF model
#     size(x) = N
#     demo = Demography{Function}(size)

#     k = length(pool.chromosomes)
#     resize!(pool.chromosomes, k + ploidy*N)
#     pool.chromosomes[k+1:k+ploidy*N] = chromosomes

#     Population{Ti,Tt,Function}(ploidy, alive, 0.0, demo)
# end

function simulate(TN::Vector{Int}, ρ=1e-8, μ=2.36e-8; indv=1, f_h=0.8, f_c=0.2)
    L = TN[1]
    N = TN[2]
    Ts = TN[3:2:end-1]
    Ns = TN[4:2:end]
    dt = 20*N + sum(Ts) - 1

    total_ind = 20*N^2 + sum(Ts .* Ns)
    Ti = total_ind > 2^31 ? Int64 : Int32
    Tl = L > 2^31 ? Int64 : Int32
    Tt = dt > 2^31 ? Int64 : Int32

    ploidy = 2
    pool = ChromosomePool{Ti,Tl,Tt}(Vector{Chromosome{Ti,Tl,Tt}}([]))
    
    chromosomes = map(1:ploidy * N) do i
        Chromosome{Ti, Tl, Tt}(i, L, Tl[], 0, 0, 0)
    end
    k = length(pool.chromosomes)
    resize!(pool.chromosomes, k + ploidy*N)
    pool.chromosomes[k+1:k+ploidy*N] = chromosomes

    alive = map(1:N) do i
        Individual{Ti}(collect((i-1)*ploidy+1:i*ploidy))
    end

    epochs(g) = sum([[N * (g <= 20*N)];[Ns[i] * (g > sum(Ts[1:i-1])+20*N && g <= sum(Ts[1:i])+20*N) for i in 1:length(Ts)]])
    demo = Demography{Function}(epochs)

    # rmap = [(0.05,Tl(1),Tl(10_000_000)),
    #         (0.2,Tl(10_000_001),Tl(11_000_000)),
    #         (0.05,Tl(11_000_000),Tl(50_000_000)),
    #         (0.3,Tl(50_000_001),Tl(52_000_000)),
    #         (0.05,Tl(52_000_001),Tl(100_000_000)),
    #         (0.3,Tl(100_000_001),Tl(105_000_000)),
    #         (0.05,Tl(105_000_000),Tl(250_000_000))]

    n_hotspots = 1000
    hot_pos = sort!(rand(1:L, n_hotspots))
    mean_width = L * f_c / n_hotspots
    widths = rand(Poisson(mean_width), n_hotspots)

    low_h = []
    up_h = []
    for i in 1:n_hotspots
        l = hot_pos[i]
        u = hot_pos[i] + widths[i]
        if u > L
            u = L
            l = L - widths[i]
        end
        push!(low_h, l)
        push!(up_h, u)
    end
    low_c = []
    up_c = []
    for i in 1:n_hotspots
        previous = i > 1 ? hot_pos[i-1] + widths[i-1] + 1 : 1
        next = hot_pos[i] - 1
        if next - previous > 0
            push!(low_c, previous)
            push!(up_c, next)
        end
    end
    if hot_pos[end] + widths[end] < L
        push!(low_c, hot_pos[end] + widths[end] + 1)
        push!(up_c, L)
    end
    pois = Poisson.(Float32.(L * ρ * [f_h, f_c]))
    ran_h = Tl(1):Tl(sum(up_h .- low_h .+ 1))
    ran_c = Tl(1):Tl(sum(up_c .- low_c .+ 1))
    ran = [ran_h, ran_c]
    rmap = RecombinationMap{Tl}(pois, ran, low_h, accumulate(+, up_h .- low_h .+ 1), low_c, accumulate(+, up_c .- low_c .+ 1))

    pop = Population{Ti,Tl,Tt,Function}(ploidy, alive, 0, demo, rmap)

    model = DiscreteGensDiploid(ρ = ρ)
    evolve!(model, pop, dt, pool)

    ibd = getIBDs(pop.alive[indv].ids[1], pop.alive[indv].ids[2], pop.time, pool)

    getIBSs(ibd, μ)
end