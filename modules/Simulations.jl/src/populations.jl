

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

# struct RecombinationMap{Tl <: Integer}
#     λ::Vector{Poisson{Float32}}
#     ran::Vector{StepRange{Tl,Tl}}
#     low_h::Vector{Tl}
#     lim_h::Vector{Tl}
#     low_c::Vector{Tl}
#     lim_c::Vector{Tl}
# end

mutable struct Population{Ti <: Integer, Tl <: Integer, Tt <: Real, Tn <: Function}
    ploidy::Ti
    alive::Vector{Individual{Ti}}
    time::Tt
    demography::Demography{Tn}
    # rmap::RecombinationMap{Tl}
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

# function map_break(breakpoint; low=[], lim=[])
#     index = searchsortedfirst(lim, breakpoint)
#     return low[index] + breakpoint - 1
# end


function simulate(TN::Vector{Int}, ρ=1e-8; indv=1)
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

    pop = Population{Ti,Tl,Tt,Function}(ploidy, alive, 0, demo)

    model = DiscreteGensDiploid(ρ = ρ)
    evolve!(model, pop, dt, pool)

    getIBDs(pop.alive[indv].ids[1], pop.alive[indv].ids[2], pop.time, pool)
end