

abstract type AbstractPopSimModel end

@kwdef struct DiscreteGensDiploid <: AbstractPopSimModel 
    ρ::Float64 = 1e-8
end

ploidy(model::DiscreteGensDiploid) = 2


randomswap(a, b) = rand(Bernoulli()) ? (a,b) : (b,a)



function evolve!(model::AbstractPopSimModel, pop::Population, dt::Real, pool::ChromosomePool)
    pop.ploidy == ploidy(model) || error("Population ploidy does not match model ploidy")
    tend = pop.time + dt

    while(pop.time < tend)
        evolve_step!(model, pop, pool)
    end
    pop
end

function evolve_step!(model::DiscreteGensDiploid, pop::Population, pool::ChromosomePool)
    Ti = typeof(pop.alive[1].ids[1])
    N = pop_size(pop)
    L = chromosome_length(pool)
    # Tl = typeof(L)
    dt = 1
    pop.time += dt
    # Tt = typeof(pop.time)

    Nnew = offspring_size(pop)

    k = length(pool.chromosomes)
    resize!(pool.chromosomes, k + 2*Nnew)

    for i in 1:Nnew
        id1 = k + 2*i - 1
        id2 = k + 2*i

        p1, p2 = sample(1:N, 2, replace = false)
        
        pid1, pid2 = randomswap(pop.alive[p1].ids[1], pop.alive[p1].ids[2])
        pool.chromosomes[id1] = crossover(Ti(id1), pid1, pid2, L , pop.time, pop.rmap)
        
        pid1, pid2 = randomswap(pop.alive[p2].ids[1], pop.alive[p2].ids[2])
        pool.chromosomes[id2] = crossover(Ti(id2), pid1, pid2, L , pop.time, pop.rmap)
    end

    pop.alive = map(1:Nnew) do i
        Individual{Ti}(collect(k+(i-1)*2+1 : k+i*2))
    end
end