

function crossover(id::Ti, pid1::Ti, pid2::Ti, L::Tl, time::Tt, rmap::RecombinationMap{Tl}) where {Ti <: Integer, Tl <: Integer, Tt <: Real}
    # breaks = n_breaks < L ? sample(1:L, n_breaks, replace = false) : collect(1:L)
    # breaks = rand(1:L, n_breaks)
    n_breaks = rand.(rmap.λ)
    breaks = rand.(rmap.ran, n_breaks)
    breaks_h = map_break.(breaks[1], low=rmap.low_h, lim=rmap.lim_h)
    breaks_c = map_break.(breaks[2], low=rmap.low_c, lim=rmap.lim_c)
    breaks = sort!(vcat(breaks_h, breaks_c))
    @assert pid1 != pid2
    Chromosome(id, L, Tl.(breaks), pid1, pid2, time)
end


function getparentat(c::Chromosome, pos::Integer; find_next::Bool=true)
    @assert 1 <= pos <= c.length
    i = searchsortedfirst(c.breakpoints, pos)
    pid = i % 2 == 1 ? c.pid1 : c.pid2
    if find_next
        nextbreak = i > length(c.breakpoints) ? c.length : c.breakpoints[i]
    else
        nextbreak = pos
    end
    pid, nextbreak
end






struct Segments{Ti <: Integer, Tl <: Integer, Tt}
    lengths::Vector{Tl}
    times::Vector{Tt}
    ids::Vector{Ti}
    type::Symbol
    para::NamedTuple
end

function Base.show(io::IO, s::Segments)
    print(io, length(s.lengths), " segments (type:", s.type, ", total length:", sum(s.lengths), ")")
end




function getIBDs(id1::Integer, id2::Integer, time::Tt, pool::ChromosomePool) where {Tt <: Real} 
    L = chromosome_length(pool)
    Tl = typeof(L)
    Ti = typeof(pool.chromosomes[1].id)
    @assert id1 != id2
    c1 = pool.chromosomes[id1]
    c2 = pool.chromosomes[id2]

    lengths = Vector{Tl}()
    ids = Vector{Ti}()
    times = Vector{Tt}()

    lpos::Tl = 1
    pos::Tl = 1
    while pos <= L
        pid1, nextbreak1 = c1.id, L
        pid2, nextbreak2 = c2.id, L

        while ((pid1 > 0) || (pid2 > 0)) && (pid1 != pid2)
            # @show pid1, pid2
            if pid1 < pid2
                pid2, nextbreak = getparentat(pool.chromosomes[pid2], pos)
                nextbreak2 = min(nextbreak2, nextbreak)
            else
                pid1, nextbreak = getparentat(pool.chromosomes[pid1], pos)
                nextbreak1 = min(nextbreak1, nextbreak)
            end
        end

        
        pid = pid1 == pid2 ? pid1 : zero(Ti)
        pos = min(nextbreak1, nextbreak2) + 1
        
        dl = pos - lpos
        dt = pid1 > 0 ? time - pool.chromosomes[pid1].time : zero(Tt)
        push!(lengths, dl)
        push!(times, dt)
        push!(ids, pid)
        lpos = pos
    end

    Segments(lengths, times, ids, :IBD, (;id1, id2, time))
end


function get_paths(row_size::Integer, s::Segments, pool::ChromosomePool; skip::Integer=0)
    Ti = Int32

    k = 0
    i = 1
    # look for a sequence of row_size consecutive segments with same ancestor
    while i <= length(s.ids) && (k < skip + 1)
        j = i
        while (j <= length(s.ids)) && (s.ids[j] == s.ids[i]) && (j-i < row_size)
            j += 1
        end
        if j-i == row_size
            k += 1
        end
        i = j
    end

    if k < skip + 1
        return nothing
    end

    # j now points to the segment right after the row
    j = j - 1
    pos = map(1:row_size) do i
        sum(s.lengths[1 : j + i - row_size - 1]) + 1
    end

    c1 = pool.chromosomes[s.para.id1]
    c2 = pool.chromosomes[s.para.id2]

    paths = Vector{Vector{Vector{Ti}}}() # pairs of paths (vectors of ids), one pair for each position
    for p in pos
        pid1 = c1.id
        pid2 = c2.id

        path1 = Vector{Ti}()
        path2 = Vector{Ti}()
        push!(path1, pid1)
        push!(path2, pid2)
        while ((pid1 > 0) || (pid2 > 0)) && (pid1 != pid2)
            if pid1 < pid2
                pid2, _ = getparentat(pool.chromosomes[pid2], p; find_next=false)
                push!(path2, pid2)
            else
                pid1, _ = getparentat(pool.chromosomes[pid1], p; find_next=false)
                push!(path1, pid1)
            end
        end
        push!(paths, [path1, path2])
    end
    push!(pos, sum(s.lengths[1:j])+1)
    return paths, pos
end


function isfullycoalesced(s::Segments)
    all(s.ids .!= 0)
end

function noncoalesced_length(s::Segments)
    sum(s.lengths[s.ids .== 0])
end


function joinIBDs(s::Segments; by_time=false, shuffle::Bool = false, seed::Int64 = 1234)
    lengths_c = copy(s.lengths)
    times_c = copy(s.times)
    ids_c = copy(s.ids)
    lengths = similar(s.lengths)
    times = similar(s.times)
    ids = similar(s.ids)
    if shuffle
        Random.seed!(seed)
        Random.shuffle!(ids_c)
        Random.shuffle!(lengths_c)
        Random.shuffle!(times_c)
        @assert any(lengths_c .!= s.lengths)
    end
    k = 1
    i = 1
    n_row = Vector{Int32}()
    # y_corr = Vector{Int32}()
    # x_corr = Vector{Int32}()
    while i <= length(ids_c)
        j = i
        lsum = zero(eltype(lengths))
        if by_time
            while (j <= length(ids_c)) && (times_c[j] == times_c[i])
                lsum += lengths_c[j]
                # if j-i > 0
                #     push!(y_corr, lengths_c[j])
                #     push!(x_corr, lengths_c[j-1])
                # end
                j += 1
            end
        else
            while (j <= length(ids_c)) && (ids_c[j] == ids_c[i])
                lsum += lengths_c[j]
                # if j-i > 0
                #     push!(y_corr, lengths_c[j])
                #     push!(x_corr, lengths_c[j-1])
                # end
                j += 1
            end
        end
        lengths[k] = lsum
        times[k] = times_c[i]
        ids[k] = ids_c[i]
        push!(n_row, j-i)
        i = j
        k += 1
    end
    resize!(lengths, k-1)
    resize!(times, k-1)
    resize!(ids, k-1)
    Segments(lengths, times, ids, :IBA, s.para), n_row     # (x_corr, y_corr)
end

function getIBSs(s::Segments, μ::Float64; allow_multiple_hits::Bool = true)
    # @assert sum(s.ids .== 0) == 0 "not coalescent segments present"
    total = sum(s.lengths)
    nmut = map(1:length(s.lengths)) do i
        l = s.lengths[i]
        t = s.times[i]
        n = rand(Poisson(2 * μ * l * t))
        min(n, l)
    end
    m = max(sum(nmut), 1)
    nl = similar(s.lengths, m)
    nt = similar(s.times, m)
    nid = similar(s.ids, m)
    k = 1
    # pos = 1
    lpos::Int32 = 0
    lpos_is_set = false
    L = 0 
    for i in 1:length(s.lengths)
        l = s.lengths[i]
        t = s.times[i]
        id = s.ids[i]
        n = nmut[i]
        if n > 0
            if allow_multiple_hits
                ps = rand(1:l, n)
                sort!(ps)
                unique!(ps)
            else
                ps = n >= l ? 
                    collect(1:l) : 
                    sort(sample(1:l, n, replace = false))
            end

            for psi in ps 
                if lpos_is_set
                    dl = psi - lpos
                    L += dl
                    nl[k] = dl
                    nt[k] = t
                    nid[k] = id
                    k += 1
                end
                lpos = psi
                lpos_is_set = true
            end
        end
        lpos -= l
    end
    @assert k <= m
    @assert allow_multiple_hits || (k == m)
    nl[k] = total - L
    nt[k] = s.times[end]
    nid[k] = s.ids[end]
    resize!(nl, k)
    resize!(nt, k)
    resize!(nid, k)
    Segments(nl, nt, nid, :IBS, merge(s.para, (;μ)))
end
    

Base.vcat(s1::Segments, s2::Segments) = Segments(
    vcat(s1.lengths, s2.lengths),
    vcat(s1.times, s2.times),
    vcat(s1.ids, s2.ids),
    s1.type,
    s1.para)



