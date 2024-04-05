function xy(h::HistogramBinnings.Histogram{T, 1, E}; mode = :density) where {T, E}
    hn = StatsBase.normalize(h; mode)
    return midpoints(h.edges[1]), hn.weights
end

function plot_demography(epochs, fits, ax; tail=1000, id="")
    old_t = maximum([sum(fit.para[3:2:end-1]) for fit in fits]) + tail
    for (fit,nepochs) in zip(fits,epochs)
        TN = fit.para[end:-1:2]
        stdTN = fit.opt.stderrors[end:-1:2]

        Polygon = matplotlib.patches.Polygon

        mean_size = []
        upp_size = []
        low_size = []
        for (n,sn) in zip(TN[1:2:end], stdTN[1:2:end])
            append!(mean_size, [n,n])
            append!(upp_size, [n+sn,n+sn])
            append!(low_size, [n-sn,n-sn])
        end

        mean_epochs = [0.]
        upp_epochs = [0.]
        low_epochs = [0.]
        for i in 1:nepochs-1
            t = sum(TN[2:2:end-1][1:i])
            st = stdTN[2:2:end-1][i]
            append!(mean_epochs, [t,t])
            if (TN[1:2:end][i] + stdTN[1:2:end][i]) > (TN[1:2:end][i+1] + stdTN[1:2:end][i+1])
                append!(upp_epochs, [t+st,t+st])
            else
                append!(upp_epochs, [t-st,t-st])
            end
            if (TN[1:2:end][i] - stdTN[1:2:end][i]) < (TN[1:2:end][i+1] - stdTN[1:2:end][i+1])
                append!(low_epochs, [t+st,t+st])
            else
                append!(low_epochs, [t-st,t-st])
            end
        end
        push!(mean_epochs, old_t)
        push!(upp_epochs, old_t)
        push!(low_epochs, old_t)

        c = "tab:" .* split("blue orange red purple olive brown cyan pink")[nepochs]

        err = Polygon(collect(zip([upp_epochs;low_epochs[end:-1:1]],[upp_size;low_size[end:-1:1]])),facecolor=c,edgecolor="none",alpha=0.5)

        ax.plot(mean_epochs, mean_size, linewidth=1, label="$nepochs epochs ($id)", color = c)
        ax.add_patch(err)
    end
end

function plot_input(TN, ax; tail=1000)
    old_t = sum(TN[3:2:end-1]) + tail
    nepochs = Int(length(TN)/2)
    TN = TN[end:-1:2]

    mean_size = []
    for n in TN[1:2:end]
        append!(mean_size, [n,n])
    end

    mean_epochs = [0.]
    for i in 1:nepochs-1
        t = sum(TN[2:2:end-1][1:i])
        append!(mean_epochs, [t,t])
    end
    push!(mean_epochs, old_t)

    ax.plot(mean_epochs, mean_size, linewidth=1, label="input demography", color = "black")
end

function tmrca(τ, TN, ρ)
    L = TN[1]
    N = TN[2]
    Ts = TN[3:2:end-1]
    Ns = TN[4:2:end]
    dt = 20*N + sum(Ts) - 1

    popSize(g) = sum([[N * (g <= 20*N)];[Ns[i] * (g > sum(Ts[1:i-1])+20*N && g <= sum(Ts[1:i])+20*N) for i in 1:length(Ts)]])

    function discrete_forward_coalescent(τ;N=N,dt=dt)
        # assuming stationary initial condition
        f(u) = 1 / (2 * popSize(u))
        return f(dt-τ) * exp(-sum([f(dt-i) for i in 0:τ]))
    end

    tmrca = (1 + 2ρ*τ*L) * discrete_forward_coalescent(τ)

    return tmrca
end

function setXrange(ax,inf,sup)
    ax.set_xlim(inf-1,sup)
    ticks = ax.get_xticks()
    ticks = ticks[(ticks .>= inf).&&(ticks .<= sup)]
    ax.set_xticks(ticks)
    ticks = ax.get_xticks(minor=true)
    ticks = ticks[(ticks .>= inf).&&(ticks .<= sup)]
    ax.set_xticks(ticks,minor=true)
    ax
end