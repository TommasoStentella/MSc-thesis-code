module EpochModel


using LinearAlgebra, Statistics
using Turing, Optim, StatsBase, Distributions, HistogramBinnings, StatsAPI
using Printf
using Logging
import DynamicPPL, ForwardDiff, Setfield

Logging.disable_logging(Logging.Warn);

include("fitting.jl")
include("mathematica-derived.jl")


export 
    hid,
	FitResult,	
	fit_epochs

# Computing

function secondderivative(f, x)
    dfdx = x -> ForwardDiff.derivative(f, x)
    ForwardDiff.derivative(dfdx, x)
end


function laplace_n(TN::Vector, s::Number)
    N = TN[2]
    y = 2 * N^2 / (1 + 2*N*s)
    # stationary solution in first epoch
    for k in 3:2:length(TN) # loop over further epochs
        T  = TN[k]
        Np = TN[k-1]
        N = TN[k+1]
        # step up or down
        gamma = N / Np
        w1 = gamma >= 1 ? gamma^2 - (gamma^2 - 1)/(2 * Np) : gamma^2
        w2 = gamma >= 1 ? (gamma^2 - 1) * Np               : zero(gamma)
        # propagate in time
        v1 = exp((-T/(2*N)) - s*T)
        v2 = (1 - v1) * ((2*N^2) / (1 + 2*N*s))
        # update value
        y = (w1 * y + w2) * v1 + v2
    end
    y
end


function laplace_n(Nv::Vector, Tv::Vector, s::Number)
	# T =        [T1, T2, ...]
	# N = [Nstat, N1, N2, ...]
    Nstat = Nv[1]
    y = 2 * Nstat^2 / (1 + 2*Nstat*s)
    # stationary solution in first epoch
    Np = Nstat
    for (T, N) in zip(Tv, Iterators.drop(Nv, 1)) # loop over further epochs
        # T  = Tv[k-2]
        # Np = Nv[k-1]
        # N =  Nv[k]
        # step up or down
        gamma = N / Np
        w1 = gamma >= 1 ? gamma^2 - (gamma^2 - 1)/(2 * Np) : gamma^2
        w2 = gamma >= 1 ? (gamma^2 - 1) * Np               : zero(gamma)
        # propagate in time
        v1 = exp((-T/(2*N)) - s*T)
        v2 = (1 - v1) * ((2*N^2) / (1 + 2*N*s))
        # update value
        y = (w1 * y + w2) * v1 + v2
        Np = N
    end
    y
end


function hid(TN::Vector, mu::Float64, r::Number)
	# TN = [L, N0, T1, N1, T2, N2, ...]
	L = TN[1]
	N = TN[end]
	(2*mu^2*L)/(N^2) * secondderivative(s -> laplace_n(TN, s), 2*mu*r)
	# prefactor        pure bliss
end

function hid(L::Number, N::Vector, T::Vector, mu::Float64, r::Number)
	# T = [T1, T2, ...]
	# N = [Nstat, N1, N2, ...]
	Nend = (length(N) == 0 ? Nstat : N[end])
	(2*mu^2*L)/(Nend^2) * secondderivative(s -> laplace_n(N, T, s), 2*mu*r)
	# prefactor           pure bliss
end

function hid_integral(TN::Vector, mu::Float64, r::Number)
    # integral of hid
	# TN = [L, N0, T1, N1, T2, N2, ...]
	L = TN[1]
	N = TN[end]
	(mu*L)/(N^2) * ForwardDiff.derivative(s -> EpochModel.laplace_n(TN, s), 2*mu*r)
	# prefactor    pure bliss
end


end # module EpochModel
