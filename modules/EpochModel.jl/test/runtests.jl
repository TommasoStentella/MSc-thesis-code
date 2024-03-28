using Test
using Distributions
using Random
using Optim
using EpochModel
using Logging
Logging.disable_logging(Logging.Warn)


@testset "Compare to Mathematica" begin
        st = map(1:1000) do i
                L = rand(1e6:1e9)
                N0 = N1 = N2 = -1.0
                while true
                        N0 = rand(1e3:1e5)
                        N1 = rand(1e3:1e5)
                        N2 = rand(1e3:1e5)
                        if (N0 < N1 > N2) || (N0 > N1 < N2)
                                break
                        end
                end
                T1 = rand(1e1:1e3)
                T2 = rand(1e1:1e3)
                mu = rand(1e-9:1e-8)
                r = rand(1:1_000_000)

                y1 = hid([L, N0], mu, r)
                y2 = hid([L, N0, T1, N1], mu, r)
                y3 = hid([L, N0, T1, N1, T2, N2], mu, r)

                y1m = EpochModel.hidm(L, N0, mu, r)
                y2m = EpochModel.hidm(L, N0, T1, N1, mu, r)
                y21m = EpochModel.hidm(L, N0, T1, N0, mu, r)
                y3m = EpochModel.hidm(L, N0, T1, N1, T2, N2, mu, r)
                y31m = EpochModel.hidm(L, N0, T1, N0, T2, N0, mu, r)

                @test abs(y1-y1m) < 1e-10
                @test abs(y1-y21m) < 1e-10 
                @test abs(y1-y31m) < 1e-10
                @test abs(y2-y2m) < 1e-10
                @test abs(y3-y3m) < 1e-10

                (;L, N0, N1, N2, T1, T2, mu, r, y1, y2, y3, y1m, y2m, y3m, 
                        d1 = y1-y1m, d2 = y2-y2m, d3 = y3-y3m,
                        d1r = abs(y1-y1m)/y1m, d2r = abs(y2-y2m)/y2m, d3r = abs(y3-y3m)/y3m)
        end
        @show maximum(abs.(getindex.(st, :d1)))
        @show maximum(abs.(getindex.(st, :d2)))
        @show maximum(abs.(getindex.(st, :d3)))
end

@testset "Fitting 1 epoch" begin
        Random.seed!(1234)
        L = 1e9
        N = 1e4
        mu = 1e-8

        r = collect(1:1_000_000)
        c = map(r->rand(Poisson(hid([L, N], mu, r))), r)

        @time f = fit_epochs_mids(r, c, mu, nepochs = 1)
        @test maximum(abs.(f.para .- [L,N])./[L,N]) < 0.01 

        # @time f = fit_epochs(r, c, mu, nepochs = 1, bin = 1000)
        # @test maximum(abs.(f.para .- [L,N])./[L,N]) < 0.01 

        # @time f = fit_epochs(r, c, mu, nepochs = 1, bin = 100)
        # @test maximum(abs.(f.para .- [L,N])./[L,N]) < 0.01 

        # @time f = fit_epochs(r, c, mu, nepochs = 1, bin = 100, space = :log)
        # @test maximum(abs.(f.para .- [L,N])./[L,N]) < 0.05 
end

@testset "Fitting 3 epochs" begin
        Random.seed!(1234)
        TN = [1e9, 1e3, 100.0, 100, 200, 1000] 
        mu = 1e-8

        r = collect(1:1_000_000)
        c = map(r->rand(Poisson(hid(TN, mu, r))), r)

        # @time f1 = fit_epochs(r, c, mu, nepochs = 1, bin = 100)
        # @show f1
        # # @time f3a = fit_epochs(r, c, mu, nepochs = 3, bin = 100)
        # @show f3a
        # @time f3b = fit_epochs(r, c, mu, nepochs = 3, bin = 1000, init = f3a.para, solver = NelderMead())
        # @show f3b
        # @time f3c = fit_epochs(r, c, mu, nepochs = 3, bin = 1000, init = f3a.para, solver = ParticleSwarm())
        # @show f3c

        # @show f3a.para_name
        @test 1== 1
end