using ZernikePolynomials
using Test

@testset "ZernikePolynomials.jl" begin

    @testset "Sequential index conversion" begin
        @test mn2Noll(0,0) == 1
        @test mn2Noll(-1,3) == 7
        @test isnan(mn2Noll(-1,0))
        @test mn2Noll(3,5) == 18

        @test mn2OSA(0,0) == 0
        @test mn2OSA(-1,3) == 7
        @test isnan(mn2OSA(-1,0))
        @test mn2OSA(3,5) == 19

        @test all([mn2OSA(OSA2mn(i)...) for i in 0:30] .== collect(0:30))
        @test all([mn2Noll(Noll2mn(i)...) for i in 1:31] .== collect(1:31))
    end

    @testset "Generation of Zernike polynomials" begin

        Z = Zernike(1,1,coord=:cartesian)
        @test Z isa Function

        Z = Zernike(1,1,coord=:polar)
        @test Z isa Function

        x = LinRange(-1,1,64)
        ϕi(i) = evaluateZernike(x,i,1.)
        r = 1:10
        N = sum(evaluateZernike(x,0,1.))
        X = round.([sum(ϕi(i) .* ϕi(j))/N for i in r, j in r])
        for i in 1:10, j in 1:10
            if i == j
                @test X[i,j] ≈ 1.0
            else
                @test X[i,j] == 0.0
            end
        end

        # coefficient estimation
        ϕ = evaluateZernike(64, [2, 3], [0.5, 0.3], index=:OSA)
        @test Zernikecoefficients(ϕ, [1, 2, 3]) ≈ [0.0, 0.5,0.3]
    end
end
