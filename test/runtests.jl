using ZernikePolynomials
using Test

@testset "ZernikePolynomials.jl" begin

    @testset "Sequential index conversion" begin
        @test Noll(NM(0,0)) === Noll(1)
        @test Noll(NM(3,-1)) === Noll(7)
        @test_throws ArgumentError NM(0,-1)
        @test Noll(NM(5,3)) === Noll(18)

        @test OSA(NM(0,0)) === OSA(0)
        @test OSA(NM(3,-1)) === OSA(7)
        @test OSA(NM(5,3)) === OSA(19)

        @test all([OSA(NM(OSA(i))) for i in 0:30] .== OSA.(0:30))
        @test all([Noll(NM(Noll(i))) for i in 1:31] .== Noll.(1:31))

        @test convert(NM, OSA(6)) === NM(OSA(6)) === NM(3, -3)
    end

    @testset "Type stability" begin
        Z = zernike(NM(4,4),coord=:cartesian)
        @test Z(0.2, 0.1) ≈ -0.002213594362117866
        @test Z(0.2f0, 0.1f0) ≈ -0.0022135945f0
        @test typeof(Z(0.2f0, 0.1f0)) == Float32
        @test typeof(Z(0.2, 0.1)) == Float64

        @test typeof(ZernikePolynomials.R(Float32, NM(1,1))(1)) == Float32
        @test typeof(ZernikePolynomials.R(Float32, NM(1,1))(1.0)) == Float64
        @test typeof(ZernikePolynomials.R(Float64, NM(1,1))(1.0)) == Float64
        @test typeof(ZernikePolynomials.R(Float64, NM(1,1))(1.0f0)) == Float64
        @test typeof(ZernikePolynomials.R(NM(1,1))(1.0)) == Float64

        @test typeof(normalization(Float32, NM(1,1))) == Float32
        @test typeof(normalization(ComplexF32, NM(1,1))) == ComplexF32
        @test typeof(normalization(NM(1,1))) == Float64
    end

    @testset "Generation of Zernike polynomials" begin

        Z = zernike(NM(1,1),coord=:cartesian)
        @test Z isa Function

        Z = zernike(NM(1,1),coord=:polar)
        @test Z isa Function

        x = LinRange(-1,1,64)
        ϕi(i) = evaluatezernike(x,OSA(i),1.)
        r = 1:10
        N = sum(evaluatezernike(x,OSA(0),1.))
        X = round.([sum(ϕi(i) .* ϕi(j))/N for i in r, j in r])
        for i in 1:10, j in 1:10
            if i == j
                @test X[i,j] ≈ 1.0
            else
                @test X[i,j] == 0.0
            end
        end

        # coefficient estimation
        ϕ = evaluatezernike(64, OSA.([2, 3]), [0.5, 0.3])
        @test zernikecoefficients(ϕ, OSA.([1, 2, 3])) ≈ [0.0, 0.5,0.3]
    end
end
