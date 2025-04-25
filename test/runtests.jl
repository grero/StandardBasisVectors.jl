using StandardBasisVectors
using Test


@testset "Basic" begin

    for T in [Float32, Float64]
        # convert a simple left-handed 2D system to a right-handed system
        Q = [-one(T) zero(T);
             zero(T) one(T)]

        V,signs,θ = StandardBasisVectors.orient_vectors(Q)

        @test signs == [-one(T), one(T)]
        @test V[:,1] ≈ [one(T),zero(T)]
        @test V[:,2] ≈ [zero(T),one(T)]

         Q = zeros(T, 3, 3)
         Q[1,1] = one(T)
         Q[2,2] = one(T)
         Q[3,3] = -one(T)
         V = standardize_basis(Q)
         @test V[:,1:2] == Q[:,1:2]
         @test V[:,3] == [zero(T), zero(T), one(T)]
     end
end
