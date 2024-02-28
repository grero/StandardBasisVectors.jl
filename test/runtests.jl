using StandardBasisVectors
using Test


@testset "Basic" begin 

    # convert a simple left-handed 2D system to a right-handed system
    Q = [-1.0 0.0;
         0.0 1.0]
    
    V,signs,θ = StandardBasisVectors.orient_vectors(Q)
    
    @test signs == [-1.0, 1.0]
    @test V[:,1] ≈ [1.0,0.0] 
    @test V[:,2] ≈ [0.0,1.0]

     Q = [1.0 0.0 0.0;
          .0 1.0 0.0;
          .0 0.0 -1.0]

     V = standardize_basis(Q)
     @test V[:,1:2] == Q[:,1:2]
     @test V[:,3] == [0.0, 0.0, 1.0]
end
