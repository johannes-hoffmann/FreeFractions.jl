@testset "FreeField" begin
    G, (a, b) = FreeField(QQ, ["a", "b"])
    for base_field in test_fields
        F, (x, y) = FreeField(base_field, ["x", "y"])
        @test F == FreeField(base_field, [:x, :y])
        @test F != G
        @test elem_type(F) == FreeFraction{elem_type(base_field)}
        @test parent_type(x) == FreeField{elem_type(base_field)}
        @test isa(x, FreeFraction)
        @test isa(x, MonomialFreeFraction)
        @test !isa(x, GenericFreeFraction)
        @test typeof(x) <: NCRingElem
        @test typeof(F) <: AbstractAlgebra.NCRing
        @test gen(F, 1) == x
        @test gens(F) == [x, y]
        @test one(F) == F(1)
        @test zero(F) == F(0)


    end
end
