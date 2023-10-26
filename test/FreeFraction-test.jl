@testset "FreeFraction" begin
    for base_field in test_fields
        @testset "String constructor over $base_field" begin
            comp_MFF = MonomialFreeFraction{elem_type(base_field)}
            @test isa(FreeFraction("x", base_field, ["x", "y"]), comp_MFF)
            comp_GFF = GenericFreeFraction{elem_type(base_field)}
            @test isa(FreeFraction("z+u", base_field, ["z", "u"]), comp_GFF)
            f = FreeFraction("z+u", base_field, ["z", "u"])
            @test parent(f) == FreeField(base_field, ["z", "u"])[1]
            @test typeof(f) == GenericFreeFraction{elem_type(base_field)}
            @test base_ring(f) == base_field
            @test nvars(f) == 2
        end
        F, (x, y) = FreeField(base_field, ["x", "y"])
        @testset "Evaluation over $base_field" begin
            f = x + y
            @test f(2, 3) == 5
            m = matrix(base_field, [1 2; 3 7])
            n = matrix(base_field, [2 3; 7 5])
            @test f(m, n) == m + n
            g = inv(x * y - y * x)
            @test_throws DomainError g(1, 1)
            @test g(m, n) == inv(m * n - n * m)
        end
    end
end
