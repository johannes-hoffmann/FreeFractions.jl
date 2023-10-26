@testset "GenericFreeFraction" begin
    for base_field in test_fields
        F, (x, y) = FreeField(base_field, ["x", "y"])
        @testset "Constructors over $base_field" begin
            row = matrix(base_field, [1 0])
            col = matrix(base_field, 2, 1, [0, 1])
            cm = identity_matrix(base_field, 2)
            vm1 = matrix(base_field, base_field.([0 -1; 0 0]))
            vm2 = zero_matrix(base_field, 2, 2)
            f = GenericFreeFraction(row, cm, [vm1, vm2], col, F)
            @test x == f
            @test f == GenericFreeFraction(x)
            @test isa(MonomialFreeFraction(x), MonomialFreeFraction)
            g = GenericFreeFraction(row, cm, [vm2, vm2], col, F)
            @test MonomialFreeFraction(g) == 0
        end
        @testset "Coercion of monomials over $base_field" begin
            m = 25 * x * y * y * x
            g_m = GenericFreeFraction(m)
            @test MonomialFreeFraction(g_m) == m
            if characteristic(base_field) == 5
                @test MonomialFreeFraction(m * inv(y)) == 0
            else
                @test_throws ErrorException MonomialFreeFraction(m * inv(y))
            end
        end
        @testset "Regularity over $base_field" begin
            @test isregular(x + y)
            @test !isregular(inv(x + y))
        end
        @testset "Arithmetic over $base_field" begin
            f = x + y
            g = y + x
            fivef = 5 * f
            @test fivef == f * 5
            @test -f == -x - y

            @test F(5) * f == fivef
            @test f * F(5) == fivef

            @test GenericFreeFraction(x) == x
            @test x == GenericFreeFraction(x)

            @test f == g
            @test f * g == x^2 + x * y + y * x + y^2
            @test f != g + 1
        end
    end
end
