import FreeFractions: isscalar, monomial_indices

@testset "MonomialFreeFraction" begin
    for base_field in test_fields
        F, (x, y) = FreeField(base_field, ["x", "y"])
        @testset "Constructors over $base_field" begin
            answer = F(42)
            @test coefficient(answer) == 42
            @test monomial_indices(answer) == Int[]
            @test parent(answer) == F
            @test typeof(answer) == MonomialFreeFraction{elem_type(base_field)}
            @test isscalar(answer)

            @test coefficient(x) == 1
            @test monomial_indices(x) == [1]
            @test parent(x) == F
            @test typeof(x) == MonomialFreeFraction{elem_type(base_field)}
            @test !isscalar(x)
            @test MonomialFreeFraction([1], F) == x
            @test MonomialFreeFraction(x) == x
        end
        @testset "Arithmetic over $base_field" begin
            @test F(1) == 1
            @test F(1) + F(2) == 3
            @test F(1) + 2 == 3
            @test 1 + F(2) == 3
            @test F(1) + base_field(2) == 3
            @test base_field(1) + F(2) == 3

            @test -F(1) == F(-1)
            @test F(1) - F(2) == -1
            @test F(1) - 2 == -1
            @test 1 - F(2) == -1

            @test F(2) * F(3) == 6
            @test F(2) * 3 == 6
            @test 2 * F(3) == 6

            @test F(2)^4 == 2^4
            @test F(2)^1 == 2
            @test F(2)^0 == one(F(2))
            @test zero(F(2))^0 == one(F(2))

            if characteristic(base_field) != 2
                half = inv(base_field(2))
                @test F(half) + half == 1
                @test half + F(half) == 1
                @test inv(F(2)) == base_field(1)//2
                @test F(2)^(-2) == base_field(1)//4
            end

            @test x != y
            @test x != F(1)
            @test F(1) != x
            @test MonomialFreeFraction(0, [1], F) == MonomialFreeFraction(0, [1, 2], F)
            @test iszero(MonomialFreeFraction(0, [1], F))

            @test x - x == 0
            @test MonomialFreeFraction(0, [1], F) - x == -x

            fivex = MonomialFreeFraction(5, [1], F)
            @test 5 * x == fivex
            @test 5 * x == x * 5
            @test 5 * x == F(5) * x
            @test x * 5 == x * F(5)
            @test MonomialFreeFraction(5, Int[], F) == 5

            @test x * y * x != x * x * y

            @test (x * 5 * y)^2 == MonomialFreeFraction(25, [1, 2, 1, 2], F)
            @test fivex^0 == one(fivex)
        end

    end
end
