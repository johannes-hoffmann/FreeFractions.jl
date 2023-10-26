@testset "RegularMinimization" begin
    for base_field in test_fields
        var_array = ["x", "y", "z", "X", "Y", "Z"]
        F6, (x, y, z, X, Y, Z) = FreeField(base_field, var_array)
        @testset "Regular minimization over $base_field" begin
            @testset "Finite factorization example over $base_field" begin
                @test x * (1 - y * x) == (1 - x * y) * x
            end
            @testset "Hua's identity over $base_field" begin
                hua = x * y * x - x + inv(inv(x) + inv(inv(y) - x))
                @test iszero(hua)
            end
            @testset "Mai example over $base_field" begin
                if characteristic(base_field) == 0
                    @test_broken dim(minimal(inv(4-x) + inv(4-x) * y * inv((4-x) - y * inv(4-x) * y) * y * inv(4-x))) == 2
                else
                    @test dim(minimal(inv(4-x) + inv(4-x) * y * inv((4-x) - y * inv(4-x) * y) * y * inv(4-x))) == 2
                end
            end
        end
    end
end
