import FreeFractions: block_upper_triangular_structure

@testset "GenericMinimization" begin
    base_field = QQ
    @testset "Block structure over $base_field" begin
        test_array = [
            (identity_matrix(base_field, 4), [1:1, 2:2, 3:3, 4:4]),
            (matrix(base_field, [1 2 3 4; 0 6 7 8; 0 1 2 3; 0 4 5 6]), [1:1, 2:4]),
            (matrix(base_field, [1 2 3 4; 5 6 7 8; 0 0 1 2; 0 0 3 4]), [1:2, 3:4]),
            (matrix(base_field, [1 2 3 4; 0 6 7 8; 0 0 1 2; 0 0 3 4]), [1:1, 2:2, 3:4]),
            (matrix(base_field, [1 2 3 4; 0 6 7 8; 0 9 1 2; 0 0 3 4]), [1:1, 2:4]),
            (matrix(base_field, [1 2 3 4; 0 6 7 8; 0 9 1 2; 5 0 3 4]), [1:4]),
            (matrix(base_field, [1 2 3 4; 5 6 7 8; 1 2 3 4; 0 0 0 5]), [1:3, 4:4]),
            (matrix(base_field, [1 2 3 4; 0 6 7 8; 0 2 3 4; 0 0 0 5]), [1:1, 2:3, 4:4]),
        ]
        for (mat, buts) in test_array
            @test buts == block_upper_triangular_structure(mat)
        end
        @test_throws DomainError block_upper_triangular_structure(zero_matrix(base_field, 4, 3))
        @test_throws ErrorException block_upper_triangular_structure(zero_matrix(base_field, 4, 4))
    end
    for base_field in test_fields
        @testset "Generic minimization over $base_field" begin
            var_array = ["x", "y", "z", "X", "Y", "Z"]
            F6, (x, y, z, X, Y, Z) = FreeField(base_field, var_array)
            @testset "Basic examples over $base_field" begin
                gen_x = GenericFreeFraction(x)
                @test iszero(inv(gen_x)-inv(x))
                @test_throws ErrorException minimal(inv(gen_x - x))
            end
            @testset "inverse of (anti-)commutator over $base_field" begin
                @test dim(minimal(inv(x * y - y * x))) == 3
                @test dim(minimal(inv(x * y + y * x))) == 3
            end
            @testset "Gardam example over $base_field" begin
                gardam = x - inv(z^(-2)*y + y^(-2)*x + x^(-1)*y*z)
                @test dim(minimal(gardam)) == 8
            end
            @testset "Schrempf example over $base_field" begin
                konrad_g = inv(x * X + inv(z * Z)) - inv(z * Y + y * Z)
                konrad = inv(z * Y + y * Z) - inv(y * Z + inv(konrad_g) + z * Y)
                if characteristic(base_field) in [2, 3]
                    @test_broken dim(minimal(konrad)) == 9
                else
                    @test dim(minimal(konrad)) == 9
                end
            end
        end
    end
end
