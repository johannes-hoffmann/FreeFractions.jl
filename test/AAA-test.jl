import FreeFractions: mat_to_col_vec, mat_to_row_vec
import FreeFractions: fill_diagonal
import FreeFractions: delete_rows_and_cols, delete_rows, delete_cols,
    delete_row_and_col, delete_row, delete_col

@testset "AAA" begin
    @testset "Conversion of matrices to single row/column matrices" begin
        m = matrix(QQ, [1 2 3; 4 5 6])
        @test mat_to_row_vec(m) == matrix(QQ, [1 2 3 4 5 6])
        @test mat_to_col_vec(m) == transpose(matrix(QQ, [1 4 2 5 3 6]))
    end
    @testset "Diagonal functionality" begin
        m = matrix(QQ, [1 2 3; 4 5 6])
        filled = fill_diagonal(m, 2)
        filled_comp = hcat(m, zero_matrix(QQ, 2, 3))
        filled_comp = vcat(filled_comp, hcat(zero_matrix(QQ, 2, 3), m))
        @test filled_comp == filled
    end
    @testset "Deletion of rows and columns" begin
        m = zero_matrix(QQ, 9, 9)
        for row in 1:9, col in 1:9
            idx = row + col - 1
            m[row, col] = idx <= 9 ? idx : idx - 9
        end
        @testset "delete_row" begin
            @test delete_row(m, 1) == m[2:9, :]
            @test delete_row(m, 2) == vcat(m[1, :], m[3:9, :])
            @test delete_row(m, 9) == m[1:8, :]
            @test_throws ErrorException delete_row(m, -3)
            @test_throws ErrorException delete_row(m, 10)
        end
        @testset "delete_col" begin
            @test delete_col(m, 1) == m[:, 2:9]
            @test delete_col(m, 2) == hcat(m[:, 1], m[:, 3:9])
            @test delete_col(m, 9) == m[:, 1:8]
            @test_throws ErrorException delete_col(m, -3)
            @test_throws ErrorException delete_col(m, 10)
        end
        @testset "delete_rows" begin
            @test delete_rows(m, Int[]) == m
            @test delete_rows(m, 1:3) == m[4:9, :]
            @test delete_rows(m, 3:5) == vcat(m[1:2, :], m[6:9, :])
            @test delete_rows(m, 7:9) == m[1:6, :]
            @test_throws ErrorException delete_rows(m, -5:3)
            @test_throws ErrorException delete_rows(m, 7:10)
            @test delete_rows(m, [8, 1, 3, 8]) == vcat(m[2, :], m[4:7, :], m[9, :])
            @test_throws ErrorException delete_rows(m, [4, 5, 10, 3])
            @test delete_rows(m, 1:9) == zero_matrix(QQ, 0, 9)
        end
        @testset "delete_cols" begin
            @test delete_cols(m, Int[]) == m
            @test delete_cols(m, 1:3) == m[:, 4:9]
            @test delete_cols(m, 3:5) == hcat(m[:, 1:2], m[:, 6:9])
            @test delete_cols(m, 7:9) == m[:, 1:6]
            @test_throws ErrorException delete_cols(m, -5:3)
            @test_throws ErrorException delete_cols(m, 7:10)
            @test delete_cols(m, [8, 1, 3, 8]) == hcat(m[:, 2], m[:, 4:7], m[:, 9])
            @test_throws ErrorException delete_cols(m, [4, 5, 10, 3])
            @test delete_cols(m, 1:9) == zero_matrix(QQ, 9, 0)
        end
        @testset "delete_row_and_col (one arg)" begin
            @test delete_row_and_col(m, 1) == m[2:9, 2:9]
            @test delete_row_and_col(m, 3)[1:2, :] == hcat(m[1:2, 1:2], m[1:2, 4:9])
            @test delete_row_and_col(m, 3)[3:8, :] == hcat(m[4:9, 1:2], m[4:9, 4:9])
            @test delete_row_and_col(m, 9) == m[1:8, 1:8]
        end
        @testset "delete_row_and_col (two args)" begin
            @test delete_row_and_col(m, 1, 3) == hcat(m[2:9, 1:2], m[2:9, 4:9])
        end
        @testset "delete_rows_and_cols (one arg)" begin
            @test delete_rows_and_cols(m, Int[]) == m
            @test delete_rows_and_cols(m, 3:5)[1:2, :] == hcat(m[1:2, 1:2], m[1:2, 6:9])
            @test delete_rows_and_cols(m, 3:5)[3:6, :] == hcat(m[6:9, 1:2], m[6:9, 6:9])
            @test delete_rows_and_cols(m, 1:9) == zero_matrix(QQ, 0, 0)
        end
        @testset "delete_rows_and_cols (two args)" begin
            @test delete_rows_and_cols(m, Int[], Int[]) == m
            @test delete_rows_and_cols(m, [2, 3, 4, 6, 8, 9], 1:4) == vcat(m[1, 5:9], m[5, 5:9], m[7, 5:9])
            @test delete_rows_and_cols(m, 1:9, 1:9) == zero_matrix(QQ, 0, 0)
        end
    end
end
