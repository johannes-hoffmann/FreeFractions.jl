import Base: getindex, setindex!
import AbstractAlgebra: nrows, ncols
import AbstractAlgebra.Generic: _checkbounds, Perm, MatrixElem

import Random: AbstractRNG, default_rng

###############################################################################
#
#   Custom random functionality
#
###############################################################################

function FFrand(rng::AbstractRNG, R, args...)
    return rand(rng, R, args...)
end

function FFrand(R, args...)
    return FFrand(default_rng(), R, args...)
end

function FFrand(Z::FlintIntegerRing, args...)
    return FFrand(Z => 1:10, args...)
end

function FFrand(Q::FlintRationalField, args...)
    return FFrand(Q => 1:10, args...)
end

###############################################################################
#
#   Extended matrix functionality
#
###############################################################################

function mat_to_col_vec(a::MatrixElem)
    new_mat = zero_matrix(base_ring(a), nrows(a) * ncols(a), 1)
    for col in 1:ncols(a)
        start = (col - 1) * nrows(a) + 1
        stop = col * nrows(a)
        new_mat[start:stop, 1] = a[:, col]
    end
    return new_mat
end

function mat_to_row_vec(a::MatrixElem)
    new_mat = zero_matrix(base_ring(a), 1, nrows(a) * ncols(a))
    for row in 1:nrows(a)
        start = (row - 1) * ncols(a) + 1
        stop = row * ncols(a)
        new_mat[1, start:stop] = a[row, :]
    end
    return new_mat
end

function matrix_unit(
    R::Ring,
    nrows::Int,
    ncols::Int,
    row::Int,
    col::Int,
    entry = 1
)
    if !_checkbounds(nrows, row)
        error("cannot create a matrix unit with $nrows rows and an entry in
            row $row")
    end
    if !_checkbounds(ncols, col)
        error("cannot create a matrix unit with $ncols columns and an entry in
            column $col")
    end
    a = zero_matrix(R, nrows, ncols)
    a[row, col] = R(entry)
    return a
end

function matrix_unit(a::MatrixElem, row::Int, col::Int, entry = 1)
    return matrix_unit(base_ring(a), nrows(a), ncols(a), row, col, entry)
end

function fill_diagonal(x::MatrixElem, blocks::Int)
    (rows, cols) = size(x)
    result = zero_matrix(base_ring(x), rows * blocks, cols * blocks)
    for k in 1:blocks
        row_range = ((k - 1) * rows + 1):(k * rows)
        col_range = ((k - 1) * cols + 1):(k * cols)
        result[row_range, col_range] = x
    end
    return result
end

function permute_rows!(x::Generic.MatrixElem, P::Generic.Perm)
    z = P * x
    for row in 1:nrows(x), col in 1:ncols(x)
        x[row, col] = z[row, col]
    end
    return x
end

function permute_rows(x::Generic.MatrixElem, P::Generic.Perm)
    return P * x
end

function permute_cols!(x::Generic.MatrixElem, P::Generic.Perm)
    z = transpose(P * transpose(x))
    for row in 1:nrows(x), col in 1:ncols(x)
        x[row, col] = z[row, col]
    end
    return x
end

function permute_cols(x::Generic.MatrixElem, P::Generic.Perm)
    return transpose(P * transpose(x))
end

function delete_rows_and_cols(
    a::MatrixElem,
    row_list::AbstractVector{Int},
    col_list::AbstractVector{Int}
)
    (nrows, ncols) = size(a)
    rows_to_delete = sort(collect(row_list))
    if !isempty(rows_to_delete)
        if first(rows_to_delete) < 1
            error("cannot delete rows with nonpositive indices")
        end
        if last(rows_to_delete) > nrows
            error("cannot delete rows with index larger than the total number of rows")
        end
        unique!(rows_to_delete)
    end
    cols_to_delete = sort(collect(col_list))
    if !isempty(cols_to_delete)
        if first(cols_to_delete) < 1
            error("cannot delete columns with nonpositive indices")
        end
        if last(cols_to_delete) > ncols
            error("cannot delete columns with index larger than the total number of rows")
        end
        unique!(cols_to_delete)
    end
    new_rows = setdiff(1:nrows, rows_to_delete)
    new_cols = setdiff(1:ncols, cols_to_delete)
    b = similar(a, length(new_rows), length(new_cols))
    row_offset = 0
    for row in eachindex(new_rows), col in eachindex(new_cols)
        b[row, col] = a[new_rows[row], new_cols[col]]
    end
    return b
end

function delete_rows_and_cols(a::MatrixElem, indices::AbstractVector{Int})
    return delete_rows_and_cols(a, indices, indices)
end

function delete_row_and_col(a::MatrixElem, row::Int, col::Int)
    return delete_rows_and_cols(a, [row], [col])
end

function delete_row_and_col(a::MatrixElem, index::Int)
    return delete_rows_and_cols(a, [index], [index])
end

function delete_rows(a::MatrixElem, rows::AbstractVector{Int})
    return delete_rows_and_cols(a, rows, Int[])
end

function delete_row(a::MatrixElem, row::Int)
    return delete_rows_and_cols(a, [row], Int[])
end

function delete_cols(a::MatrixElem, cols::AbstractVector{Int})
    return delete_rows_and_cols(a, Int[], cols)
end

function delete_col(a::MatrixElem, col::Int)
    return delete_rows_and_cols(a, Int[], [col])
end
