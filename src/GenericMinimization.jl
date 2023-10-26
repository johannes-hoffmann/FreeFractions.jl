###############################################################################
#
#   Block structure detection
#
###############################################################################

function block_upper_triangular_structure(a::MatElem)
    check_square(a)
    dim = nrows(a)
    if dim == 0
        return Int[]
    end
    for idx in 1:dim
        if iszero_row(a, idx) || iszero_column(a, idx)
            error("invalid free fraction")
        end
    end
    bds = [1:dim]
    current_column = 1
    while true
        if a[dim, current_column] != 0
            return bds
        end
        last_start_index = first(last(bds))
        next_estimate = last_start_index + 1
        for row in dim:-1:(last_start_index + 1)
            if a[row, current_column] != 0
                next_estimate = max(next_estimate, row + 1)
            end
        end
        for col in (current_column + 1):dim
            if col == next_estimate
                bds[lastindex(bds)] = last_start_index:(next_estimate - 1)
                push!(bds, next_estimate:dim)
                current_column = next_estimate
                break
            end
            for row in dim:-1:next_estimate
                if a[row, col] != 0
                    next_estimate = row + 1
                    break
                end
            end
            if next_estimate == dim + 1
                return bds
            end
        end
    end
    return bds
end

function block_upper_triangular_structure(a::FreeFraction)
    return block_upper_triangular_structure(representation_matrix(a))
end

function number_of_pivot_blocks(a::FreeFraction)
    return length(block_upper_triangular_structure(a))
end

###############################################################################
#
#   Pivot refinement
#
###############################################################################

function pivot_refinement(a::FreeFraction; kwargs...)
    return pivot_refinement!(deepcopy(a); kwargs...)
end

function pivot_refinement!(a::FreeFraction; kwargs...)
    buts = block_upper_triangular_structure(a)
    for block in buts
        pivot_refinement!(a, block; kwargs...)
    end
    return a
end

function pivot_refinement!(
    a::FreeFraction,
    block::UnitRange{Int};
    max_permutation_size::Int = 6
)
    if length(block) == 1
        return a
    end
    if length(block) <= max_permutation_size
        G = SymmetricGroup(dim(a))
        H = SymmetricGroup(length(block))
        iota = emb(G, collect(block))
        counter = 1
        for p in elements!(H)
            if iota(p)[1] == 1 # the first column has to be fixed
                permute_cols!(a, iota(p))
                pivot_refinement_row!(a, block)
                if !(block in block_upper_triangular_structure(a))
                    return a
                end
                permute_cols!(a, inv(iota(p)))
            end
            permute_rows!(a, iota(p))
            pivot_refinement_column!(a, block)
            if !(block in block_upper_triangular_structure(a))
                return a
            end
            permute_rows!(a, inv(iota(p)))
            counter += 1
        end
        return a
    else
        pivot_refinement_row!(a, block)
        if !(block in block_upper_triangular_structure(a))
            return a
        end
        pivot_refinement_column!(a, block)
        if !(block in block_upper_triangular_structure(a))
            return a
        end
    end
end

function pivot_refinement_row!(a::FreeFraction, block::UnitRange{Int})
    len = length(block)
    for col in first(block):(last(block) - 1)
        eq_mat = constant_matrix(a)[block, first(block):col]
        for mat in variable_matrices(a)
            eq_mat = hcat(eq_mat, mat[block, first(block):col])
        end
        (nul, ker) = left_kernel(eq_mat)
        if nul >= last(block) - col
            T = ker
            for pos in 1:len
                T_new = vcat(matrix_unit(base_ring(a), 1, len, 1, pos), T)
                r = rank(T_new)
                if r > rank(T)
                    T = T_new
                    if r == len
                        break
                    end
                end
            end
            if rank(T) < len
                error("could not extend basis during pivot_refinement_row!")
            end
            col_block = first(block):dim(a)
            temp = constant_matrix(a)[block, col_block]
            constant_matrix(a)[block, col_block] = T * temp
            for mat in variable_matrices(a)
                mat[block, col_block] = T * mat[block, col_block]
            end
            column_vector(a)[block, :] = T * column_vector(a)[block, :]
            pivot_refinement!(a, first(block):col)
            pivot_refinement!(a, (col + 1):last(block))
            break
        end
    end
    return a
end

function pivot_refinement_column!(a::FreeFraction, block::UnitRange{Int})
    len = length(block)
    for row in (first(block) + 1):last(block)
        eq_mat = constant_matrix(a)[row:last(block), block]
        for mat in variable_matrices(a)
            eq_mat = vcat(eq_mat, mat[row:last(block), block])
        end
        (nul, ker) = right_kernel(eq_mat)
        if nul >= row - first(block)
            start_pos = 1
            T = ker
            if first(block) == 1
                start_pos += 1
                non_zero_pos = 0
                for col in 1:ncols(T)
                    if T[1, col] != 0
                        non_zero_pos = col
                        break
                    end
                end
                if non_zero_pos == 0
                    continue
                elseif non_zero_pos > 1
                    swap_cols!(T, 1, non_zero_pos)
                end
                for col in 2:ncols(T)
                    if T[1, col] != 0
                        add_column!(T, -T[1, col], 1, col)
                    end
                end
            end
            for pos in start_pos:len
                T_new = hcat(T, matrix_unit(base_ring(a), len, 1, pos, 1))
                r = rank(T_new)
                if r > rank(T)
                    T = T_new
                    if r == len
                        break
                    end
                end
            end
            if rank(T) < len
                error("could not extend basis during pivot_refinement_column!")
            end
            constant_matrix(a)[1:last(block), block] *= T
            for mat in variable_matrices(a)
                mat[1:last(block), block] *= T
            end
            pivot_refinement!(a, first(block):(row - 1))
            pivot_refinement!(a, row:last(block))
            break
        end
    end
    return a
end

###############################################################################
#
#   Minimization
#
###############################################################################

function extend(a::GenericFreeFraction)
    R = base_ring(a)
    new_dim = dim(a) + 1
    new_row = matrix_unit(R, 1, new_dim, 1, 1)
    new_col = vcat(zero_matrix(R, 1, 1), column_vector(a))
    new_const = block_diagonal_matrix([identity_matrix(R, 1), constant_matrix(a)])
    new_const[1, 2] = R(-1)
    var_mats = variable_matrices(a)
    new_var = map(mat -> block_diagonal_matrix([zero_matrix(R, 1, 1), mat]), var_mats)
    return GenericFreeFraction(new_row, new_const, new_var, new_col, parent(a))
end

function contract(a::GenericFreeFraction)
    if any(mat -> mat[1, :] != 0, variable_matrices(a))
        error("cannot contract - non-scalar entry in first row")
    end
    if constant_matrix(a)[1, 2] == 0
        for col in 3:dim(a)
            if constant_matrix(a)[1, col] != 0
                swap_cols!(a, 2, col)
                break
            end
        end
    end
    if constant_matrix(a)[1, 2] == 0
        error("cannot contract - no non-zero entry found")
    end
    multiply_column!(a, -inv(constant_matrix(a)[1, 2]), 2)
    for col in 3:dim(a)
        if constant_matrix(a)[1, col] != 0
            add_column!(a, constant_matrix(a)[1, col], 2, col)
        end
    end
    b = delete_row_and_col(a, 1, 1)
    row_vector(b)[1, 1] = 1
    return b
end

function solve_minimization_equations(
    a::FreeFraction,
    k_block::UnitRange{<:Int},
    left::Bool
)
    R = base_ring(a)
    n = dim(a)
    if left
        if !(1 <= first(k_block) <= last(k_block) < n)
            error("wrong block parameter for solving left equations: $k_block")
        end
        T_block = (last(k_block) + 1):n
        U_block = k_block
    else
        if !(1 < first(k_block) <= last(k_block) <= n)
            error("wrong block parameter for solving right equations: $k_block")
        end
        T_block = k_block
        U_block = 1:(first(k_block) - 1)
    end
    T_len = length(T_block)
    U_len = length(U_block)
    upper_half = fill_diagonal(constant_matrix(a)[T_block, T_block], U_len)
    for mat in variable_matrices(a)
        temp = mat[T_block, T_block]
        upper_half = hcat(upper_half, fill_diagonal(temp, U_len))
    end
    if left
        upper_right = fill_diagonal(column_vector(a)[T_block, 1], U_len)
        upper_half = hcat(upper_half, upper_right)
    end
    if first(U_block) == 1
        mod_col_block = 2:last(U_block)
    else
        mod_col_block = U_block
    end
    cmat = transpose(constant_matrix(a)[U_block, mod_col_block])
    lower_half = kronecker_product(cmat, identity_matrix(R, T_len))
    for mat in variable_matrices(a)
        vmat = transpose(mat[U_block, mod_col_block])
        new_col = kronecker_product(vmat, identity_matrix(R, T_len))
        lower_half = hcat(lower_half, new_col)
    end
    if left
        lower_half = hcat(lower_half, zero_matrix(R, nrows(lower_half), U_len))
    end
    lhs = vcat(upper_half, lower_half)
    rhs = zero_matrix(R, 1, 0)
    for mat in matrix_pencil(a)
        for row in U_block
            rhs = hcat(rhs, -mat[row, T_block])
        end
    end
    if left
        rhs = hcat(rhs, -transpose(column_vector(a)[U_block, 1]))
    end
    T = zero_matrix(R, U_len, T_len)
    U = zero_matrix(R, U_len, T_len)
    x = nothing
    try
        x = solve_left(lhs, rhs)
    catch
        return (false, T, U)
    end
    for row in 1:U_len, col in 1:T_len
        T[row, col] = x[1, T_len * (row - 1) + col]
    end
    if first(U_block) == 1
        for row in 1:(U_len - 1), col in 1:T_len
            U[row + 1, col] = x[1, T_len * (row - 1) + col + U_len * T_len]
        end
    else
        for row in 1:U_len, col in 1:T_len
            U[row, col] = x[1, T_len * (row - 1) + col + U_len * T_len]
        end
    end
    return (true, T, U)
end

function minimal_generic(
    a::GenericFreeFraction;
    warn::Bool = true,
    min_inv::Bool = true
)
    b = _minimal_generic(a; warn = false)
    if min_inv && dim(b) > 1
        inv_b = _minimal_generic(inv(b); warn = false)
        inv_inv_b = _minimal_generic(inv(inv_b); warn = false)
        if dim(inv_inv_b) < dim(b)
            b = minimal_generic(inv_inv_b; warn = false, min_inv)
        end
    end
    if dim(b) > 1 && warn
        @warn "minimization of non-regular free fractions is not guaranteed to
            be optimal"
    end
    return b
end

function _minimal_generic(a::GenericFreeFraction; warn::Bool = true)
    b = pivot_refinement(a)
    R = base_ring(a)
    k = 2
    while k <= number_of_pivot_blocks(b) # catches case of only one pivot block
        buts = block_upper_triangular_structure(b)
        m = length(buts)
        kprime = m + 1 - k
        kprime_block = buts[kprime]
        # check left
        (solvable, T, U) = solve_minimization_equations(b, kprime_block, true)
        if solvable
            if kprime == 1
                return parent(a)(0)
            end
            P = identity_matrix(R, dim(b))
            Q = identity_matrix(R, dim(b))
            col_block = (last(kprime_block) + 1):dim(b)
            P[kprime_block, col_block] = T
            Q[kprime_block, col_block] = U
            b = P * b * Q
            b = delete_rows_and_cols(b, kprime_block)
            if k > max(2, (m + 1) // 2)
                k = k - 1
            end
            continue
        end
        if kprime == 1 # check left extended in special case
            kprime_block_ext = kprime_block .+ 1
            b_ext = extend(b)
            (solvable, T, U) = solve_minimization_equations(b_ext, kprime_block_ext, true)
            if solvable
                P = identity_matrix(R, dim(b_ext))
                Q = identity_matrix(R, dim(b_ext))
                col_block = (last(kprime_block_ext) + 1):dim(b_ext)
                P[kprime_block_ext, col_block] = T
                Q[kprime_block_ext, col_block] = U
                b_ext = P * b_ext * Q
                b_ext = delete_rows_and_cols(b_ext, kprime_block_ext)
                b = contract(b_ext)
                if k > max(2, (m + 1) // 2)
                    k = k - 1
                end
                continue
            end
        end
        # check right
        k_block = buts[k]
        (solvable, T, U) = solve_minimization_equations(b, k_block, false)
        if solvable
            P = identity_matrix(R, dim(b))
            Q = identity_matrix(R, dim(b))
            row_block = 1:(first(k_block) - 1)
            P[row_block, k_block] = T
            Q[row_block, k_block] = U
            b = P * b * Q
            b = delete_rows_and_cols(b, k_block)
            if k > max(2, (m + 1) // 2)
                k = k - 1
            end
            continue
        end
        k += 1
    end
    block = last(block_upper_triangular_structure(b))
    d = dim(b)
    if block == d:d
        if column_vector(b)[d, 1] == 0
            b = delete_rows_and_cols(b, block)
        end
    end
    if all(x -> x == 0, variable_matrices(b))
        b = MonomialFreeFraction(b)
    end
    if dim(b) > 1 && warn
        @warn "minimization of non-regular free fractions is not guaranteed to
            be optimal"
    end
    return b
    # TODO: return Pb such that P*v_b=[0...0 lambda_prime]
end
