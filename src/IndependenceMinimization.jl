###############################################################################
#
#   Checking family independence (two-sided)
#
###############################################################################

function family_independence_at_point(
    a::FreeFraction{T},
    left::Bool,
    args::MatElem{T}...
) where T
    mat_size = nrows(args[1])
    d = dim(a)
    g = nvars(a)
    family = try
        if left
            evaluate_left_family(a, args...)
        else
            evaluate_right_family(a, args...)
        end
    catch
        return (-1, zero_matrix(base_ring(a), 0, 0))
    end
    new_mat = zero_matrix(base_ring(a), 0, mat_size * mat_size)
    for new_row in 1:d
        start = (new_row - 1) * mat_size + 1
        stop = new_row * mat_size
        if left
            l_range = start:stop
            r_range = 1:mat_size
        else
            l_range = 1:mat_size
            r_range = start:stop
        end
        temp = mat_to_row_vec(family[l_range, r_range])
        new_mat = vcat(new_mat, temp)
    end
    (nullity, ker) = left_kernel(new_mat)
    if left
        return (nullity, ker)
    else
        return (nullity, transpose(ker))
    end
end

function family_independence(
    a::FreeFraction,
    left::Bool;
    max_size::Int = 10
)
    if max_size < 1
        error("max_size has to be strictly positive")
    end
    R = base_ring(a)
    candidate_space = identity_matrix(R, dim(a))
    for mat_size in 2:max_size
        test_mats = []
        for var_index in 1:nvars(a)
            push!(test_mats, matrix(R, FFrand(R, mat_size, mat_size)))
        end
        (nullity, ker) = family_independence_at_point(a, left, test_mats...)
        if nullity == 0
            return (true, ker)
        elseif nullity > 0
            old_rank = rank(candidate_space)
            if left
                candidate_space = dirty_row_intersect(candidate_space, ker)
            else
                candidate_space = dirty_col_intersect(candidate_space, ker)
            end
            if rank(candidate_space) == 0
                return (true, candidate_space)
            end
        end
    end
    return (false, candidate_space)
end

###############################################################################
#
#   Checking family independence (one-sided)
#
###############################################################################

function left_family_independence_at_point(
    a::FreeFraction{T},
    args::MatElem{T}...
) where T
    return family_independence_at_point(a, true, args...)
end

function right_family_independence_at_point(
    a::FreeFraction{T},
    args::MatElem{T}...
) where T
    return family_independence_at_point(a, false, args...)
end

function left_family_independence(a::FreeFraction)
    return family_independence(a, true)
end

function right_family_independence(a::FreeFraction)
    return family_independence(a, false)
end

###############################################################################
#
#   Ad-hoc intersection of vector spaces
#
###############################################################################

function dirty_col_intersect(a::MatElem{T}, b::MatElem{T}) where T
    if nrows(a) != nrows(b)
        error("number of rows must be equal")
    end
    (nullity, ker) = nullspace(hcat(a, -b))
    return a * ker[1:ncols(a), :]
end

function dirty_row_intersect(a::MatElem{T}, b::MatElem{T}) where T
    if ncols(a) != ncols(b)
        error("number of columns must be equal")
    end
    (nullity, ker) = left_kernel(vcat(a, -b))
    return ker[:, 1:nrows(a)] * a
end

###############################################################################
#
#   Eliminate last entry of family (two-sided)
#
###############################################################################

function eliminate_last_family_entry(
    a::GenericFreeFraction{T},
    c::MatElem{T},
    left::Bool
) where T
    n = dim(a)
    if (left && size(c) != (1, n)) || (!left && size(c) != (n, 1))
        error("wrong dimensions: (size, n, left) is ($(size(c)), $n, $left)")
    end
    if iszero(c)
        error("cannot eliminate with a zero vector")
    end
    last_non_zero_index = 0
    for index in n:-1:1
        if (left && c[1, index] != 0) || (!left && c[index, 1] != 0)
            last_non_zero_index = index
            break
        end
    end
    if left
        swap = swap_cols
        multiply! = multiply_column!
        add! = add_column!
        mod = identity
    else
        swap = swap_rows
        multiply! = multiply_row!
        add! = add_row!
        mod = transpose
    end
    new_c = mod(swap(c, last_non_zero_index, n))
    new_a = swap(a, last_non_zero_index, n)
    multiply!(new_a, inv(new_c[1, n]), n)
    for index in 1:(n - 1)
        add!(new_a, -new_c[1, index], n, index)
    end
    return new_a
end

###############################################################################
#
#   Eliminate last entry of family (one-sided)
#
###############################################################################

function eliminate_last_left_family_entry(
    a::GenericFreeFraction{T},
    c::MatElem{T}
) where T
    return eliminate_last_family_entry(a, c, true)
end

function eliminate_last_right_family_entry(
    a::GenericFreeFraction{T},
    c::MatElem{T}
) where T
    return eliminate_last_family_entry(a, c, false)
end

###############################################################################
#
#   Minimization (two-sided)
#
###############################################################################

function guess_smaller_representation_normalized(
    a::GenericFreeFraction,
    left::Bool
)
    n = dim(a)
    for index in 1:n
        if left
            candidate = admissible!(delete_row_and_col(a, index, n))
        else
            candidate = delete_row_and_col(a, n, index)
        end
        if probabilistic_fullness_of_representation_matrix(candidate)
            return (true, candidate)
        end
    end
    return (false, a)
end

function guess_smaller_representation(
    a::GenericFreeFraction,
    left::Bool;
    max_size::Int = 10
)
    (flag, ker) = family_independence(a, left; max_size)
    if flag
        return a
    end
    if left
        mod = identity
    else
        mod = transpose
    end
    nullity = size(mod(ker), 1)
    for index in 1:nullity
        temp = eliminate_last_family_entry(a, mod(mod(ker)[index, :]), left)
        (flag, cand) = guess_smaller_representation_normalized(temp, left)
        if flag
            return cand
        end
    end
    return a
end

function minimal_independence(a::FreeFraction; max_size::Int = 10)
    n = dim(a)
    b = deepcopy(a)
    new_n = 0
    while n != new_n
        n = dim(b)
        b = guess_smaller_representation(b, true; max_size)
        b = guess_smaller_representation(b, false; max_size)
        new_n = dim(b)
    end
    return b
end

###############################################################################
#
#   Minimization (one-sided)
#
###############################################################################

function guess_smaller_representation_left_normalized(a::GenericFreeFraction)
    return guess_smaller_representation_normalized(a, true)
end

function guess_smaller_representation_right_normalized(a::GenericFreeFraction)
    return guess_smaller_representation_normalized(a, false)
end

function guess_smaller_representation_left(a::GenericFreeFraction)
    return guess_smaller_representation(a, true)
end

function guess_smaller_representation_right(a::GenericFreeFraction)
    return guess_smaller_representation(a, false)
end
