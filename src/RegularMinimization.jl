function compute_controllability_space(a::GenericFreeFraction)
    d = dim(a)
    v = column_vector(a)
    if d == 1
        return v
    end
    if iszero(v)
        return zero_matrix(base_ring(a), d, 0)
    end
    b = v
    s = v
    base_mat = v
    old_size = ncols(b)
    while true
        m = zero_matrix(base_ring(a), d, 0)
        for varindex in 1:nvars(a), col in 1:ncols(s)
            m = hcat(m, variable_matrices(a)[varindex] * s[:, col])
        end
        s = zero_matrix(base_ring(a), d, 0)
        for col in 1:ncols(m)
            test_mat = hcat(base_mat, m[:, col])
            if rank(test_mat) == ncols(test_mat)
                s = hcat(s, m[:, col])
                b = hcat(b, m[:, col])
                base_mat = test_mat
                if ncols(b) == d
                    return b
                end
            end
        end
        if ncols(b) == old_size
            return b
        end
        old_size = ncols(b)
    end
end

function recover_controllability_space(a::GenericFreeFraction)
    v = column_vector(a)
    d = dim(a)
    if iszero(v) || d == 1
        return v
    end
    c_mat = transpose(row_vector(a))
    b_mat = v
    sb = v
    sc = c_mat
    base_mat = v
    old_size = ncols(b_mat)
    var_mats = variable_matrices(a)
    while true
        mb = zero_matrix(base_ring(a), d, 0)
        mc = zero_matrix(base_ring(a), d, 0)
        for varindex in 1:nvars(a), col in 1:ncols(sb)
            mb = hcat(mb, var_mats[varindex] * sb[:, col])
            mc = hcat(mc, transpose(var_mats[varindex]) * sc[:, col])
        end
        sb = zero_matrix(base_ring(a), d, 0)
        sc = zero_matrix(base_ring(a), d, 0)
        for col in 1:ncols(mb)
            test_mat = hcat(base_mat, mb[:, col])
            if rank(test_mat) == ncols(test_mat)
                sb = hcat(sb, mb[:, col])
                b_mat = hcat(b_mat, mb[:, col])
                sc = hcat(sc, mc[:, col])
                c_mat = hcat(c_mat, mc[:, col])
                base_mat = test_mat
                if ncols(b_mat) == d
                    return (b_mat, c_mat)
                end
            end
        end
        if ncols(b_mat) == old_size
            return (b_mat, c_mat)
        end
        old_size = ncols(b_mat)
    end
end

function compute_complement(b::MatElem{T}) where T<:FieldElem
    n = nrows(b)
    c = similar(b, n, 0)
    if ncols(b) == n
        return c
    end
    e = identity_matrix(base_ring(b), n)
    if iszero(b)
        return e
    end
    base_mat = b
    for col in 1:n
        test_mat = hcat(base_mat, e[:, col])
        if rank(test_mat) == ncols(test_mat)
            c = hcat(c, e[:, col])
            base_mat = test_mat
            if ncols(base_mat) == n
                @assert ncols(b) + ncols(c) == n
                @assert rank(hcat(b, c)) == n
                return c
            end
        end
    end
end

function cutdown(a::GenericFreeFraction, n::Integer, offset::Integer)
    start = 1 + offset
    stop = n + offset
    g = nvars(a)
    if iszero(row_vector(a)) || iszero(column_vector(a))
        return parent(a)(0)
    end
    return GenericFreeFraction(
        row_vector(a)[:, start:stop],
        constant_matrix(a)[start:stop, start:stop],
        map(mat -> mat[start:stop, start:stop], variable_matrices(a)),
        column_vector(a)[start:stop, :],
        parent(a)
    )
end

function minimal_regular(a::GenericFreeFraction)
    if !isregular(a)
        error("minimal_regular only applies to regular free fractions")
    end
    MDR = monic_descriptor_realization(a)
    b = compute_controllability_space(MDR)
    c = compute_complement(b)
    P = hcat(b, c)
    MDR = inv(P) * MDR * P
    MDR = cutdown(MDR, ncols(b), 0)
    if dim(MDR) == 0
        return parent(a)(0)
    elseif dim(MDR) == 1
        if all(x -> x == 0, variable_matrices(MDR))
            return MonomialFreeFraction(MDR)
        end
        return MDR
    end
    MDR = transpose(MDR)
    b = compute_controllability_space(MDR)
    c = compute_complement(b)
    P = hcat(c, b)
    MDR = inv(P) * MDR * P
    MDR = cutdown(MDR, ncols(b), ncols(c))
    if dim(MDR) == 0
        return parent(a)(0)
    elseif dim(MDR) == 1
        if all(x -> x == 0, variable_matrices(MDR))
            return MonomialFreeFraction(MDR)
        end
        return MDR
    end
    MDR = transpose(MDR)
    return MDR
    # return clean_rhs(MDR) TODO
    # TODO: beautify
end

function minimal_regular(
    a::GenericFreeFraction{T},
    points::Union{T, Int, Rational{Int}}...
) where T
    try
        a(points...)
    catch e
        if isa(e, DomainError)
            error("free fraction is not regular in the given point")
        else
            rethrow(e)
        end
    end
    if all(x -> x == 0, points)
        @debug "points is zero"
        return minimal_regular(a)
    end
    b = forward_shift(a, points...)
    b = minimal_regular(b)
    return backward_shift(b, points...)
end
function minimal_regular(
    a::FreeFraction{T},
    args::Array{<:Union{T, Int, Rational{Int}}, 1}
) where T
    return minimal_regular(a, args...)
end
