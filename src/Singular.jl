###############################################################################
#
#   Gröbner-based fullness test (Cohn/Reutenauer)
#
###############################################################################

function is_nonzero(a::FreeFraction)
    n = dim(a) + 1
    var_names = String[]
    for row in 1:n
        for col in 1:n
            push!(var_names, "a_$(row)_$(col)")
        end
    end
    for row in 1:n
        for col in 1:n
            push!(var_names, "b_$(row)_$(col)")
        end
    end
    SR, varSR = Singular.PolynomialRing(base_ring(a), var_names)
    A = matrix(SR, n, n, varSR[1:(n * n)])
    B = matrix(SR, n, n, varSR[(n * n + 1):(2 * n * n)])
    (detA, detB) = (det(A), det(B))
    L_mats = typeof(A)[]
    L = zero_matrix(SR, n, n)
    L[1, 2:n] = change_base_ring(SR, row_vector(a))
    L[2:n, 1] = change_base_ring(SR, column_vector(a))
    L[2:n, 2:n] = change_base_ring(SR, constant_matrix(a))
    push!(L_mats, L)
    for mat in variable_matrices(a)
        L = zero_matrix(SR, 1, n)
        lowerHalf = hcat(zero_matrix(SR, n - 1, 1), change_base_ring(SR, mat))
        L = vcat(L, lowerHalf)
        push!(L_mats, L)
    end
    @debug "building ALB_mats..."
    ALB_mats = typeof(A)[]
    for L in L_mats
        display(L)
        push!(ALB_mats, A * L * B)
    end
    @debug "beginning ideal tests..."
    for k in 1:n
        @debug "  k is now $k"
        Ik = Ideal(SR, SR(0))
        for ALB_mat in ALB_mats
            for row in 1:k, col in k:n
                Ik += Ideal(SR, ALB_mat[row, col])
            end
        end
        display(Ik)
        Ik = std(Ik)
        Ik = std(Ik + Ideal(SR, detA - 1))
        Ik = std(Ik + Ideal(SR, detB - 1))
        if Singular.reduce(SR(1), Ik) != SR(0)
            return false
        end
    end
    return true
end
