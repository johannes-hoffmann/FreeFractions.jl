# given A, B, C as vectors of matrices,
# find matrices X and Y solving
#     AX + YB = C   (A, B, C are called linking constraints)
# under the optional assumptions
#     XF = G
#     YM = N
#     PX = Q
#     SY = T
# given by vectors of matrices F, G, M, N, P, Q, S, T
function simultaneous_matrix_solve(
    input::Dict{Symbol, <:Vector{<:MatElem{T}}}
) where T
    keys = [:A, :B, :C, :F, :G, :M, :N, :P, :Q, :S, :T]
    if !all(x -> haskey(input, x), [:A, :B, :C])
        error("all linking constraints have to be defined")
    end
    K = base_ring(first(input[:A]))
    for key in keys
        if haskey(input, key)
            first_size = size(first(input[key]))
            for mat in input[key]
                if base_ring(mat) != K
                    error("key $key: all input matrices must share the base
                        field $K")
                end
                if size(mat) != first_size
                    error("all matrices within an array must share the same
                        size")
                end
            end
        else
            setdiff!(keys, [key])
        end
    end
    possible_key_groups = [[:A, :B, :C], [:F, :G], [:M, :N], [:P, :Q], [:S, :T]]
    key_groups = [[:A, :B, :C], [:F, :G], [:M, :N], [:P, :Q], [:S, :T]] # dirty, but works
    for key_group in possible_key_groups
        @debug "current group:" key_group
        if all(x->haskey(input, x), key_group)
            len = length(input[first(key_group)])
            for i in 2:length(key_group)
                if length(input[key_group[i]]) != len
                    error("key group $key_group: all input arrays must have
                        the same length within the key group")
                end
            end
        else
            setdiff!(key_groups, [key_group])
        end
    end
    @debug "final group:" key_groups
    # sizes
    α = ncols(first(input[:A])) # α_indicators = [:A, :G, :P]
    β = ncols(first(input[:B])) # β_indicators = [:B, :C, :F, :Q]
    γ = nrows(first(input[:A])) # γ_indicators = [:A, :C, :N, :S]
    δ = nrows(first(input[:B])) # δ_indicators = [:B, :M, :T]
    check_list = [
        (:F, β, nrows),
        (:G, α, nrows),
        (:M, δ, nrows),
        (:N, γ, nrows),
        (:P, α, ncols),
        (:Q, β, ncols),
        (:S, γ, ncols),
        (:T, δ, ncols)
    ]
    for (key, val, func) in check_list
        if haskey(input, key) && val != func(first(input[key]))
            error("wrong dimensions in matrix $key")
        end
    end
    optional_conditions = [
        [:F, :G, ncols],
        [:M, :N, ncols],
        [:P, :Q, nrows],
        [:S, :T, nrows]
    ]
    for (l, r, func) in optional_conditions
        if [l, r] in key_groups
            if func(first(input[l])) != func(first(input[r]))
                error("$func must be the same for $l and $r")
            end
        end
    end
    # create system with linking constraint
    upper = zero_matrix(K, α * β, 0)
    lower = zero_matrix(K, γ * δ, 0)
    rhs = zero_matrix(K, 1, 0)
    for mat in input[:A]
        temp = kronecker_product(transpose(mat), identity_matrix(K, β))
        upper = hcat(upper, temp)
    end
    for mat in input[:B]
        lower = hcat(lower, fill_diagonal(mat, γ))
    end
    @assert ncols(upper) == ncols(lower)
    for mat in input[:C]
        for row in 1:nrows(mat)
            rhs = hcat(rhs, mat[row, :])
        end
    end
    if [:F, :G] in key_groups # if F,G present, add column constraint for X
        for mat in input[:F]
            upper = hcat(upper, fill_diagonal(mat, α))
        end
        col_number = α * ncols(first(input[:F])) * length(input[:F])
        lower = hcat(lower, zero_matrix(K, γ * δ, col_number))
        for mat in input[:G]
            for row in 1:nrows(mat)
                rhs = hcat(rhs, mat[row, :])
            end
        end
        @assert ncols(upper) == ncols(lower)
    end
    if [:M, :N] in key_groups # if M,N present, add column constraint for Y
        col_number = γ * ncols(first(input[:M])) * length(input[:M])
        upper = hcat(upper, zero_matrix(K, α * β, col_number))
        for mat in input[:M]
            lower = hcat(lower, fill_diagonal(mat, γ))
        end
        for mat in input[:N]
            for row in 1:nrows(mat)
                rhs = hcat(rhs, mat[row, :])
            end
        end
        @assert ncols(upper) == ncols(lower)
    end
    if [:P, :Q] in key_groups # if P,Q present, add row constraint for X
        for mat in input[:P]
            temp = kronecker_product(transpose(mat), identity_matrix(K, β))
            upper = hcat(upper, temp)
        end
        col_number = β * nrows(first(input[:P])) * length(input[:P])
        lower = hcat(lower, zero_matrix(K, γ * δ, col_number))
        for mat in input[:Q]
            for row in 1:nrows(mat)
                rhs = hcat(rhs, mat[row, :])
            end
        end
        @assert ncols(upper) == ncols(lower)
    end
    if [:S, :T] in key_groups # if S,T present, add row constraint for Y
        col_number = δ * nrows(first(input[:S])) * length(input[:S])
        upper = hcat(upper, zero_matrix, α * β, col_number)
        for mat in input[:S]
            temp = kronecker_product(transpose(mat), identity_matrix(K, δ))
            lower = hcat(lower, temp)
        end
        for mat in input[:T]
            for row in 1:nrows(mat)
                rhs = hcat(rhs, mat[row, :])
            end
        end
        @assert ncols(upper) == ncols(lower)
    end
    lhs = vcat(upper, lower)
    # solve system
    X = zero_matrix(K, α, β)
    Y = zero_matrix(K, γ, δ)
    x = nothing
    kernel_dimension = left_kernel(lhs)[1]
    if kernel_dimension > 0
        @warn "left kernel dimension: $kernel_dimension"
    else
        @debug "left kernel dimension: $kernel_dimension"
    end
    try
        x = solve_left(lhs, rhs)
    catch
        return (false, X, Y)
    end
    for row in 1:α, col in 1:β
        X[row, col] = x[1, β * (row - 1) + col]
    end
    for row in 1:γ, col in 1:δ
        Y[row, col] = x[1, δ * (row - 1) + col + γ * δ]
    end
    return (true, X, Y)
end

function sms_test(f::FreeFraction, g::FreeFraction)
    input = Dict([(:A, -matrix_pencil(f))])
    input[:B] = matrix_pencil(g)
    input[:C] = fill(zero_matrix(base_ring(f), dim(f), dim(g)), nvars(f) + 1)
    input[:P] = [row_vector(f)]
    input[:Q] = [row_vector(g)]
    input[:M] = [column_vector(g)]
    input[:N] = [column_vector(f)]
    (flag, U, T) = simultaneous_matrix_solve(input)
    if !flag
        error("no solution found")
    end
    u_cond = row_vector(f) * U == row_vector(g)
    t_cond = transpose(T * column_vector(g) == column_vector(f))
    l_cond = matrix_pencil(T * g) == matrix_pencil(f * U)
    e_cond = T == U
    if !(u_cond && t_cond && l_cond)
        @error "wrong solution"
    end
    if e_cond
        return T
    end
    return (T, U)
end
