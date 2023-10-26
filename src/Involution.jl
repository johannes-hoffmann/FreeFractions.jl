###############################################################################
#
#   Involution.jl : Involutions for free fields/fractions
#
###############################################################################

struct Involution{T} <: Function where T<:FieldElem
    var_involution::MatElem{T}
    base_involution::Function
    nvars::Int

    function Involution(var_involution::MatElem{T}, base_involution::Function) where T<:FieldElem
        if !is_square(var_involution)
            error("matrix defining involutive behavior on variables must be square")
        end
        if !is_invertible(var_involution)
            error("matrix defining involutive behavior on variables must be invertible")
        end
        new{T}(var_involution, base_involution, size(var_involution, 1))
    end
end

###############################################################################
#
#   Constructors
#
###############################################################################

function selfadjoint_involution(base_field::Field, n::Int, base_involution::Function = x->x)
    return Involution(identity_matrix(base_field, n), base_involution)
end

###############################################################################
#
#   Field getters
#
###############################################################################

function base_involution(star::Involution)
    return star.base_involution
end

function var_involution(star::Involution)
    return star.var_involution
end

function nvars(star::Involution)
    return star.nvars
end

###############################################################################
#
#   Applying the involution
#
###############################################################################

function (star::Involution{T})(M::MatElem{T}) where T<:FieldElem
    return matrix(base_ring(M), base_involution(star).(transpose(M)))
end

function (star::Involution{T})(a::MonomialFreeFraction{T}) where T<:FieldElem
    return star(GenericFreeFraction(a)) # for now
end

function (star::Involution{T})(a::GenericFreeFraction{T}) where T<:FieldElem
    if nvars(star) != nvars(a)
        error("dimension mismatch: involution has size $(nvars(star)), but the free fraction has $(nvars(a)) variables")
    end
    row_vec = star(column_vector(a))
    col_vec = star(row_vector(a))
    const_mat = star(constant_matrix(a))
    var_mats = dense_matrix_type(T)[]
    for k in 1:nvars(a)
        var_mat = zero_matrix(base_ring(a), dim(a), dim(a))
        for j in 1:nvars(a)
            var_mat += star(variable_matrices(a)[j]) * var_involution(star)[j, k]
        end
        push!(var_mats, var_mat)
    end
    b = GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, parent(a))
    return b
end

###############################################################################
#
#   Self-adjoint minimal representations
#
###############################################################################

function symmetric_representation(a::FreeFraction{T}) where T<:FieldElem
    star = Involution(identity_matrix(base_ring(a), nvars(a)), x->x)
    return selfadjoint_representation(a, star)
end

function selfadjoint_representation(a::FreeFraction{T}, star::Involution{T}) where T<:FieldElem
    (min_flag, b) = minimal_candidate(a)
    input = Dict([(:A, -matrix_pencil(b))])
    input[:B] = matrix_pencil(star(b))
    input[:C] = fill(zero_matrix(base_ring(b), dim(b), dim(b)), nvars(b) + 1)
    input[:M] = [star(row_vector(b))]
    input[:N] = [column_vector(b)]
    (solve_flag, X, Y) = simultaneous_matrix_solve(input)
    if !solve_flag
        if min_flag
            error("no self-adjoint representation found - check if input represents a self-adjoint free fraction")
        end
        error("no self-adjoint representation found, but minimality could not be guaranteed")
    end
    return inv(Y) * b
end