###############################################################################
#
#   MonomialFreeFraction.jl : implementation monomial of free fractions
#
###############################################################################

struct MonomialFreeFraction{T} <: FreeFraction{T}
    coefficient::T
    monomial_indices::Array{<:Int, 1}
    parent::FreeField{T}

    function MonomialFreeFraction(
        coefficient::T,
        monomial_indices::Array{<:Int, 1},
        parent::FreeField{T}
    ) where T
        for index in 1:length(monomial_indices)
            if !(monomial_indices[index] in 1:nvars(parent))
                error("indices in monomial_indices must be in the range of 1
                    to $(nvars(parent)), but index $index contains
                    $(monomial_indices[index]).")
            end
        end
        new{T}(coefficient, monomial_indices, parent)
    end
end

###############################################################################
#
#   Constructors
#
###############################################################################

function MonomialFreeFraction(
    x,
    monomial_indices::Array{<:Int, 1},
    parent::FreeField
)
    return MonomialFreeFraction(base_ring(parent)(x), monomial_indices, parent)
end

function MonomialFreeFraction(
    monomial_indices::Array{<:Int, 1},
    parent::FreeField
)
    return MonomialFreeFraction(1, monomial_indices, parent)
end

function MonomialFreeFraction(x, parent::FreeField)
    return MonomialFreeFraction(x, Int[], parent)
end

function MonomialFreeFraction(a::MonomialFreeFraction)
    return a
end

###############################################################################
#
#   Data type
#
###############################################################################

function promote_rule(
    ::Type{MonomialFreeFraction{T}},
    ::Type{MonomialFreeFraction{T}}
) where T <: FieldElement
    return MonomialFreeFraction{T}
end

function promote_rule(
    ::Type{MonomialFreeFraction{T}},
    ::Type{U}
) where {T <: FieldElement, U <: RingElement}
    promote_rule(T, U) == T ? MonomialFreeFraction{T} : Union{}
end

###############################################################################
#
#   Field getters
#
###############################################################################

function coefficient(a::MonomialFreeFraction)
    return a.coefficient
end

function monomial_indices(a::MonomialFreeFraction)
    return a.monomial_indices
end

function parent(a::MonomialFreeFraction)
    return a.parent
end

###############################################################################
#
#   Hashing
#
###############################################################################

function hash(a::MonomialFreeFraction, h::UInt)
    b = 0xedee5746b948beab%UInt
    b = xor(b, hash(coefficient(a), h))
    b = xor(b, hash(monomial_indices(a), h))
    return xor(b, hash(parent(a), h))
end

###############################################################################
#
#   Copying
#
###############################################################################

function copy(a::MonomialFreeFraction)
    return MonomialFreeFraction(coefficient(a), monomial_indices(a), parent(a))
end

function deepcopy_internal(a::MonomialFreeFraction, dict::IdDict)
    return MonomialFreeFraction(
        deepcopy(coefficient(a)),
        deepcopy(monomial_indices(a)),
        parent(a)
    )
end

###############################################################################
#
#   Information extraction
#
###############################################################################

function monomial_length(a::MonomialFreeFraction)
    return length(monomial_indices(a))
end

function isscalar(a::MonomialFreeFraction)
    return iszero(monomial_length(a))
end

function isregular(::MonomialFreeFraction)
    return true
end

###############################################################################
#
#   Linearization information and manipulation
#
###############################################################################

function row_vector(a::MonomialFreeFraction)
    if iszero(dim(a))
        return zero_matrix(base_ring(a), 1, 0)
    end
    return matrix_unit(base_ring(a), 1, dim(a), 1, 1)
end

function constant_matrix(a::MonomialFreeFraction)
    return identity_matrix(base_ring(a), dim(a))
end

function variable_matrices(a::MonomialFreeFraction{T}) where T
    d = dim(a)
    var_mats = dense_matrix_type(T)[]
    for _ in 1:nvars(a)
        push!(var_mats, zero_matrix(base_ring(a), d, d))
    end
    for varindex in 1:monomial_length(a)
        mat_index = monomial_indices(a)[varindex]
        var_mats[mat_index][varindex, varindex + 1] = base_ring(a)(-1)
    end
    return var_mats
end

function column_vector(a::MonomialFreeFraction)
    if iszero(dim(a))
        return zero_matrix(base_ring(a), 0, 1)
    end
    return matrix_unit(base_ring(a), dim(a), 1, dim(a), 1, coefficient(a))
end

function dim(a::MonomialFreeFraction)
    if iszero(a) && isscalar(a)
        return 0
    end
    return monomial_length(a) + 1
end

function isminimal(a::MonomialFreeFraction)
    if iszero(a) && dim(a) > 0
        return false
    end
    return true
end

function monic_descriptor_realization(a::MonomialFreeFraction)
    return GenericFreeFraction(a)
end

function transpose(a::MonomialFreeFraction)
    return transpose(GenericFreeFraction(a))
end

function admissible!(a::MonomialFreeFraction)
    return a
end

function admissible(a::MonomialFreeFraction)
    return deepcopy(a)
end

function clean_rhs!(a::MonomialFreeFraction)
    return a
end

function clean_rhs(a::MonomialFreeFraction)
    return deepcopy(a)
end

function minimal_candidate(a::MonomialFreeFraction)
    return (true, minimal(a))
end

function minimal(a::MonomialFreeFraction)
    if iszero(a)
        return parent(a)(0)
    end
    return a
end

function minimal_generic(a::MonomialFreeFraction; kwargs...)
    return minimal(a)
end

function minimal_regular(a::MonomialFreeFraction)
    return minimal(a)
end

function minimal_independence(a::MonomialFreeFraction; kwargs...)
    return minimal(a)
end

###############################################################################
#
#   Display functions
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", a::MonomialFreeFraction)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::MonomialFreeFraction)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function expressify(a::MonomialFreeFraction; context = nothing)
    prod = Expr(:call, :*, expressify(coefficient(a); context))
    for idx in monomial_indices(a)
        push!(prod.args, gen_symbol(parent(a), idx))
    end
    return prod
end

function needs_parentheses(a::MonomialFreeFraction)
    return needs_parentheses(coefficient(a))
end

function displayed_with_minus_in_front(a::MonomialFreeFraction)
    return displayed_with_minus_in_front(coefficient(a))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::MonomialFreeFraction)
    return MonomialFreeFraction(-coefficient(a), monomial_indices(a),
        parent(a)
    )
end

function inv(a::MonomialFreeFraction)
    if iszero(a)
        error("division by zero")
    end
    if isscalar(a)
        return MonomialFreeFraction(inv(coefficient(a)), parent(a))
    end
    if intermediate_minimization(parent(a))
        d = dim(a) - 1
        R = base_ring(a)
        row_vec = matrix_unit(R, 1, d, 1, 1)
        const_mat = zero_matrix(R, d, d)
        for col in 2:ncols(const_mat)
            const_mat[col - 1, col] = -R(1)
        end
        var_mats = [zero_matrix(R, d, d)]
        for _ in 2:nvars(a)
            push!(var_mats, zero_matrix(R, d, d))
        end
        for index in 1:monomial_length(a)
            j = d - index + 1
            if index == monomial_length(a)
                var_mats[monomial_indices(a)[index]][j, j] = R(coefficient(a))
            else
                var_mats[monomial_indices(a)[index]][j, j] = R(1)
            end
        end
        col_vec = matrix_unit(R, d, 1, d, 1)
        F = parent(a)
        return GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, F)
    end
    return inv(GenericFreeFraction(a))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(a::MonomialFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    F = determine_parent(a, b)
    if iszero(a)
        if intermediate_minimization(F)
            return minimal(b)
        end
        return b
    end
    if iszero(b)
        if intermediate_minimization(F)
            return minimal(a)
        end
        return a
    end
    if monomial_indices(a) == monomial_indices(b)
        new_coeff = coefficient(a) + coefficient(b)
        if new_coeff == 0
            return parent(a)(0)
        end
        return MonomialFreeFraction(new_coeff, monomial_indices(a), F)
    end
    return GenericFreeFraction(a) + GenericFreeFraction(b)
end

function *(a::MonomialFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    F = determine_parent(a, b)
    if iszero(a) || iszero(b)
        return F(0)
    end
    newCoefficient = coefficient(a) * coefficient(b)
    newMonomialIndices = vcat(monomial_indices(a), monomial_indices(b))
    return MonomialFreeFraction(newCoefficient, newMonomialIndices, F)
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

function +(t::T, a::MonomialFreeFraction{T}) where T<:RingElem # ambiguity
    if iszero(t)
        return a
    elseif iszero(a) || isscalar(a)
        return MonomialFreeFraction(t + coefficient(a), parent(a))
    else
        c = coefficient(a)
        b = GenericFreeFraction(a)
        constant_matrix(b)[1, dim(a)] = -t
        len = monomial_length(a)
        variable_matrices(b)[monomial_indices(a)[len]][len, len + 1] *= c
        column_vector(b)[len + 1, 1] = 1
        return b
    end
end

function +(a::MonomialFreeFraction{T}, t::T) where T<:RingElem # ambiguity
    return t + a
end

function *(t::T, a::MonomialFreeFraction{T}) where T<:RingElem # ambiguity
    if iszero(t)
        return parent(a)(0)
    end
    new_coeff = t * coefficient(a)
    return MonomialFreeFraction(new_coeff, monomial_indices(a), parent(a))
end

function *(a::MonomialFreeFraction{T}, t::T) where T<:RingElem # ambiguity
    if iszero(t)
        return parent(a)(0)
    end
    new_coeff = coefficient(a) * t
    return MonomialFreeFraction(new_coeff, monomial_indices(a), parent(a))
end

###############################################################################
#
#   Vector space base changing
#
###############################################################################

function *(m::MatElem{T}, a::MonomialFreeFraction{T}) where T
    return m * GenericFreeFraction(a)
end

function *(a::MonomialFreeFraction{T}, m::MatElem{T}) where T
    return GenericFreeFraction(a) * m
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::MonomialFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    if parent(a) != parent(b)
        return false
    end
    if coefficient(a) != coefficient(b)
        return false
    elseif iszero(a) # in this case, both a and b are zero
        return true
    end
    return monomial_indices(a) == monomial_indices(b)
end

# function ==(a::MonomialFreeFraction{T}, n::T) where T
#     if isscalar(a) && coefficient(a) == n
#         return true
#     else
#         return false
#     end
# end

function isone(a::MonomialFreeFraction)
    return isscalar(a) && isone(coefficient(a))
end

function iszero(a::MonomialFreeFraction)
    return iszero(coefficient(a))
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function _evaluate(
    a::MonomialFreeFraction{T},
    mode::Symbol,
    args::MatElem{T}...
) where T
    if mode == :left || mode == :right || mode == :matrix
        return _evaluate(GenericFreeFraction(a), mode, args...)
    end
    result = diagonal_matrix(coefficient(a), nrows(args[1]))
    for varindex in monomial_indices(a)
        result *= args[varindex]
    end
    return result
end

function _evaluate(
    a::MonomialFreeFraction{T},
    mode::Symbol,
    args::Union{T, Int, Rational{Int}}...
) where T
    if mode == :left || mode == :right || mode == :matrix
        return _evaluate(GenericFreeFraction(a), mode, args...)
    end
    result = coefficient(a)
    for varindex in monomial_indices(a)
        result *= parent(a)(args[varindex])
    end
    return result
end

###############################################################################
#
#   Linear argument shift
#
###############################################################################

function forward_shift(
        a::MonomialFreeFraction{T},
        points::Union{T, Int, Rational{Int}}...
    ) where T
    return forward_shift(GenericFreeFraction(a), points...)
end
