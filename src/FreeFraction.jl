###############################################################################
#
#   Freefraction.jl : Elements of free (skew) fields
#
###############################################################################

# abstract type FreeFraction{T<:FieldElem} <: NCRingElem end

include("MonomialFreeFraction.jl")
include("GenericFreeFraction.jl")

###############################################################################
#
#   Constructors
#
###############################################################################

function FreeFraction(
    s::AbstractString,
    input_field::Field,
    input_vars::Array{String, 1};
    minimize::Bool = true
)
    F, F_vars = FreeField(input_field, input_vars, minimize)
    vardict = Dict{String, MonomialFreeFraction{elem_type(input_field)}}()
    for varindex in 1:nvars(F)
        vardict[input_vars[varindex]] = F_vars[varindex]
    end
    return F(eval(replaceVariables(Meta.parse(s), vardict)))
end

function replaceVariables(ex::Union{Symbol,Expr,Number}, d::Dict{String,<:Any})
    if typeof(ex) == Symbol
        if String(ex) in keys(d)
            return d[String(ex)]
        else
            return ex
        end
    end
    if typeof(ex) == Expr # recursively replace in ex.args
        return Expr(ex.head, map(x -> replaceVariables(x, d), ex.args)...)
    end
    return ex
end

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

function parent_type(::Type{<:FreeFraction{T}}) where T
    return FreeField{T}
end

function elem_type(::Type{FreeField{T}}) where T
    return FreeFraction{T}
end

function isexact_type(::Type{T}) where {S<:FieldElem, T<:FreeFraction{S}}
    return isexact_type(S)
end

function show_minus_one(::Type{T}) where {S<:FieldElem, T<:FreeFraction{S}}
    return show_minus_one(S)
end

function check_parent(a::FreeFraction, b::FreeFraction, throw::Bool = true)
    c = parent(a) == parent(b)
    if throw && !c
        error("Incompatible free fields in operation")
    end
    return c
end

###############################################################################
#
#   Information extraction
#
###############################################################################

function base_ring(a::FreeFraction)
    return base_ring(parent(a))
end

function nvars(a::FreeFraction)
    return nvars(parent(a))
end

function display_ring(a::FreeFraction)
    return display_ring(parent(a))
end

function matrix_pencil(a::FreeFraction)
    return vcat([constant_matrix(a)], variable_matrices(a))
end

function extended_matrix_pencil(a::FreeFraction)
    R = base_ring(a)
    A0 = hcat(zero_matrix(R, 1, 1), row_vector(a))
    A0 = vcat(A0, hcat(column_vector(a), constant_matrix(a)))
    pencil = [A0]
    for mat in variable_matrices(a)
        push!(pencil, block_diagonal_matrix([zero_matrix(R, 1, 1), mat]))
    end
    return pencil
end

function extended_matrix_pencil_tensor(a::FreeFraction)
    pencil = extended_matrix_pencil(a)
    R = base_ring(a)
    tensor = zeros(R, dim(a) + 1, dim(a) + 1, nvars(a) + 1)
    for varindex in eachindex(pencil)
        tensor[:, :, varindex] .= pencil[varindex]
    end
    return tensor
end

function determine_parent(a::FreeFraction{T}, b::FreeFraction{T}) where T
    if !check_parent(a, b, false)
        error("could not determine a common parent")
    end
    if !intermediate_minimization(parent(a))
        return parent(a)
    end
    return parent(b)
end

###############################################################################
#
#   Display functions
#
###############################################################################

function needs_parentheses(a::FreeFraction)
    return true
end

function displayed_with_minus_in_front(a::FreeFraction)
    return false
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function isunit(a::FreeFraction)
    return !iszero(a)
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function -(a::FreeFraction{T}, b::FreeFraction{T}) where T
    return a + (-b)
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

function +(n::Union{Integer, Rational}, a::FreeFraction)
    return base_ring(a)(n) + a
end

function +(a::FreeFraction, n::Union{Integer, Rational})
    return a + base_ring(a)(n)
end

function -(a::FreeFraction{T}, n::Union{Integer, Rational}) where T
    return a + (-n)
end

function -(n::Union{Integer, Rational}, a::FreeFraction{T}) where T
    return n + (-a)
end

function *(n::Union{Integer, Rational}, a::FreeFraction)
    return base_ring(a)(n) * a
end

function *(a::FreeFraction, n::Union{Integer, Rational})
    return a * base_ring(a)(n)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::FreeFraction, n::Integer)
    if n < 0
        return inv(a^(-n))
    elseif n == 0
        return one(a)
    elseif n == 1
        return a
    else # n is at least 2
        b = a^div(n, 2)
        b *= b
        if isodd(n)
            b *= a
        end
        return b
    end
end

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function add!(
    c::FreeFraction{T},
    a::FreeFraction{T},
    b::FreeFraction{T}
) where T
    c = a + b
    return c
end

function addeq!(c::FreeFraction{T}, a::FreeFraction{T}) where T
    c += a
    return c
end

function addmul!(
    c::FreeFraction{T},
    a::FreeFraction{T},
    b::FreeFraction{T},
    t::FreeFraction{T}
) where T
    c += a * b
    return c
end

function mul!(
    c::FreeFraction{T},
    a::FreeFraction{T},
    b::FreeFraction{T}
) where T
    c = a * b
    return c
end

function zero!(c::FreeFraction)
    return zero(parent(c))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact_left(a::FreeFraction, b::FreeFraction)
    check_parent(a, b)
    if iszero(b)
        if iszero(a)
            error("no unique quotient exists")
        else
            throw(ErrorException("division by zero"))
        end
    end
    return inv(b) * a
end

function divexact_right(a::FreeFraction, b::FreeFraction)
    check_parent(a, b)
    if iszero(b)
        if iszero(a)
            error("no unique quotient exists")
        else
            throw(ErrorException("division by zero"))
        end
    end
    return a * inv(b)
end

###############################################################################
#
#   Comparison
#
###############################################################################

###############################################################################
#
#   Adhoc comparison
#
###############################################################################

function ==(a::FreeFraction, n::Union{Integer, Rational})
    return a == base_ring(a)(n)
end

function ==(n::Union{Integer, Rational}, a::FreeFraction)
    return a == base_ring(a)(n)
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate_left_family(a::FreeFraction, args...)
    return _evaluate(a, :left, args...)
end

function evaluate_right_family(a::FreeFraction, args...)
    return _evaluate(a, :right, args...)
end

function evaluate_matrix(a::FreeFraction, args...)
    return _evaluate(a, :matrix, args...)
end

function evaluate(a::FreeFraction, args...)
    return _evaluate(a, :both, args...)
end

function (a::FreeFraction)(args...)
    return evaluate(a, args...)
end

function (a::FreeFraction)(args::Array{<:Any, 1})
    return a(args...)
end

function evaluate_random(
    a::FreeFraction;
    mat_size::Int = 5
)
    test_mats = []
    R = base_ring(a)
    for var_index in 1:nvars(a)
        push!(test_mats, matrix(R, FFrand(R, mat_size, mat_size)))
    end
    return a(test_mats...)
end

###############################################################################
#
#   Linear argument shift
#
###############################################################################

function forward_shift(a::FreeFraction, args::Array{<:Any, 1})
    return forward_shift(a, args...)
end

function backward_shift(a::FreeFraction, b...)
    return forward_shift(a, map(x -> -x, b)...)
end
