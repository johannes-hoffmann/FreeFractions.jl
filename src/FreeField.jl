###############################################################################
#
#   FreeField.jl : Free (skew) fields
#
###############################################################################

struct FreeField{T} <: NCRing where T<:FieldElem
    base_ring::Field
    variable_symbols::Array{Symbol, 1} # variables
    nvars::Int # number of variables
    intermediate_minimization::Bool

    function FreeField(
        base_ring::Field,
        variable_symbols::Array{Symbol, 1},
        intermediate_minimization::Bool = true
    )
        T = elem_type(base_ring)
        n = length(variable_symbols)
        new{T}(base_ring, variable_symbols, n, intermediate_minimization)
    end
end

###############################################################################
#
#   Constructors
#
###############################################################################

function FreeField(
    R::Field,
    input_variables::Array{String, 1},
    minimize::Bool = true
)
    symbol_variables = [Symbol(x) for x in input_variables]
    new_free_field = FreeField(R, symbol_variables, minimize)
    return tuple(new_free_field, gens(new_free_field))
end

###############################################################################
#
#   Field getters
#
###############################################################################

function base_ring(F::FreeField)
    return F.base_ring
end
function symbols(F::FreeField)
    return F.variable_symbols
end
function nvars(F::FreeField)
    return F.nvars
end
function intermediate_minimization(F::FreeField)
    return F.intermediate_minimization
end

###############################################################################
#
#   Information extraction
#
###############################################################################

function characteristic(F::FreeField)
    return characteristic(base_ring(F))
end

function gen_symbol(F::FreeField, varindex::Int)
    return symbols(F)[varindex]
end

function gen_string(F::FreeField, varindex::Int)
    return String(gen_symbol(F, varindex))
end

function gen(F::FreeField, varindex::Int)
    return MonomialFreeFraction(base_ring(F)(1), [varindex], F)
end

function gens(F::FreeField)
    return [gen(F, varindex) for varindex in 1:nvars(F)]
end

function display_ring(F::FreeField)
    R, vars = PolynomialRing(base_ring(F), [String(x) for x in symbols(F)])
    return R
end

###############################################################################
#
#   Display functions
#
###############################################################################

function show(io::IO, F::FreeField)
    output = "Free Field in $(String(gen_symbol(F, 1)))"
    for varindex in 2:nvars(F)
        output *= ", $(String(gen_symbol(F, varindex)))"
    end
    output *= " over $(base_ring(F))"
    print(io, output)
    return
end

###############################################################################
#
#   Create basic objects
#
###############################################################################

function (F::FreeField{T})(t::T) where T
    return MonomialFreeFraction(t, F)
end

function (F::FreeField)(n::Integer)
    return F(base_ring(F)(n))
end

function (F::FreeField)(n::Number)
    return F(base_ring(F)(n))
end

function (F::FreeField)()
    return F(0)
end

function (F::FreeField{T})(a::FreeFraction{T}) where T
    if parent(a) != F
        error("Unable to coerce free fraction")
    end
    return a
end

function one(F::FreeField)
    return F(1)
end

function zero(F::FreeField)
    return F(0)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(F::FreeField{T}, G::FreeField{T}) where T
    return base_ring(F) == base_ring(G) && symbols(F) == symbols(G)
end
