module FreeFractions

using Nemo
import AbstractAlgebra

try
    import Singular: Ideal, std
    include("Singular.jl")
    export is_nonzero
catch e
    m = "Singular not available, additional functionality based on Gröbner
        technology deactivated."
    @warn m
    @debug "Original exception was:\n" * e.msg
end

import Base: +, -, *, ^, ==
import Base: adjoint, convert, copy, deepcopy_internal, hash, inv, isequal,
    isone, iszero, one, setindex!, show, transpose, zero

import AbstractAlgebra: base_ring, dim, elem_type, evaluate, gen, gens, nvars,
    parent, parent_type, promote_rule,
    Field, NCRing, Ring,
    add!, addeq!, mul!, zero!,
    add_column!, add_column, add_row!, add_row,
    multiply_column!, multiply_column, multiply_row!, multiply_row,
    swap_cols!, swap_cols, swap_rows!, swap_rows,
    expressify,
    check_square

import AbstractAlgebra.Generic: characteristic, check_parent,
    elements!, emb, isexact_type, isunit, solve_left

export FreeField, FreeFraction, GenericFreeFraction, MonomialFreeFraction
export constant_matrix, dim, variable_matrices, row_vector, column_vector
export matrix_pencil, extended_matrix_pencil, extended_matrix_pencil_tensor
export gen, gens, characteristic, nvars, parent
export forward_shift, backward_shift, scalar, coefficient
export evaluate_left_family, evaluate_right_family, evaluate_matrix, evaluate
export isregular, isminimal, isunit
export minimal, minimal_generic, minimal_regular, minimal_independence
export pivot_refinement!, pivot_refinement, block_upper_triangular_structure,
    issymmetric, ispersymmetric, clean_rhs, clean_rhs!

abstract type FreeFraction{T<:FieldElem} <: NCRingElem end

include("AAA.jl")
include("FreeField.jl")
include("FreeFraction.jl")

include("Involution.jl")

include("RegularMinimization.jl")
include("GenericMinimization.jl")
include("IndependenceMinimization.jl")


include("MatrixSolve.jl")

end # module
