###############################################################################
#
#   GenericFreeFraction.jl : generic implementation of free fractions
#
###############################################################################

struct GenericFreeFraction{T} <: FreeFraction{T}
    row_vector::MatElem{T}
    constant_matrix::MatElem{T}
    variable_matrices::Array{K,1} where K<:MatElem{T}
    column_vector::MatElem{T}
    parent::FreeField{T}

    function GenericFreeFraction(
        row_vector::MatElem{T},
        constant_matrix::MatElem{T},
        variable_matrices::Array{K,1},
        column_vector::MatElem{T},
        parent::FreeField{T}
    ) where K<:MatElem{T} where T<:FieldElem
        if nrows(row_vector) != 1
            error("row_vector must be a matrix with exactly one row")
        end
        if ncols(column_vector) != 1
            error("column_vector must be a matrix with exactly one column")
        end
        if !issquare(constant_matrix)
            error("constant_matrix has to be a quadratic matrix")
        end
        if ncols(constant_matrix) != ncols(row_vector)
            e = "row_vector has to have the same number of columns as "
            e *= "constant_matrix"
            error(e)
        end
        if nrows(constant_matrix) != nrows(column_vector)
            e = "column_vector has to have the same number of rows as "
            e *= "constant_matrix"
            error(e)
        end
        if length(variable_matrices) != nvars(parent)
            e = "the number of matrices in variable_matrices must match the "
            e *= "number of variables in parent"
            error(e)
        end
        for mat in variable_matrices
            if size(mat) != size(constant_matrix)
                e = "all matrices in variable_matrices have to be of the "
                e *= "same size as constant_matrix"
                error(e)
            end
        end
        new{T}(row_vector, constant_matrix, variable_matrices, column_vector,
            parent)
    end
end

###############################################################################
#
#   Constructors
#
###############################################################################

function GenericFreeFraction(
    row_vec::MatElem{T},
    var_mats::Array{K,1},
    col_vec::MatElem{T},
    parent::FreeField{T}
) where K<:MatElem{T} where T<:FieldElem
    const_mat = var_mats[1]
    var_mats = deleteat!(var_mats,1)
    return GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, parent)
end

function GenericFreeFraction(a::MonomialFreeFraction)
    return GenericFreeFraction(row_vector(a), constant_matrix(a),
        variable_matrices(a), column_vector(a), parent(a))
end

###############################################################################
#
#   Additional constructors for other linearization types
#
###############################################################################

function MonomialFreeFraction(a::GenericFreeFraction)
    b = a
    if intermediate_minimization(parent(a))
        b = minimal(b)
    end
    if all(x -> x == 0, variable_matrices(b))
        if rank(constant_matrix(b)) < dim(b)
            error("invalid free fraction")
        end
        t = (row_vector(b) * inv(constant_matrix(b)) * column_vector(b))[1, 1]
        return parent(b)(t)
    end
    if intermediate_minimization(parent(a))
        indices = Int[]
        rem = b
        for pos in 1:(dim(b) - 1)
            p = rand(SymmetricGroup(nvars(b)))
            for index in 1:nvars(b)
                cand = inv(gen(parent(b), p[index])) * rem
                if dim(cand) == dim(rem) - 1
                    push!(indices, p[index])
                    rem = cand
                    break
                end
            end
            if length(indices) != pos
                break
            end
        end
        if isa(rem, MonomialFreeFraction)
            return MonomialFreeFraction(indices, parent(b)) * rem
        end
    end
    error("could not coerce input into a MonomialFreeFraction")
end

###############################################################################
#
#   Data type
#
###############################################################################

function promote_rule(
    ::Type{GenericFreeFraction{T}},
    ::Type{GenericFreeFraction{T}}
) where T <: FieldElement
    return GenericFreeFraction{T}
end

function promote_rule(
    ::Type{GenericFreeFraction{T}},
    ::Type{U}
) where {T <: FieldElement, U <: RingElement}
    promote_rule(T, U) == T ? GenericFreeFraction{T} : Union{}
end

###############################################################################
#
#   Field getters
#
###############################################################################

function row_vector(a::GenericFreeFraction)
    return a.row_vector
end

function constant_matrix(a::GenericFreeFraction)
    return a.constant_matrix
end

function variable_matrices(a::GenericFreeFraction)
    return a.variable_matrices
end

function column_vector(a::GenericFreeFraction)
    return a.column_vector
end

function parent(a::GenericFreeFraction)
    return a.parent
end

###############################################################################
#
#   Hashing
#
###############################################################################

function hash(a::GenericFreeFraction, h::UInt)
    b = 0xcc83cb1dc1ec3cf2%UInt
    b = xor(b, hash(row_vector(a), h))
    b = xor(b, hash(constant_matrix(a), h))
    b = xor(b, hash(variable_matrices(a), h))
    b = xor(b, hash(column_vector(a), h))
    return xor(b, hash(parent(a), h))
end

###############################################################################
#
#   Copying
#
###############################################################################

function copy(a::GenericFreeFraction)
    return GenericFreeFraction(
        row_vector(a),
        constant_matrix(a),
        variable_matrices(a),
        column_vector(a),
        parent(a)
    )
end
function deepcopy_internal(a::GenericFreeFraction, dict::IdDict)
    return GenericFreeFraction(
        deepcopy(row_vector(a)),
        deepcopy(constant_matrix(a)),
        deepcopy(variable_matrices(a)),
        deepcopy(column_vector(a)),
        parent(a)
    )
end

###############################################################################
#
#   Information extraction
#
###############################################################################

function isregular(a::GenericFreeFraction)
    return dim(a) == rank(constant_matrix(a))
end

###############################################################################
#
#   Linearization information and manipulation
#
###############################################################################

function representation_matrix(a::FreeFraction)
    R = display_ring(parent(a))
    rep_mat = change_base_ring(R, constant_matrix(a))
    vars = gens(R)
    var_mats = variable_matrices(a)
    for varindex in 1:nvars(a)
        rep_mat += vars[varindex] * change_base_ring(R, var_mats[varindex])
    end
    return rep_mat
end

function dim(a::GenericFreeFraction)
    return size(constant_matrix(a), 1)
end

function isminimal(a::GenericFreeFraction; kwargs...)
    if dim(a) == 0
        return true
    elseif dim(a) == 1
        if iszero(a)
            return false
        end
        return true
    end
    try
        if isminimal_by_family_independence(a)
            return true
        end
    catch
        # in case of failure, continue with other tests
    end
    (flag, candidate) = minimal_candidate(a; kwargs...)
    if dim(candidate) < dim(a)
        return false
    end # now a and candidate have the same dimension
    if flag
        return true
    end # no information available
    return missing
end

function isminimal_by_family_independence(a::GenericFreeFraction)
    flag = left_family_independence(a)[1]
    return flag && right_family_independence(a)[1]
end

function monic_descriptor_realization(a::GenericFreeFraction)
    if !isregular(a)
        error("monic descriptor realizations exist only for regular free
            fractions")
    end
    const_mat_inv = inv(constant_matrix(a))
    return GenericFreeFraction(
        row_vector(a),
        identity_matrix(base_ring(a), dim(a)),
        map(mat -> const_mat_inv * mat, variable_matrices(a)),
        const_mat_inv * column_vector(a),
        parent(a)
    )
end

function transpose(a::GenericFreeFraction)
    return GenericFreeFraction(
        transpose(column_vector(a)),
        transpose(constant_matrix(a)),
        map(mat -> transpose(mat), variable_matrices(a)),
        transpose(row_vector(a)),
        parent(a)
    )
end

function admissible!(a::GenericFreeFraction)
    if dim(a) == 0
        return a
    end
    if row_vector(a) == matrix_unit(base_ring(a), 1, dim(a), 1, 1)
        return a
    end
    if row_vector(a) == 0
        row_vector(a)[1, 1] = base_ring(a)(1)
        column_vector(a)[:, :] = zero_matrix(base_ring(a), dim(a), 1)
        return a # parent(a)(0)
    end
    if row_vector(a)[1, 1] == 0
        for col in 2:dim(a)
            if row_vector(a)[1, col] != 0 # at least one will be non-zero
                swap_cols!(a, 1, col)
                break
            end
        end
    end
    multiply_column!(a, inv(row_vector(a)[1, 1]), 1)
    for col in 2:dim(a)
        if row_vector(a)[1, col] != 0
            add_column!(a, -row_vector(a)[1, col], 1, col)
        end
    end
    return a
end

function admissible(a::GenericFreeFraction)
    b = deepcopy(a)
    return admissible!(b)
end

function clean_rhs!(a::GenericFreeFraction)
    d = dim(a)
    if column_vector(a) == matrix_unit(base_ring(a), d, 1, d, 1)
        return a
    end
    if column_vector(a) == 0
        return a
    end
    if column_vector(a)[d, 1] == 0 # at least one will be non-zero
        for row in (d - 1):-1:1
            if column_vector(a)[row, 1] != 0
                swap_rows!(a, row, d)
                break
            end
        end
    end
    multiply_row!(a, inv(column_vector(a)[d, 1]), d)
    for row in 1:(d - 1)
        if column_vector(a)[row, 1] != 0
            add_row!(a, -column_vector(a)[row, 1], d, row)
        end
    end
    return a
end

function clean_rhs(a::GenericFreeFraction)
    b = deepcopy(a)
    return clean_rhs!(b)
end

# docs: no testing for correctness after independence
function minimal(
    a::GenericFreeFraction;
    warn::Bool = true,
    independence::Bool = true,
    verify::Bool = false,
    kwargs...
)
    try
        if isminimal_by_family_independence(a)
            return a
        end
    catch
        # in case of failure, continue with other tests
    end
    (flag, candidate) = minimal_candidate(a; independence, kwargs...)
    if warn && !flag
        @warn "minimality of the returned candidate could not be guaranteed"
    end
    return candidate
end

function minimal_candidate(
    a::GenericFreeFraction;
    regular::Bool = false,
    generic::Bool = true,
    independence::Bool = true,
    independence_max_size::Int = 10,
    kwargs...
)
    if iszero(row_vector(a)) || iszero(column_vector(a))
        return (true, parent(a)(0))
    end
    if dim(a) == 1 # TODO corner case zero!
        return (true, a) # TODO deepcopy?
    end
    if regular && isregular(a)
        return (true, minimal_regular(a))
    end
    b = deepcopy(a)
    if generic
        b = minimal_generic(a; warn = false, kwargs...)
        if regular && isregular(b)
            return minimal_candidate(b)
        end
        if number_of_pivot_blocks(b) == dim(b)
            return (true, b)
        end
    end
    if independence
        try
            b = minimal_independence(b, max_size = independence_max_size)
            if isminimal_by_family_independence(b)
                return (true, b)
            end
        catch
            # in case of failure, continue with other tests
        end
    end
    if dim(b) < dim(a)
        return minimal_candidate(b; regular, generic, independence, kwargs...)
    end
    return (false, b)
end

###############################################################################
#
#   Display functions
#
###############################################################################

function show(io::IO, a::GenericFreeFraction)
    output = "Free fraction over $(parent(a)) "
    output *= "represented in dimension $(dim(a)) by "
    output *= "u = $(row_vector(a)), C = $(constant_matrix(a)), "
    output *= "Q = $(variable_matrices(a)), v = $(column_vector(a))"
    print(io, output)
end

function show(io::IO, ::MIME"text/plain", a::GenericFreeFraction)
    d = dim(a)
    col_vec = column_vector(a)
    output = "Free fraction in $(String(gen_symbol(parent(a), 1)))"
    for varindex in 2:nvars(a)
        output *= ", $(String(gen_symbol(parent(a), varindex)))"
    end
    output *= " over $(base_ring(a)) represented in dimension $d by\n"
    print(io, output)
    if d == 0
        print(io, "[ ] [ ]⁻¹ [ ]")
        return
    end
    output = ""
    # build ustring
    ustring = "["
    for col in 1:d
        ustring *= " $(row_vector(a)[1, col])"
    end
    ustring *= " ] "
    left_margin = length(ustring)
    # calculate column widths
    rep_mat = representation_matrix(a)
    col_widths = zeros(Integer, d)
    for row in 1:d, col in 1:d
        col_widths[col] = max(col_widths[col], length(string(rep_mat[row, col])))
    end
    #calculate width of the column_vector
    col_vec_width = 0
    for row in 1:d
        col_vec_width = max(col_vec_width, length(string(col_vec[row, 1])))
    end
    # build output
    if d == 1
        output = ustring * "[ $(rep_mat[1,1]) ]⁻¹ [ $(col_vec[1, 1]) ]"
        print(io, output)
        return
    end
    for row in 1:d
        output = ""
        if row == 1 # before matrix entries
            output = lpad("[", left_margin + 1)
        elseif row == d
            output = ustring * "["
        else
            output = lpad("[", left_margin + 1)
        end
        for col in 1:d # matrix entries
            output *= lpad(rep_mat[row, col], col_widths[col] + 1)
        end
        if row == 1 # after matrix entries
            output *= " ]⁻¹ [ " * lpad(col_vec[row, 1], col_vec_width) * " ]\n"
        elseif row == d
            output *= " ]   [ " * lpad(col_vec[row, 1], col_vec_width) * " ]"
        else
            output *= " ]   [ " * lpad(col_vec[row, 1], col_vec_width) * " ]\n"
        end
        print(io, output)
    end
end

###############################################################################
#
#   Unary operations
#
###############################################################################

function -(a::GenericFreeFraction)
    return parent(a)(-1) * a
end

function inv(a::GenericFreeFraction{T}) where T
    # if a = (u, Q, v), then inv(a) is
    #          [ -v Q ]  [ 0 ]
    # [ 1 0 ], [  0 u ], [ 1 ]
    d = dim(a)
    R = base_ring(a)
    row_vec = matrix_unit(R, 1, d + 1, 1, 1)
    col_vec = matrix_unit(R, d + 1, 1, d + 1, 1)
    upper = hcat(-column_vector(a), constant_matrix(a))
    lower = hcat(zero_matrix(R, 1, 1), row_vector(a))
    const_mat = vcat(upper, lower)
    var_mats = dense_matrix_type(T)[]
    for mat in variable_matrices(a)
        upper = hcat(zero_matrix(R, d, 1), mat)
        var_mat = vcat(upper, zero_matrix(R, 1, d + 1))
        push!(var_mats, var_mat)
    end
    F = parent(a)
    res = GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, F)
    if intermediate_minimization(F)
        res = minimal(res; warn = false, independence = false, min_inv = false)
    end
    return res
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(a::GenericFreeFraction{T}, b::GenericFreeFraction{T}) where T
    # if a = (a.u, a.Q, a.v) and similar for b, then a+b is
    #            [ a.Q -a.Q*transpose(a.u)*b.u ]  [ a.v ]
    # [ a.u 0 ], [  0            b.Q           ], [ b.v ]
    F = determine_parent(a, b)
    R = base_ring(a)
    row_vec = hcat(row_vector(a), zero_matrix(R, 1, dim(b)))
    col_vec = vcat(column_vector(a), column_vector(b))
    upper_right = transpose(row_vector(a)) * row_vector(b)
    # build new const_mat
    upper = hcat(constant_matrix(a), -constant_matrix(a) * upper_right)
    lower = hcat(zero_matrix(R, dim(b), dim(a)), constant_matrix(b))
    const_mat = vcat(upper, lower)
    # build new var_mats
    var_mats = dense_matrix_type(T)[]
    (var_matsA, var_matsB) = (variable_matrices(a), variable_matrices(b))
    for varindex in 1:nvars(a)
        upper = hcat(var_matsA[varindex], -var_matsA[varindex] * upper_right)
        lower = hcat(zero_matrix(R, dim(b), dim(a)), var_matsB[varindex])
        push!(var_mats, vcat(upper, lower))
    end
    res = GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, F)
    if intermediate_minimization(F)
        res = minimal(res; warn = false, independence = false)
    end
    return res
end

function *(a::GenericFreeFraction{T}, b::GenericFreeFraction{T}) where T
    # if a = (a.u, a.Q, a.v) and similar for b, then a*b is
    #            [ a.Q -a.v*b.u ]  [  0  ]
    # [ a.u 0 ], [  0     b.Q   ], [ b.v ]
    F = determine_parent(a, b)
    R = base_ring(a)
    row_vec = hcat(row_vector(a), zero_matrix(R, 1, dim(b)))
    col_vec = vcat(zero_matrix(R, dim(a), 1), column_vector(b))
    upper = hcat(constant_matrix(a), -column_vector(a) * row_vector(b))
    lower = hcat(zero_matrix(R, dim(b), dim(a)), constant_matrix(b))
    const_mat = vcat(upper, lower)
    var_mats = dense_matrix_type(T)[]
    (var_matsA, var_matsB) = (variable_matrices(a), variable_matrices(b))
    for varindex in 1:nvars(a)
        push!(var_mats, block_diagonal_matrix([var_matsA[varindex], var_matsB[varindex]]))
    end
    res = GenericFreeFraction(row_vec, const_mat, var_mats, col_vec, F)
    if intermediate_minimization(F)
        res = minimal(res; warn = false, independence = false)
    end
    return res
end

###############################################################################
#
#   Ad hoc binary operations
#
###############################################################################

function +(a::MonomialFreeFraction{T}, b::GenericFreeFraction{T}) where T
    return GenericFreeFraction(a) + b
end

function +(a::GenericFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    return a + GenericFreeFraction(b)
end

function *(t::T, a::GenericFreeFraction{T}) where T<:RingElem # ambiguity
    if iszero(t)
        return parent(a)(0)
    else
        return GenericFreeFraction(row_vector(a), constant_matrix(a),
            variable_matrices(a), t * column_vector(a), parent(a))
    end
end

function *(a::GenericFreeFraction{T}, t::T) where T<:RingElem # ambiguity
    if iszero(t)
        return parent(a)(0)
    else
        return GenericFreeFraction(row_vector(a), constant_matrix(a),
            variable_matrices(a), column_vector(a) * t, parent(a))
    end
end

function *(a::MonomialFreeFraction{T}, b::GenericFreeFraction{T}) where T
    if isscalar(a)
        return coefficient(a) * b
    end
    return GenericFreeFraction(a) * b
end

function *(a::GenericFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    if isscalar(b)
        return a * coefficient(b)
    end
    return a * GenericFreeFraction(b)
end

###############################################################################
#
#   Vector space base changing
#
###############################################################################

function *(m::MatElem{T}, a::GenericFreeFraction{T}) where T
    if base_ring(m) != base_ring(a)
        error("incompatible base rings")
    end
    if !issquare(m)
        error("only available for square matrices")
    end
    if ncols(m) != dim(a)
        throw(DimensionMismatch("number of columns of the matrix does not match
            dimension of the free fraction representation"))
    end
    return GenericFreeFraction(
        row_vector(a),
        m * constant_matrix(a),
        map(mat -> m * mat, variable_matrices(a)),
        m * column_vector(a),
        parent(a)
    )
end

function *(a::GenericFreeFraction{T}, m::MatElem{T}) where T
    if base_ring(a) != base_ring(m)
        error("incompatible base rings")
    end
    if !issquare(m)
        error("only available for square matrices")
    end
    if nrows(m) != dim(a)
        throw(DimensionMismatch("number of rows of the matrix does not match
            dimension of the free fraction representation"))
    end
    return GenericFreeFraction(
        row_vector(a) * m,
        constant_matrix(a) * m,
        map(mat -> mat * m, variable_matrices(a)),
        column_vector(a),
        parent(a)
    )
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::GenericFreeFraction{T}, b::GenericFreeFraction{T}; kwargs...) where T
    if parent(a) != parent(b)
        return false
    end
    same_lin = true # check if a and b have the exact same linearization
    for prop in [row_vector, column_vector, constant_matrix, variable_matrices]
        same_lin = same_lin && prop(a) == prop(b)
        if !same_lin
            break
        end
    end
    if same_lin
        return true
    end
    return iszero(a - b; kwargs...)
end

function isone(a::GenericFreeFraction)
    return iszero(a - 1)
end

function iszero(a::GenericFreeFraction; warn::Bool = true, kwargs...)
    if iszero(row_vector(a)) || iszero(column_vector(a))
        return true
    end
    (flag, candidate) = minimal_candidate(a; kwargs...)
    if dim(candidate) == 0
        return true
    end
    if flag && dim(candidate) > 0
        return false
    end
    if warn
        @warn "since minimization is not guaranteed to work optimally for
            non-regular free fractions, some linearizations of zero might
            not be recognized as zero"
    end
    return false
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

function ==(a::GenericFreeFraction{T}, b::MonomialFreeFraction{T}) where T
    return a == GenericFreeFraction(b)
end

function ==(a::MonomialFreeFraction{T}, b::GenericFreeFraction{T}) where T
    return b == a
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function _evaluate(
    a::GenericFreeFraction{T},
    mode::Symbol,
    args::MatElem{T}...,
) where T
    if length(args) != nvars(a)
        error("wrong number of arguments: need $(nvars(a)) matrices,
            but got $(length(args))")
    end
    mat_size = nrows(args[1])
    for mat in args
        if size(mat) != (mat_size, mat_size)
            error("matrix dimensions of arguments not compatible")
        end
    end
    id = identity_matrix(base_ring(a), mat_size)
    res = kronecker_product(constant_matrix(a), id)
    for var_index in 1:nvars(a)
        var_mat = variable_matrices(a)[var_index]
        res += kronecker_product(var_mat, args[var_index])
    end
    if mode == :matrix
        return res
    end
    if rank(res) < dim(a) * mat_size
        e_msg = "free fraction cannot be evaluated at this matrix point"
        if mode == :left || mode == :right
            e_msg = "$(mode) family of " * e_msg
        end
        throw(DomainError(args, e_msg))
    end
    res = inv(res)
    if mode != :left
        new_row = kronecker_product(row_vector(a), id)
        res = new_row * res
    end
    if mode != :right
        new_col = kronecker_product(column_vector(a), id)
        res = res * new_col
    end
    return res
end

function _evaluate(
    a::GenericFreeFraction{T},
    mode::Symbol,
    args::Union{T, Int}...
) where T
    if length(args) != nvars(a)
        error("wrong number of arguments: need $(nvars(a)) scalars,
            but got $(length(args))")
    end
    res = constant_matrix(a)
    for var_index in 1:nvars(a)
        res += args[var_index] * variable_matrices(a)[var_index]
    end
    if mode == :matrix
        return res
    end
    if rank(res) < dim(a)
        e_msg = "free fraction cannot be evaluated at this point"
        if mode == :left || mode == :right
            e_msg = "$(mode) family of " * e_msg
        end
        throw(DomainError(args, e_msg))
    end
    res = inv(res)
    if mode != :left
        res = row_vector(a) * res
    end
    if mode != :right
        res = res * column_vector(a)
    end
    if mode != :left && mode != :right
        return res[1, 1]
    end
    return res
end

###############################################################################
#
#   Probabilistic fullness of representation matrix
#
###############################################################################

function probabilistic_fullness_of_representation_matrix(
    a::GenericFreeFraction;
    max_size::Int = 5
)
    for mat_size in 1:max_size
        test_mats = []
        R = base_ring(a)
        for _ in 1:nvars(a)
            push!(test_mats, matrix(R, FFrand(R, mat_size, mat_size)))
        end
        m = evaluate_matrix(a, test_mats...)
        if rank(m) == dim(a) * mat_size
            return true
        end
    end
    return false
end

###############################################################################
#
#   Linear argument shift
#
###############################################################################

function forward_shift(
        a::GenericFreeFraction{T},
        points::Union{T, Int, Rational{Int}}...
    ) where T
    if length(points) != nvars(a)
        error("expected $(nvars(a)) shift arguments, but received
            $(length(points))")
    end
    d = dim(a)
    b = deepcopy(a)
    R = base_ring(b)
    cm = constant_matrix(b)
    for var_index in 1:nvars(b)
        mat = variable_matrices(b)[var_index]
        for row in 1:d, col in 1:d
            if mat[row, col] != 0
                cm[row, col] += mat[row, col] * R(points[var_index])
            end
        end
    end
    return b
end

###############################################################################
#
#   Deleting rows and columns
#
###############################################################################

function delete_rows_and_cols(
    a::GenericFreeFraction,
    row_list::AbstractVector{Int},
    col_list::AbstractVector{Int}
)
    rows = unique!(sort!(collect(row_list)))
    cols = unique!(sort!(collect(col_list)))
    new_col = delete_rows(column_vector(a), rows)
    new_row = delete_cols(row_vector(a), cols)
    new_const = delete_rows_and_cols(constant_matrix(a), rows, cols)
    new_var = map(mat -> delete_rows_and_cols(mat, rows, cols),
        variable_matrices(a))
    return GenericFreeFraction(new_row, new_const, new_var, new_col, parent(a))
end

function delete_rows_and_cols(
    a::GenericFreeFraction,
    indices::AbstractVector{Int}
)
    return delete_rows_and_cols(a, indices, indices)
end

function delete_row_and_col(a::GenericFreeFraction, row::Int, col::Int)
    return delete_rows_and_cols(a, [row], [col])
end

function delete_row_and_col(a::GenericFreeFraction, index::Int)
    return delete_rows_and_cols(a, [index], [index])
end

function delete_rows(a::GenericFreeFraction, rows::AbstractVector{Int})
    return delete_rows_and_cols(a, rows, Int[])
end

function delete_row(a::GenericFreeFraction, row::Int)
    return delete_rows_and_cols(a, [row], Int[])
end

function delete_cols(a::GenericFreeFraction, cols::AbstractVector{Int})
    return delete_rows_and_cols(a, Int[], cols)
end

function delete_col(a::GenericFreeFraction, col::Int)
    return delete_rows_and_cols(a, Int[], [col])
end

###############################################################################
#
#   Elementary matrix transformations
#
###############################################################################

function add_column!(
    a::GenericFreeFraction,
    scalar::RingElement,
    start_col::Int,
    target_col::Int,
    rows::UnitRange{Int} = 1:dim(a)
)
    s = base_ring(a)(scalar)
    add_column!(row_vector(a), s, start_col, target_col, 1:1)
    add_column!(constant_matrix(a), s, start_col, target_col, rows)
    add_column!.(variable_matrices(a), s, start_col, target_col, Ref(rows))
    return a
end

function add_column(
    a::GenericFreeFraction,
    scalar::RingElement,
    start_col::Int,
    target_col::Int,
    rows::UnitRange{Int} = 1:dim(a)
)
    b = deepcopy(a)
    return add_column!(b, scalar, start_col, target_col, rows)
end

function add_row!(
    a::GenericFreeFraction,
    scalar::RingElement,
    start_row::Int,
    target_row::Int,
    cols::UnitRange{Int} = 1:dim(a)
)
    s = base_ring(a)(scalar)
    add_row!(column_vector(a), s, start_row, target_row, 1:1)
    add_row!(constant_matrix(a), s, start_row, target_row, cols)
    add_row!.(variable_matrices(a), s, start_row, target_row, Ref(cols))
    return a
end

function add_row(
    a::GenericFreeFraction,
    scalar::RingElement,
    start_row::Int,
    target_row::Int,
    cols::UnitRange{Int} = 1:dim(a)
)
    b = deepcopy(a)
    return add_row!(b, scalar, start_row, target_row, cols)
end

function multiply_column!(
    a::GenericFreeFraction,
    scalar::RingElement,
    col::Int,
    rows::UnitRange{Int} = 1:dim(a)
)
    s = base_ring(a)(scalar)
    multiply_column!(row_vector(a), s, col, 1:1)
    multiply_column!(constant_matrix(a), s, col, rows)
    multiply_column!.(variable_matrices(a), s, col, Ref(rows))
    return a
end

function multiply_column(
    a::GenericFreeFraction,
    scalar::RingElement,
    col::Int,
    rows::UnitRange{Int} = 1:dim(a)
)
    b = deepcopy(a)
    return multiply_column!(b, scalar, col, rows)
end

function multiply_row!(
    a::GenericFreeFraction,
    scalar::RingElement,
    row::Int,
    cols::UnitRange{Int} = 1:dim(a)
)
    s = base_ring(a)(scalar)
    multiply_row!(column_vector(a), s, row, 1:1)
    multiply_row!(constant_matrix(a), s, row, cols)
    multiply_row!.(variable_matrices(a), s, row, Ref(cols))
    return a
end

function multiply_row(
    a::GenericFreeFraction,
    scalar::RingElement,
    row::Int,
    cols::UnitRange{Int} = 1:dim(a)
)
    b = deepcopy(a)
    return multiply_row!(b, scalar, row, cols)
end

function swap_cols!(a::GenericFreeFraction, i::Int, j::Int)
    swap_cols!(row_vector(a), i, j)
    swap_cols!(constant_matrix(a), i, j)
    swap_cols!.(variable_matrices(a), i, j)
    return a
end

function swap_cols(a::GenericFreeFraction, i::Int, j::Int)
    b = deepcopy(a)
    return swap_cols!(b, i, j)
end

function swap_rows!(a::GenericFreeFraction, i::Int, j::Int)
    swap_rows!(column_vector(a), i, j)
    swap_rows!(constant_matrix(a), i, j)
    swap_rows!.(variable_matrices(a), i, j)
    return a
end

function swap_rows(a::GenericFreeFraction, i::Int, j::Int)
    b = deepcopy(a)
    return swap_rows!(b, i, j)
end

###############################################################################
#
#   General row/column permutations
#
###############################################################################

function permute_cols!(a::GenericFreeFraction, p::Perm)
    permute_cols!(row_vector(a), p)
    permute_cols!(constant_matrix(a), p)
    permute_cols!.(variable_matrices(a), Ref(p))
    return a
end

function permute_rows!(a::GenericFreeFraction, p::Perm)
    permute_rows!(column_vector(a), p)
    permute_rows!(constant_matrix(a), p)
    permute_rows!.(variable_matrices(a), Ref(p))
    return(a)
end
