# FreeFractions.jl

## About

FreeFractions is a Julia package to be used with the [Nemo](https://github.com/Nemocas/Nemo.jl) package.
It provides an implementation of free skew fields according to the [Noncommutative ring Interface](https://nemocas.github.io/AbstractAlgebra.jl/latest/ncrings/) of Nemo's backend [AbstractAlgebra](https://github.com/Nemocas/AbstractAlgebra.jl).

This package was originally developed within the [Collaborative Research Center TRR 195](https://www.computeralgebra.de/sfb/), which is funded by the Deutsche Forschungsgemeinschaft (DFG).

Work is still ongoing, so everything might change at some point.
In particular, this documentation is a work in progress.

## Installation and setup

In the package manager (accessible in the REPL via `]`):

```julia
(@v1.9) pkg> add https://github.com/johannes-hoffmann/FreeFractions.jl
```

Now leave the package manager (via backspace) and load the package (to effectively use it, you will also need to load Nemo):

```julia
julia> using FreeFractions; using Nemo
```

Note: FreeFractions includes some functionality that is only available if the package [Singular](https://github.com/oscar-system/Singular.jl) is present.
If Singular is not present, a warning to this effect will be shown.

## Quick start guide if you only want a reasonably small representation

To get a representation of $(x\cdot y+y\cdot x)^{-1}$ (as a non-commutative rational function in variables $x$ and $y$ with rational coeffcients), you can do the following:

```julia
julia> F, (x, y) = FreeField(QQ, ["x", "y"]);
julia> f = inv(x*y+y*x)
Free fraction in x, y over Rational field represented in dimension 3 by
          [ -y 0 1 ]⁻¹ [ 0 ]
          [ -x 1 0 ]   [ 0 ]
[ 1 0 0 ] [  0 y x ]   [ 1 ]
```

The representation follows the convention $f=u\cdot A^{-1}\cdot v$, to get the extended matrix pencil of 
$
\begin{pmatrix}
0&u\\
v&A
\end{pmatrix}
$
you can use either `extended_matrix_pencil(f)` for a vector of matrices or `extended_matrix_pencil_tensor(f)` for a three-dimensional array.

You can use any Nemo `Field` instead of `QQ` as long as its arithmetic is exact.
For (reasonable) real or complex numbers, you can use `CalciumField`:
```julia
julia> C = CalciumField();
julia> F, (x, y) = FreeField(C, ["x", "y"]);
julia> f = sqrt(C(2)) * C(1im) * (x*y-y*x)
Free fraction in x, y over Exact complex field represented in dimension 4 by
            [ 1 -1*x -1*y    0 ]⁻¹ [                                                              0 ]
            [ 0    1    0    y ]   [                                                              0 ]
            [ 0    0    1 -1*x ]   [                                                              0 ]
[ 1 0 0 0 ] [ 0    0    0    1 ]   [ -1.41421*I {-a*b where a = 1.41421 [a^2-2=0], b = I [b^2+1=0]} ]
```
Here you get a representation of $\mathrm{i}\sqrt{2}(xy-yx)$. You can use `ComplexF64.(extended_matrix_pencil_tensor(f))` to get a tensor in Julia's float arithmetic if you are interested in doing further numeric computations .

## Detailed usage

### The parent object: `FreeField`

Create a `FreeField`, the parent object of a `FreeFraction`, over a Nemo `Field` (e.g. the rational numbers `QQ`) in some non-commuting variables (e.g. `x` and `y`):

```julia
julia> F, (x, y) = FreeField(QQ, ["x", "y"]);

julia> F
Free Field in x, y over Rational Field

julia> typeof(F)
FreeField{QQFieldElem}

julia> base_ring(F)
Rational Field

julia> nvars(F)
2

julia> gen(F, 1)
x

julia> gens(F)
2-element Array{MonomialFreeFraction{QQFieldElem},1}:
 x
 y

julia> characteristic(F)
0

julia> F(42)
42
```

### The element object: `FreeFraction`

A rational expression of a free fraction belonging to a `FreeField` can be constructed easily:

```julia
julia> F, (x, y) = FreeField(QQ, ["x", "y"]);

julia> f = x * y * 7 * x
7*x*y*x

julia> g = inv(x * y + y * x)
Free fraction in x, y over Rational Field represented in dimension 3 by
          [ y  0  1 ]⁻¹ [ 0 ]
          [ x  1  0 ]   [ 0 ]
[ 1 0 0 ] [ 0 -y -x ]   [ 1 ]
```

Why do `f` and `g` look different?
The abstract type `FreeFraction` has two concrete subtypes: `MonomialFreeFraction` and `GenericFreeFraction`.
Monomials (including scalars) are represented by the former, everything else by the latter.
Representations of monomials can be converted back and forth.
In most cases, the user will not have to think about this distinction.

```julia
julia> typeof(f), typeof(g)
(MonomialFreeFraction{fmpq}, GenericFreeFraction{fmpq})

julia> h = GenericFreeFraction(f)
Free fraction in x, y over Rational Field represented in dimension 4 by
            [ 1 -x  0  0 ]⁻¹ [ 0 ]
            [ 0  1 -y  0 ]   [ 0 ]
            [ 0  0  1 -x ]   [ 0 ]
[ 1 0 0 0 ] [ 0  0  0  1 ]   [ 7 ]

julia> MonomialFreeFraction(h)
7*x*y*x
```

Alternatively, a `FreeFraction` can be constructed without explicitly creating the parent `FreeField` (which will still be constructed in the background):

```julia
julia> m = FreeFraction("inv(x * y + y * x)", QQ, ["x", "y"])
Free fraction in x, y over Rational Field represented in dimension 3 by
          [ y  0  1 ]⁻¹ [ 0 ]
          [ x  1  0 ]   [ 0 ]
[ 1 0 0 ] [ 0 -y -x ]   [ 1 ]

julia> parent(m) == F
true
```

Further information can be gained:

```julia
julia> dim(f)
3

julia> isregular(f)
false

julia> iszero(f)
false

julia> isunit(f)
true
```

If the creation `FreeFraction`s seems slow this is due to the fact that after performing an addition, multiplication, or inversion of `FreeFraction`s, the resulting linear representation will be automatically minimized<sup>[1](#min_note)</sup>.
Minimality can be checked automatically, and the minimization process can be called manually as well:

```julia
julia> isminimal(g)
true

julia> minimal(g)
Free fraction in x, y over Rational Field represented in dimension 3 by
          [ y  0  1 ]⁻¹ [ 0 ]
          [ x  1  0 ]   [ 0 ]
[ 1 0 0 ] [ 0 -y -x ]   [ 1 ]
```

You can try to evaluate any `FreeFraction` at tuples of scalars or square matrices of the same sizes.
If this is not possible for the current representation, an error will be thrown.

```julia
julia> f = inv(x * y + y * x)
Free fraction in x, y over Rational Field represented in dimension 3 by
          [ y  0  1 ]⁻¹ [ 0 ]
          [ x  1  0 ]   [ 0 ]
[ 1 0 0 ] [ 0 -y -x ]   [ 1 ]

julia> f(1, 2)
1//4

julia> m1 = matrix(QQ, [1 2; 3 4])
[1  2]
[3  4]

julia> m2 = matrix(QQ, [5 6; 7 8])
[5  6]
[7  8]

julia> f(m1, m2)
[ -6//7   1//2]
[37//56  -3//8]

julia> f(0, 0)
ERROR: DomainError with (0, 0):
free fraction cannot be evaluated at this point
```

---

<a name="min_note">1</a>: To suppress automatic intermediate minimization, the `FreeField` constructor accepts an optional Boolean parameter:

```julia
julia> F, (x, y) = FreeField(QQ, ["x", "y"], false)
(Free Field in x, y over Rational Field, MonomialFreeFraction{fmpq}[x, y])

julia> f = inv(x * y + y * x)
Free fraction in x, y over Rational Field represented in dimension 7 by
                  [  0 1 -x  0 -1  0  0 ]⁻¹ [ 0 ]
                  [  0 0  1 -y  0  0  0 ]   [ 0 ]
                  [ -1 0  0  1  0  0  0 ]   [ 0 ]
                  [  0 0  0  0  1 -y  0 ]   [ 0 ]
                  [  0 0  0  0  0  1 -x ]   [ 0 ]
                  [ -1 0  0  0  0  0  1 ]   [ 0 ]
[ 1 0 0 0 0 0 0 ] [  0 1  0  0  0  0  0 ]   [ 1 ]

julia> FreeFractions.intermediate_minimization(F)
false
```

Note that two `FreeField`s that only differ in this `intermediate_minimization` parameter will be considered as equal by `==` since they do not differ mathematically.
In the process of adding or multiplying two `FreeFraction`s, intermediate minimization steps will only be carried out if both parent objects indicate this, i.e. their `intermediate_minimization` parameter is set to `true`.
