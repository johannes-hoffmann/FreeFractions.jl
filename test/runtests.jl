using FreeFractions
using Nemo
import AbstractAlgebra
using Test

test_fields = [QQ, GF(2), GF(3), GF(5)]

include("FreeFractions-test.jl")
