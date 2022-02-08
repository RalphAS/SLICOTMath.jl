using SLICOTMath
const SLICOT=SLICOTMath
using Test
using LinearAlgebra
using LinearAlgebra: BlasInt
const BlasBool = BlasInt

LSAME(x,y) = uppercase(x) == uppercase(y)

# helpers for getting data from the examples
parsex(T,s) = parse(T,s)
# CAUTION: this relies on absence of spaces in complex constants
function parsex(::Type{ComplexF64}, s)
  if s == "0"
    return zero(ComplexF64)
  end
  s1 = split(s,",")
  if length(s1) == 1
    return parse(Float64,s1[1])+0.0im
  end
  if s1[2][1] in ['+','-']
    impart = s1[2][1:end-1] * "im"
  else
    impart = '+' * s1[2][1:end-1] * "im"
  end
  s2 = s1[1][2:end] * impart
  parse(ComplexF64, s2)
end
function parsex(::Type{Bool}, s)
    return s[2] in ('T','t')
end

const quiet = Ref(true)
const tdir = joinpath(@__DIR__,"SLICOT-Reference")

if !isdir(joinpath(tdir,"examples"))
    @error """Package tests require that 'test/SLICOT-Reference' hold the input data.
This can be prepared by installing a copy of the SLICOT distribution there,
or making a link to one.
"""
    # force an early exit
    @testset "installation" begin
        @test false
    end
end

files = readdir(@__DIR__)
funclist = [replace(s, "_test.jl" => "") for s in files if (length(s)>7) && (s[end-7:end] == "_test.jl")]

for func2test in funclist
    testcode = lowercase(func2test) * "_test.jl"
    include(testcode)

    rname=Symbol("run_"*func2test)
    eval(:(runme=$rname))

    @testset "$func2test" begin
        datfile = joinpath(tdir, "examples", uppercase(func2test) * ".dat")

        if quiet[]
            runme(datfile, devnull)
        else
            runme(datfile)
        end
    end # testset
end # for

