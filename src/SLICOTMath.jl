module SLICOTMath
using SLICOT_jll
using LinearAlgebra
using LinearAlgebra: BlasInt
using DocStringExtensions

function chkargsok(ret::BlasInt)
    if ret < 0
        throw(ArgumentError("invalid argument #$(-ret) to SLICOT call"))
    end
end
include("slicot_m.jl")
include("m_docs.jl")


end
