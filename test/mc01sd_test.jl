# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01sd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 10
    P = Array{Float64,1}(undef, DPMAX+1)
    MANT = Array{Float64,1}(undef, DPMAX+1)
    E = Array{BlasInt,1}(undef, DPMAX+1)
    IWORK = Array{BlasInt,1}(undef, DPMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP = parse(BlasInt, vs[1])
if ( DP<=-1 || DP>DPMAX )
else

    vs = String[]
    _isz = DP+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P[1:_isz] .= parsex.(Float64, vs)
    
# interp call 1

    S, T, INFO = SLICOT.mc01sd!(DP, P, MANT, E)
    println(io, "S = $S")
    println(io, "T = $T")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
    BETA = 2.0
end # if
end # if
    close(f)
end # run_X()
