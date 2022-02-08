# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01md(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 20
    P = Array{Float64,1}(undef, DPMAX+1)
    Q = Array{Float64,1}(undef, DPMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP = parse(BlasInt, vs[1])
    ALPHA = parse(Float64, replace(vs[2],'D'=>'E'))
    K = parse(BlasInt, vs[3])
if ( DP<=-1 || DP>DPMAX )
else

    vs = String[]
    _isz = DP+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P[1:_isz] .= parsex.(Float64, vs)
    
# interp call 1

    INFO = SLICOT.mc01md!(DP, ALPHA, K, P, Q)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
end # if
end # if
    close(f)
end # run_X()
