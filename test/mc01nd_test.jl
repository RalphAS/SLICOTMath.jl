# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01nd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 20
    P = Array{Float64,1}(undef, DPMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP = parse(BlasInt, vs[1])
    XR = parse(Float64, replace(vs[2],'D'=>'E'))
    XI = parse(Float64, replace(vs[3],'D'=>'E'))
if ( DP<=-1 || DP>DPMAX )
else

    vs = String[]
    _isz = DP+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P[1:_isz] .= parsex.(Float64, vs)
    
# interp call 1

    VR, VI, INFO = SLICOT.mc01nd!(DP, XR, XI, P)
    println(io, "VR = $VR")
    println(io, "VI = $VI")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
end # if
end # if
    close(f)
end # run_X()
