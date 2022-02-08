# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01rd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DP1MAX = 10
    DP2MAX = 10
    DP3MAX = 10
    LENP3 = max(DP1MAX+DP2MAX,DP3MAX)+1
    P1 = Array{Float64,1}(undef, DP1MAX+1)
    P2 = Array{Float64,1}(undef, DP2MAX+1)
    P3 = Array{Float64,1}(undef, DP1MAX+DP2MAX+DP3MAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP1 = parse(BlasInt, vs[1])
if ( DP1<=-2 || DP1>DP1MAX )
else

    vs = String[]
    _isz = DP1+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P1[1:_isz] .= parsex.(Float64, vs)
    
    vs = split(readline(f))
    DP2 = parse(BlasInt, vs[1])
if ( DP2<=-2 || DP2>DP2MAX )
else

    vs = String[]
    _isz = DP2+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P2[1:_isz] .= parsex.(Float64, vs)
    
    vs = split(readline(f))
    DP3 = parse(BlasInt, vs[1])
if ( DP3<=-2 || DP3>DP3MAX )
else

    vs = String[]
    _isz = DP3+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P3[1:_isz] .= parsex.(Float64, vs)
    
end # if
    vs = split(readline(f))
    ALPHA = parse(Float64, replace(vs[1],'D'=>'E'))
# interp call 1

    DP3, INFO = SLICOT.mc01rd!(DP1, DP2, DP3, ALPHA, P1, P2, P3)
    println(io, "DP3 = $DP3")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( DP3>=0 )
end # if
end # if
end # if
end # if
    close(f)
end # run_X()
