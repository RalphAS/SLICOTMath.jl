# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc01td(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 10
    P = Array{Float64,1}(undef, DPMAX+1)
    DWORK = Array{Float64,1}(undef, 2*DPMAX+2)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP = parse(BlasInt, vs[1])
    DICO =  vs[2][1]
    if ( DP<=-1 || DP>DPMAX )
        @error "Illegal DP=$DP"
    end

    DPP = DP

    vs = String[]
    _isz = DP+1
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P[1:_isz] .= parsex.(Float64, vs)
    close(f)

# interp call 1

    DP, STABLE, NZ, INFO, IWARN = SLICOT.mc01td!(DICO, DP, P)
    @test INFO == 0
    INFO == 0 || return
    println(io, "DP = $DP")
    println(io, "STABLE = $STABLE")
    println(io, "NZ = $NZ")

    if ( IWARN!=0 )
        @warn "mc01td returned IWARN=$IWARN"
        return
    end # if
    stable = Bool(STABLE)
    if ( stable )
        println(io, "P is stable")
    else
        println(io, "P is unstable")
        println(io, "unstable zeros: $NZ")
    end
end # run_mc01td()
