# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc01wd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 10
    P = Array{Float64,1}(undef, DPMAX+1)
    Q = Array{Float64,1}(undef, DPMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DP = parse(BlasInt, vs[1])
    if ( DP<=-1 || DP>DPMAX )
        @error "Illegal DP=$DP"
    end

    vs = String[]
    _isz = DP+1
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    P[1:_isz] .= parsex.(Float64, vs)

    vs = split(readline(f))
    U1 = parse(Float64, replace(vs[1],'D'=>'E'))
    U2 = parse(Float64, replace(vs[2],'D'=>'E'))
    close(f)
    # interp call 1

    INFO = SLICOT.mc01wd!(DP, P, U1, U2, Q)
    @test INFO == 0
    INFO == 0 || return

    println(io, "Q (power,coefft):")
    show(io, "text/plain", hcat(0:DP-2,Q[3:DP+1]))
    println(io)
    println(io, "R (power, coefft):")
    show(io, "text/plain", hcat(0:1,[Q[1]+Q[2]*U2,Q[2]]))
    println(io)
end # run_mc01wd()
