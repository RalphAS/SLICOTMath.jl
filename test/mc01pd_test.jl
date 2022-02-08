# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc01pd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 10
    REZ = Array{Float64,1}(undef, KMAX)
    IMZ = Array{Float64,1}(undef, KMAX)
    P = Array{Float64,1}(undef, KMAX+1)
    DWORK = Array{Float64,1}(undef, KMAX+1)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    if ( K<0 || K>KMAX )
        @error "Illegal K=$K"
    end

    vs = String[]
    _isz,_jsz = (K,2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
        REZ[i],IMZ[i] = parsex.(Float64, vs[2i-1:2i])
    end
    close(f)
# interp call 1

    INFO = SLICOT.mc01pd!(K, REZ, IMZ, P)
    @test INFO == 0
    INFO == 0 || return


    println(io, "P:")
    show(io, "text/plain", P[1:K+1])
    println(io)

end # run_mc01pd()
