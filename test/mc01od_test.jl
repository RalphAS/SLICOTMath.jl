# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc01od(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 10
    REZ = Array{Float64,1}(undef, KMAX)
    IMZ = Array{Float64,1}(undef, KMAX)
    REP = Array{Float64,1}(undef, KMAX+1)
    IMP = Array{Float64,1}(undef, KMAX+1)
    DWORK = Array{Float64,1}(undef, 2*KMAX+2)
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

    INFO = SLICOT.mc01od!(K, REZ, IMZ, REP, IMP)
    @test INFO == 0
    INFO == 0 || return

# unable to translate write loop:
#  write , REP(I+1), IMP(I+1), I = 0,K 
# interp output 1
    println(io, "P:")
    show(io, "text/plain", REP[1:K+1]+im*IMP[1:K+1])
    println(io)

end # run_mc01od()
