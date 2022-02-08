# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb05nd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDEX = NMAX
    LDEXIN = NMAX
    LDWORK = NMAX*( NMAX+1 )
    A = Array{Float64,2}(undef, LDA,NMAX)
    EX = Array{Float64,2}(undef, LDEX,NMAX)
    EXINT = Array{Float64,2}(undef, LDEXIN,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    DELTA = parse(Float64, replace(vs[2],'D'=>'E'))
    TOL = parse(Float64, replace(vs[3],'D'=>'E'))
    if ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
    # interp call 1

    INFO = SLICOT.mb05nd!(N, DELTA, A, EX, EXINT, TOL, LDWORK)
    @test INFO == 0
    INFO == 0 || return

# interp output 1
    println(io, "EX:")
    _nc = N
    _nr = N
    show(io, "text/plain", EX[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "EXINT:")
    _nc = N
    _nr = N
    show(io, "text/plain", EXINT[1:_nr,1:_nc])
    println(io)

end # run_mb05nd()
