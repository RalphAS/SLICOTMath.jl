# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03od(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    NMAX = 10
    MMAX = 10
    LDA = NMAX
    LDTAU = min(MMAX,NMAX)
    LDWORK = 3*NMAX + 1
    A = Array{Float64,2}(undef, LDA,NMAX)
    JPVT = Array{BlasInt,1}(undef, NMAX)
    TAU = Array{Float64,1}(undef, LDTAU)
    SVAL = Array{Float64,1}(undef, 3)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    JOBQR =  vs[3][1]
    RCOND = parse(Float64, replace(vs[4],'D'=>'E'))
    SVLMAX = parse(Float64, replace(vs[5],'D'=>'E'))
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end
    if ( M<0 || M>MMAX )
        @error "Illegal M=$M"
    end

    vs = String[]
    _isz,_jsz = (M,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
    # interp call 1

    JPVT[1:N] .= 0
    RANK, INFO = SLICOT.mb03od!(JOBQR, M, N, A, JPVT, RCOND, SVLMAX, TAU, SVAL)
    @test INFO == 0
    println(io, "RANK = $RANK")

# interp output 1
    println(io, "JPVT:")
    _nr = N
    show(io, "text/plain", JPVT[1:_nr])
    println(io)

# interp output 2
    println(io, "SVAL:")
    _nr = 3
    show(io, "text/plain", SVAL[1:_nr])
    println(io)

end # run_mb03od()
