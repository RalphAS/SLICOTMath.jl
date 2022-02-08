# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03rd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 10
    LDA = NMAX
    LDX = NMAX
    LDWORK = 3*NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    X = Array{Float64,2}(undef, LDX,NMAX)
    BLSIZE = Array{BlasInt,1}(undef, NMAX)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    BWORK = Array{Bool,1}(undef, NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    PMAX = parse(Float64, replace(vs[2],'D'=>'E'))
    TOL = parse(Float64, replace(vs[3],'D'=>'E'))
    JOBX =  vs[4][1]
    SORT =  vs[5][1]
    if ( N<0 || N>NMAX )
        @error "Invalid N=$N"
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
# interp call 1

    #FOREIGN.dgees!( 'Vectors', 'Not sorted', SELECT, N, A, LDA, SDIM, WR, WI, X, LDX, DWORK, LDWORK, BWORK, INFO )
    Aschur = schur(A[1:N,1:N])
    A[1:N,1:N] .= Aschur.T
    X[1:N,1:N] .= Aschur.Z

# interp call 2

    NBLCKS, INFO = SLICOT.mb03rd!(JOBX, SORT, N, PMAX, A, X, BLSIZE, WR, WI, TOL)
    println(io, "NBLCKS = $NBLCKS")
    @test INFO == 0
    INFO == 0 || return

    if ( INFO!=0 )
        @warn "mb03rd returns INFO=$INFO"
    end
# interp output 1
    println(io, "BLSIZE:")
    _nr = NBLCKS
    show(io, "text/plain", BLSIZE[1:_nr])
    println(io)

# interp output 2
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "X:")
    _nc = N
    _nr = N
    show(io, "text/plain", X[1:_nr,1:_nc])
    println(io)

    close(f)
end # run_mb03rd()
