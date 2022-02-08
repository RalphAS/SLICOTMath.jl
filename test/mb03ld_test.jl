# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03ld(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDA = NMAX÷2
    LDB = NMAX÷2
    LDDE = NMAX÷2
    LDFG = NMAX÷2
    LDQ = 2*NMAX
    LDWORK = 8*NMAX*NMAX + max( 8*NMAX + 32, NMAX÷2 + 168, 272 )
    LIWORK = max( 32, NMAX + 12, NMAX*2 + 3 )
    A = Array{Float64,2}(undef, LDA,NMAX÷2)
    DE = Array{Float64,2}(undef, LDDE,NMAX÷2+1)
    B = Array{Float64,2}(undef, LDB,NMAX÷2)
    FG = Array{Float64,2}(undef, LDFG,NMAX÷2+1)
    Q = Array{Float64,2}(undef, LDQ,2*NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX÷2)
    ALPHAI = Array{Float64,1}(undef, NMAX÷2)
    BETA = Array{Float64,1}(undef, NMAX÷2)
    BWORK = Array{BlasBool,1}(undef, NMAX÷2)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    COMPQ =  vs[1][1]
    ORTH =  vs[2][1]
    N = parse(BlasInt, vs[3])
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end

    M = N÷2

    vs = String[]
    _isz,_jsz = (M,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       DE[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       FG[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    NEIG, INFO = SLICOT.mb03ld!(COMPQ, ORTH, N, A, DE, B, FG, Q, ALPHAR, ALPHAI, BETA, LIWORK)
    @test INFO == 0
    INFO == 0 || return
    println(io, "NEIG = $NEIG")

# interp output 1
    println(io, "A:")
    _nc = M
    _nr = M
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "DE:")
    _nc = M+1
    _nr = M
    show(io, "text/plain", DE[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "B:")
    _nc = M
    _nr = M
    show(io, "text/plain", B[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "FG:")
    show(io, "text/plain", FG[1:M,2:M+1])
    println(io)

# interp output 5
    println(io, "ALPHAR:")
    _nr = M
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

# interp output 6
    println(io, "ALPHAI:")
    _nr = M
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

# interp output 7
    println(io, "BETA:")
    _nr = M
    show(io, "text/plain", BETA[1:_nr])
    println(io)

if ( LSAME( COMPQ, 'C' ) && NEIG>0 )
# interp output 8
    println(io, "Q:")
    _nc = NEIG
    _nr = N
    show(io, "text/plain", Q[1:_nr,1:_nc])
    println(io)

end # if

end # run_mb03ld()
