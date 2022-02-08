# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02nd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LMAX = 20
    LDC = max( MMAX, NMAX+LMAX )
    LDX = NMAX
    LENGQ = 2*min(MMAX,NMAX+LMAX)-1
    LIWORK = NMAX+2*LMAX
    LDWORK = max(2, max( MMAX, NMAX+LMAX ) + 2*min( MMAX, NMAX+LMAX ), min( MMAX, NMAX+LMAX ) + max( ( NMAX+LMAX )*( NMAX+LMAX-1 )÷2, MMAX*( NMAX+LMAX-( MMAX-1 )÷2 ) ) + max( 6*(NMAX+LMAX)-5, LMAX*LMAX + max( NMAX+LMAX, 3*LMAX ) ) )
    LBWORK = NMAX+LMAX
    C = Array{Float64,2}(undef, LDC,NMAX+LMAX)
    X = Array{Float64,2}(undef, LDX,LMAX)
    Q = Array{Float64,1}(undef, LENGQ)
    INUL = Array{BlasBool,1}(undef, NMAX+LMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    BWORK = Array{BlasBool,1}(undef, LBWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    L = parse(BlasInt, vs[3])
    RANK = parse(BlasInt, vs[4])
    THETA = parse(Float64, replace(vs[5],'D'=>'E'))
    TOL = parse(Float64, replace(vs[6],'D'=>'E'))
    RELTOL = parse(Float64, replace(vs[7],'D'=>'E'))
if ( M<0 || M>MMAX )
elseif ( N<0 || N>NMAX )
elseif ( L<0 || L>LMAX )
elseif ( RANK>min( MMAX, NMAX ) )
elseif ( RANK<0 && THETA<ZERO )
else

    vs = String[]
    _isz,_jsz = (M,N+L)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       C[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    RANK1 = RANK
    THETA1 = THETA
# interp call 1

    RANK, THETA, INFO, IWARN = SLICOT.mb02nd!(M, N, L, RANK, THETA, C, X, Q, INUL, TOL, RELTOL)
    println(io, "RANK = $RANK")
    println(io, "THETA = $THETA")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

if ( INFO!=0 )
else
if ( IWARN!=0 )
else
end # if
    MINMNL = min( M, N+L )
    LOOP = MINMNL - 1
    println(io, "Q:")
    show(io, "text/plain", Bidiagonal(Q[1:MINMNL],Q[MINMNL+1:MINMNL+LOOP],'U'))
    println(io)
    println(io, "X:")
    show(io, "text/plain", X[1:N,1:L])
    println(io)

# interp output 1
    println(io, "C:")
    _nc = N+L
    _nr = max( M, N + L )
    show(io, "text/plain", C[1:_nr,1:_nc])
    println(io)

    inul = convert(Vector{Bool}, INUL[1:N+L])
    println(io, "INUL:")
    show(io, "text/plain", inul[1:N+L])
    println(io)
end # if
end # if
    close(f)
end # run_mb02nd()
