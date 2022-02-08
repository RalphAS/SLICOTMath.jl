# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04xd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LDA = MMAX
    LDU = MMAX
    LDV = NMAX
    MAXMN = max( MMAX, NMAX )
    MNMIN = min( MMAX, NMAX )
    LENGQ = 2*MNMIN-1
    LDWORK = max( 2*NMAX, NMAX*( NMAX+1 )÷2 ) + max( 2*MNMIN + MAXMN, 8*MNMIN - 5 )
    A = Array{Float64,2}(undef, LDA,NMAX)
    U = Array{Float64,2}(undef, LDU,MMAX)
    V = Array{Float64,2}(undef, LDV,NMAX)
    Q = Array{Float64,1}(undef, LENGQ)
    INUL = Array{BlasBool,1}(undef, MAXMN)
    DWORK = Array{Float64,1}(undef, LDWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    RANK = parse(BlasInt, vs[3])
    THETA = parse(Float64, replace(vs[4],'D'=>'E'))
    TOL = parse(Float64, replace(vs[5],'D'=>'E'))
    RELTOL = parse(Float64, replace(vs[6],'D'=>'E'))
    JOBU =  vs[7][1]
    JOBV =  vs[8][1]
if ( M<0 || M>MMAX )
elseif ( N<0 || N>NMAX )
elseif ( RANK>MNMIN )
elseif ( RANK<0 && THETA<ZERO )
else

    vs = String[]
    _isz,_jsz = (M,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    RANK1 = RANK
    THETA1 = THETA
# interp call 1
    INUL .= 0
    RANK, THETA, INFO, IWARN = SLICOT.mb04xd!(JOBU, JOBV, M, N, RANK, THETA, A, U, V, Q, INUL, TOL, RELTOL)
    println(io, "RANK = $RANK")
    println(io, "THETA = $THETA")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")
    inul = convert(Vector{Bool},INUL)

if ( INFO!=0 )
else
if ( IWARN!=0 )
else
end # if
    LJOBUA = LSAME( JOBU, 'A' )
    LJOBUS = LSAME( JOBU, 'S' )
    LJOBVA = LSAME( JOBV, 'A' )
    LJOBVS = LSAME( JOBV, 'S' )
    WANTU = LJOBUA||LJOBUS
    WANTV = LJOBVA||LJOBVS
    MINMN = min( M, N )
    LOOP = MINMN - 1
    println(io, "Q:")
    show(io, "text/plain", Bidiagonal(Q[1:MINMN],Q[MINMN+1:MINMN+LOOP],'U'))
    println(io)
if ( WANTU )
    NCOLU = LJOBUS ? MINMN : M
# interp output 1
    println(io, "U:")
    _nc = NCOLU
    _nr = M
    show(io, "text/plain", U[1:_nr,1:_nc])
    println(io)
    println(io, "INUL:")
    show(io, "text/plain", inul[1:NCOLU])
    println(io)

end # if
if ( WANTV )
    NCOLV = LJOBVS ? MINMN : N
# interp output 2
    println(io, "V:")
    _nc = NCOLV
    _nr = N
    show(io, "text/plain", V[1:_nr,1:_nc])
    println(io)
    println(io, "INUL:")
    show(io, "text/plain", inul[1:NCOLV])
    println(io)

end # if
end # if
end # if
    close(f)
end # run_mb04xd()
