# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02md(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LMAX = 20
    LDC = max( MMAX,NMAX+LMAX )
    LDX = NMAX
    LDWORK = MMAX*(NMAX+LMAX) + max( 3*min(MMAX,NMAX+LMAX) + max(MMAX,NMAX+LMAX), 5*min(MMAX,NMAX+LMAX), 3*LMAX )
    LIWORK = LMAX
    LENGS = min( MMAX, NMAX+LMAX )
    C = Array{Float64,2}(undef, LDC,NMAX+LMAX)
    S = Array{Float64,1}(undef, LENGS)
    X = Array{Float64,2}(undef, LDX,LMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    L = parse(BlasInt, vs[3])
    JOB =  vs[4][1]
if ( LSAME( JOB, 'R' ) )
    vs = split(readline(f))
    TOL = parse(Float64, replace(vs[1],'D'=>'E'))
elseif ( LSAME( JOB, 'T' ) )
    vs = split(readline(f))
    RANK = parse(BlasInt, vs[1])
    SDEV = parse(Float64, replace(vs[2],'D'=>'E'))
    TOL = SDEV
elseif ( LSAME( JOB, 'N' ) )
    vs = split(readline(f))
    RANK = parse(BlasInt, vs[1])
    TOL = parse(Float64, replace(vs[2],'D'=>'E'))
else
    RANK = 0
    vs = split(readline(f))
    SDEV = parse(Float64, replace(vs[1],'D'=>'E'))
    TOL = SDEV
end # if
if ( M<0 || M>MMAX )
elseif ( N<0 || N>NMAX )
elseif ( L<0 || L>LMAX )
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
# interp call 1

    RANK, INFO, IWARN = SLICOT.mb02md!(JOB, M, N, L, RANK, C, S, X, TOL)
    println(io, "RANK = $RANK")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

if ( INFO!=0 )
else
if ( IWARN!=0 )
else
end # if
# unable to translate write statement:
#   in do block [('40', 'J', 'L'), ('20', 'I', 'N')]
    #   write X(I,J)
    println(io, "X:")
    _nr = N
    _nc = L
    show(io, "text/plain", X[1:_nr,1:_nc])
    println(io)
# interp output 1
    println(io, "S:")
    _nr = min( M, N+L )
    show(io, "text/plain", S[1:_nr])
    println(io)

end # if
end # if
    close(f)
end # run_mb02md()
