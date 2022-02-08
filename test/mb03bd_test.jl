# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03bd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 6
    NMAX = 50
    LDA1 = NMAX
    LDA2 = NMAX
    LDQ1 = NMAX
    LDQ2 = NMAX
    LDWORK = KMAX + max( 2*NMAX, 8*KMAX )
    LIWORK = 2*KMAX + NMAX
    QIND = Array{BlasInt,1}(undef, KMAX)
    S = Array{BlasInt,1}(undef, KMAX)
    A = Array{Float64,3}(undef, LDA1,LDA2,KMAX)
    Q = Array{Float64,3}(undef, LDQ1,LDQ2,KMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX)
    ALPHAI = Array{Float64,1}(undef, NMAX)
    BETA = Array{Float64,1}(undef, NMAX)
    SCAL = Array{BlasInt,1}(undef, NMAX)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    DEFL =  vs[2][1]
    COMPQ =  vs[3][1]
    K = parse(BlasInt, vs[4])
    N = parse(BlasInt, vs[5])
    H = parse(BlasInt, vs[6])
    ILO = parse(BlasInt, vs[7])
    IHI = parse(BlasInt, vs[8])
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end

    vs = String[]
    _isz = K
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    S[1:_isz] .= parsex.(BlasInt, vs)


    vs = String[]
    _isz,_jsz,_ksz = (N,N,K)
    while length(vs) < _isz*_jsz*_ksz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for k in 1:_ksz
      for i in 1:_isz
        _i0 = (i-1)*_jsz + (k-1)*_jsz*_isz
        A[i,1:_jsz,k] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
      end
    end
    if ( LSAME( COMPQ, 'P' ) )

    vs = String[]
    _isz = K 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    QIND[1:_isz] .= parsex.(BlasInt, vs)

    end # if
    close(f)
# interp call 1

    INFO, IWARN = SLICOT.mb03bd!(JOB, DEFL, COMPQ, QIND, K, N, H, ILO, IHI, S, A, Q, ALPHAR, ALPHAI, BETA, SCAL, LIWORK, LDWORK)
    @test INFO == 0
    println(io, "IWARN = $IWARN")

    if ( LSAME( JOB, 'S' ) || LSAME( JOB, 'T' ) )
# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    _nk = K
    show(io, "text/plain", A[1:_nr,1:_nc,1:_nk])
    println(io)

    end # if
    if ( LSAME( COMPQ, 'U' ) || LSAME( COMPQ, 'I' ) )
# interp output 2
    println(io, "Q:")
    _nc = N
    _nr = N
    _nk = K
    show(io, "text/plain", Q[1:_nr,1:_nc,1:_nk])
    println(io)

    elseif ( LSAME( COMPQ, 'P' ) )
        for L in 1:K
            if ( QIND[L]>0 )
                println(io, "factor ",QIND[L])
#   write QIND( L )
# unable to translate write loop:
#  write ( I, J, QIND( L ) ), J = 1, N 
# interp output 3
                show(io, "text/plain", Q[1:N,1:N,QIND[L]])
                println(io)

            end # if
        end # for
    end # if
# interp output 4
    println(io, "ALPHAR:")
    _nr = N
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

# interp output 5
    println(io, "ALPHAI:")
    _nr = N
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

# interp output 6
    println(io, "BETA:")
    _nr = N
    show(io, "text/plain", BETA[1:_nr])
    println(io)

# interp output 7
    println(io, "SCAL:")
    _nr = N
    show(io, "text/plain", SCAL[1:_nr])
    println(io)
end # run_mb03bd()
