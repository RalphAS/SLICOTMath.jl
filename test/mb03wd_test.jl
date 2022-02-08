# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03wd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    PMAX = 20
    LDA1 = NMAX
    LDA2 = NMAX
    LDTAU = NMAX-1
    LDZ1 = NMAX
    LDZ2 = NMAX
    LDZTA = NMAX
    LDWORK = max( NMAX, NMAX + PMAX - 2 )
    ZERO = 0.0e0
    ONE = 1.0e0
    A = Array{Float64,3}(undef, LDA1,LDA2,PMAX)
    TAU = Array{Float64,2}(undef, LDTAU,PMAX)
    Z = Array{Float64,3}(undef, LDZ1,LDZ2,PMAX)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
    AS = Array{Float64,3}(undef, LDA1,LDA2,PMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    ZTA = Array{Float64,2}(undef, LDZTA,NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    P = parse(BlasInt, vs[2])
    ILO = parse(BlasInt, vs[3])
    IHI = parse(BlasInt, vs[4])
    ILOZ = parse(BlasInt, vs[5])
    IHIZ = parse(BlasInt, vs[6])
    JOB =  vs[7][1]
    COMPZ =  vs[8][1]
    if ( N<0 || N>min( LDA1, LDA2 ) )
        @error "Illegal N=$N"
    end
    if ( P<=0 || P>PMAX )
        @error "Illegal P=$P"
    end

    _isz,_jsz,_ksz = (N,N,P)
    for k in 1:_ksz
      vs = String[]
      while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
      end
      for i in 1:_isz
        _i0 = (i-1)*_jsz
        A[i,1:_jsz,k] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
      end
    end
    close(f)
# interp call 1

    #FOREIGN.dlacpy!( 'F', N, N, A(1,1,K), LDA1, AS(1,1,K), LDA1 )
    AS[1:N,1:N,1:P] .= A[1:N,1:N,1:P]

# interp call 2

    INFO = SLICOT.mb03vd!(N, P, ILO, IHI, A, TAU)

    if ( INFO!=0 )
        @warn "mb03vd returns info=$INFO"
        return
    end
    if ( LSAME( COMPZ, 'V' ) )
# interp call 3

    #FOREIGN.dlacpy!( 'L', N, N, A(1,1,K), LDA1, Z(1,1,K), LDZ1 )
        for k in 1:P
            Z[1:N,1:N,k] .= tril(A[1:N,1:N,k])
        end

# interp call 4

        INFO = SLICOT.mb03vy!(N, P, ILO, IHI, Z, TAU)
        @test INFO == 0
    INFO == 0 || return

        if ( INFO > 0 )
            @error "mb03vy returns info=$INFO"
        end
        # deviate from source here; continue even if COMPZ != 'V'
    end
# interp call 5

    INFO = SLICOT.mb03wd!(JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, A, Z, WR, WI, LDWORK)

    if ( INFO > 0 )
        @warn "mb03wd returns info=$INFO"
        println(io, "valid eigvals:")
        imin = max(ILO, INFO+1)
        show(io, "text/plain", WR[imin:N] + im * WI[imin:N])
        println(io)
    end # if

# interp call 6

    INFO = SLICOT.mb03wx!(ILO-1, P, A, WR, WI)

    if  IHI<N
# interp call 7

        INFO = SLICOT.mb03wx!(N-IHI, P, view(A,IHI+1:N, IHI+1:N, 1:P))

    end
# unable to translate write statement:
#   in do block [('40', 'I', 'N')]
#   write WR(I), WI(I)
# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    _nk = P
    show(io, "text/plain", A[1:_nr,1:_nc,1:_nk])
    println(io)

# interp output 2
    println(io, "Z:")
    _nc = N
    _nr = N
    _nk = P
    show(io, "text/plain", Z[1:_nr,1:_nc,1:_nk])
    println(io)

    SSQ = ZERO
    for K in 1:P
        KP1 = mod(K,P) + 1
# interp call 8

    # FOREIGN.dgemm!( 'T', 'N', N, N, N, ONE, Z(1,1,K), LDZ1, AS(1,1,K), LDA1, ZERO, ZTA, LDZTA )
        zta = Z[1:N,1:N,K]' * AS[1:N,1:N,K]
# interp call 9

    # FOREIGN.dgemm!( 'N', 'N', N, N, N, ONE, ZTA, LDZTA, Z(1,1,KP1), LDZ1, -ONE, A(1,1,K), LDA1 )
        A[1:N,1:N,K] .-= zta * Z[1:N,1:N,KP1]
        # SSQ = DLAPY2( SSQ, DLANGE( 'Frobenius', N, N, A(1,1,K), LDA1, DWORK ) )
        SSQ = hypot(SSQ, norm(A[1:N,1:N,K]))
    end
    println(io, "decomposition residual: ",SSQ)
end # run_mb03wd()
