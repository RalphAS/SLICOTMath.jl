# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03kd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 6
    NMAX = 50
    LDA1 = NMAX
    LDA2 = NMAX
    LDQ1 = NMAX
    LDQ2 = NMAX
    LDWORK = max( 42*KMAX + NMAX, 80*KMAX - 48, KMAX + max( 2*NMAX, 8*KMAX ) )
    LIWORK = max( 4*KMAX, 2*KMAX + NMAX )
    HUND = 1.0e2
    ZERO = 0.0e0
    QIND = Array{BlasInt,1}(undef, KMAX)
    S = Array{BlasInt,1}(undef, KMAX)
    A = Array{Float64,3}(undef, LDA1,LDA2,KMAX)
    Q = Array{Float64,3}(undef, LDQ1,LDQ2,KMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX)
    ALPHAI = Array{Float64,1}(undef, NMAX)
    BETA = Array{Float64,1}(undef, NMAX)
    SCAL = Array{BlasInt,1}(undef, NMAX)
    ND = Array{BlasInt,1}(undef, KMAX)
    NI = Array{BlasInt,1}(undef, KMAX)
    SELECT = Array{Bool,1}(undef, NMAX)
    T = Array{Float64,1}(undef, NMAX*NMAX*KMAX)
    IXT = Array{BlasInt,1}(undef, KMAX)
    QK = Array{Float64,1}(undef, NMAX*NMAX*KMAX)
    IXQ = Array{BlasInt,1}(undef, KMAX)
# WARNING: desperate attempt to initialize M
    M = 0
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    LDQ = Array{BlasInt,1}(undef, KMAX)
    LDT = Array{BlasInt,1}(undef, KMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    DEFL =  vs[2][1]
    COMPQ =  vs[3][1]
    STRONG =  vs[4][1]
    K = parse(BlasInt, vs[5])
    N = parse(BlasInt, vs[6])
    H = parse(BlasInt, vs[7])
    ILO = parse(BlasInt, vs[8])
    IHI = parse(BlasInt, vs[9])
    if ( N<0 || N>NMAX )
        @error "illegal N"
    end
    TOL = HUND

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

    if ( LSAME( COMPQ, 'U' ) )
        vs = String[]
        _isz,_jsz,_ksz = (N,N,K)
        while length(vs) < _isz*_jsz*_ksz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for k in 1:_ksz
            for i in 1:_isz
                _i0 = (i-1)*_jsz + (k-1)*_jsz*_isz
                Q[i,1:_jsz,k] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
            end
        end
    end

    if ( LSAME( COMPQ, 'P' ) )

        vs = String[]
        _isz = K
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        QIND[1:_isz] .= parsex.(BlasInt, vs)
        vs = String[]
        _isz,_jsz,_ksz = (N,N,K)
        while length(vs) < _isz*_jsz*_ksz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for k in 1:_ksz
            if QIND[k] > 0
                for i in 1:_isz
                    _i0 = (i-1)*_jsz + (k-1)*_jsz*_isz
                    Q[i,1:_jsz,QIND[k]] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
                end
            end
        end
    end # if
    close(f)
    if LSAME(JOB,'E')
        JOB = 'S'
    end
    println(io, "JOB = $JOB")
    println(io, "COMPQ = $COMPQ")
# interp call 1

    INFO, IWARN = SLICOT.mb03bd!(JOB, DEFL, COMPQ, QIND, K, N, H, ILO, IHI, S, A, Q, ALPHAR, ALPHAI, BETA, SCAL, LIWORK, LDWORK)
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

    if ( IWARN!=0 )
        @warn "mb03bd returned IWARN=$IWARN"
    end
    verbose = false
    if verbose
        ## added for diagnostic
        if ( LSAME( JOB, 'S' ) || LSAME( JOB, 'T' ) )
            # interp output 1
            println(io, "after Schur (HessIdx=$H) A:")
            _nc = N
            _nr = N
            _nk = K
            show(io, "text/plain", A[1:_nr,1:_nc,1:_nk])
            println(io)
        end # if
        if ( LSAME( COMPQ, 'U' ) || LSAME( COMPQ, 'I' ) )
            println(io, "Q:")
            _nc = N
            _nr = N
            _nk = K
            show(io, "text/plain", Q[1:_nr,1:_nc,1:_nk])
            println(io)
        elseif ( LSAME( COMPQ, 'P' ) )
            println(io, "TODO: show Q for COMPQ='P'")
            for L in 1:K
            #if ( QIND[L]>0 )
                #println(io, "Q:")
                #show(io, "text/plain", Q)
                #println(io)
                #nd # if
            end # for
        end # if
    end ## diagnostic

    # prepare data for calling MB03KD with screwy data structures and factor ordering
    for L in 1:K
        ND[L] = max(1,N)
        NI[L] = 0
        LDT[L] = max(1,N)
        IXT[L] = (L-1)*LDT[L]*N + 1
        LDQ[L] = max(1,N)
        IXQ[L] = IXT[L]
        if ( L<= K÷2 )
            i = S[ K - L + 1 ]
            S[K-L+1] = S[L]
            S[L] = i
        end # if
    end
    for L in 1:K
        # interp call 2
        #    FOREIGN.dlacpy!( 'Full', N, N, A( 1, 1, K-L+1 ), LDA1, T( IXT( L ) ), LDT( L ) )
        ir1 = IXT[L]
        for ic in 1:N
            T[ir1:ir1+N-1] .= A[1:N,ic,K-L+1]
            ir1 += LDT[L]
        end
    end

    if verbose
        println(io, "stored T:")
        for k in 1:K
            println(io, "K=$k")
            Ttmp = reshape(T[IXT[k]:IXT[k]+LDT[k]*N-1],N,N)
            show(io, "text/plain", Ttmp)
            println(io)
        end
    end # diagnostic

    if ( LSAME( COMPQ, 'U' ) || LSAME( COMPQ, 'I' ) )
        COMPQ = 'U'

        for L in 1:K
            # interp call 3
            # FOREIGN.dlacpy!( 'Full', N, N, Q( 1, 1, K-L+1 ), LDQ1, QK( IXQ( L ) ), LDQ( L ) )
            ir1 = IXQ[L]
            for ic in 1:N
                QK[ir1:ir1+N-1] .= Q[1:N,ic,K-L+1]
                ir1 += LDQ[L]
            end
        end
    elseif ( LSAME( COMPQ, 'P' ) )
        COMPQ = 'W'
        for L in 1:K
            if QIND[L] < 0
                QIND[L] = 2
            end
            P = QIND[L]
            if P != 0
                ir1 = IXQ[P]
                Qk[ir1:ir1+N-1] .= Q[1:N,ic,K-P+1]
                ir1 += LDQ[P]
            end
        end
    end # if

    if verbose    # diag
        println(io, "stored Q:")
        for k in 1:K
            println(io, "K=$k")
            Ttmp = reshape(QK[IXQ[k]:IXQ[k]+LDQ[k]*N-1],N,N)
            show(io, "text/plain", Ttmp)
            println(io)
        end
    end
    SELECT .= false
    for i in 1:N
        SELECT[i] = ALPHAR[i] < 0
    end
        # interp output 1
    println(io, "ALPHAR:")
    _nr = N
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

    # interp output 2
    println(io, "ALPHAI:")
    _nr = N
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

    # interp output 3
    println(io, "BETA:")
    _nr = N
    show(io, "text/plain", BETA[1:_nr])
    println(io)

    # interp output 4
    println(io, "SCAL:")
    _nr = N
    show(io, "text/plain", SCAL[1:_nr])
    println(io)

    # interp call 4
    select = convert(Vector{BlasBool}, SELECT)
    println(io, "select:"); show(io, "text/plain", select[1:K]); println()
    M, INFO = SLICOT.mb03kd!(COMPQ, QIND, STRONG, K, N, H, ND, NI, S, select, T, LDT, IXT, QK, LDQ, IXQ, TOL)
    @test INFO == 0
    if INFO != 0
        println(io, "mb03kd returned INFO = $INFO")
        return
    end

    if ( INFO!=0 )
        #@error "mb03kd returned INFO=$INFO"
    end
    Tf = zeros(N,N,K)
    for L in 1:K
        P = K - L + 1
        ir1 = IXT[P]
        for ic in 1:N
            Tf[1:N,ic,L] .= T[ir1:ir1+N-1]
            ir1 += LDT[P]
        end
    end
    # unable to translate write loop:
    #  write ( IXT( P ) + I - 1 + ( J - 1 )*LDT( P ) ), J = 1, N 
    # interp output 5
    println(io, "T:")
    show(io, "text/plain", Tf)
    println(io)

    if ( LSAME( COMPQ, 'U' ) || LSAME( COMPQ, 'I' ) )
        Qf = zeros(N,N,K)
        for L in 1:K
            P = K - L + 1
            ir1 = IXQ[P]
            for ic in 1:N
                Qf[1:N,ic,L] .= QK[ir1:ir1+N-1]
                ir1 += LDQ[P]
            end
        end
        # unable to translate write loop:
        #  write ( IXQ( P ) + I - 1 + ( J - 1 )*LDQ( P ) ), J = 1, N
        # interp output 6
        println(io, "QK:")
        show(io, "text/plain", Qf)
        println(io)

    elseif ( LSAME( COMPQ, 'W' ) )
        for L in 1:K
            if ( QIND[L]>0 )
                println(io, "factor $(QIND[L]):")
                Qf = zeros(N,N)
                P = K - QIND[L] + 1
                ir1 = IXQ[P]
                for ic in 1:N
                    Qf[1:N,ic] .= QK[ir1:ir1+N-1]
                    ir1 += LDQ[P]
                end
                show(io, "text/plain", Qf)
                println(io)
            end
        end
            # unable to translate write loop:
            #  write ( IXQ( P ) + I - 1 + ( J - 1 )*LDQ( P ) ), J = 1, N 
            # interp output 7
    end # if
end # module
