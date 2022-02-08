# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03zd(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    NMAX = 200
    LDG = NMAX
    LDRES = 2*NMAX
    LDS = NMAX
    LDT = NMAX
    LDU1 = NMAX
    LDU2 = NMAX
    LDUS = 2*NMAX
    LDUU = 2*NMAX
    LDV1 = NMAX
    LDV2 = NMAX
    LDWORK = 3*NMAX*NMAX + 7*NMAX
    SELECT = Array{BlasBool,1}(undef, NMAX)
    SCALE = Array{Float64,1}(undef, NMAX)
    S = Array{Float64,2}(undef, LDS,NMAX)
    T = Array{Float64,2}(undef, LDT,NMAX)
    G = Array{Float64,2}(undef, LDG,NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX)
    U2 = Array{Float64,2}(undef, LDU2,NMAX)
    V1 = Array{Float64,2}(undef, LDV1,NMAX)
    V2 = Array{Float64,2}(undef, LDV2,NMAX)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
    US = Array{Float64,2}(undef, LDUS,2*NMAX)
    UU = Array{Float64,2}(undef, LDUU,2*NMAX)
    IWORK = Array{BlasInt,1}(undef, 2*NMAX)
    LWORK = Array{BlasBool,1}(undef, 2*NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RES = Array{Float64,2}(undef, LDRES,NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    ILO = parse(BlasInt, vs[2])
    WHICH =  vs[3][1]
    METH =  vs[4][1]
    STAB =  vs[5][1]
    BALANC =  vs[6][1]
    ORTBAL =  vs[7][1]
    if ( N<=0 || N>NMAX )
        @error "illegal N=$N"
    end
    if LSAME(WHICH,'S')
        vs = String[]
        _isz = N
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
       SELECT[1:_isz] .= parsex.(Bool, vs[1:_isz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       S[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    if LSAME(WHICH,'A') && LSAME(METH,'L')
        vs = String[]
        _isz,_jsz = (N,N)
        while length(vs) < _isz*_jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for i in 1:_isz
            _i0 = (i-1)*_jsz
            G[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
        end
    end

    if LSAME(BALANC,'P') || LSAME(BALANC,'S') || LSAME(BALANC,'B')
        vs = String[]
        _isz = N
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
       SCALE[1:_isz] .= parsex.(Float64, vs[1:_isz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U1[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U2[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       V1[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       V2[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    M, INFO = SLICOT.mb03zd!(WHICH, METH, STAB, BALANC, ORTBAL, SELECT, N, 2*N, ILO, SCALE, S, T, G, U1, U2, V1, V2, WR, WI, US, UU, IWORK)
    println(io, "M = $M")
    @test INFO == 0
    INFO == 0 || return

    println(io, "eigvals:")
    show(io, "text/plain", WR[1:N] + im *WI[1:N])
    println(io)

    if ( LSAME( STAB, 'S' )||LSAME( STAB, 'B' ) )
# interp output 1
        println(io, "US:")
        _nc = M
        _nr = 2*N
        show(io, "text/plain", US[1:_nr,1:_nc])
        println(io)

        if ( LSAME( ORTBAL, 'B' )||LSAME( BALANC, 'N' )|| LSAME( BALANC, 'P' ) )
# interp call 2

            # FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, 2*N, ONE, US, LDUS, US, LDUS, ZERO, RES, LDRES )
            BLAS.gemm!('T','N',ONE,US[1:2N,1:M],US[1:2N,1:M],ZERO,view(RES,1:M,1:M))
            resid = norm(RES[1:M,1:M] - I)
            @test resid < 1e-12
            println(io, "orth. resid norm of US: $resid")
        end # if
# interp call 3

        # FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, N, ONE, US, LDUS, US(N+1,1), LDUS, ZERO, RES, LDRES )
        BLAS.gemm!('T','N',ONE,US[1:N,1:M],US[N+1:2N,1:M],ZERO,view(RES,1:M,1:M))

# interp call 4

        # FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, N, -ONE, US(N+1,1), LDUS, US, LDUS, ONE, RES, LDRES )
        BLAS.gemm!('T','N',-ONE,US[N+1:2N,1:M],US[1:N,1:M],ONE,view(RES,1:M,1:M))

        resid = norm(RES[1:M,1:M])
        @test resid < 1e-12
        println(io, "sympl. resid norm of US: $resid")
    end
    if ( LSAME( STAB, 'U' )||LSAME( STAB, 'B' ) )
        # interp output 2
        println(io, "UU:")
        _nc = M
        _nr = 2*N
        show(io, "text/plain", UU[1:_nr,1:_nc])
        println(io)

        if ( LSAME( ORTBAL, 'B' )||LSAME( BALANC, 'N' )|| LSAME( BALANC, 'P' ) )
# interp call 5

            #FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, 2*N, ONE, UU, LDUU, UU, LDUU, ZERO, RES, LDRES )
            BLAS.gemm!('T','N',ONE,UU[1:2N,1:M],UU[1:2N,1:M],ZERO,view(RES,1:M,1:M))
            resid = norm(RES[1:M,1:M] - I)
            @test resid < 1e-12
            println(io, "orth. resid norm of UU: $resid")
        end # if
# interp call 6

        # FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, N, ONE, UU, LDUU, UU(N+1,1), LDUU, ZERO, RES, LDRES )
        BLAS.gemm!('T','N',ONE,UU[1:N,1:M],UU[N+1:2N,1:M],ZERO,view(RES,1:M,1:M))

# interp call 7

        #FOREIGN.dgemm!( 'Transpose', 'No Transpose', M, M, N, -ONE, UU(N+1,1), LDUU, UU, LDUU, ONE, RES, LDRES )
        BLAS.gemm!('T','N',-ONE,UU[N+1:2N,1:M],UU[1:N,1:M],ONE,view(RES,1:M,1:M))

        resid = norm(RES[1:M,1:M])
        @test resid < 1e-12
        println(io, "sympl. resid norm of UU: $resid")

    end # if
end # run_mb03zd()
