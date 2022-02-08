# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02kd(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    KMAX = 20
    LMAX = 20
    MMAX = 20
    NMAX = 20
    RMAX = 20
    LDB  = LMAX*NMAX
    LDC  = KMAX*MMAX
    LDTC = MMAX*KMAX
    LDTR = KMAX
    LDWORK = 2*( KMAX*LMAX + KMAX*RMAX + LMAX*RMAX + 1 )*( MMAX + NMAX )
    TC = Array{Float64,2}(undef, LDTC,LMAX)
    TR = Array{Float64,2}(undef, LDTR,NMAX*LMAX)
    B = Array{Float64,2}(undef, LDB,RMAX)
    C = Array{Float64,2}(undef, LDC,RMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    L = parse(BlasInt, vs[2])
    M = parse(BlasInt, vs[3])
    N = parse(BlasInt, vs[4])
    R = parse(BlasInt, vs[5])
    LDBLK =  vs[6][1]
    TRANS =  vs[7][1]
    if ( K<=0 || K>KMAX )
        @error "Illegal K=$K"
    elseif ( L<=0 || L>LMAX )
        @error "Illegal L=$L"
    elseif ( M<=0 || M>MMAX )
        @error "Illegal M=$M"
    elseif ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
    elseif ( R<=0 || R>RMAX )
        @error "Illegal R=*R"
    end
    if ( LSAME( LDBLK, 'R' ) )

        vs = String[]
        _isz,_jsz = ((M-1)*K,L)
        while length(vs) < _isz*_jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for i in 1:_isz
            _i0 = (i-1)*_jsz
            TC[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
        end
        vs = String[]
        _isz,_jsz = (K,N*L)
        while length(vs) < _isz*_jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for i in 1:_isz
            _i0 = (i-1)*_jsz
            TR[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
        end
    else
        vs = String[]
        _isz,_jsz = (M*K,L)
        while length(vs) < _isz*_jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for i in 1:_isz
            _i0 = (i-1)*_jsz
            TC[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
        end
        vs = String[]
        _isz,_jsz = (K,(N-1)*L)
        while length(vs) < _isz*_jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        for i in 1:_isz
            _i0 = (i-1)*_jsz
            TR[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
        end

    end # if
    if ( LSAME( TRANS, 'N' ) )

    vs = String[]
    _isz,_jsz = (N*L,R)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    else

    vs = String[]
    _isz,_jsz = (M*K,R)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    end # if
    close(f)
    ALPHA = ONE
    BETA = ZERO
# interp call 1

    INFO = SLICOT.mb02kd!(LDBLK, TRANS, K, L, M, N, R, ALPHA, BETA, TC, TR, B, C)
    @test INFO == 0
    INFO == 0 || return
    if ( LSAME( TRANS, 'N' ) )
        # interp output 1
    println(io, "C:")
    _nc = R
    _nr = M*K
    show(io, "text/plain", C[1:_nr,1:_nc])
    println(io)

    else
# interp output 2
    println(io, "C:")
    _nc = R
    _nr = N*L
    show(io, "text/plain", C[1:_nr,1:_nc])
    println(io)

    end # if
end # run_mb02kd()
