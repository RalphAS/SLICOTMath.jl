# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02id(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX  = 20
    LMAX  = 20
    MMAX = 20
    NMAX = 20
    RBMAX = 20
    RCMAX = 20
    LDB  = KMAX*MMAX
    LDC  = KMAX*MMAX
    LDTC = MMAX*KMAX
    LDTR = KMAX
    LDWORK = 2*NMAX*LMAX*( LMAX + KMAX ) + ( 6 + NMAX )*LMAX + MMAX*KMAX*( LMAX + 1 ) + RBMAX + RCMAX
    TC = Array{Float64,2}(undef, LDTC,LMAX)
    TR = Array{Float64,2}(undef, LDTR,NMAX*LMAX)
    B = Array{Float64,2}(undef, LDB,RBMAX)
    C = Array{Float64,2}(undef, LDC,RCMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    L = parse(BlasInt, vs[2])
    M = parse(BlasInt, vs[3])
    N = parse(BlasInt, vs[4])
    RB = parse(BlasInt, vs[5])
    RC = parse(BlasInt, vs[6])
    JOB =  vs[7][1]
if ( K<=0 || K>KMAX )
elseif ( L<=0 || L>LMAX )
elseif ( M<=0 || M>MMAX )
elseif ( N<=0 || N>NMAX )
elseif ( ( LSAME( JOB, 'O' ) || LSAME( JOB, 'A' ) ) && ( ( RB<=0 ) || ( RB>RBMAX ) ) )
elseif ( ( LSAME( JOB, 'U' ) || LSAME( JOB, 'A' ) ) && ( ( RC<=0 ) || ( RC>RCMAX ) ) )
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
    _isz,_jsz = (K,N*L-L)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       TR[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
if ( LSAME( JOB, 'O' ) || LSAME( JOB, 'A' ) )

    vs = String[]
    _isz,_jsz = (M*K,RB)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
if ( LSAME( JOB, 'U' ) || LSAME( JOB, 'A' ) )

    vs = String[]
    _isz,_jsz = (N*L,RC)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       C[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
# interp call 1

    INFO = SLICOT.mb02id!(JOB, K, L, M, N, RB, RC, TC, TR, B, C)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( JOB, 'O' ) || LSAME( JOB, 'A' ) )
# interp output 1
    println(io,"B:")
    _nc = RB
    _nr = N*L
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

end # if
if ( LSAME( JOB, 'U' ) || LSAME( JOB, 'A' ) )
# interp output 2
    println(io,"C:")
    _nc = RC
    _nr = M*K
    show(io,"text/plain",C[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
    close(f)
end # run_X()
