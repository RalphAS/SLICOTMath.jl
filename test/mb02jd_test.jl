# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02jd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 10
    LMAX = 10
    MMAX = 20
    NMAX = 20
    LDR  = NMAX*LMAX
    LDQ  = MMAX*KMAX
    LDTC = MMAX*KMAX
    LDTR = KMAX
    LDWORK = ( MMAX*KMAX + NMAX*LMAX ) *( LMAX + 2*KMAX ) + 6*LMAX + MMAX*KMAX + NMAX*LMAX
    TC = Array{Float64,2}(undef, LDTC,LMAX)
    TR = Array{Float64,2}(undef, LDTR,NMAX*LMAX)
    Q = Array{Float64,2}(undef, LDQ,NMAX*LMAX)
    R = Array{Float64,2}(undef, LDR,NMAX*LMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    L = parse(BlasInt, vs[2])
    M = parse(BlasInt, vs[3])
    N = parse(BlasInt, vs[4])
    JOB =  vs[5][1]
if ( K<=0 || K>KMAX )
elseif ( L<=0 || L>LMAX )
elseif ( M<=0 || M>MMAX )
elseif ( N<=0 || N>NMAX )
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
    S = (min(M*K,N*L)+L-1) ÷ L
# interp call 1

    INFO = SLICOT.mb02jd!(JOB, K, L, M, N, 0, S, TC, TR, Q, R)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( JOB, 'Q' ) )
# interp output 1
    println(io,"Q:")
    _nc = min( N*L, M*K )
    _nr = M*K
    show(io,"text/plain",Q[1:_nr,1:_nc])
    println(io,)

end # if
# interp output 2
    println(io,"R:")
    _nc = min( N*L, M*K )
    _nr = N*L
    show(io,"text/plain",R[1:_nr,1:_nc])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
