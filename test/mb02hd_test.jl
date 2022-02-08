# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02hd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 20
    LMAX  = 20
    MMAX = 20
    MLMAX = 10
    NMAX = 20
    NUMAX = 10
    LDRB = ( MLMAX + NUMAX + 1 )*LMAX
    LDTC = ( MLMAX + 1 )*KMAX
    LDTR = KMAX
    LDWORK = LDRB*LMAX + ( 2*NUMAX + 1 )*LMAX*KMAX + 2*LDRB*( KMAX + LMAX ) + LDRB + 6*LMAX
    TC = Array{Float64,2}(undef, LDTC,LMAX)
    TR = Array{Float64,2}(undef, LDTR,NMAX*LMAX)
    RB = Array{Float64,2}(undef, LDRB,NMAX*LMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    L = parse(BlasInt, vs[2])
    M = parse(BlasInt, vs[3])
    ML = parse(BlasInt, vs[4])
    N = parse(BlasInt, vs[5])
    NU = parse(BlasInt, vs[6])
    TRIU =  vs[7][1]
if ( K<0 || K>KMAX )
elseif ( L<0 || L>LMAX )
elseif ( M<=0 || M>MMAX )
elseif ( ML<0 || ML>MLMAX )
elseif ( N<=0 || N>NMAX )
elseif ( NU<0 || NU>NUMAX )
else

    vs = String[]
    _isz,_jsz = (ML*K+K,L)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       TC[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (K,NU*L)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       TR[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    S = (min(M*K,N*L)+L-1) ÷ L
# interp call 1

    INFO = SLICOT.mb02hd!(TRIU, K, L, M, ML, N, NU, 0, S, TC, TR, RB)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
    LENR = ( ML + NU + 1 )*L
    LENR = min( LENR, N*L )
# interp output 1
    println(io,"RB:")
    _nc = min( N*L, M*K )
    _nr = LENR
    show(io,"text/plain",RB[1:_nr,1:_nc])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
