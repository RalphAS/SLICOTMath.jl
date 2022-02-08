# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02gd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 20
    NMAX = 20
    NLMAX = 20
    LDRB = ( NLMAX + 1 )*KMAX
    LDT = KMAX*NMAX
    LDWORK = ( NLMAX + 1 )*KMAX*KMAX + ( 3 + NLMAX )*KMAX
    T = Array{Float64,2}(undef, LDT,NMAX*KMAX)
    RB = Array{Float64,2}(undef, LDRB,NMAX*KMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    K = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    NL = parse(BlasInt, vs[3])
    TRIU =  vs[4][1]
    TYPET = 'R'
    M = ( NL + 1 )*K
if ( N<=0 || N>NMAX )
elseif ( NL<=0 || NL>NLMAX )
elseif ( K<=0 || K>KMAX )
else

    vs = String[]
    _isz,_jsz = (K,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    INFO = SLICOT.mb02gd!(TYPET, TRIU, K, N, NL, 0, N, T, RB)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( TRIU, 'T' ) )
    SIZR = NL*K + 1
else
    SIZR = ( NL + 1 )*K
end # if
# interp output 1
    println(io,"RB:")
    _nc = N*K
    _nr = SIZR
    show(io,"text/plain",RB[1:_nr,1:_nc])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
