# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02ed(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 20
    NMAX = 20
    LDB = KMAX*NMAX
    LDT = KMAX*NMAX
    LDWORK = NMAX*KMAX*KMAX + ( NMAX+2 )*KMAX
    T = Array{Float64,2}(undef, LDT,KMAX*NMAX)
    B = Array{Float64,2}(undef, LDB,KMAX*NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    K = parse(BlasInt, vs[2])
    NRHS = parse(BlasInt, vs[3])
    TYPET =  vs[4][1]
    M = N*K
if ( N<=0 || N>NMAX )
else
if ( K<=0 || K>KMAX )
else
if ( NRHS<=0 || NRHS>KMAX*NMAX )
else
if ( LSAME( TYPET, 'R' ) )

    vs = String[]
    _isz,_jsz = (K,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
else

    vs = String[]
    _isz,_jsz = (M,K)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
if ( LSAME( TYPET, 'R' ) )

    vs = String[]
    _isz,_jsz = (NRHS,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
else

    vs = String[]
    _isz,_jsz = (M,NRHS)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
# interp call 1

    INFO = SLICOT.mb02ed!(TYPET, K, N, NRHS, T, B, LDWORK)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( TYPET, 'R' ) )
# interp output 1
    println(io,"B:")
    _nc = M
    _nr = NRHS
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

else
# interp output 2
    println(io,"B:")
    _nc = NRHS
    _nr = M
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
end # if
end # if
    close(f)
end # run_X()
