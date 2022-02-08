# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02vd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LDA = NMAX
    LDB = MMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    IPIV = Array{BlasInt,1}(undef, NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    TRANS =  vs[3][1]
if ( N<0 || N>NMAX )
else

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
if ( M<0 || M>MMAX )
else

    vs = String[]
    _isz,_jsz = (M,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    INFO = SLICOT.mb02vd!(TRANS, M, N, A, IPIV, B)
    @test INFO == 0
    INFO == 0 || return

if ( INFO==0 )
# interp output 1
    println(io,"B:")
    _nc = N
    _nr = M
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

else
end # if
end # if
end # if
    close(f)
end # run_X()
