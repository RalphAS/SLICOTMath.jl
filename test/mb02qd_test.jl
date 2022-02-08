# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb02qd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    MMAX = 20
    NRHSMX = 20
    LDA = MMAX
    LDB = max( MMAX, NMAX )
    LDWORK = max(   min( MMAX, NMAX) + 3*NMAX + 1, 2*min( MMAX, NMAX) + NRHSMX )
    A = Array{Float64,2}(undef, LDA,NMAX)
    B = Array{Float64,2}(undef, LDB,NRHSMX)
    Y = Array{Float64,1}(undef, NMAX*NRHSMX)
    JPVT = Array{BlasInt,1}(undef, NMAX)
    SVAL = Array{Float64,1}(undef, 3)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    NRHS = parse(BlasInt, vs[3])
    RCOND = parse(Float64, replace(vs[4],'D'=>'E'))
    SVLMAX = parse(Float64, replace(vs[5],'D'=>'E'))
    JOB =  vs[6][1]
    INIPER =  vs[7][1]
if ( M<0 || M>MMAX )
else
if ( N<0 || N>NMAX )
else
if ( NRHS<0 || NRHS>NRHSMX )
else

    vs = String[]
    _isz,_jsz = (M,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,NRHS)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    RANK, INFO = SLICOT.mb02qd!(JOB, INIPER, M, N, NRHS, RCOND, SVLMAX, A, B, Y, JPVT, SVAL, LDWORK)
    println(io, "RANK = $RANK")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"B:")
    _nc = NRHS
    _nr = N
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
end # if
    close(f)
end # run_X()
