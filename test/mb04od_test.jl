# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb04od(datfile, io=stdout)
    ZERO  = 0.0e0
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    PMAX = 20
    LDA = PMAX
    LDB = NMAX
    LDC = PMAX
    LDR = NMAX
    LDWORK = max( NMAX-1,MMAX )
    R = Array{Float64,2}(undef, LDR,NMAX)
    A = Array{Float64,2}(undef, LDA,NMAX)
    B = Array{Float64,2}(undef, LDB,MMAX)
    C = Array{Float64,2}(undef, LDC,MMAX)
    TAU = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    M = parse(BlasInt, vs[2])
    P = parse(BlasInt, vs[3])
    UPLO =  vs[4][1]
if ( N<0 || N>NMAX )
else
if ( M<0 || M>MMAX )
else
if ( P<0 || P>PMAX )
else

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       R[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (P,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (P,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       C[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    SLICOT.mb04od!(UPLO, N, M, P, R, A, B, C, TAU)


    triu!(R)
# interp output 1
    println(io,"R:")
    _nc = N
    _nr = N
    show(io,"text/plain",R[1:_nr,1:_nc])
    println(io,)

if ( M>0 )
# interp output 2
    println(io,"B:")
    _nc = M
    _nr = N
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

if ( P>0 )
# interp output 3
    println(io,"C:")
    _nc = M
    _nr = P
    show(io,"text/plain",C[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
end # if
end # if
    close(f)
end # run_X()
