# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb04ud(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LDA = MMAX
    LDE = MMAX
    LDQ = MMAX
    LDZ = NMAX
    LDWORK = max( NMAX,MMAX )
    A = Array{Float64,2}(undef, LDA,NMAX)
    E = Array{Float64,2}(undef, LDE,NMAX)
    Q = Array{Float64,2}(undef, LDQ,MMAX)
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    ISTAIR = Array{BlasInt,1}(undef, MMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    TOL = parse(Float64, replace(vs[3],'D'=>'E'))
if ( M<0 || M>MMAX )
elseif ( N<0 || N>NMAX )
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
    _isz,_jsz = (M,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       E[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    JOBQ = 'N'
    JOBZ = 'N'
# interp call 1

    RANKE, INFO = SLICOT.mb04ud!(JOBQ, JOBZ, M, N, A, E, Q, Z, ISTAIR, TOL)
    println(io, "RANKE = $RANKE")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"A:")
    _nc = N
    _nr = M
    show(io,"text/plain",A[1:_nr,1:_nc])
    println(io,)

# interp output 2
    println(io,"E:")
    _nc = N
    _nr = M
    show(io,"text/plain",E[1:_nr,1:_nc])
    println(io,)

# interp output 3
    println(io,"ISTAIR:")
    _nr = M
    show(io,"text/plain",ISTAIR[1:_nr])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
