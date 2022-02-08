# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb04fd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDA  = NMAX÷2
    LDB = NMAX÷2
    LDDE = NMAX÷2
    LDFG = NMAX÷2
    LDQ = NMAX
    LDWORK = 3*NMAX*NMAX÷4
    A = Array{Float64,2}(undef, LDA,NMAX÷2)
    DE = Array{Float64,2}(undef, LDDE,NMAX÷2+1)
    B = Array{Float64,2}(undef, LDB,NMAX÷2)
    FG = Array{Float64,2}(undef, LDFG,NMAX÷2+1)
    Q = Array{Float64,2}(undef, LDQ,NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX÷2)
    ALPHAI = Array{Float64,1}(undef, NMAX÷2)
    BETA = Array{Float64,1}(undef, NMAX÷2)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, NMAX÷2+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    COMPQ =  vs[2][1]
    N = parse(BlasInt, vs[3])

    vs = String[]
    _isz,_jsz = (N÷2,N÷2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N÷2,N÷2+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       DE[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N÷2,N÷2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N÷2,N÷2+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       FG[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
if ( N<0 || N>NMAX || mod( N, 2 )!=0 )
else
# interp call 1

    INFO = SLICOT.mb04fd!(JOB, COMPQ, N, A, DE, B, FG, Q, ALPHAR, ALPHAI, BETA)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( JOB, 'T' ) )
# interp output 1
    println(io,"A:")
    _nc = N÷2
    _nr = N÷2
    show(io,"text/plain",A[1:_nr,1:_nc])
    println(io,)

end # if
# interp output 2
    println(io,"DE:")
    _nc = N÷2+1
    _nr = N÷2
    show(io,"text/plain",DE[1:_nr,1:_nc])
    println(io,)

if ( LSAME( JOB, 'T' ) )
# interp output 3
    println(io,"B:")
    _nc = N÷2
    _nr = N÷2
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

end # if
# interp output 4
    println(io,"FG:")
    _nc = N÷2+1
    _nr = N÷2
    show(io,"text/plain",FG[1:_nr,1:_nc])
    println(io,)

if ( !LSAME( COMPQ, 'N' ) )
# interp output 5
    println(io,"Q:")
    _nc = N
    _nr = N
    show(io,"text/plain",Q[1:_nr,1:_nc])
    println(io,)

end # if
# interp output 6
    println(io,"ALPHAR:")
    _nr = N÷2
    show(io,"text/plain",ALPHAR[1:_nr])
    println(io,)

# interp output 7
    println(io,"ALPHAI:")
    _nr = N÷2
    show(io,"text/plain",ALPHAI[1:_nr])
    println(io,)

# interp output 8
    println(io,"BETA:")
    _nr = N÷2
    show(io,"text/plain",BETA[1:_nr])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
