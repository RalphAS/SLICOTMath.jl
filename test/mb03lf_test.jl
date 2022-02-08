# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb03lf(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDB = NMAX÷2
    LDFG = NMAX÷2
    LDQ = 2*NMAX
    LDU = NMAX
    LDZ  = NMAX
    LDWORK = 10*NMAX*NMAX + max( NMAX*NMAX + max( NMAX÷2 + 252, 432 ), max( 8*NMAX +  48, 171 ) )
    LIWORK = max( NMAX + 18, NMAX÷2 + 48, 5*NMAX÷2 + 1 )
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX÷2)
    FG = Array{Float64,2}(undef, LDFG,NMAX÷2+1)
    Q = Array{Float64,2}(undef, LDQ,2*NMAX)
    U = Array{Float64,2}(undef, LDU,2*NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX÷2)
    ALPHAI = Array{Float64,1}(undef, NMAX÷2)
    BETA = Array{Float64,1}(undef, NMAX÷2)
    BWORK = Array{BlasBool,1}(undef, NMAX÷2)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    COMPQ =  vs[1][1]
    COMPU =  vs[2][1]
    ORTH =  vs[3][1]
    N = parse(BlasInt, vs[4])
if ( N<0 || N>NMAX || mod( N, 2 )!=0 )
else
    M = N÷2

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       Z[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       FG[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    NEIG, INFO, IWARN = SLICOT.mb03lf!(COMPQ, COMPU, ORTH, N, Z, B, FG, Q, U, ALPHAR, ALPHAI, BETA, LIWORK)
    println(io, "NEIG = $NEIG")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

if ( INFO!=0 )
else
# interp output 1
    println(io,"Z:")
    _nc = N
    _nr = N
    show(io,"text/plain",Z[1:_nr,1:_nc])
    println(io,)

# interp output 2
    println(io,"ALPHAR:")
    _nr = M
    show(io,"text/plain",ALPHAR[1:_nr])
    println(io,)

# interp output 3
    println(io,"ALPHAI:")
    _nr = M
    show(io,"text/plain",ALPHAI[1:_nr])
    println(io,)

# interp output 4
    println(io,"BETA:")
    _nr = M
    show(io,"text/plain",BETA[1:_nr])
    println(io,)

if ( LSAME( COMPQ, 'C' ) && NEIG>0 )
# interp output 5
    println(io,"Q:")
    _nc = NEIG
    _nr = N
    show(io,"text/plain",Q[1:_nr,1:_nc])
    println(io,)

end # if
if ( LSAME( COMPU, 'C' ) && NEIG>0 )
# interp output 6
    println(io,"U:")
    _nc = NEIG
    _nr = N
    show(io,"text/plain",U[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
    close(f)
end # run_X()
