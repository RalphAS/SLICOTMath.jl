# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb04ed(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 60
    LDB = NMAX÷2
    LDFG = NMAX÷2
    LDQ = NMAX
    LDU1 = NMAX÷2
    LDU2 = NMAX÷2
    LDWORK = 3*NMAX^2÷2 + max( NMAX, 24 ) + 3
    LDZ = NMAX
    LIWORK = NMAX + 9
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX÷2)
    FG = Array{Float64,2}(undef, LDFG,NMAX÷2+1)
    Q = Array{Float64,2}(undef, LDQ,NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX÷2)
    U2 = Array{Float64,2}(undef, LDU2,NMAX÷2)
    ALPHAR = Array{Float64,1}(undef, NMAX÷2)
    ALPHAI = Array{Float64,1}(undef, NMAX÷2)
    BETA = Array{Float64,1}(undef, NMAX÷2)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    COMPQ =  vs[2][1]
    COMPU =  vs[3][1]
    N = parse(BlasInt, vs[4])

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
if ( LSAME( COMPU, 'U' ) )

    vs = String[]
    _isz,_jsz = (N÷2,N÷2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U1[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N÷2,N÷2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U2[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
if ( N<0 || N>NMAX || mod( N, 2 )!=0 )
else
# interp call 1

    INFO = SLICOT.mb04ed!(JOB, COMPQ, COMPU, N, Z, B, FG, Q, U1, U2, ALPHAR, ALPHAI, BETA, LIWORK)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"Z:")
    _nc = N
    _nr = N
    show(io,"text/plain",Z[1:_nr,1:_nc])
    println(io,)

# interp output 2
    println(io,"B:")
    _nc = N÷2
    _nr = N÷2
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

# interp output 3
    println(io,"FG:")
    _nc = N÷2+1
    _nr = N÷2
    show(io,"text/plain",FG[1:_nr,1:_nc])
    println(io,)

# interp output 4
    println(io,"ALPHAR:")
    _nr = N÷2
    show(io,"text/plain",ALPHAR[1:_nr])
    println(io,)

# interp output 5
    println(io,"ALPHAI:")
    _nr = N÷2
    show(io,"text/plain",ALPHAI[1:_nr])
    println(io,)

# interp output 6
    println(io,"BETA:")
    _nr = N÷2
    show(io,"text/plain",BETA[1:_nr])
    println(io,)

# interp output 7
    println(io,"Q:")
    _nc = N
    _nr = N
    show(io,"text/plain",Q[1:_nr,1:_nc])
    println(io,)

if ( !LSAME( COMPU, 'N' ) )
# interp output 8
    println(io,"U1:")
    _nc = N÷2
    _nr = N÷2
    show(io,"text/plain",U1[1:_nr,1:_nc])
    println(io,)

# interp output 9
    println(io,"U2:")
    _nc = N÷2
    _nr = N÷2
    show(io,"text/plain",U2[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
    close(f)
end # run_X()
