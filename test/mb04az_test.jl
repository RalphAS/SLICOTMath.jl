# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04az(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDB = NMAX
    LDC =   NMAX
    LDD = NMAX
    LDFG = NMAX
    LDQ = 2*NMAX
    LDU = NMAX
    LDWORK = 18*NMAX*NMAX + NMAX + max( 2*NMAX, 24 ) + 3
    LDZ = NMAX
    LIWORK = 2*NMAX + 9
    LZWORK = 8*NMAX + 28
    Z = Array{ComplexF64,2}(undef, LDZ,NMAX)
    B = Array{ComplexF64,2}(undef, LDB,NMAX)
    FG = Array{ComplexF64,2}(undef, LDFG,NMAX)
    D = Array{ComplexF64,2}(undef, LDD,NMAX)
    C = Array{ComplexF64,2}(undef, LDC,NMAX)
    Q = Array{ComplexF64,2}(undef, LDQ,2*NMAX)
    U = Array{ComplexF64,2}(undef, LDU,2*NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX)
    ALPHAI = Array{Float64,1}(undef, NMAX)
    BETA = Array{Float64,1}(undef, NMAX)
    BWORK = Array{BlasBool,1}(undef, NMAX)
    ZWORK = Array{ComplexF64,1}(undef, LZWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    COMPQ =  vs[2][1]
    COMPU =  vs[3][1]
    N = parse(BlasInt, vs[4])
    if ( N<0 || N>NMAX || mod( N, 2 )!=0 )
        @error "Illegal N=$N"
    end
    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       Z[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end
    vs = String[]
    _isz,_jsz = (N÷2,N÷2)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end
    vs = String[]
    _isz,_jsz = (N÷2,N÷2+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       FG[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end

# interp call 1
    close(f)

    INFO = SLICOT.mb04az!(JOB, COMPQ, COMPU, N, Z, B, FG, D, C, Q, U, ALPHAR, ALPHAI, BETA, LIWORK, BWORK)
    @test INFO == 0
    INFO == 0 || return


    M = N÷2
if ( LSAME( JOB, 'T' ) )
# interp output 1
    println(io, "Z:")
    _nc = N
    _nr = N
    show(io, "text/plain", Z[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "B:")
    _nc = N
    _nr = N
    show(io, "text/plain", UpperTriangular(B[1:_nr,1:_nc]))
    println(io)

# interp output 3
    println(io, "FG:")
    _nc = N
    _nr = N
    show(io, "text/plain", UpperTriangular(FG[1:_nr,1:_nc]))
    println(io)

# interp output 4
    println(io, "D:")
    _nc = N
    _nr = N
    show(io, "text/plain", D[1:_nr,1:_nc])
    println(io)

# interp output 5
    println(io, "C:")
    _nc = N
    _nr = N
    show(io, "text/plain", C[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( COMPQ, 'C' ) )
# interp output 6
    println(io, "Q:")
    _nc = 2*N
    _nr = 2*N
    show(io, "text/plain", Q[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( COMPU, 'C' ) )
# interp output 7
    println(io, "U:")
    _nc = 2*N
    _nr = N
    show(io, "text/plain", U[1:_nr,1:_nc])
    println(io)

end # if
# interp output 8
    println(io, "ALPHAR:")
    _nr = N
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

# interp output 9
    println(io, "ALPHAI:")
    _nr = N
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

# interp output 10
    println(io, "BETA:")
    _nr = N
    show(io, "text/plain", BETA[1:_nr])
    println(io)

end # run_mb04az()
