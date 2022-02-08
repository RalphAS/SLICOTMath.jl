# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04ad(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDH = NMAX
    LDQ1  = NMAX
    LDQ2 = NMAX
    LDT = NMAX
    LDU11 = NMAX÷2
    LDU12 = NMAX÷2
    LDU21  = NMAX÷2
    LDU22 = NMAX÷2
    LDWORK = 3*NMAX*NMAX + max( NMAX, 48 ) + 6
    LDZ = NMAX
    LIWORK = NMAX + 18
    ZERO = 0.0e0
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    H = Array{Float64,2}(undef, LDH,NMAX)
    Q1 = Array{Float64,2}(undef, LDQ1,NMAX)
    Q2 = Array{Float64,2}(undef, LDQ2,NMAX)
    U11 = Array{Float64,2}(undef, LDU11,NMAX÷2)
    U12 = Array{Float64,2}(undef, LDU12,NMAX÷2)
    U21 = Array{Float64,2}(undef, LDU21,NMAX÷2)
    U22 = Array{Float64,2}(undef, LDU22,NMAX÷2)
    T = Array{Float64,2}(undef, LDT,NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX÷2)
    ALPHAI = Array{Float64,1}(undef, NMAX÷2)
    BETA = Array{Float64,1}(undef, NMAX÷2)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    COMPQ1 =  vs[2][1]
    COMPQ2 =  vs[3][1]
    COMPU1 =  vs[4][1]
    COMPU2 =  vs[5][1]
    N = parse(BlasInt, vs[6])
    if ( N<0 || N>NMAX )
        @error "illegal N=$N"
    end

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
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       H[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    INFO = SLICOT.mb04ad!(JOB, COMPQ1, COMPQ2, COMPU1, COMPU2, N, Z, H, Q1, Q2, U11, U12, U21, U22, T, ALPHAR, ALPHAI, BETA, LIWORK)
    @test INFO == 0
    if ( INFO!=0 )
        @warn "mb04ad returns info=$INFO"
        return
    end
    M = N÷2
# interp call 2

    #FOREIGN.dlaset!( 'Full', M, M, ZERO, ZERO, Z( M+1, 1 ), LDZ )
    Z[M+1:M,1:M] .= ZERO

# interp output 1
    println(io, "T:")
    _nc = N
    _nr = N
    show(io, "text/plain", T[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "Z:")
    _nc = N
    _nr = N
    show(io, "text/plain", Z[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "H:")
    _nc = N
    _nr = N
    show(io, "text/plain", H[1:_nr,1:_nc])
    println(io)

if ( LSAME( COMPQ1, 'I' ) )
# interp output 4
    println(io, "Q1:")
    _nc = N
    _nr = N
    show(io, "text/plain", Q1[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( COMPQ2, 'I' ) )
# interp output 5
    println(io, "Q2:")
    _nc = N
    _nr = N
    show(io, "text/plain", Q2[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( COMPU1, 'I' ) )
# interp output 6
    println(io, "U11:")
    _nc = M
    _nr = M
    show(io, "text/plain", U11[1:_nr,1:_nc])
    println(io)

# interp output 7
    println(io, "U12:")
    _nc = M
    _nr = M
    show(io, "text/plain", U12[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( COMPU2, 'I' ) )
# interp output 8
    println(io, "U21:")
    _nc = M
    _nr = M
    show(io, "text/plain", U21[1:_nr,1:_nc])
    println(io)

# interp output 9
    println(io, "U22:")
    _nc = M
    _nr = M
    show(io, "text/plain", U22[1:_nr,1:_nc])
    println(io)

end # if
# interp output 10
    println(io, "ALPHAR:")
    _nr = M
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

# interp output 11
    println(io, "ALPHAI:")
    _nr = M
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

# interp output 12
    println(io, "BETA:")
    _nr = M
    show(io, "text/plain", BETA[1:_nr])
    println(io)

end # run_mb04ad()
