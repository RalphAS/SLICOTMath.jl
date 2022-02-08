# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04bd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 50
    LDA = NMAX÷2
    LDB = NMAX÷2
    LDC1 = NMAX÷2
    LDC2 = NMAX÷2
    LDDE = NMAX÷2
    LDF = NMAX÷2
    LDQ1 = NMAX
    LDQ2 = NMAX
    LDVW = NMAX÷2
    LDWORK = 2*NMAX*NMAX + max( 4*NMAX, 36 )
    LIWORK = max( NMAX + 12, 2*NMAX + 3 )
    A = Array{Float64,2}(undef, LDA,NMAX÷2)
    DE = Array{Float64,2}(undef, LDDE,NMAX÷2+1)
    C1 = Array{Float64,2}(undef, LDC1,NMAX÷2)
    VW = Array{Float64,2}(undef, LDVW,NMAX÷2+1)
    Q1 = Array{Float64,2}(undef, LDQ1,NMAX)
    Q2 = Array{Float64,2}(undef, LDQ2,NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX÷2)
    F = Array{Float64,2}(undef, LDF,NMAX÷2)
    C2 = Array{Float64,2}(undef, LDC2,NMAX÷2)
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
    N = parse(BlasInt, vs[4])
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end

    M = N÷2

    vs = String[]
    _isz,_jsz = (M,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       DE[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       C1[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (M,M+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       VW[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    INFO = SLICOT.mb04bd!(JOB, COMPQ1, COMPQ2, N, A, DE, C1, VW, Q1, Q2, B, F, C2, ALPHAR, ALPHAI, BETA, LIWORK, LDWORK)
    @test INFO == 0
    INFO == 0 || return

# interp output 1
    println(io, "A:")
    _nc = M
    _nr = M
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "DE:")
    show(io, "text/plain", DE[1:M,2:M+1])
    println(io)

# interp output 3
    println(io, "B:")
    _nc = M
    _nr = M
    show(io, "text/plain", B[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "F:")
    _nc = M
    _nr = M
    show(io, "text/plain", UpperTriangular(F[1:_nr,1:_nc]))
    println(io)

# interp output 5
    println(io, "C1:")
    _nc = M
    _nr = M
    show(io, "text/plain", C1[1:_nr,1:_nc])
    println(io)

# interp output 6
    println(io, "C2:")
    _nc = M
    _nr = M
    show(io, "text/plain", C2[1:_nr,1:_nc])
    println(io)

# interp output 7
    println(io, "VW:")
    show(io, "text/plain", VW[1:M,2:M+1])
    println(io)

# interp output 8
    println(io, "ALPHAR:")
    _nr = M
    show(io, "text/plain", ALPHAR[1:_nr])
    println(io)

# interp output 9
    println(io, "ALPHAI:")
    _nr = M
    show(io, "text/plain", ALPHAI[1:_nr])
    println(io)

# interp output 10
    println(io, "BETA:")
    _nr = M
    show(io, "text/plain", BETA[1:_nr])
    println(io)

if ( !LSAME( COMPQ1, 'N' ) )
# interp output 11
    println(io, "Q1:")
    _nc = N
    _nr = N
    show(io, "text/plain", Q1[1:_nr,1:_nc])
    println(io)

end # if
if ( !LSAME( COMPQ2, 'N' ) )
# interp output 12
    println(io, "Q2:")
    _nc = N
    _nr = N
    show(io, "text/plain", Q2[1:_nr,1:_nc])
    println(io)
end
end # run_mb04bd()
