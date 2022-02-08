# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04vd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    LDA  = MMAX
    LDE = MMAX
    LDQ = MMAX
    LDZ = NMAX
    LINUK = max( NMAX,MMAX+1 )
    LIWORK = NMAX
    LDWORK = max( NMAX,MMAX )
    ZERO = 0.0e0
    ONE = 1.0e0
    A = Array{Float64,2}(undef, LDA,NMAX)
    E = Array{Float64,2}(undef, LDE,NMAX)
    Q = Array{Float64,2}(undef, LDQ,MMAX)
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    ISTAIR = Array{BlasInt,1}(undef, MMAX)
# WARNING: desperate attempt to initialize NBLCKS
    NBLCKS = 0
# WARNING: desperate attempt to initialize NBLCKI
    NBLCKI = 0
    IMUK = Array{BlasInt,1}(undef, LINUK)
    INUK = Array{BlasInt,1}(undef, LINUK)
    IMUK0 = Array{BlasInt,1}(undef, NMAX)
    MNEI = Array{BlasInt,1}(undef, 3)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    TOL = parse(Float64, replace(vs[3],'D'=>'E'))
    MODE =  vs[4][1]
    if ( M<0 || M>MMAX )
        @error "Illegal M=$M"
    end
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end

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
    close(f)
    JOBQ = 'I'
    JOBZ = 'I'
# interp call 1

    RANKE, INFO = SLICOT.mb04ud!(JOBQ, JOBZ, M, N, A, E, Q, Z, ISTAIR, TOL)
    println(io, "RANKE = $RANKE")
    @test INFO == 0
    INFO == 0 || return

    JOBQ = 'U'
    JOBZ = 'U'

# interp call 2

    NBLCKS, NBLCKI, INFO = SLICOT.mb04vd!(MODE, JOBQ, JOBZ, M, N, RANKE, A, E, Q, Z, ISTAIR,
                                          IMUK, INUK, IMUK0, MNEI, TOL)
    @test INFO == 0
    INFO == 0 || return

# interp output 1
    println(io, "Q:")
    _nc = M
    _nr = M
    show(io, "text/plain", Q[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "E:")
    _nc = N
    _nr = M
    show(io, "text/plain", E[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "A:")
    _nc = N
    _nr = M
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "Z:")
    _nc = N
    _nr = N
    show(io, "text/plain", Z[1:_nr,1:_nc])
    println(io)

if ( ! LSAME( MODE, 'S' ) )
# interp output 5
    println(io, "IMUK:")
    _nr = NBLCKS
    show(io, "text/plain", IMUK[1:_nr])
    println(io)

# interp output 6
    println(io, "INUK:")
    _nr = NBLCKS
    show(io, "text/plain", INUK[1:_nr])
    println(io)

else
# interp output 7
    println(io, "IMUK:")
    _nr = NBLCKS
    show(io, "text/plain", IMUK[1:_nr])
    println(io)

# interp output 8
    println(io, "INUK:")
    _nr = NBLCKS
    show(io, "text/plain", INUK[1:_nr])
    println(io)

# interp output 9
    println(io, "IMUK0:")
    _nr = NBLCKI
    show(io, "text/plain", IMUK0[1:_nr])
    println(io)

# interp output 10
    println(io, "MNEI:")
    _nr = 3
    show(io, "text/plain", MNEI[1:_nr])
    println(io)
end # if
end # run_mb04vd()
