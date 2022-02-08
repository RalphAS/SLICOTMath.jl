# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04dd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 100
    LDA = NMAX
    LDQG = NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    SCALE = Array{Float64,1}(undef, NMAX)
    DUMMY = Array{Float64,1}(undef, 1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    JOB =  vs[2][1]
    if ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       QG[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
    # interp call 1

    ILO, INFO = SLICOT.mb04dd!(JOB, N, A, QG, SCALE)
    @test INFO == 0
    INFO == 0 || return
    println(io, "ILO = $ILO")

    # interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

    # interp output 2
    println(io, "QG:")
    _nc = N+1
    _nr = N
    show(io, "text/plain", QG[1:_nr,1:_nc])
    println(io)

    if ( ILO>1 )
        subd = hypot(norm(tril(view(A,2:N,1:ILO-1))),
                    norm(tril(view(QG,1:N,1:ILO-1))))
        println(io, "subdiagonal norm: $subd")
    end # if
end # run_mb04dd()
