# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04dy(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDQG = NMAX
    LDWORK = NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    D = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    JOBSCL =  vs[2][1]
    if ( N<0 || N>NMAX )
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
    _isz,_jsz = (N,N)
    while length(vs) < _isz*(_jsz+1)÷2
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    _i0 = 0
    for j in 1:_jsz
       QG[j,j+1:_jsz+1] .= parsex.(Float64, vs[_i0+1:_i0+_jsz-j+1])
       _i0 += _jsz-j+1
    end
    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*(_jsz+1)÷2
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    _i0 = 0
    for j in 1:_jsz
       QG[j:_jsz,j] .= parsex.(Float64, vs[_i0+1:_i0+_jsz-j+1])
       _i0 += _jsz-j+1
    end
    close(f)

    # interp call 1

    INFO = SLICOT.mb04dy!(JOBSCL, N, A, QG, D)
    @test INFO == 0
    INFO != 0 && return

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

    if ( LSAME( JOBSCL, 'S' ) )
        # interp output 3
        println(io, "D:")
        _nr = N
        show(io, "text/plain", D[1:_nr])
        println(io)

    elseif ( LSAME( JOBSCL, '1' ) || LSAME( JOBSCL, 'O' ) )
        println(io, "tau = ",D[1])
    end # if

end # run_mb04dy()
