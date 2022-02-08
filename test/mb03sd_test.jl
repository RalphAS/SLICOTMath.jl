# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03sd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDQG = NMAX
    LDWORK = NMAX*( NMAX+1 )
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
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

    INFO = SLICOT.mb03sd!(JOBSCL, N, A, QG, WR, WI, LDWORK)
    @test INFO == 0
    INFO != 0 && return

    println(io, "eigvals:")
    show(io, "text/plain", vcat(WR[1:N] + im*WI[1:N],
                 -WR[N:-1:1] - im*WI[N:-1:1]))
    println(io)

end # run_mb03sd()
