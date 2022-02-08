# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb05md(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDV = NMAX
    LDY = NMAX
    LDWORK = 4*NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    V = Array{Float64,2}(undef, LDV,NMAX)
    Y = Array{Float64,2}(undef, LDY,NMAX)
    VALR = Array{Float64,1}(undef, NMAX)
    VALI = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, NMAX)
    f = open(datfile,"r")
    readline(f)
    BALANC = 'N'
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    DELTA = parse(Float64, replace(vs[2],'D'=>'E'))
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
    close(f)
    # interp call 1

    INFO = SLICOT.mb05md!(BALANC, N, DELTA, A, V, Y, VALR, VALI, LDWORK)
    @test INFO == 0
    INFO == 0 || return

# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# unable to translate write loop:
#  write (I), VALI(I), I = 1,N 
# interp output 2
    println(io, "eigvals:")
    show(io, "text/plain", VALR[1:N]+im*VALI[1:N])
    println(io)

# interp output 3
    println(io, "V:")
    _nc = N
    _nr = N
    show(io, "text/plain", V[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "Y:")
    _nc = N
    _nr = N
    show(io, "text/plain", Y[1:_nr,1:_nc])
    println(io)

end # run_mb05md()
