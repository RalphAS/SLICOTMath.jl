# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03qg(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 10
    LDA = NMAX
    LDE = NMAX
    LDU = NMAX
    LDV = NMAX
    LDWORK = 8*NMAX + 16
    A = Array{Float64,2}(undef, LDA,NMAX)
    E = Array{Float64,2}(undef, LDE,NMAX)
    U = Array{Float64,2}(undef, LDU,NMAX)
    V = Array{Float64,2}(undef, LDV,NMAX)
    BETA = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    WI = Array{Float64,1}(undef, NMAX)
    WR = Array{Float64,1}(undef, NMAX)
    BWORK = Array{Bool,1}(undef, NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    NLOW = parse(BlasInt, vs[2])
    NSUP = parse(BlasInt, vs[3])
    ALPHA = parse(Float64, replace(vs[4],'D'=>'E'))
    DICO =  vs[5][1]
    STDOM =  vs[6][1]
    JOBU =  vs[7][1]
    JOBV =  vs[8][1]
    if ( N<0 || N>NMAX )
        @error "invalid N=$N"
    else

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
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       E[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    # FOREIGN.dgges!( 'Vectors', 'Vectors', 'Not sorted', DELCTG, N, A, LDA, E, LDE, NDIM, WR, WI, BETA, U, LDU, V, LDV, DWORK, LDWORK, BWORK, INFO )
    Aschur = schur(A[1:N,1:N], E[1:N,1:N])
    WR[1:N] .= real.(Aschur.alpha)
        WI[1:N] .= imag.(Aschur.alpha)
        U[1:N,1:N] .= Aschur.Q
        V[1:N,1:N] .= Aschur.Z
        A[1:N,1:N] .= Aschur.S
        E[1:N,1:N] .= Aschur.T

    println(io, "orig T:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)
    println(io, "eigvals:")
    show(io, "text/plain", (WR[1:N] .+ im*WI[1:N]) ./ BETA[1:N])
    println(io)

# interp call 2

    NDIM, INFO = SLICOT.mb03qg!(DICO, STDOM, JOBU, JOBV, N, NLOW, NSUP, ALPHA, A, E, U, V)
    println(io, "NDIM = $NDIM")
    @test INFO == 0

        if ( INFO!=0 )
            println(io, "INFO = $INFO")
            return
        end
# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "E:")
    _nc = N
    _nr = N
    show(io, "text/plain", E[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "U:")
    _nc = N
    _nr = N
    show(io, "text/plain", U[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "V:")
    _nc = N
    _nr = N
    show(io, "text/plain", V[1:_nr,1:_nc])
    println(io)


end # if

    close(f)
end # run_mb03qg()
