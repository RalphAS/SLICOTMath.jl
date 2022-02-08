# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03md(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    NMAX = 20
    Q = Array{Float64,1}(undef, NMAX)
    E = Array{Float64,1}(undef, NMAX-1)
    Q2 = Array{Float64,1}(undef, NMAX)
    E2 = Array{Float64,1}(undef, NMAX-1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    THETA = parse(Float64, replace(vs[2],'D'=>'E'))
    L = parse(BlasInt, vs[3])
    TOL = parse(Float64, replace(vs[4],'D'=>'E'))
    RELTOL = parse(Float64, replace(vs[5],'D'=>'E'))
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end
    if ( L<0 || L>N )
        @error "Illegal L=$L"
    end

    vs = String[]
    _isz = N
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    Q[1:_isz] .= parsex.(Float64, vs)

    vs = String[]
    _isz = N-1
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    E[1:_isz] .= parsex.(Float64, vs)
    close(f)

    println(io, "J:")
    show(io, "text/plain", Bidiagonal(Q[1:N],E[1:N-1],:U))
    println(io)
    Q2[N] = Q[N]^2
    PIVMIN = Q2[N]
    for i in 1:N-1
        Q2[i] = Q[i]^2
        E2[i] = E[i]^2
        PIVMIN = max( PIVMIN, Q2[i], E2[i] )
    end
    SAFMIN = floatmin(1.0)
    PIVMIN = max( PIVMIN*SAFMIN, SAFMIN )
    TOL = max( TOL, ZERO )
    if RELTOL <= 0
        RELTOL = 2.0 * eps(1.0)
    end
# interp call 1

    L, THETA, INFO, IWARN = SLICOT.mb03md!(N, L, THETA, Q, E, Q2, E2, PIVMIN, TOL, RELTOL)
    @test INFO == 0
    println(io, "L = $L")
    println(io, "THETA = $THETA")
    println(io, "IWARN = $IWARN")

end # run_mb03md()
