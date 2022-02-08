# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04gd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    NMAX = 10
    MMAX = 10
    LDA = MMAX
    LDTAU = min(MMAX,NMAX)
    LDWORK = 3*MMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    JPVT = Array{BlasInt,1}(undef, MMAX)
    TAU = Array{Float64,1}(undef, LDTAU)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
    end
    if ( M<0 || M>MMAX )
        @error "Illegal M=$M"
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
    _isz = M
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    JPVT[1:_isz] .= parsex.(BlasInt, vs)
    close(f)

# interp call 1

    INFO = SLICOT.mb04gd!(M, N, A, JPVT, TAU)
    @test INFO == 0
    INFO == 0 || return

    if ( INFO!=0 )
        @warn "mb04gd returns info=$INFO"
        return
    end

# interp output 1
    println(io, "JPVT:")
    _nr = M
    show(io, "text/plain", JPVT[1:_nr])
    println(io)

if ( M>=N )
    if  N>1 
# interp call 2

        #FOREIGN.dlaset!( 'Lower', N-1, N-1, ZERO, ZERO, A(M-N+2,1), LDA )
        triu!(view(A,M-N+1:M,1:N-1))

    end
else
# interp call 3

    #FOREIGN.dlaset!( 'Full', M, N-M-1, ZERO, ZERO, A, LDA )
    A[1:M,1:N-M-1] .= ZERO

# interp call 4

    #FOREIGN.dlaset!( 'Lower', M, M, ZERO, ZERO, A(1,N-M), LDA )
    triu!(view(A,1:M,N-M:N-1),1)

end # if
# interp output 2
    println(io, "A:")
    _nc = N
    _nr = M
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

end # run_mb04gd()
