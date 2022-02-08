# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03xp(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    NMAX = 200
    LDA = NMAX
    LDB = NMAX
    LDQ = NMAX
    LDRES = NMAX
    LDWORK = NMAX
    LDZ = NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    B = Array{Float64,2}(undef, LDA,NMAX)
    Q = Array{Float64,2}(undef, LDQ,NMAX)
    Z = Array{Float64,2}(undef, LDZ,NMAX)
    ALPHAR = Array{Float64,1}(undef, NMAX)
    ALPHAI = Array{Float64,1}(undef, NMAX)
    BETA = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RES = Array{Float64,2}(undef, LDRES,3*NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    ILO = parse(BlasInt, vs[2])
    IHI = parse(BlasInt, vs[3])
    if ( N<=0 || N>NMAX )
        @error "illegal N=$N"
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
# interp call 1

    #FOREIGN.dlacpy!( 'All', N, N, A, LDA, RES(1,N+1), LDRES )
    RES[1:N,N+1:2N] .= A[1:N,1:N]


    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 2

    #FOREIGN.dlacpy!( 'All', N, N, B, LDB, RES(1,2*N+1), LDRES )
    RES[1:N,2*N+1:3N] .= B[1:N,1:N]

# interp call 3

    INFO = SLICOT.mb03xp!('S', 'I', 'I', N, ILO, IHI, A, B, Q, Z, ALPHAR, ALPHAI, BETA, LDWORK)
    @test INFO == 0
    INFO == 0 || return

# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp call 4

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, RES(1,N+1), LDRES, Z, LDZ, ZERO, RES, LDRES )
    BLAS.gemm!('N','N',ONE,RES[1:N,N+1:2N],Z[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 5

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, -ONE, Q, LDQ, A, LDA, ONE, RES, LDRES )
    BLAS.gemm!('N','N',-ONE,Q[1:N,1:N],A[1:N,1:N],ONE,view(RES,1:N,1:N))
    resid = norm(RES[1:N,1:N])
    println(io, "Frobenius norm of A Z - Q S: $resid")
    @test resid < 1e-12
    println(io, "B:")
    _nc = N
    _nr = N
    show(io, "text/plain", B[1:_nr,1:_nc])
    println(io)

# interp call 6

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, RES(1,2*N+1), LDRES, Q, LDQ, ZERO, RES, LDRES )
    BLAS.gemm!('N','N',ONE,RES[1:N,2*N+1:3N],Q[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 7

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, -ONE, Z, LDZ, B, LDB, ONE, RES, LDRES )
    BLAS.gemm!('N','N',-ONE,Z[1:N,1:N],B[1:N,1:N],ONE,view(RES,1:N,1:N))
    resid = norm(RES[1:N,1:N])
    println(io, "Frobenius norm of B Q - Z T:  $resid")
    @test resid < 1e-12
    println(io, "Q:")
    _nc = N
    _nr = N
    show(io, "text/plain", Q[1:_nr,1:_nc])
    println(io)

# interp call 8

    #FOREIGN.dgemm!( 'Transpose', 'No Transpose', N, N, N, ONE, Q, LDQ, Q, LDQ, ONE, RES, LDRES )
    BLAS.gemm!('T','N',ONE,Q[1:N,1:N],Q[1:N,1:N],ZERO,view(RES,1:N,1:N))
    resid =  norm(RES[1:N,1:N] - I)
    println(io, "orth. residual norm of Q: $resid")# interp output 4
    @test resid < 1e-12
    println(io, "Z:")
    _nc = N
    _nr = N
    show(io, "text/plain", Z[1:_nr,1:_nc])
    println(io)

# interp call 9

    #FOREIGN.dgemm!( 'Transpose', 'No Transpose', N, N, N, ONE, Z, LDZ, Z, LDZ, ONE, RES, LDRES )
    BLAS.gemm!('T','N',ONE,Z[1:N,1:N],Z[1:N,1:N],ZERO,view(RES,1:N,1:N))
    resid = norm(RES[1:N,1:N] - I)
    println(io, "orth. residual norm of Z: $resid")# unable to translate write statement:
    @test resid < 1e-12
#   in do block [('70', 'I', 'N')]
    #   write ALPHAR(I), ALPHAI(I), BETA(I)
    println(io, "ALPHA:")
    show(io, "text/plain", ALPHAR[1:N]+im*ALPHAI[1:N])
    println(io)
    println(io, "BETA:")
    show(io, "text/plain", BETA[1:N])
    println(io)

end # run_mb03xp()
