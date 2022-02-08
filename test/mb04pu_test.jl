# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04pu(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    TWO = 2.0e0
    NIN = 5
    NOUT = 6
    NMAX = 100
    LDA  = NMAX
    LDQG = NMAX
    LDRES  = NMAX
    LDU1 = NMAX
    LDU2 = NMAX
    LDWORK = 2*NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    CS = Array{Float64,1}(undef, 2*NMAX)
    TAU = Array{Float64,1}(undef, NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX)
    U2 = Array{Float64,2}(undef, LDU2,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RES = Array{Float64,2}(undef, LDRES,3*NMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    if ( N<=0 || N>NMAX )
        @error "invalid N=$N"
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
    _isz,_jsz = (N,N+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       QG[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 2

    #FOREIGN.dlacpy!( 'All', N, N+1, QG, LDQG, RES(1,2*N+1), LDRES )
    RES[1:N,2*N+1:3N+1] .= QG[1:N,1:N+1]

    close(f)
# interp call 3

    INFO = SLICOT.mb04pu!(N, 1, A, QG, CS, TAU, LDWORK)
    @test INFO == 0
    INFO == 0 || return

    # yes, this is from the original
    # INFO = 0
    # if ( INFO!=0 )
    # else
# interp call 4

    #FOREIGN.dlacpy!( 'Lower', N, N, A, LDA, U1, LDU1 )
    U1[1:N,1:N] .= tril(A[1:N,1:N])

# interp call 5

    #FOREIGN.dlacpy!( 'Lower', N, N, QG, LDQG, U2, LDU2 )
    U2[1:N,1:N] .= tril(QG[1:N,1:N])

# interp call 6

    INFO = SLICOT.mb04wp!(N, 1, U1, U2, CS, TAU)
    @test INFO == 0
    INFO == 0 || return

    if  N>2
# interp call 7

        # FOREIGN.dlaset!( 'Lower', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
        triu!(view(A,2:N,1:N-2))

    end
    if  N>1
# interp call 8

        # FOREIGN.dlaset!( 'Lower', N-1, N-1, ZERO, ZERO, QG(2,1), LDQG )
        triu!(view(QG,1:N,1:N-1))

    end
# interp output 1
    println(io, "U1:")
    _nc = N
    _nr = N
    show(io, "text/plain", U1[1:_nr,1:_nc])
    println(io)
    println(io, "U2:")
    _nc = N
    _nr = N
    show(io, "text/plain", U2[1:_nr,1:_nc])
    println(io)

# unable to translate write statement:
#  write MA02JD( .FALSE., .FALSE., N, U1, LDU1, U2, LDU2, RES, LDRES )
    resnorm = SLICOT.ma02jd!(false,false,N,U1,U2)
    println(io, "U orth. residual norm = $resnorm")
    @test resnorm < 1e-12
# interp output 2
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "QG:")
    _nc = N+1
    _nr = N
    show(io, "text/plain", QG[1:_nr,1:_nc])
    println(io)

# interp call 9

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U1, LDU1, A, LDA, ZERO, RES, LDRES )
    BLAS.gemm!('N','N',ONE,U1[1:N,1:N],A[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 10

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, RES, LDRES, U1, LDU1, ONE, RES(1,N+1), LDRES )
    BLAS.gemm!('N','T',-ONE,RES[1:N,1:N],U1[1:N,1:N],ONE,view(RES,1:N,N+1:2N))

# interp call 11

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, ONE, U2, LDU2, A, LDA, ZERO, RES, LDRES )
    BLAS.gemm!('N','T',ONE,U2[1:N,1:N],A[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 12

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, ONE, RES, LDRES, U2, LDU2, ONE, RES(1,N+1), LDRES )
    BLAS.gemm!('N','T',ONE,RES[1:N,1:N],U2[1:N,1:N],ONE,view(RES,1:N,N+1:2N))

# interp call 13

    #FOREIGN.dsymm!( 'Right', 'Upper', N, N, ONE, QG(1,2), LDQG, U1, LDU1, ZERO, RES, LDRES )
    BLAS.symm!('R','U',ONE,QG[1:N,2:N+1],U1[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 14

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, RES, LDRES, U2, LDU2, ONE, RES(1,N+1), LDRES )
    BLAS.gemm!('N','T',-ONE,RES[1:N,1:N],U2[1:N,1:N],ONE,view(RES,1:N,N+1:2N))

    RES[1:N,1:N] .= U2[1:N,1:N]
# interp call 15

    for i in 1:N
        # FOREIGN.dscal!( N, QG(I,I), RES(1,I), 1 )
        RES[1:N,i] .*= QG[i,i]
    end

# interp call 16

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, RES, LDRES, U1, LDU1, ONE, RES(1,N+1), LDRES )
    BLAS.gemm!('N','T',-ONE,RES[1:N,1:N],U1[1:N,1:N],ONE,view(RES,1:N,N+1:2N))

# interp call 17

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U2, LDU2, A, LDA, ZERO, RES, LDRES )
    BLAS.gemm!('N','N',ONE,U2[1:N,1:N],A[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 18

    #FOREIGN.dsyr2k!( 'Lower', 'No Transpose', N, N, ONE, RES, LDRES, U1, LDU1, ONE, RES(1,2*N+1), LDRES )
    BLAS.syr2k!('L','N',ONE,RES[1:N,1:N],U1[1:N,1:N],ONE,view(RES,1:N,2*N+1:3N))

# interp call 19
    #FOREIGN.dscal!( N, ONE/TWO, QG(1,2), LDQG+1 )
    for i in 1:N
        QG[i,i+1] *= 1/2
    end

    RES[1:N,1:N] .= U2[1:N,1:N]
# interp call 20

    #FOREIGN.dtrmm!(  'Right', 'Upper' , 'No Transpose', 'Not unit', N, N, ONE, QG(1,2), LDQG, RES, LDRES )
    BLAS.trmm!('R','U','N','N',ONE,QG[1:N,2:N+1],view(RES,1:N,1:N))

# interp call 21

    #FOREIGN.dsyr2k!( 'Lower', 'No Transpose', N, N, ONE, RES, LDRES, U2, LDU2, ONE, RES(1,2*N+1), LDRES )
    BLAS.syr2k!('L','N',ONE,RES[1:N,1:N],U2[1:N,1:N],ONE,view(RES,1:N,2N+1:3N))

# interp call 22

    for i in 1:N
    #FOREIGN.dsyr!( 'Lower', N, -QG(I,I), U1(1,I), 1, RES(1,2*N+1), LDRES )
        BLAS.syr!('L',-QG[i,i],U1[1:N,i],view(RES,1:N,2*N+1:3N))
    end

# interp call 23

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U1, LDU1, A, LDA, ZERO, RES, LDRES )
    BLAS.gemm!('N','N',ONE,U1[1:N,1:N],A[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 24

    #FOREIGN.dsyr2k!( 'Upper', 'No Transpose', N, N, ONE, RES, LDRES, U2, LDU2, ONE, RES(1,2*N+2), LDRES )
    BLAS.syr2k!('U','N',ONE,RES[1:N,1:N],U2[1:N,1:N],ONE,view(RES,1:N,2*N+2:3*N+1))

    RES[1:N,1:N] .= U1[1:N,1:N]
# interp call 25

    #FOREIGN.dtrmm!(  'Right', 'Upper' , 'No Transpose', 'Not unit', N, N, ONE, QG(1,2), LDQG, RES, LDRES )
    BLAS.trmm!('R','U','N','N',ONE,QG[1:N,2:N+1],view(RES,1:N,1:N))

# interp call 26

    #FOREIGN.dsyr2k!( 'Upper', 'No Transpose', N, N, -ONE, RES, LDRES, U1, LDU1, ONE, RES(1,2*N+2), LDRES )
    BLAS.syr2k!('U','N',-ONE,RES[1:N,1:N],U1[1:N,1:N],ONE,view(RES,1:N,2*N+2:3*N+1))

# interp call 27
    for i in 1:N
        #FOREIGN.dsyr!( 'Upper', N, QG(I,I), U2(1,I), 1, RES(1,2*N+2), LDRES )
        BLAS.syr!('U',QG[i,i],U2[1:N,i],view(RES,1:N,2*N+2:3*N+1))
    end

# unable to translate write statement:
#  write MA02ID( 'Hamiltonian', 'Frobenius', N, RES(1,N+1), LDRES, RES(1,2*N+1), LDRES, DWORK )
    resnorm = SLICOT.ma02id!('H','F',N,view(RES,1:N,N+1:2*N),view(RES,1:N,2*N+1:3*N),DWORK)
    println(io, "residual norm = $resnorm")
    @test resnorm < 1e-12
end # run_mb04pu()
