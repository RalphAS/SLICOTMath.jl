# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03xz(datfile, io=stdout)
    ZER = 0.0e0
    ZERO = 0.0+0.0im
    ONE = 1.0+0.0im
    NIN = 5
    NOUT = 6
    NMAX = 100
    LDA   = 2*NMAX
    LDAE  = 2*NMAX
    LDAS  = NMAX
    LDGE  = 2*NMAX
    LDQG  = 2*NMAX
    LDQGE = 2*NMAX
    LDQGS  = NMAX
    LDRES = 2*NMAX
    LDU1  = 2*NMAX
    LDU2   = 2*NMAX
    LDWORK = 20*NMAX*NMAX + 12*NMAX + 2
    LZWORK = 12*NMAX - 2
    A = Array{ComplexF64,2}(undef, LDA,2*NMAX)
    QG = Array{ComplexF64,2}(undef, LDQG,2*NMAX+1)
    U1 = Array{ComplexF64,2}(undef, LDU1,2*NMAX)
    U2 = Array{ComplexF64,2}(undef, LDU2,2*NMAX)
    WR = Array{Float64,1}(undef, 2*NMAX)
    WI = Array{Float64,1}(undef, 2*NMAX)
    SCALE = Array{Float64,1}(undef, NMAX)
    BWORK = Array{BlasBool,1}(undef, 2*NMAX)
    AS = Array{ComplexF64,2}(undef, LDAS,NMAX)
    QGS = Array{ComplexF64,2}(undef, LDQGS,NMAX+1)
    DWORK = Array{Float64,1}(undef, LDWORK)
    QGE = Array{ComplexF64,2}(undef, LDQGE,2*NMAX+1)
    GE = Array{ComplexF64,2}(undef, LDGE,2*NMAX)
    AE = Array{ComplexF64,2}(undef, LDAE,2*NMAX)
    RES = Array{ComplexF64,2}(undef, LDRES,2*NMAX)
    ZWORK = Array{ComplexF64,1}(undef, LZWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    BALANC =  vs[2][1]
    JOB =  vs[3][1]
    JOBU =  vs[4][1]
    if ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
    end

    M = 2*N
    A[1:M,1:M] .= ZERO
    QG[1:M,1:2M] .= ZERO

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       A[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end
    AS[1:N,1:N] .= A[1:N,1:N]

    vs = String[]
    _isz,_jsz = (N,N+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       QG[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    #FOREIGN.zlacpy!( 'All', N, N+1, QG, LDQG, QGS, LDQGS )
    QGS[1:N,1:N+1] .= QG[1:N,1:N+1]

# interp call 2

    ILO, INFO = SLICOT.mb03xz!(BALANC, JOB, JOBU, N, A, QG, U1, U2, WR, WI, SCALE, BWORK)
    println(io, "ILO = $ILO")
    @test INFO == 0
    INFO == 0 || return

    println(io, "eigvals:")
    show(io, "text/plain", WR[1:M] + im*WI[1:M])
    println(io)
    if ( LSAME( JOB, 'S' )||LSAME( JOB, 'G' ) )
# interp output 1
    println(io, "A:")
    _nc = M
    _nr = M
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

    end # if
    if ( LSAME( JOB, 'G' ) )
# interp output 2
    println(io, "QG:")
    _nc = M
    _nr = M
    show(io, "text/plain", QG[1:_nr,1:_nc])
    println(io)

    end # if
    if ( LSAME( JOB, 'G' )&&LSAME( JOBU, 'U' ) )
# interp call 3

        ILO, INFO = SLICOT.mb04dz!(BALANC, N, AS, QGS, DWORK)
        println(io, "ILO = $ILO")
        @test INFO == 0
    INFO == 0 || return

# interp call 4

        #FOREIGN.zlaset!( 'Lower', M-1, M-1, ZERO, ZERO, A(2,1), LDA )
        triu!(view(A,1:M,1:M-1))

# interp call 5

        SLICOT.ma02ez!('U', 'C', 'N', M, QG)

        # compute Ae,Ge,Qe
        AE[1:N,1:N] .= im * imag.(AS[1:N,1:N])
        AE[N+1:2N,1:N] .= -im * real.(AS[1:N,1:N])
        for j in 1:N
            AE[1:N,j+N] .= -AE[N+1:2N,j]
            AE[N+1:2N,j+N] .= AE[1:N,j]
        end
        QGE[1:N,1:N+1] .= im * imag.(QGS[1:N,1:N+1])
        for j in 1:N+1
            QGE[N+j:2N,j] .= -im * real.(QGS[j:N,j])
        end
        for j in 1:N
            QGE[1:j,j+N+1] .= -QGE[1:j,j+1]
            QGE[N+1:2N,j+N+1] .= QGE[1:N,j+1]
        end
        QGE[N+1:2N,N+1] .= QGE[1:N,1]

# interp call 7

        SLICOT.ma02ez!('L', 'T', 'N', N, view(QGE,N+1:2N,1:N+1))

# interp call 8

        SLICOT.ma02ez!('U', 'T', 'N', N, view(QGE,1:N,N+2:2N+2))

# interp call 9

        # FOREIGN.zlacpy!( 'Upper', M, M, QGE(1,2), LDQGE, GE, LDGE )
        GE[1:M,1:M] .= triu(QGE[1:M,2:M+1])

# interp call 10

        SLICOT.ma02ez!('U', 'T', 'S', M, GE)

# interp call 11

        SLICOT.ma02ez!('L', 'T', 'S', M, QGE)

# interp call 12

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, AE, LDAE, U1, LDU1, ZERO, RES, LDRES )
        BLAS.gemm!('N','N',true,AE[1:M,1:M],U1[1:M,1:M],false,view(RES,1:M,1:M))

# interp call 13

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, -ONE, GE, LDGE, U2, LDU2, ONE, RES, LDRES )
        BLAS.gemm!('N','N',-ONE,GE[1:M,1:M],U2[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 14

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, -ONE, U1, LDU1, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','N',-ONE,U1[1:M,1:M],A[1:M,1:M],true,view(RES,1:M,1:M))

        # TEMP = ZLANGE( 'Frobenius', M, M, RES, LDRES, DWORK )
        TEMP = norm(RES[1:M,1:M])
# interp call 15

        #FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, AE, LDAE, U2, LDU2, ZERO, RES, LDRES )
        BLAS.gemm!('N','N',true,AE[1:M,1:M],U2[1:M,1:M],false,view(RES,1:M,1:M))

# interp call 16

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, GE, LDGE, U1, LDU1, ONE, RES, LDRES )
        BLAS.gemm!('N','N',true,GE[1:M,1:M],U1[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 17

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, -ONE, U1, LDU1, QG, LDQG, ONE, RES, LDRES )
        BLAS.gemm!('N','N',-ONE,U1[1:M,1:M],QG[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 18

        # FOREIGN.zgemm!( 'No Transpose', 'Conj Transpose', M, M, M, ONE, U2, LDU2, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','C',true,U2[1:M,1:M],A[1:M,1:M],true,view(RES,1:M,1:M))

        # TEMP = DLAPY2( TEMP, ZLANGE( 'Frobenius', M, M, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP,norm(RES[1:M,1:M]))
# interp call 19

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, true, QGE, LDQGE, U1, LDU1, ZERO, RES, LDRES )
        BLAS.gemm!('N','N',true,QGE[1:M,1:M],U1[1:M,1:M],false,view(RES,1:M,1:M))

# interp call 20

        # FOREIGN.zgemm!( 'Transpose', 'No Transpose', M, M, M, -ONE, AE, LDAE, U2, LDU2, ONE, RES, LDRES )
        BLAS.gemm!('T','N',-ONE,AE[1:M,1:M],U2[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 21

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, U2, LDU2, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','N',true,U2[1:M,1:M],A[1:M,1:M],true,view(RES,1:M,1:M))

        # TEMP = DLAPY2( TEMP, ZLANGE( 'Frobenius', M, M, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP,norm(RES[1:M,1:M]))
# interp call 22

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, QGE, LDQGE, U2, LDU2, ZERO, RES, LDRES )
        BLAS.gemm!('N','N',true,QGE[1:M,1:M],U2[1:M,1:M],false,view(RES,1:M,1:M))

# interp call 23

        # FOREIGN.zgemm!( 'Transpose', 'No Transpose', M, M, M, ONE, AE, LDAE, U1, LDU1, ONE, RES, LDRES )
        BLAS.gemm!('T','N',true,AE[1:M,1:M],U1[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 24

        # FOREIGN.zgemm!( 'No Transpose', 'No Transpose', M, M, M, ONE, U2, LDU2, QG, LDQG, ONE, RES, LDRES )
        BLAS.gemm!('N','N',true,U2[1:M,1:M],QG[1:M,1:M],true,view(RES,1:M,1:M))

# interp call 25

        # FOREIGN.zgemm!( 'No Transpose', 'Conj Transpose', M, M, M, ONE, U1, LDU1, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','C',true,U1[1:M,1:M],A[1:M,1:M],true,view(RES,1:M,1:M))

        #TEMP = DLAPY2( TEMP, ZLANGE( 'Frobenius', M, M, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP,norm(RES[1:M,1:M]))
        println(io, "residual norm = $TEMP")
    end # if
    if ( !LSAME( JOB, 'E' )&&LSAME( JOBU, 'U' ) )
        # interp output 3
        println(io, "U1:")
        _nc = M
        _nr = M
        show(io, "text/plain", U1[1:_nr,1:_nc])
        println(io)
        println(io, "U2:")
        _nc = M
        _nr = M
        show(io, "text/plain", U1[2:_nr,1:_nc])
        println(io)

        # unable to translate write statement:
        #  write MA02JZ( .FALSE., .FALSE., M, U1, LDU1, U2, LDU2, RES, LDRES )
        resid = SLICOT.ma02jz!(false,false,M,U1,U2)
        println(io, "U residual = $resid")
    end # if
    if ( LSAME( BALANC, 'S' )||LSAME( BALANC, 'B' ) )
        println(io, "SCALE:")
        show(io, "text/plain", SCALE[1:N])
        println(io)
    end # if
end # run_mb03xz()
