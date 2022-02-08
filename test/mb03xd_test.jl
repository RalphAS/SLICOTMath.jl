# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03xd(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    NMAX = 100
    LDA = NMAX
    LDQG = NMAX
    LDRES = NMAX
    LDT = NMAX
    LDU1 = NMAX
    LDU2 = NMAX
    LDV1 = NMAX
    LDV2 = NMAX
    LDWORK = 3*NMAX*NMAX + 7*NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    T = Array{Float64,2}(undef, LDT,NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX)
    U2 = Array{Float64,2}(undef, LDU2,NMAX)
    V1 = Array{Float64,2}(undef, LDV1,NMAX)
    V2 = Array{Float64,2}(undef, LDV2,NMAX)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
    SCALE = Array{Float64,1}(undef, NMAX)
    RES = Array{Float64,2}(undef, LDRES,3*NMAX+1)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    BALANC =  vs[2][1]
    JOB =  vs[3][1]
    JOBU =  vs[4][1]
    JOBV =  vs[5][1]
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

    # FOREIGN.dlacpy!( 'All', N, N, A, LDA, RES(1,N+1), LDRES )
    RES[1:N,N+1:2*N] .= A[1:N,1:N]

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
# interp call 2

    #FOREIGN.dlacpy!( 'All', N, N+1, QG, LDQG, RES(1,2*N+1), LDRES )
    RES[1:N,2*N+1:3*N+1] .= QG[1:N,1:N+1]

# interp call 3

    ILO, INFO = SLICOT.mb03xd!(BALANC, JOB, JOBU, JOBV, N, A, QG, T, U1, U2, V1, V2, WR, WI, SCALE)
    println(io, "ILO = $ILO")
    @test INFO == 0
    INFO == 0 || return

    if ( LSAME( JOB, 'S' )||LSAME( JOB, 'G' ) )
# interp output 1
        println(io, "A:")
        _nc = N
        _nr = N
        show(io, "text/plain", A[1:_nr,1:_nc])
        println(io)

# interp output 2
        println(io, "T:")
        _nc = N
        _nr = N
        show(io, "text/plain", T[1:_nr,1:_nc])
        println(io)

    end # if
    if ( LSAME( JOB, 'G' ) )
        println(io, "QG:")
        show(io, "text/plain", QG[1:N,2:N+1])
        println(io)

    end # if
    if ( LSAME( JOB, 'G' )&&LSAME( JOBU, 'U' )&& LSAME( JOBV, 'V' ) )
# interp call 4

        ILO, INFO = SLICOT.mb04dd!(BALANC, N, view(RES,1:N,N+1:2*N),
                                   view(RES,1:N,2*N+1:3*N), DWORK)

# interp call 5

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, RES(1,N+1), LDRES, V1, LDV1, ZERO, RES, LDRES )
        RES[1:N,1:N] .= RES[1:N,N+1:2*N] * V1[1:N,1:N]

# interp call 6

        #FOREIGN.dsymm!( 'Left', 'Upper', N, N, -ONE, RES(1,2*N+2), LDRES, V2, LDV2, ONE, RES, LDRES )
        BLAS.symm!('L','U',-ONE,view(RES,1:N,2*N+2:3*N+1),view(V2,1:N,1:N),ONE,view(RES,1:N,1:N))

# interp call 7

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, -ONE, U1, LDU1, T, LDT, ONE, RES, LDRES )
        BLAS.gemm!('N','N',-ONE,U1[1:N,1:N],T[1:N,1:N],ONE,view(RES,1:N,1:N))

        # TEMP = DLANGE( 'Frobenius', N, N, RES, LDRES, DWORK )
        TEMP = norm(RES[1:N,1:N])
# interp call 8

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, RES(1,N+1), LDRES, V2, LDV2, ZERO, RES, LDRES )
        BLAS.gemm!('N','N',ONE,RES[1:N,N+1:2*N],V2[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 9

        #FOREIGN.dsymm!( 'Left', 'Upper', N, N, ONE, RES(1,2*N+2), LDRES, V1, LDV1, ONE, RES, LDRES )
        BLAS.symm!('L','U',ONE,RES[1:N,2*N+2:3*N+1],V1[1:N,1:N],ONE,view(RES,1:N,1:N))

# interp call 10

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, -ONE, U1, LDU1, QG(1,2), LDQG, ONE, RES, LDRES )
        BLAS.gemm!('N','N',-ONE,U1[1:N,1:N],QG[1:N,2:N+1],ONE,view(RES,1:N,1:N))

# interp call 11

        # FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, U2, LDU2, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','T',-ONE,U2[1:N,1:N],A[1:N,1:N],ONE,view(RES,1:N,1:N))

        # TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP, norm(RES[1:N,1:N]))
# interp call 12

        # FOREIGN.dsymm!( 'Left', 'Lower', N, N, ONE, RES(1,2*N+1), LDRES, V1, LDV1, ZERO, RES, LDRES )
        BLAS.symm!('L','L',ONE,RES[1:N,2*N+1:3*N],V1[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 13

        # FOREIGN.dgemm!( 'Transpose', 'No Transpose', N, N, N, ONE, RES(1,N+1), LDRES, V2, LDV2, ONE, RES, LDRES )
        BLAS.gemm!('T','N',ONE,RES[1:N,N+1:2*N],V2[1:N,1:N],ONE,view(RES,1:N,1:N))

# interp call 14

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U2, LDU2, T, LDT, ONE, RES, LDRES )
        BLAS.gemm!('N','N',ONE,U2[1:N,1:N],T[1:N,1:N],ONE,view(RES,1:N,1:N))

        # TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP, norm(RES[1:N,1:N]))
# interp call 15

        # FOREIGN.dsymm!( 'Left', 'Lower', N, N, ONE, RES(1,2*N+1), LDRES, V2, LDV2, ZERO, RES, LDRES )
        BLAS.symm!('L','L',ONE,RES[1:N,2*N+1:3*N],V2[1:N,1:N],ZERO,view(RES,1:N,1:N))

# interp call 16

        # FOREIGN.dgemm!( 'Transpose', 'No Transpose', N, N, N, -ONE, RES(1,N+1), LDRES, V1, LDV1, ONE, RES, LDRES )
        BLAS.gemm!('T','N',-ONE,RES[1:N,N+1:2*N],V1[1:N,1:N],ONE,view(RES,1:N,1:N))

# interp call 17

        # FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U2, LDU2, QG(1,2), LDQG, ONE, RES, LDRES )
        BLAS.gemm!('N','N',ONE,U2[1:N,1:N],QG[1:N,2:N+1],ONE,view(RES,1:N,1:N))

# interp call 18

        # FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, U1, LDU1, A, LDA, ONE, RES, LDRES )
        BLAS.gemm!('N','T',-ONE,U1[1:N,1:N],A[1:N,1:N],ONE,view(RES,1:N,1:N))

        # TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES, LDRES, DWORK ) )
        TEMP = hypot(TEMP, norm(RES[1:N,1:N]))
        println(io, "decomposition residual norm = $TEMP")
    end # if
    if ( LSAME( JOBU, 'U' ) )
        # interp output 4
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
        resid = SLICOT.ma02jd!(false,false,N,U1,U2)
        println(io, "U residual norm = $resid")
    end # if
    if ( LSAME( JOBV, 'V' ) )
        # interp output 5
        println(io, "V1:")
        _nc = N
        _nr = N
        show(io, "text/plain", V1[1:_nr,1:_nc])
        println(io)
        println(io, "V2:")
        _nc = N
        _nr = N
        show(io, "text/plain", V2[1:_nr,1:_nc])
        println(io)

# unable to translate write statement:
#  write MA02JD( .FALSE., .FALSE., N, V1, LDV1, V2, LDV2, RES, LDRES )
        resid = SLICOT.ma02jd!(false,false,N,V1,V2)
        println(io, "V residual norm = $resid")
    end # if
    if ( LSAME( BALANC, 'S' )||LSAME( BALANC, 'B' ) )
        println(io, "SCALE:")
        show(io, "text/plain", SCALE[1:N])
        println(io)
    end # if
end # run_mb03xd()
