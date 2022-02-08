# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04tb(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    NBMAX = 64
    NMAX = 421
    LDA = NMAX
    LDB = NMAX
    LDG = NMAX
    LDQ = NMAX
    LDRES = NMAX
    LDU1 = NMAX
    LDU2 = NMAX
    LDV1 = NMAX
    LDV2 = NMAX
    LDWORK = NBMAX*( 16*NMAX + 1 )
    A = Array{Float64,2}(undef, LDA,NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX)
    G = Array{Float64,2}(undef, LDG,NMAX)
    Q = Array{Float64,2}(undef, LDQ,NMAX)
    CSL = Array{Float64,1}(undef, 2*NMAX)
    CSR = Array{Float64,1}(undef, 2*NMAX)
    TAUL = Array{Float64,1}(undef, NMAX)
    TAUR = Array{Float64,1}(undef, NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX)
    U2 = Array{Float64,2}(undef, LDU2,NMAX)
    V1 = Array{Float64,2}(undef, LDV1,NMAX)
    V2 = Array{Float64,2}(undef, LDV2,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RES = Array{Float64,2}(undef, LDRES,5*NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    TRANA =  vs[2][1]
    TRANB =  vs[3][1]
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
    RES[1:N,1:N] .= A[1:N,1:N]

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    #FOREIGN.dlacpy!( 'All', N, N, B, LDB, RES(1,N+1), LDRES )
    RES[1:N,N+1:2N] .= B[1:N,1:N]


    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       G[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 2

    #FOREIGN.dlacpy!( 'All', N, N, G, LDG, RES(1,2*N+1), LDRES )
    RES[1:N,2*N+1:3N] .= G[1:N,1:N]


    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       Q[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 3

    #FOREIGN.dlacpy!( 'All', N, N, Q, LDQ, RES(1,3*N+1), LDRES )
    RES[1:N,3*N+1:4N] .= Q[1:N,1:N]
    close(f)

# interp call 4

    INFO = SLICOT.mb04tb!(TRANA, TRANB, N, 1, A, B, G, Q, CSL, CSR, TAUL, TAUR)
    @test INFO == 0
    INFO == 0 || return

    U1[1:N,1:N] .= A[1:N,1:N]
    U2[1:N,1:N] .= Q[1:N,1:N]
# interp call 5

    INFO = SLICOT.mb04wr!('U', TRANA, N, 1, U1, U2, CSL, TAUL)
    @test INFO == 0
    INFO == 0 || return

    V2[1:N,1:N] .= Q[1:N,1:N]
    V1[1:N,1:N] .= B[1:N,1:N]
# interp call 6

    INFO = SLICOT.mb04wr!('V', TRANB, N, 1, V1, V2, CSR, TAUR)
    @test INFO == 0
    INFO == 0 || return

    if ( LSAME( TRANA, 'N' ) )
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
        #  write MA02JD( .FALSE., .FALSE., N, U1, LDU1, U2, LDU2, RES(1,4*N+1), LDRES )
        resnorm = SLICOT.ma02jd!(false,false,N,U1,U2)
        println(io, "U orth. residual norm = $resnorm")
        @test resnorm < 1e-12
    else
        # interp output 2
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
        resnorm = SLICOT.ma02jd!(true,false,N,U1,U2)
        println(io, "U orth. residual norm = $resnorm")
        @test resnorm < 1e-12

        # unable to translate write statement:
        #  write MA02JD( .TRUE., .FALSE., N, U1, LDU1, U2, LDU2, RES(1,4*N+1), LDRES )
    end # if
# interp call 7

    #FOREIGN.dlaset!( 'All', N, N, ZERO, ZERO, Q, LDQ )
    Q[1:N,1:N] .= ZERO

    if ( LSAME( TRANA, 'N' ) )
# interp call 8

        #FOREIGN.dlaset!( 'Lower', N-1, N-1, ZERO, ZERO, A(2,1), LDA )
        triu!(view(A,1:N,1:N-1))

        # interp output 3
        println(io, "A:")
        _nc = N
        _nr = N
        show(io, "text/plain", A[1:_nr,1:_nc])
        println(io)
        println(io, "G:")
        show(io, "text/plain", G[1:_nr,1:_nc])
        println(io)

    else
# interp call 9

        #FOREIGN.dlaset!( 'Upper', N-1, N-1, ZERO, ZERO, A(1,2), LDA )
        tril!(view(A,1:N-1,1:N))

        # interp output 4
        println(io, "A:")
        _nc = N
        _nr = N
        show(io, "text/plain", A[1:_nr,1:_nc])
        println(io)
        println(io, "G:")
        show(io, "text/plain", G[1:_nr,1:_nc])
        println(io)

    end # if
    if ( LSAME( TRANB, 'N' ) )
        if ( N>1 )
# interp call 10

            #FOREIGN.dlaset!( 'Upper', N-2, N-2, ZERO, ZERO, B(1,3), LDB )
            tril!(view(B,1:N-1,2:N))

        end # if
# interp output 5
        _nc = N
        _nr = N
        println(io, "Q:")
        show(io, "text/plain", Q[1:_nr,1:_nc])
        println(io)
        println(io, "B:")
        show(io, "text/plain", B[1:_nr,1:_nc])
        println(io)

    else
        if ( N>1 )
            triu!(view(B,2:N-1,1:N-2))
        end # if
        _nc = N
        _nr = N
        println(io, "Q:")
        show(io, "text/plain", Q[1:_nr,1:_nc])
        println(io)
        println(io, "B:")
        show(io, "text/plain", B[1:_nr,1:_nc])
        println(io)
# interp output 6

    end # if
    if ( LSAME( TRANB, 'N' ) )
        TRANV1 = 'T'
    else
        TRANV1 = 'N'
    end # if
# interp call 11

    #FOREIGN.dgemm!( TRANA, TRANV1, N, N, N, ONE, RES, LDRES, V1, LDV1, ZERO, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANA,TRANV1,ONE,RES[1:N,1:N],V1[1:N,1:N],ZERO,view(RES,1:N,4*N+1:5N))

# interp call 12

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, -ONE, RES(1,2*N+1), LDRES, V2, LDV2, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N','T',-ONE,RES[1:N,2*N+1:3N],V2[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

# interp call 13

    #FOREIGN.dgemm!( TRANA, TRANA, N, N, N, -ONE, U1, LDU1, A, LDA, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANA,TRANA,-ONE,U1[1:N,1:N],A[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

    #TEMP = DLANGE( 'Frobenius', N, N, RES(1,4*N+1), LDRES, DWORK )
    TEMP = norm(RES[1:N,4*N+1:5N])
# interp call 14

    #FOREIGN.dgemm!( TRANA, 'Transpose', N, N, N, ONE, RES, LDRES, V2, LDV2, ZERO, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANA,'T',ONE,RES[1:N,1:N],V2[1:N,1:N],ZERO,view(RES,1:N,4*N+1:5N))

# interp call 15

    #FOREIGN.dgemm!( 'No Transpose', TRANV1, N, N, N, ONE, RES(1,2*N+1), LDRES, V1, LDV1, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N',TRANV1,ONE,RES[1:N,2*N+1:3N],V1[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))
# interp call 16

    #FOREIGN.dgemm!( TRANA, 'No Transpose', N, N, N, -ONE, U1, LDU1, G, LDG, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANA,'N',-ONE,U1[1:N,1:N],G[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

# interp call 17

    #FOREIGN.dgemm!( 'No Transpose', TRANB, N, N, N, -ONE, U2, LDU2, B, LDB, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N',TRANB,-ONE,U2[1:N,1:N],B[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

    #TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES(1,4*N+1), LDRES, DWORK ) )
    TEMP = hypot(TEMP,norm(RES[1:N,4*N+1:5N]))
# interp call 18

    #FOREIGN.dgemm!( 'No Transpose', TRANV1, N, N, N, ONE, RES(1,3*N+1), LDRES, V1, LDV1, ZERO, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N',TRANV1,ONE,RES[1:N,3*N+1:4N],V1[1:N,1:N],ZERO,view(RES,1:N,4*N+1:5N))

# interp call 19

    #FOREIGN.dgemm!( TRANB, 'Transpose', N, N, N, -ONE, RES(1,N+1), LDRES, V2, LDV2, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANB,'T',-ONE,RES[1:N,N+1:2N], V2[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

# interp call 20

    #FOREIGN.dgemm!( 'No Transpose', TRANA, N, N, N, ONE, U2, LDU2, A, LDA, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N',TRANA,ONE,U2[1:N,1:N],A[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

    #TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES(1,4*N+1), LDRES, DWORK ) )
    TEMP = hypot(TEMP,norm(RES[1:N,4*N+1:5N]))
# interp call 21

    #FOREIGN.dgemm!( 'No Transpose', 'Transpose', N, N, N, ONE, RES(1,3*N+1), LDRES, V2, LDV2, ZERO, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N','T',ONE,RES[1:N,3*N+1:4N],V2[1:N,1:N],ZERO,view(RES,1:N,4*N+1:5N))

# interp call 22

    #FOREIGN.dgemm!( TRANB, TRANV1, N, N, N, ONE, RES(1,N+1), LDRES, V1, LDV1, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANB,TRANV1,ONE,RES[1:N,N+1:2N],V1[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

# interp call 23

    #FOREIGN.dgemm!( 'No Transpose', 'No Transpose', N, N, N, ONE, U2, LDU2, G, LDG, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!('N','N',ONE,U2[1:N,1:N],G[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

# interp call 24

    #FOREIGN.dgemm!( TRANA, TRANB, N, N, N, -ONE, U1, LDU1, B, LDB, ONE, RES(1,4*N+1), LDRES )
    BLAS.gemm!(TRANA,TRANB,-ONE,U1[1:N,1:N],B[1:N,1:N],ONE,view(RES,1:N,4*N+1:5N))

    #TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES(1,4*N+1), LDRES, DWORK ) )
    TEMP = hypot(TEMP,norm(RES[1:N,4*N+1:5N]))
    println(io, "residual H V - U R norm: $TEMP")
    @test TEMP < 1e-12
    if ( LSAME( TRANB, 'N' ) )
        # interp output 7
        _nc = N
        _nr = N
        println(io, "V1:")
        show(io, "text/plain", V1[1:_nr,1:_nc])
        println(io)
        println(io, "V2:")
        show(io, "text/plain", V2[1:_nr,1:_nc])
        println(io)

        # unable to translate write statement:
        #  write MA02JD( .TRUE., .TRUE., N, V1, LDV1, V2, LDV2, RES(1,4*N+1), LDRES )
        resnorm = SLICOT.ma02jd!(true,true,N,V1,V2)
        println(io, "V orth. residual norm = $resnorm")
        @test resnorm < 1e-12
    else
# interp output 8
        _nc = N
        _nr = N
        println(io, "V1:")
        show(io, "text/plain", V1[1:_nr,1:_nc])
        println(io)
        println(io, "V2:")
        show(io, "text/plain", V1[1:_nr,1:_nc])
        println(io)

        # unable to translate write statement:
        #  write MA02JD( .FALSE., .TRUE., N, V1, LDV1, V2, LDV2, RES(1,4*N+1), LDRES )
        resnorm = SLICOT.ma02jd!(false,true,N,V1,V2)
        println(io, "V orth. residual norm = $resnorm")
        @test resnorm < 1e-12
    end # if
end # run_mb04tb()
