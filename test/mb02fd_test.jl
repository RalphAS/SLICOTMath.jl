# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02fd(datfile, io=stdout)
    ZERO = 0.0e0
    ONE = 1.0e0
    NIN = 5
    NOUT = 6
    ITMAX = 10
    KMAX = 20
    NMAX = 20
    LDR = NMAX*KMAX
    LDT = KMAX
    LDWORK = ( NMAX + 1 )*KMAX
    S = Array{BlasInt,1}(undef, ITMAX)
    tdim = NMAX*KMAX
    #T = Array{Float64,2}(undef, LDT,NMAX*KMAX)
    T = fill(zero(Float64), LDT,NMAX*KMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    #R = Array{Float64,2}(undef, LDR,NMAX*KMAX)
    R = fill(zero(Float64),LDR,NMAX*KMAX)
    V = Array{Float64,1}(undef, NMAX*KMAX)
    W = Array{Float64,1}(undef, NMAX*KMAX)
    Z = Array{Float64,1}(undef, NMAX*KMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    K = parse(BlasInt, vs[2])
    IT = parse(BlasInt, vs[3])
    TYPET = 'R'
    M = N*K
    local NNRM
    if ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
    elseif ( K<=0 || K>KMAX )
        @error "Illegal K=$K"
    elseif ( IT<=0 || IT>ITMAX )
        @error "Illegal IT=$IT"
    end

    vs = String[]
    _isz = IT
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    S[1:_isz] .= parsex.(BlasInt, vs)


    vs = String[]
    _isz,_jsz = (K,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
    P = 0
    POS = 1
    for SCIT in 1:IT # do 90
# interp call 1

        #Rtmp = R[POS:LDR,POS:tdim]
        INFO = SLICOT.mb02fd!(TYPET, K, N, P, S[SCIT], view(T,1:LDT,POS:tdim),
                              view(R,POS:LDR,POS:tdim), LDWORK)
        #R[POS:LDR,POS:tdim] .= Rtmp

        @test INFO == 0
        if ( INFO!=0 )
            return
        end # if
        S1 = S[SCIT] + P
        if ( S1==0 )
            LEN = N*K
            # interp call 2

            #FOREIGN.dlaset!( 'All', LEN, 1, ONE, ONE, V, 1 )
            V[1:LEN] .= ONE

            # interp call 3
            for PIT in 1:5 # do 30
                for i in 1:N
                    #FOREIGN.dgemv!( 'NoTranspose', K, LEN-(I-1)*K, ONE, T, LDT, V((I-1)*K+1), 1, ZERO, W((I-1)*K+1), 1 )
                    #W[(i-1)*K+1:i*K] .= T[1:K,1:LEN-(i-1)*K] * V[(i-1)*K+1:LEN]
                    _nc = LEN-(i-1)*K
                    BLAS.gemv!('N',ONE,T[1:K,1:_nc],V[(i-1)*K+1:LEN],ZERO,view(W,(i-1)*K+1:i*K))
                end

                # interp call 4

                for i in 1:N-1
                    #FOREIGN.dgemv!( 'Transpose', K, (N-I)*K, ONE, T(1,K+1), LDT, V((I-1)*K+1), 1, ONE, W(I*K+1), 1 )
                    #W[i*K+1:N*K] .= T[1:K,K+1:(N-i+1)*K]' * V[(i-1)*K+1:i*K]
                    BLAS.gemv!('T',ONE,T[1:K,K+1:(N-i+1)*K],V[(i-1)*K+1:i*K],ONE,view(W,i*K+1:N*K))
                end
                # interp call 5

                #FOREIGN.dcopy!( LEN, W, 1, V, 1 )
                V[1:LEN] .= W[1:LEN]
            end # for PIT
            #NNRM = DNRM2( LEN, V, 1 )
            NNRM = norm(V[1:LEN])
            # interp call 6

            #FOREIGN.dscal!( LEN, ONE/NNRM, V, 1 )
            V[1:LEN] .*= ONE / NNRM

        else
            LEN = ( N - S1 )*K
            # interp call 7

            #FOREIGN.dlaset!( 'All', LEN, 1, ONE, ONE, V, 1 )
            V[1:LEN] .= ONE

            for PIT in 1:5 # do 80
                POSR = ( S1 - 1 )*K + 1
                for i in 1:N-S1 # do 40
                    # interp call 8

                    #FOREIGN.dgemv!( 'NoTranspose', K, LEN-(I-1)*K, ONE, T(1,POSR+K), LDT, V((I-1)*K+1), 1, ZERO, W((I-1)*K+1), 1 )
                    _nc = LEN-(i-1)*K
                    #W[(i-1)*K+1:i*K] .= T[1:K,POSR+K:POSR+K+_nc-1] * V[(i-1)*K+1:LEN]
                    BLAS.gemv!('N',ONE,T[1:K,POSR+K:POSR+K+_nc-1],V[(i-1)*K+1:LEN],ZERO,
                               view(W,(i-1)*K+1:i*K))
                end

                for i in 1:N-S1 # do 50
                    # interp call 9

                    #FOREIGN.dtrmv!( 'U', 'N', 'N', K, R[POSR,POSR], LDR, V[(I-1)*K+1], 1 )
                    #V[(i-1)*K+1:i*K] .= UpperTriangular(R[POSR:POSR+K-1,POSR:POSR+K-1]) * V[(i-1)*K+1:i*K]
                    BLAS.trmv!('U','N','N',R[POSR:POSR+K-1,POSR:POSR+K-1],view(V,(i-1)*K+1:i*K))

                    # interp call 10

                    #FOREIGN.dgemv!( 'N', K, LEN-I*K, ONE, R[POSR,POSR+K], LDR, V[I*K+1], 1, ONE, V[(I-1)*K+1], 1 )
                    #V[(i-1)*K+1:i*K] .= R[POSR:POSR+K-1,POSR+K:POSR+K+LEN-i*K-1] * V[i*K+1:LEN]
                    BLAS.gemv!('N',ONE,R[POSR:POSR+K-1,POSR+K:POSR+K+LEN-i*K-1],
                               V[i*K+1:LEN],ONE,view(V,(i-1)*K+1:i*K))
                end
                # interp call 11

                #FOREIGN.dlaset!( 'All', LEN, 1, ZERO, ZERO, Z, 1 )
                Z[1:LEN] .= ZERO

                for i in 1:N-S1 # do 60
                    # interp call 12

                    #FOREIGN.dgemv!( 'T', K, LEN-I*K, ONE, R[POSR,POSR+K], LDR, V[(I-1)*K+1], 1, ONE, Z[I*K+1], 1 )
                    #Z[i*K+1:LEN] .+= R[POSR:POSR+K-1,POSR:POSR+LEN-i*K-1]' * V[(i-1)*K+1:i*K]
                    BLAS.gemv!('T',ONE,R[POSR:POSR+K-1,POSR:POSR+LEN-i*K-1],V[(i-1)*K+1:i*K],
                               ONE,view(Z,i*K+1:LEN))

                    # interp call 13

                    #FOREIGN.dtrmv!( 'U', 'T', 'N', K, R[POSR,POSR], LDR, V[(I-1)*K+1], 1 )
                    #V[(i-1)*K+1:i*K] .= UpperTriangular(R[POSR:POSR+K-1,POSR:POSR+K-1])' * V[(i-1)*K+1:i*K]
                    BLAS.trmv!('U','T','N',R[POSR:POSR+K-1,POSR:POSR+K-1],
                               view(V,(i-1)*K+1:i*K))
                    # interp call 14

                    #FOREIGN.daxpy!( K, ONE, V((I-1)*K+1), 1, Z((I-1)*K+1), 1 )
                    Z[(i-1)*K+1:i*K] .+= V[(i-1)*K+1:i*K]
                end

                # interp call 15

                # FOREIGN.dlaset!( 'All', LEN, 1, ZERO, ZERO, V, 1 )
                V[1:LEN] .= ZERO

                for i in 1:N-S1 # do 70
                    # interp call 16

                    #FOREIGN.dgemv!( 'T', K, LEN-(I-1)*K, ONE, T(1,POSR+K), LDT, W((I-1)*K+1), 1, ONE, V((I-1)*K+1), 1 )
                    #V[(i-1)*K+1:LEN] .+= T[1:K,POSR+K:POSR+K+LEN-(i-1)*K-1]' * W[(i-1)*K+1:i*K]
                    BLAS.gemv!('T',ONE,T[1:K,POSR+K:POSR+K+LEN-(i-1)*K-1],W[(i-1)*K+1:i*K],
                               ONE,view(V,(i-1)*K+1:LEN))
                end
                # interp call 17

                #FOREIGN.daxpy!( LEN, -ONE, Z, 1, V, 1 )
                #V[1:LEN] .-= Z[1:LEN]
                BLAS.axpy!(-ONE,Z[1:LEN],view(V,1:LEN))

                # NNRM = DNRM2( LEN, V, 1 )
                NNRM = norm(V[1:LEN])
                # interp call 18

                # FOREIGN.dscal!( LEN, -ONE/NNRM, V, 1 )
                V[1:LEN] .*= (-ONE / NNRM)
            end # for / do 80
            POS = ( S1 - 1 )*K + 1
            P = S1
        end # if
        println(io, "rows: $(P*K); norm of Schur complement: $NNRM")

    end # for / do 90
    # interp output 1
    println(io, "R:")
    _nc = M
    _nr = P*K
    show(io, "text/plain", R[1:_nr,1:_nc])
    println(io)


end # run_mb02fd()
