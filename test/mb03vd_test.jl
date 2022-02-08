# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03vd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    PMAX = 20
    LDA1 = NMAX
    LDA2 = NMAX
    LDQ1 = NMAX
    LDQ2 = NMAX
    LDTAU = NMAX-1
    LDWORK = NMAX
    ZERO = 0.0e0
    ONE = 1.0e0
    A = Array{Float64,3}(undef, LDA1,LDA2,PMAX)
    TAU = Array{Float64,2}(undef, LDTAU,PMAX)
    Q = Array{Float64,3}(undef, LDQ1,LDQ2,PMAX)
    AS = Array{Float64,3}(undef, LDA1,LDA2,PMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    QTA = Array{Float64,2}(undef, LDQ1,NMAX)
    SSQ = 0.0
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    P = parse(BlasInt, vs[2])
    ILO = parse(BlasInt, vs[3])
    IHI = parse(BlasInt, vs[4])
    if ( N<0 || N>min( LDA1, LDA2 ) )
        @error "illegal N=$N"
    end
    if ( P<=0 || P>PMAX )
        @error "illegal P=$P"
    end

    _isz,_jsz,_ksz = (N,N,P)
    for k in 1:_ksz
      vs = String[]
      while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
      end
      for i in 1:_isz
        _i0 = (i-1)*_jsz
        A[i,1:_jsz,k] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
      end
    end
# interp call 1

    #FOREIGN.dlacpy!( 'F', N, N, A(1,1,K), LDA1, AS(1,1,K), LDA1 )
    AS = copy(A)

# interp call 2

    INFO = SLICOT.mb03vd!(N, P, ILO, IHI, A, TAU)
    @test INFO == 0
    INFO == 0 || return

    if ( INFO!=0 )
        @warn "mb03vd returns info=$INFO"
    end
# interp call 3

  for K in 1:P
    # FOREIGN.dlacpy!( 'L', N, N, A(1,1,K), LDA1, Q(1,1,K), LDQ1 )
        Q[1:N,1:N,K] .= tril(A[1:N,1:N,K])

      if ( N>1 )
          if ( N>2 && K==1 )
# interp call 4

              # FOREIGN.dlaset!( 'L', N-2, N-2, ZERO, ZERO, A(3,1,K), LDA1 )
              triu!(view(A,1:N,1:N,K),-1)

          elseif ( K>1 )
# interp call 5

              # FOREIGN.dlaset!( 'L', N-1, N-1, ZERO, ZERO, A(2,1,K), LDA1 )
              triu!(view(A,1:N,1:N,K))
          end # if
      end # if
  end # for K

    println(io, "A:")
    _nc = N
    _nr = N
    _nk = P
    show(io, "text/plain", A[1:_nr,1:_nc,1:_nk])
    println(io)

# interp call 6

    INFO = SLICOT.mb03vy!(N, P, ILO, IHI, Q, TAU)
    @test INFO == 0
    INFO == 0 || return
    println(io, "Q:")
    _nc = N
    _nr = N
    _nk = P
    show(io, "text/plain", Q[1:_nr,1:_nc,1:_nk])
    println(io)

    if ( INFO!=0 )
        @warn "mb03vy returns info=$INFO"
    else
    SSQ = ZERO
    for K in 1:P
        KP1 = mod(K,P)+1
# interp call 7

    #FOREIGN.dgemm!( 'T', 'N', N, N, N, ONE, Q(1,1,K), LDQ1, AS(1,1,K), LDA1, ZERO, QTA, LDQ1 )
        qta = Q[1:N,1:N,K]' * AS[1:N,1:N,K]
# interp call 8

    #FOREIGN.dgemm!( 'N', 'N', N, N, N, ONE, QTA, LDQ1, Q(1,1,KP1), LDQ1, -ONE, A(1,1,K), LDA1 )
        mul!(view(A,1:N,1:N,K),qta,view(Q,1:N,1:N,KP1),1,-1)
        #SSQ = DLAPY2( SSQ, DLANGE( 'Frobenius', N, N, A(1,1,K), LDA1, DWORK ) )
        SSQ = hypot(SSQ, norm(view(A,1:N,1:N,K)))
    end
    println(io, "norm(Q'*A*Q - Aout) = ", SSQ)
    end
    close(f)
end # run_mb03vd()
