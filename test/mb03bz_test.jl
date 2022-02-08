# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb03bz(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 6
    NMAX = 50
    LDA1 = NMAX
    LDA2 = NMAX
    LDQ1 = NMAX
    LDQ2 = NMAX
    LDWORK = NMAX
    LZWORK = NMAX
    S = Array{BlasInt,1}(undef, KMAX)
    A = Array{ComplexF64,3}(undef, LDA1,LDA2,KMAX)
    Q = Array{ComplexF64,3}(undef, LDQ1,LDQ2,KMAX)
    ALPHA = Array{ComplexF64,1}(undef, NMAX)
    BETA = Array{ComplexF64,1}(undef, NMAX)
    SCAL = Array{BlasInt,1}(undef, NMAX)
    ZWORK = Array{ComplexF64,1}(undef, LZWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    JOB =  vs[1][1]
    COMPQ =  vs[2][1]
    K = parse(BlasInt, vs[3])
    N = parse(BlasInt, vs[4])
    ILO = parse(BlasInt, vs[5])
    IHI = parse(BlasInt, vs[6])
if ( N<0 || N>NMAX )
else

    vs = String[]
    _isz = K 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    S[1:_isz] .= parsex.(BlasInt, vs)
    

    vs = String[]
    _isz,_jsz,_ksz = (N,N,K)
    while length(vs) < _isz*_jsz*_ksz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for k in 1:_ksz
      for i in 1:_isz
        _i0 = (i-1)*_jsz + (k-1)*_jsz*_isz
        A[i,1:_jsz,k] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
      end
    end
# interp call 1

    INFO = SLICOT.mb03bz!(JOB, COMPQ, K, N, ILO, IHI, S, A, Q, ALPHA, BETA, SCAL, LDWORK, LZWORK)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( JOB, 'S' ) )
# interp output 1
    println(io,"A:")
    _nc = N
    _nr = N
    _nk = K
    show(io,"text/plain",A[1:_nr,1:_nc,1:_nk])
    println(io,)

end # if
if ( !LSAME( COMPQ, 'N' ) )
# interp output 2
    println(io,"Q:")
    _nc = N
    _nr = N
    _nk = K
    show(io,"text/plain",Q[1:_nr,1:_nc,1:_nk])
    println(io,)

end # if
# interp output 3
    println(io,"ALPHA:")
    _nr = N
    show(io,"text/plain",ALPHA[1:_nr])
    println(io,)

# interp output 4
    println(io,"BETA:")
    _nr = N
    show(io,"text/plain",BETA[1:_nr])
    println(io,)

# interp output 5
    println(io,"SCAL:")
    _nr = N
    show(io,"text/plain",SCAL[1:_nr])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
