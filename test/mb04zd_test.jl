# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04zd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDQG = NMAX
    LDU = NMAX
    LDWORK = ( NMAX+NMAX )*( NMAX+NMAX+1 )
    ZERO = 0.0e0
    ONE = 1.0e0
    A = Array{Float64,2}(undef, LDA,NMAX)
    QG = Array{Float64,2}(undef, LDQG,NMAX+1)
    U = Array{Float64,2}(undef, LDU,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    COMPU =  vs[2][1]
    if ( N<0 || N>NMAX )
        @error "Illegal N=$N"
        return
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
    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*(_jsz+1)÷2
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    _i0 = 0
    for j in 1:_jsz
       QG[j,j+1:_jsz+1] .= parsex.(Float64, vs[_i0+1:_i0+_jsz-j+1])
       _i0 += _jsz-j+1
    end
    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*(_jsz+1)÷2
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    _i0 = 0
    for j in 1:_jsz
       QG[j:_jsz,j] .= parsex.(Float64, vs[_i0+1:_i0+_jsz-j+1])
       _i0 += _jsz-j+1
    end
    close(f)
    # interp call 1

    INFO = SLICOT.mb04zd!(COMPU, N, A, QG, U)

    @test INFO == 0
    INFO != 0 && return

# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "QG:")
    _nc = N
    _nr = N+1
    show(io, "text/plain", QG[1:_nr,1:_nc])
    println(io)

    Ht = zeros(2N,2N)
    for i in 1:N
        Ht[:,i] .= vcat(A[i,1:N],QG[1:i-1,i+1],QG[i,i+1:N+1])
        Ht[:,N+i] .= vcat(QG[i,1:i-1],QG[i:N,i],-A[1:N,i])
    end
    H = Matrix(Ht')
    println(io, "sq. reduced Hamiltonian:")
    show(io, "text/plain", H)
    println(io)
    # the original does this w/ BLAS in bits, but why bother?
    Hsq = H*H
    println(io, "squared sq. reduced Hamiltonian:")
    show(io, "text/plain", Hsq)
    println(io)
    llquadres = norm(Hsq[N+1:N,1:N])
    @test llquadres < 1e-12
    A2_hess_err =  norm(tril(Hsq[1:N,1:N],-2))
    @test A2_hess_err < 1e-12

end # run_mb04zd()
