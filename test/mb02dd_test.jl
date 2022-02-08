# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02dd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    KMAX = 20
    MMAX = 20
    NMAX = 20
    LDG = KMAX*( MMAX + NMAX )
    LDL = KMAX*( MMAX + NMAX )
    LDR = KMAX*( MMAX + NMAX )
    LDT = KMAX*( MMAX + NMAX )
    LDWORK = ( MMAX + NMAX - 1 )*KMAX
    LCS = 3*LDWORK
    tdim = KMAX*(MMAX+NMAX)
    T = Array{Float64,2}(undef, LDT,KMAX*(MMAX+NMAX))
    G = Array{Float64,2}(undef, LDG,KMAX*(MMAX+NMAX))
    R = Array{Float64,2}(undef, LDR,KMAX*(MMAX+NMAX))
    L = Array{Float64,2}(undef, LDL,KMAX*(MMAX+NMAX))
    CS = Array{Float64,1}(undef, LCS)
    DWORK = Array{Float64,1}(undef, LDWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    K = parse(BlasInt, vs[2])
    M = parse(BlasInt, vs[3])
    JOB =  vs[4][1]
    TYPET =  vs[5][1]
    S = ( N + M )*K
if ( N<=0 || N>NMAX )
else
if ( K<=0 || K>KMAX )
else
if ( M<=0 || M>MMAX )
else
if ( LSAME( TYPET, 'R' ) )

    vs = String[]
    _isz,_jsz = (K,S)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
else

    vs = String[]
    _isz,_jsz = (S,K)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
# interp call 1

    INFO = SLICOT.mb02cd!(JOB, TYPET, K, N, T, G, R, L, CS, LCS, LDWORK)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io, "R:")
    _nc = N*K
    _nr = N*K
    show(io, "text/plain", R[1:_nr,1:_nc])
    println(io)

if ( LSAME( JOB, 'R' ) || LSAME( JOB, 'A' ) )
if ( LSAME( TYPET, 'R' ) )
# interp output 2
    println(io, "G:")
    _nc = N*K
    _nr = 2*K
    show(io, "text/plain", G[1:_nr,1:_nc])
    println(io)

else
# interp output 3
    println(io, "G:")
    _nc = 2*K
    _nr = N*K
    show(io, "text/plain", G[1:_nr,1:_nc])
    println(io)

end # if
end # if
if ( LSAME( JOB, 'A' ) )
# interp output 4
    println(io, "L:")
    _nc = N*K
    _nr = N*K
    show(io, "text/plain", L[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( TYPET, 'R' ) )
# interp call 2

    # FOREIGN.dlacpy!( 'All', N*K, K, R(1,(N-1)*K+1), LDR, R(K+1,N*K+1), LDR )
    for ir in 1:N*K
        for ic in 1:K
            R[K+ir,N*K+ic] = R[ir,(N-1)*K+ic]
        end
    end

# interp call 3

    INFO = SLICOT.mb02dd!(JOB, TYPET, K, M, N, view(T,1:K,N*K+1:tdim), T, G,
                          view(R,:,N*K+1:(N*M)*K), view(L,N*K+1:LDL,:), CS, LDWORK)
    @test INFO == 0
    INFO == 0 || return

else
# interp call 4

    #FOREIGN.dlacpy!( 'All', K, N*K, R((N-1)*K+1,1), LDR, R(N*K+1,K+1), LDR )
    for ir in 1:K
        for ic in 1:N*K
            R[N*K+ir,K+ic] = R[(N-1)*K+ir,ic]
        end
    end

# interp call 5

    INFO = SLICOT.mb02dd!(JOB, TYPET, K, M, N, view(T,N*K+1:(N+M)*K,:), T, G, view(R,N*K+1:(N+M)*K,:), view(L,:, N*K+1:(N+M)*K), CS, LDWORK)
    @test INFO == 0
    INFO == 0 || return

end # if
if ( INFO!=0 )
else
# interp output 5
    println(io, "R:")
    _nc = S
    _nr = S
    show(io, "text/plain", R[1:_nr,1:_nc])
    println(io)

if ( LSAME( JOB, 'R' ) || LSAME( JOB, 'A' ) )
if ( LSAME( TYPET, 'R' ) )
# interp output 6
    println(io, "G:")
    _nc = S
    _nr = 2*K
    show(io, "text/plain", G[1:_nr,1:_nc])
    println(io)

else
# interp output 7
    println(io, "G:")
    _nc = 2*K
    _nr = S
    show(io, "text/plain", G[1:_nr,1:_nc])
    println(io)

end # if
end # if
if ( LSAME( JOB, 'A' ) )
# interp output 8
    println(io, "L:")
    _nc = S
    _nr = S
    show(io, "text/plain", L[1:_nr,1:_nc])
    println(io)

end # if
end # if
end # if
end # if
end # if
end # if
    close(f)
end # run_mb02dd()
