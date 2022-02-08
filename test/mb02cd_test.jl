# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02cd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    KMAX = 20
    NMAX = 20
    LDG = 2*KMAX
    LDL = NMAX*KMAX
    LDR = NMAX*KMAX
    LDT = KMAX
    LDWORK = ( NMAX - 1 )*KMAX
    LCS = 3*LDWORK
    T = Array{Float64,2}(undef, LDT,NMAX*KMAX)
    G = Array{Float64,2}(undef, LDG,NMAX*KMAX)
    R = Array{Float64,2}(undef, LDR,NMAX*KMAX)
    L = Array{Float64,2}(undef, LDL,NMAX*KMAX)
    CS = Array{Float64,1}(undef, LCS)
    DWORK = Array{Float64,1}(undef, LDWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    K = parse(BlasInt, vs[2])
    JOB =  vs[3][1]
    TYPET = 'R'
    M = N*K
if ( N<=0 || N>NMAX )
else
if ( K<=0 || K>KMAX )
else

    vs = String[]
    _isz,_jsz = (K,M)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       T[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    INFO = SLICOT.mb02cd!(JOB, TYPET, K, N, T, G, R, L, CS, LCS, LDWORK)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( JOB, 'G' ) || LSAME( JOB, 'A' ) || LSAME( JOB, 'L' ) || LSAME( JOB, 'R' ) )
# interp call 2

    #FOREIGN.dlaset!( 'Full', K, K, ZERO, ZERO, G(K+1,1), LDG )
    G[K+1:2*K,1:K] .= ZERO

# interp output 1
    println(io, "G:")
    _nc = M
    _nr = 2*K
    show(io, "text/plain", G[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( JOB, 'L' ) || LSAME( JOB, 'A' ) )
# interp output 2
    println(io, "L:")
    _nc = M
    _nr = M
    show(io, "text/plain", L[1:_nr,1:_nc])
    println(io)

end # if
if ( LSAME( JOB, 'R' ) || LSAME( JOB, 'A' ) || LSAME( JOB, 'O' ) )
# interp output 3
    println(io, "R:")
    _nc = M
    _nr = M
    show(io, "text/plain", R[1:_nr,1:_nc])
    println(io)

end # if
end # if
end # if
end # if
    close(f)
end # run_mb02cd()
