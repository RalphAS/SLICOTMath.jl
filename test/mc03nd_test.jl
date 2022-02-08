# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc03nd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DPMAX = 2
    MPMAX = 5
    NPMAX = 4
    LDP1 = MPMAX
    LDP2 = NPMAX
    LDNULL = NPMAX
    LDKER1 = NPMAX
    LDKER2 = NPMAX
    M = DPMAX*MPMAX
    N = ( DPMAX-1 )*MPMAX+NPMAX
    LIWORK = M+2*max( N,M+1 )+N
    LDWORK = M*N^2+2*M*N+2*N^2
    P = Array{Float64,3}(undef, LDP1,LDP2,DPMAX+1)
    GAM = Array{BlasInt,1}(undef, M+1)
    NULLSP = Array{Float64,2}(undef, LDNULL,(M+1))
    KER = Array{Float64,3}(undef, LDKER1,LDKER2,M+1)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
    DWORK = Array{Float64,1}(undef, LDWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    MP = parse(BlasInt, vs[1])
    NP = parse(BlasInt, vs[2])
    DP = parse(BlasInt, vs[3])
    TOL = parse(Float64, replace(vs[4],'D'=>'E'))
if ( MP<0 || MP>MPMAX )
elseif ( NP<0 || NP>NPMAX )
elseif ( DP<=0 || DP>DPMAX )
else

    _isz,_jsz,_ksz = (MP,NP,DP + 1)
    for k in 1:_ksz
      for i in 1:_isz
        vs = String[]
        while length(vs) < _jsz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        P[i,1:_jsz,k] .= parsex.(Float64, vs)
      end
    end
# interp call 1

    DK, INFO = SLICOT.mc03nd!(MP, NP, DP, P, GAM, NULLSP, KER, TOL, IWORK, LDWORK)
    println(io, "DK = $DK")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
elseif ( DK<0 )
else
    NK = 0
    M1 = 0
    for i in 1:DK+1
        NK = NK + GAM[i]
        M1 = M1 + GAM[i]*i
    end
    
# interp output 1
    println(io, "NULLSP:")
    _nc = M1
    _nr = NP
    show(io, "text/plain", NULLSP[1:_nr,1:_nc])
    println(io)
    println(io, "KER:")
    for i in 1:NP
        for j in 1:NK
            println(io, "elt($i,$j: ",KER[i,j,1:DK+1])
        end
    end
end # if
end # if
    close(f)
end # run_mc03nd()
