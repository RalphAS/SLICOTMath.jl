# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc03md(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    CP1MAX = 10
    CP2MAX = 10
    DP1MAX = 10
    DP2MAX = 10
    DP3MAX = 20
    RP1MAX = 10
    LDP11 = RP1MAX
    LDP12 = CP1MAX
    LDP21 = CP1MAX
    LDP22 = CP2MAX
    LDP31 = RP1MAX
    LDP32 = CP2MAX
    P1 = Array{Float64,3}(undef, LDP11,LDP12,DP1MAX+1)
    P2 = Array{Float64,3}(undef, LDP21,LDP22,DP2MAX+1)
    P3 = Array{Float64,3}(undef, LDP31,LDP32,DP3MAX+1)
    DWORK = Array{Float64,1}(undef, CP1MAX)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    RP1 = parse(BlasInt, vs[1])
    CP1 = parse(BlasInt, vs[2])
    CP2 = parse(BlasInt, vs[3])
if ( RP1<0 || RP1>RP1MAX )
elseif ( CP1<0 || CP1>CP1MAX )
elseif ( CP2<0 || CP2>CP2MAX )
else
    vs = split(readline(f))
    DP1 = parse(BlasInt, vs[1])
if ( DP1<=-2 || DP1>DP1MAX )
else

    _isz,_jsz,_ksz = (RP1,CP1,DP1 + 1)
    for k in 1:_ksz
      for j in 1:_jsz
        vs = String[]
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        P1[1:_isz,j,k] .= parsex.(Float64, vs)
      end
    end
    vs = split(readline(f))
    DP2 = parse(BlasInt, vs[1])
if ( DP2<=-2 || DP2>DP2MAX )
else

    _isz,_jsz,_ksz = (CP1,CP2,DP2 + 1)
    for k in 1:_ksz
      for j in 1:_jsz
        vs = String[]
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        P2[1:_isz,j,k] .= parsex.(Float64, vs)
      end
    end
    vs = split(readline(f))
    DP3 = parse(BlasInt, vs[1])
if ( DP3<=-2 || DP3>DP3MAX )
else

    _isz,_jsz,_ksz = (RP1,CP2,DP3 + 1)
    for k in 1:_ksz
      for j in 1:_jsz
        vs = String[]
        while length(vs) < _isz
            append!(vs, replace.(split(readline(f)),'D'=>'E'))
        end
        P3[1:_isz,j,k] .= parsex.(Float64, vs)
      end
    end
    vs = split(readline(f))
    ALPHA = parse(Float64, replace(vs[1],'D'=>'E'))
# interp call 1

    DP3, INFO = SLICOT.mc03md!(RP1, CP1, CP2, DP1, DP2, DP3, ALPHA, P1, P2, P3)
    println(io, "DP3 = $DP3")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
    for i in 1:RP1
        for j in 1:CP2
            println(io, "elt($i,$j): ",P3[i,j,1:DP3+1])
        end
    end
end # if
end # if
end # if
end # if
end # if
    close(f)
end # run_mc03md()
