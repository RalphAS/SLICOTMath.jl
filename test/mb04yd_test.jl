# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb04yd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    MMAX = 20
    NMAX = 20
    MNMIN = min( MMAX, NMAX )
    LDU = MMAX
    LDV = NMAX
    LDWORK = 6*MNMIN - 5
    Q = Array{Float64,1}(undef, MNMIN)
    E = Array{Float64,1}(undef, MNMIN-1)
    U = Array{Float64,2}(undef, LDU,MNMIN)
    V = Array{Float64,2}(undef, LDV,MNMIN)
    INUL = Array{BlasBool,1}(undef, MNMIN)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    M = parse(BlasInt, vs[1])
    N = parse(BlasInt, vs[2])
    THETA = parse(Float64, replace(vs[3],'D'=>'E'))
    RANK = parse(BlasInt, vs[4])
    TOL = parse(Float64, replace(vs[5],'D'=>'E'))
    RELTOL = parse(Float64, replace(vs[6],'D'=>'E'))
    JOBU =  vs[7][1]
    JOBV =  vs[8][1]
    MINMN = min( M, N )
if ( M<0 || M>MMAX )
elseif ( N<0 || N>NMAX )
elseif ( RANK>MINMN )
elseif ( RANK<0 && THETA<ZERO )
else

    vs = String[]
    _isz = MINMN 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    Q[1:_isz] .= parsex.(Float64, vs)
    

    vs = String[]
    _isz = MINMN-1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    E[1:_isz] .= parsex.(Float64, vs)
    
    RANK1 = RANK
    LJOBUU = LSAME( JOBU, 'U' )
    LJOBVU = LSAME( JOBV, 'U' )
    # CHECKME: partly translated condition
    if  LJOBUU 

    vs = String[]
    _isz,_jsz = (M,MINMN)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    end
    # CHECKME: partly translated condition
    if  LJOBVU 

    vs = String[]
    _isz,_jsz = (N,MINMN)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       V[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    end

    INUL .= false
    # CHECKME: partly translated condition
    if  LJOBUU||LJOBVU 

    vs = String[]
    _isz = MINMN 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    INUL[1:_isz] .= parsex.(Bool, vs)
    
    end
# interp call 1

    RANK, THETA, INFO, IWARN = SLICOT.mb04yd!(JOBU, JOBV, M, N, RANK, THETA, Q, E, U, V, INUL, TOL, RELTOL, LDWORK)
    println(io, "RANK = $RANK")
    println(io, "THETA = $THETA")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

if ( INFO!=0 )
else
if ( IWARN!=0 )
end # if
if ( !LSAME( JOBV, 'N' ) )
# interp output 1
    println(io,"V:")
    _nc = MINMN
    _nr = N
    show(io,"text/plain",V[1:_nr,1:_nc])
    println(io,)

end # if
if ( !LSAME( JOBU, 'N' ) )
# interp output 2
    println(io,"U:")
    _nc = MINMN
    _nr = M
    show(io,"text/plain",U[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
    close(f)
end # run_X()
