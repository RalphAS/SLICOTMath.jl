# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb03td(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    NMAX = 100
    LDA  = NMAX
    LDG  = NMAX
    LDRES  = NMAX
    LDU1 = NMAX
    LDU2 = NMAX
    LDWORK = 8*NMAX
    SELECT = Array{Bool,1}(undef, NMAX)
    LOWER = Array{Bool,1}(undef, NMAX)
    A = Array{Float64,2}(undef, LDA,NMAX)
    G = Array{Float64,2}(undef, LDG,NMAX)
    U1 = Array{Float64,2}(undef, LDU1,NMAX)
    U2 = Array{Float64,2}(undef, LDU2,NMAX)
    WR = Array{Float64,1}(undef, NMAX)
    WI = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RES = Array{Float64,2}(undef, LDRES,NMAX)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    TYP =  vs[2][1]
    COMPU =  vs[3][1]
if ( N<=0 || N>NMAX )
else

    vs = String[]
    _isz = N 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    SELECT[1:_isz] .= parsex.(Bool, vs)


    vs = String[]
    _isz = N 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    LOWER[1:_isz] .= parsex.(Bool, vs)


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
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       G[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
if ( LSAME( COMPU, 'U' ) )

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U1[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       U2[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
end # if
# interp call 1
    select = convert(Vector{BlasBool},SELECT)
    lower = convert(Vector{BlasBool},LOWER)
    M, INFO = SLICOT.mb03td!(TYP, COMPU, select, lower, N, A, G, U1, U2, WR, WI, LDWORK)
    println(io, "M = $M")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
if ( LSAME( COMPU, 'U' ) )
# interp output 1
    println(io, "U1:")
    _nc = N
    _nr = N
    show(io, "text/plain", U1[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "U2:")
    _nc = N
    _nr = N
    show(io, "text/plain", U2[1:_nr,1:_nc])
    println(io)

# unable to translate write statement:
#  write MA02JD( .FALSE., .FALSE., N, U1, LDU1, U2, LDU2, RES, LDRES )
end # if
# interp output 3
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

    Gfull = zeros(N,N)
if ( LSAME( TYP, 'S' ) )
    # interp output 4
    for i in 1:N
        for j in 1:i-1
            Gfull[i,j] = -G[j,i]
        end
        for j in i+1:N
            Gfull[i,j] = G[i,j]
        end
    end
else
    for i in 1:N
        for j in 1:i-1
            Gfull[i,j] = G[j,i]
        end
        for j in i:N
            Gfull[i,j] = G[i,j]
        end
    end
end # if
    println(io, "G:")
    _nc = N
    _nr = N
    show(io, "text/plain", Gfull[1:_nr,1:_nc])
    println(io)

end # if
end # if
    close(f)
end # run_mb03td()
