# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb04dp(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 10
    LDA  = NMAX
    LDC = NMAX
    LDDE = NMAX
    LDVW = NMAX
    A = Array{Float64,2}(undef, LDA,NMAX)
    DE = Array{Float64,2}(undef, LDDE,NMAX+1)
    C = Array{Float64,2}(undef, LDC,NMAX)
    VW = Array{Float64,2}(undef, LDVW,NMAX+1)
    LSCALE = Array{Float64,1}(undef, NMAX)
    RSCALE = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, 8*NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    JOB =  vs[2][1]
    THRESH = parse(Float64, replace(vs[3],'D'=>'E'))
    if ( N<=0 || N>NMAX )
        @error "Illegal N=$N"
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
    _isz,_jsz = (N,N+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       DE[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       C[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N+1)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       VW[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    ILO, INFO, IWARN = SLICOT.mb04dp!(JOB, N, THRESH, A, DE, C, VW, LSCALE, RSCALE, DWORK)
    @test INFO == 0
    INFO == 0 || return
    println(io, "ILO = $ILO")

# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "DE:")
    _nc = N+1
    _nr = N
    show(io, "text/plain", DE[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "C:")
    _nc = N
    _nr = N
    show(io, "text/plain", C[1:_nr,1:_nc])
    println(io)

# interp output 4
    println(io, "VW:")
    _nc = N+1
    _nr = N
    show(io, "text/plain", VW[1:_nr,1:_nc])
    println(io)

# interp output 5
    println(io, "LSCALE:")
    _nr = N
    show(io, "text/plain", LSCALE[1:_nr])
    println(io)

# interp output 6
    println(io, "RSCALE:")
    _nr = N
    show(io, "text/plain", RSCALE[1:_nr])
    println(io)

if ( LSAME( JOB, 'S' ) || LSAME( JOB, 'B' ) )
if ( !( THRESH==-2 || THRESH==-4 ) )
    println(io, "initial norms: ", DWORK[1:2])
    println(io, "initial norms: ", DWORK[3:4])
    println(io, "final threshold: ", DWORK[5])
else
    println(io, "IWARN = $IWARN")
end # if
end # if

end # run_mb04dp()
