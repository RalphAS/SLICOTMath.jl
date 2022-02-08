# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb4dlz(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 10
    LDA = NMAX
    LDB = NMAX
    A = Array{ComplexF64,2}(undef, LDA,NMAX)
    B = Array{ComplexF64,2}(undef, LDB,NMAX)
    DWORK = Array{Float64,1}(undef, 8*NMAX)
    LSCALE = Array{Float64,1}(undef, NMAX)
    RSCALE = Array{Float64,1}(undef, NMAX)
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
       A[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(ComplexF64, vs[_i0+1:_i0+_jsz])
    end
    close(f)
# interp call 1

    ILO, IHI, INFO, IWARN = SLICOT.mb4dlz!( JOB, N, THRESH, A, B,
                                     LSCALE, RSCALE, DWORK)
    @test INFO == 0
    INFO == 0 || return
    println(io, "ILO = $ILO")
    println(io, "IHI = $IHI")

# interp output 1
    println(io, "A:")
    _nc = N
    _nr = N
    show(io, "text/plain", A[1:_nr,1:_nc])
    println(io)

# interp output 2
    println(io, "B:")
    _nc = N
    _nr = N
    show(io, "text/plain", B[1:_nr,1:_nc])
    println(io)

# interp output 3
    println(io, "LSCALE:")
    _nr = N
    show(io, "text/plain", LSCALE[1:_nr])
    println(io)

# interp output 4
    println(io, "RSCALE:")
    _nr = N
    show(io, "text/plain", RSCALE[1:_nr])
    println(io)

    if ( LSAME( JOB, 'S' ) || LSAME( JOB, 'B' ) )
        if ( !( THRESH==-2 || THRESH==-4 ) )
            # interp output 5
            println(io, "initial 1-norms:")
            _nr = 2
            show(io, "text/plain", DWORK[1:_nr])
            println(io)

            println(io, "final 1-norms:")
            show(io, "text/plain", DWORK[3:4])
            println(io)

# interp output 7
            println(io, "final threshold:")
            show(io, "text/plain", DWORK[5])
            println(io)
        else
            @warn "mb4dlz returned iwarn=$IWARN"
        end
    end
end # run_mb4dlz()
