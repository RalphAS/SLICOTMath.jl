# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb01td(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    LDB = NMAX
    LDWORK = NMAX-1
    A = Array{Float64,2}(undef, LDA,NMAX)
    B = Array{Float64,2}(undef, LDB,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
if ( N<0 || N>NMAX )
else

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
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
# interp call 1

    INFO = SLICOT.mb01td!(N, A, B)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"B:")
    _nc = N
    _nr = N
    show(io,"text/plain",B[1:_nr,1:_nc])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
