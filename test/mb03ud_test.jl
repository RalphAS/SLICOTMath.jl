# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb03ud(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 10
    LDA = NMAX
    LDQ = NMAX
    LDWORK = max( 1, 5*NMAX )
    A = Array{Float64,2}(undef, LDA,NMAX)
    Q = Array{Float64,2}(undef, LDQ,NMAX)
    SV = Array{Float64,1}(undef, NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    JOBQ =  vs[2][1]
    JOBP =  vs[3][1]
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
# interp call 1

    INFO = SLICOT.mb03ud!(JOBQ, JOBP, N, A, Q, SV)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"SV:")
    _nr = N
    show(io,"text/plain",SV[1:_nr])
    println(io,)

if ( LSAME( JOBP, 'V' ) )
# interp output 2
    println(io,"A:")
    _nc = N
    _nr = N
    show(io,"text/plain",A[1:_nr,1:_nc])
    println(io,)

end # if
if ( LSAME( JOBQ, 'V' ) )
# interp output 3
    println(io,"Q:")
    _nc = N
    _nr = N
    show(io,"text/plain",Q[1:_nr,1:_nc])
    println(io,)

end # if
end # if
end # if
    close(f)
end # run_X()
