# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mb05od(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NMAX = 20
    LDA = NMAX
    NDIAG = 9
    LDWORK = NMAX*( 2*NMAX+NDIAG+1 )+NDIAG
    A = Array{Float64,2}(undef, LDA,NMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, NMAX)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    DELTA = parse(Float64, replace(vs[2],'D'=>'E'))
    BALANC =  vs[3][1]
if ( N<=0 || N>NMAX )
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

    MDIG, IDIG, INFO, IWARN = SLICOT.mb05od!(BALANC, N, NDIAG, DELTA, A, LDWORK)
    println(io, "MDIG = $MDIG")
    println(io, "IDIG = $IDIG")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

if ( INFO!=0 )
else
# interp output 1
    println(io,"A:")
    _nc = N
    _nr = N
    show(io,"text/plain",A[1:_nr,1:_nc])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
