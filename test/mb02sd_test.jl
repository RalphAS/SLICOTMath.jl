# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mb02sd(datfile, io=stdout)
    ZERO = 0.0e0
    NIN = 5
    NOUT = 6
    NMAX = 20
    NRHMAX = 20
    LDB = NMAX
    LDH = NMAX
    LDWORK = 3*NMAX
    LIWORK = NMAX
    H = Array{Float64,2}(undef, LDH,NMAX)
    IPIV = Array{BlasInt,1}(undef, NMAX)
    B = Array{Float64,2}(undef, LDB,NRHMAX)
    DWORK = Array{Float64,1}(undef, LDWORK)
    IWORK = Array{BlasInt,1}(undef, LIWORK)
f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    N = parse(BlasInt, vs[1])
    NRHS = parse(BlasInt, vs[2])
    NORM =  vs[3][1]
    TRANS =  vs[4][1]
if ( N<0 || N>NMAX )
else

    vs = String[]
    _isz,_jsz = (N,N)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       H[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
if ( NRHS<0 || NRHS>NRHMAX )
else

    vs = String[]
    _isz,_jsz = (N,NRHS)
    while length(vs) < _isz*_jsz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    for i in 1:_isz
       _i0 = (i-1)*_jsz
       B[i,1:_jsz] .= parsex.(Float64, vs[_i0+1:_i0+_jsz])
    end
    if  N>2 
# interp call 1

    #FOREIGN.dlaset!( 'Lower', N-2, N-2, ZERO, ZERO, H(3,1), LDH )
    triu!(H, -1)

    end
# interp call 2

    INFO = SLICOT.mb02sd!(N, H, IPIV)
    @test INFO == 0
    INFO == 0 || return

    #HNORM = DLANHS( NORM, N, H, LDH, DWORK )
    HNORM = ccall((LinearAlgebra.BLAS.@blasfunc(dlanhs_), LinearAlgebra.LAPACK.liblapack),
                  Float64, (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                            Ptr{Float64}), NORM, N, H, LDH, DWORK)
# interp call 3

    RCOND, INFO = SLICOT.mb02td!(NORM, N, HNORM, H, IPIV)
    println(io, "RCOND = $RCOND")
    @test INFO == 0
    INFO == 0 || return

if ( INFO==0 && RCOND>eps(0.9) )
# interp call 4

    INFO = SLICOT.mb02rd!(TRANS, N, NRHS, H, IPIV, B)
    @test INFO == 0
    INFO == 0 || return

else
end # if
# interp output 1
    println(io, "B:")
    _nc = NRHS
    _nr = N
    show(io, "text/plain", B[1:_nr,1:_nc])
    println(io)

end # if
end # if
    close(f)
end # run_mb02sd()
