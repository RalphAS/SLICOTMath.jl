# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01xd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    NEV = 3
    LDWORK = 42
    ZERO = 0.0e0
    EVR = Array{Float64,1}(undef, NEV)
    EVI = Array{Float64,1}(undef, NEV)
    EVQ = Array{Float64,1}(undef, NEV)
    DWORK = Array{Float64,1}(undef, LDWORK)
    RT = Array{Float64,1}(undef, 2)
    RTS = Array{ComplexF64,1}(undef, 3)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    ALPHA = parse(Float64, replace(vs[1],'D'=>'E'))
    BETA = parse(Float64, replace(vs[2],'D'=>'E'))
    GAMMA = parse(Float64, replace(vs[3],'D'=>'E'))
    DELTA = parse(Float64, replace(vs[4],'D'=>'E'))
# interp call 1

    INFO = SLICOT.mc01xd!(ALPHA, BETA, GAMMA, DELTA, EVR, EVI, EVQ)
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
# interp output 1
    println(io,"EVR:")
    _nr = 3
    show(io,"text/plain",EVR[1:_nr])
    println(io,)

# interp output 2
    println(io,"EVI:")
    _nr = 3
    show(io,"text/plain",EVI[1:_nr])
    println(io,)

# interp output 3
    println(io,"EVQ:")
    _nr = 3
    show(io,"text/plain",EVQ[1:_nr])
    println(io,)

if ( true )
    RTS = (EVR + im * EVI) ./ EVQ
# interp output 4
    println(io,"RTS:")
    _nr = 3
    show(io,"text/plain",RTS[1:_nr])
    println(io,)

end # if
end # if
    close(f)
end # run_X()
