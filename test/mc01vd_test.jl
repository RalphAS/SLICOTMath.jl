# Julia code
# Copyright (c) 2022 the SLICOTMath.jl developers
# Portions extracted from SLICOT-Reference distribution:
# Copyright (c) 2002-2020 NICONET e.V.

function run_mc01vd(datfile, io=stdout)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    A = parse(Float64, replace(vs[1],'D'=>'E'))
    B = parse(Float64, replace(vs[2],'D'=>'E'))
    C = parse(Float64, replace(vs[3],'D'=>'E'))
# interp call 1

    Z1RE, Z1IM, Z2RE, Z2IM, INFO = SLICOT.mc01vd!(A, B, C)
    println(io, "Z1RE = $Z1RE")
    println(io, "Z1IM = $Z1IM")
    println(io, "Z2RE = $Z2RE")
    println(io, "Z2IM = $Z2IM")
    @test INFO == 0
    INFO == 0 || return

if ( INFO!=0 )
else
end # if
    close(f)
end # run_X()
