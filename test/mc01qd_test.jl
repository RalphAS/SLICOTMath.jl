# Portions translated from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.
function run_mc01qd(datfile, io=stdout)
    NIN = 5
    NOUT = 6
    DAMAX = 10
    DBMAX = 10
    A = Array{Float64,1}(undef, DAMAX+1)
    B = Array{Float64,1}(undef, DBMAX+1)
    RQ = Array{Float64,1}(undef, DAMAX+1)
    f = open(datfile,"r")
    readline(f)
    vs = split(readline(f))
    DA = parse(BlasInt, vs[1])
    if ( DA<=-2 || DA>DAMAX )
        @error "Illegal DA=$DA"
    end

    vs = String[]
    _isz = DA+1 
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    A[1:_isz] .= parsex.(Float64, vs)

    vs = split(readline(f))
    DB = parse(BlasInt, vs[1])
    DBB = DB
    if ( DB<=-1 || DB>DBMAX )
        @error "Illegal DB=$DB"
    end

    vs = String[]
    _isz = DB+1
    while length(vs) < _isz
        append!(vs, replace.(split(readline(f)),'D'=>'E'))
    end
    B[1:_isz] .= parsex.(Float64, vs)
    close(f)
# interp call 1

    DB, INFO, IWARN = SLICOT.mc01qd!(DA, DB, A, B, RQ)
    println(io, "DB = $DB")
    @test INFO == 0
    INFO == 0 || return
    println(io, "IWARN = $IWARN")

    if ( IWARN!=0 )
        println(io, "degree reduced from $DBB to $DB")
    end # if
    DQ = DA - DB
    DR = DB - 1
    IMAX = DQ
    if DR > IMAX
        IMAX = DR
    end
    println(io, "polynomials (power Q-coefft R-coefft):")
    for i in 0:IMAX
        if ( i<=DQ && i <=DR )
            println(io, "$i $(RQ[DB+i+1]) $(RQ[i+1])")
        elseif ( i<=DQ )
            println(io, "$i $(RQ[DB+i+1])")
        else
            println(io, "$i $(RQ[i+1])")
        end # if
    end

end # run_mc01qd()
