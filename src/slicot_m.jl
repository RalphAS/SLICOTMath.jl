const BlasBool = BlasInt

# automatically built wrappers for SLICOT


"""
$(TYPEDSIGNATURES)
returns (yr,  yi)
"""
function ma01ad!(xr::Number, xi::Number)

    yr = Ref{Float64}()
    yi = Ref{Float64}()

    ccall((:ma01ad_, libslicot), Cvoid, (Ref{Float64}, Ref{Float64},
            Ptr{Float64}, Ptr{Float64}), xr, xi, yr, yi)

    return yr[], yi[]
end


"""
$(TYPEDSIGNATURES)
returns (alpha, beta,  scal)
"""
function ma01bd!(base::Number, lgbas::Number, k::Integer,
    s::AbstractVector{BlasInt}, a::AbstractVector{Float64},
    inca::Integer)

    alpha = Ref{Float64}()
    beta = Ref{Float64}()
    scal = Ref{BlasInt}()

    ccall((:ma01bd_, libslicot), Cvoid, (Ref{Float64}, Ref{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}), base, lgbas,
            k, s, a, inca, alpha, beta, scal)

    return alpha[], beta[], scal[]
end


"""
$(TYPEDSIGNATURES)
returns (alpha, beta,  scal)
"""
function ma01bz!(base::Number, k::Integer, s::AbstractVector{BlasInt},
    a::AbstractVector{ComplexF64}, inca::Integer)

    alpha = Ref{ComplexF64}()
    beta = Ref{ComplexF64}()
    scal = Ref{BlasInt}()

    ccall((:ma01bz_, libslicot), Cvoid, (Ref{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{BlasInt}), base,
            k, s, a, inca, alpha, beta, scal)

    return alpha[], beta[], scal[]
end


"""
$(TYPEDSIGNATURES)
"""
function ma01cd!(a::Number, ia::Integer, b::Number, ib::Integer)


    jlres = ccall((:ma01cd_, libslicot), BlasInt, (Ref{Float64},
            Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}), a, ia, b, ib)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02ad!(job::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))

    ccall((:ma02ad_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Clong), job, m, n, a, lda, b, ldb, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02bd!(side::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))

    ccall((:ma02bd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), side,
            m, n, a, lda, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02bz!(side::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))

    ccall((:ma02bz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong),
            side, m, n, a, lda, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02cd!(n::Integer, kl::Integer, ku::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))

    ccall((:ma02cd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, kl, ku, a,
            lda)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02cz!(n::Integer, kl::Integer, ku::Integer,
    a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))

    ccall((:ma02cz_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, kl, ku,
            a, lda)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02dd!(job::AbstractChar, uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, ap::AbstractVector{Float64})

    lda = max(1,stride(a,2))

    ccall((:ma02dd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Clong, Clong), job, uplo, n, a, lda, ap, 1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02ed!(uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))

    ccall((:ma02ed_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02es!(uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))

    ccall((:ma02es_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02ez!(uplo::AbstractChar, trans::AbstractChar,
    skew::AbstractChar, n::Integer, a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))

    ccall((:ma02ez_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Clong, Clong, Clong), uplo, trans, skew, n, a, lda, 1,
            1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (x1, c, s,  info)
"""
function ma02fd!(ini_x1::Number, x2::Number)

    x1 = Ref{Float64}(ini_x1)
    c = Ref{Float64}()
    s = Ref{Float64}()
    info = Ref{BlasInt}()

    ccall((:ma02fd_, libslicot), Cvoid, (Ptr{Float64}, Ref{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}), x1, x2, c, s,
            info)
    chkargsok(info[])

    return x1[], c[], s[], info[]
end


"""
$(TYPEDSIGNATURES)
"""
function ma02gd!(n::Integer, a::AbstractMatrix{Float64}, k1::Integer,
    k2::Integer, ipiv::AbstractVector{BlasInt}, incx::Integer)

    lda = max(1,stride(a,2))

    ccall((:ma02gd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}), n, a, lda, k1, k2, ipiv, incx)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02gz!(n::Integer, a::AbstractMatrix{ComplexF64},
    k1::Integer, k2::Integer, ipiv::AbstractVector{BlasInt},
    incx::Integer)

    lda = max(1,stride(a,2))

    ccall((:ma02gz_, libslicot), Cvoid, (Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), n, a, lda,
            k1, k2, ipiv, incx)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02hd!(job::AbstractChar, m::Integer, n::Integer,
    diag::Number, a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))

    jlres = ccall((:ma02hd_, libslicot), BlasBool, (Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Clong), job, m, n, diag, a, lda, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02hz!(job::AbstractChar, m::Integer, n::Integer,
    diag::Complex, a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))

    jlres = ccall((:ma02hz_, libslicot), BlasBool, (Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
            Ptr{ComplexF64}, Ref{BlasInt}, Clong), job, m, n, diag,
            a, lda, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02id!(typ::AbstractChar, norm::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))

    jlres = ccall((:ma02id_, libslicot), Float64, (Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong, Clong),
            typ, norm, n, a, lda, qg, ldqg, dwork, 1, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02iz!(typ::AbstractChar, norm::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, qg::AbstractMatrix{ComplexF64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))

    jlres = ccall((:ma02iz_, libslicot), Float64, (Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Clong,
            Clong), typ, norm, n, a, lda, qg, ldqg, dwork, 1, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02jd!(ltran1::Bool, ltran2::Bool, n::Integer,
    q1::AbstractMatrix{Float64}, q2::AbstractMatrix{Float64})

    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldres = max(1,n)
    res = Matrix{Float64}(undef, ldres,n)

    jlres = ccall((:ma02jd_, libslicot), Float64, (Ref{BlasBool},
            Ref{BlasBool}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
            ltran1, ltran2, n, q1, ldq1, q2, ldq2, res, ldres)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02jz!(ltran1::Bool, ltran2::Bool, n::Integer,
    q1::AbstractMatrix{ComplexF64}, q2::AbstractMatrix{ComplexF64})

    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldres = max(1,n)
    res = Matrix{ComplexF64}(undef, ldres,n)

    jlres = ccall((:ma02jz_, libslicot), Float64, (Ref{BlasBool},
            Ref{BlasBool}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}), ltran1, ltran2, n, q1, ldq1, q2, ldq2,
            res, ldres)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02md!(norm::AbstractChar, uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))

    jlres = ccall((:ma02md_, libslicot), Float64, (Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Clong, Clong), norm, uplo, n, a, lda,
            dwork, 1, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02mz!(norm::AbstractChar, uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))

    jlres = ccall((:ma02mz_, libslicot), Float64, (Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Clong, Clong), norm, uplo, n, a, lda,
            dwork, 1, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02nz!(uplo::AbstractChar, trans::AbstractChar,
    skew::AbstractChar, n::Integer, k::Integer, l::Integer,
    a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))

    ccall((:ma02nz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong),
            uplo, trans, skew, n, k, l, a, lda, 1, 1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function ma02od!(skew::AbstractChar, m::Integer,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))

    jlres = ccall((:ma02od_, libslicot), BlasInt, (Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Clong), skew, m, a, lda, de, ldde, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
"""
function ma02oz!(skew::AbstractChar, m::Integer,
    a::AbstractMatrix{ComplexF64}, de::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))

    jlres = ccall((:ma02oz_, libslicot), BlasInt, (Ref{UInt8},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Clong), skew, m, a, lda,
            de, ldde, 1)

    return jlres
end


"""
$(TYPEDSIGNATURES)
returns (nzr,  nzc)
"""
function ma02pd!(m::Integer, n::Integer, a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    nzr = Ref{BlasInt}()
    nzc = Ref{BlasInt}()

    ccall((:ma02pd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
            m, n, a, lda, nzr, nzc)

    return nzr[], nzc[]
end


"""
$(TYPEDSIGNATURES)
returns (nzr,  nzc)
"""
function ma02pz!(m::Integer, n::Integer,
    a::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))
    nzr = Ref{BlasInt}()
    nzc = Ref{BlasInt}()

    ccall((:ma02pz_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}), m, n, a, lda, nzr, nzc)

    return nzr[], nzc[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01kd!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    k::Integer, alpha::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, beta::Number,
    c::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()

    ccall((:mb01kd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc,
            info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ld!(uplo::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, alpha::Number, beta::Number,
    r::AbstractMatrix{Float64}, a::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, ldwork::Integer)

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01ld_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), uplo, trans, m, n, alpha,
            beta, r, ldr, a, lda, x, ldx, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01md!(uplo::AbstractChar, n::Integer, alpha::Number,
    a::AbstractMatrix{Float64}, x::AbstractVector{Float64},
    incx::Integer, beta::Number, y::AbstractVector{Float64},
    incy::Integer)

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01md_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Clong), uplo, n, alpha, a, lda, x, incx, beta, y, incy,
            1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01nd!(uplo::AbstractChar, n::Integer, alpha::Number,
    x::AbstractVector{Float64}, incx::Integer,
    y::AbstractVector{Float64}, incy::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01nd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo,
            n, alpha, x, incx, y, incy, a, lda, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01oc!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, x::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()

    ccall((:mb01oc_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), uplo, trans,
            n, alpha, beta, r, ldr, h, ldh, x, ldx, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01od!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, x::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, ldwork::Integer)

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    ldx = max(1,stride(x,2))
    lde = max(1,stride(e,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01od_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), uplo, trans,
            n, alpha, beta, r, ldr, h, ldh, x, ldx, e, lde, dwork,
            ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01oe!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, e::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    lde = max(1,stride(e,2))
    info = Ref{BlasInt}()

    ccall((:mb01oe_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Clong, Clong), uplo, trans, n, alpha,
            beta, r, ldr, h, ldh, e, lde, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01oh!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, a::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01oh_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Clong, Clong), uplo, trans, n, alpha,
            beta, r, ldr, h, ldh, a, lda, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01oo!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    h::AbstractMatrix{Float64}, x::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, p::AbstractMatrix{Float64})

    ldh = max(1,stride(h,2))
    ldx = max(1,stride(x,2))
    lde = max(1,stride(e,2))
    ldp = max(1,stride(p,2))
    info = Ref{BlasInt}()

    ccall((:mb01oo_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), uplo, trans,
            n, h, ldh, x, ldx, e, lde, p, ldp, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01os!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    h::AbstractMatrix{Float64}, x::AbstractMatrix{Float64},
    p::AbstractMatrix{Float64})

    ldh = max(1,stride(h,2))
    ldx = max(1,stride(x,2))
    ldp = max(1,stride(p,2))
    info = Ref{BlasInt}()

    ccall((:mb01os_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), uplo, trans, n, h, ldh, x, ldx, p, ldp,
            info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ot!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, t::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    lde = max(1,stride(e,2))
    ldt = max(1,stride(t,2))
    info = Ref{BlasInt}()

    ccall((:mb01ot_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Clong, Clong), uplo, trans, n, alpha,
            beta, r, ldr, e, lde, t, ldt, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01pd!(scun::AbstractChar, type::AbstractChar, m::Integer,
    n::Integer, kl::Integer, ku::Integer, anrm::Number, nbl::Integer,
    nrows::Integer, a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01pd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), scun, type,
            m, n, kl, ku, anrm, nbl, nrows, a, lda, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01qd!(type::AbstractChar, m::Integer, n::Integer,
    kl::Integer, ku::Integer, cfrom::Number, cto::Number,
    nbl::Integer, nrows::Integer, a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01qd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong), type, m, n, kl, ku,
            cfrom, cto, nbl, nrows, a, lda, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rb!(side::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, m::Integer, n::Integer, alpha::Number,
    beta::Number, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb01rb_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), side, uplo, trans, m, n, alpha,
            beta, r, ldr, a, lda, b, ldb, info, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rd!(uplo::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, alpha::Number, beta::Number,
    r::AbstractMatrix{Float64}, a::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, ldwork::Integer)

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01rd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), uplo, trans, m, n, alpha,
            beta, r, ldr, a, lda, x, ldx, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rh!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, x::AbstractMatrix{Float64},
    ldwork::Integer)

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01rh_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), uplo, trans, n, alpha, beta, r, ldr, h,
            ldh, x, ldx, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rt!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, x::AbstractMatrix{Float64},
    ldwork::Integer)

    ldr = max(1,stride(r,2))
    lde = max(1,stride(e,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01rt_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), uplo, trans, n, alpha, beta, r, ldr, e,
            lde, x, ldx, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ru!(uplo::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, alpha::Number, beta::Number,
    r::AbstractMatrix{Float64}, a::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, ldwork::Integer)

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldx = max(1,stride(x,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01ru_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), uplo, trans, m, n, alpha,
            beta, r, ldr, a, lda, x, ldx, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rw!(uplo::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, a::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb01rw_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong, Clong), uplo, trans, m, n, a, lda, z, ldz, dwork,
            info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01rx!(side::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, m::Integer, n::Integer, alpha::Number,
    beta::Number, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb01rx_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), side, uplo, trans, m, n, alpha,
            beta, r, ldr, a, lda, b, ldb, info, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ry!(side::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, m::Integer, alpha::Number, beta::Number,
    r::AbstractMatrix{Float64}, h::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, dwork::AbstractVector{Float64})

    ldr = max(1,stride(r,2))
    ldh = max(1,stride(h,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb01ry_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong, Clong, Clong), side, uplo, trans, m, alpha, beta,
            r, ldr, h, ldh, b, ldb, dwork, info, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb01sd!(jobs::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64}, r::AbstractVector{Float64},
    c::AbstractVector{Float64})

    lda = max(1,stride(a,2))

    ccall((:mb01sd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Clong), jobs, m, n, a, lda, r, c, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb01ss!(jobs::AbstractChar, uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, d::AbstractVector{Float64})

    lda = max(1,stride(a,2))

    ccall((:mb01ss_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Clong, Clong), jobs, uplo, n, a, lda, d, 1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01td!(n::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n-1)

    ccall((:mb01td_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}), n, a, lda, b, ldb, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ud!(side::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, alpha::Number, h::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64})

    ldh = max(1,stride(h,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb01ud_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), side, trans,
            m, n, alpha, h, ldh, a, lda, b, ldb, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01uw!(side::AbstractChar, trans::AbstractChar, m::Integer,
    n::Integer, alpha::Number, h::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, ldwork::Integer)

    ldh = max(1,stride(h,2))
    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb01uw_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), side, trans,
            m, n, alpha, h, ldh, a, lda, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01ux!(side::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, m::Integer, n::Integer, alpha::Number,
    t::AbstractMatrix{Float64}, a::AbstractMatrix{Float64})

    ldt = max(1,stride(t,2))
    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb01ux_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong), side, uplo, trans, m, n, alpha, t, ldt, a, lda,
            dwork, ldwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (mc, nc,  info)
"""
function mb01vd!(trana::AbstractChar, tranb::AbstractChar,
    ma::Integer, na::Integer, mb::Integer, nb::Integer, alpha::Number,
    beta::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, c::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    mc = Ref{BlasInt}()
    nc = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb01vd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            trana, tranb, ma, na, mb, nb, alpha, beta, a, lda, b,
            ldb, c, ldc, mc, nc, info, 1, 1)
    chkargsok(info[])

    return mc[], nc[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01wd!(dico::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, hess::AbstractChar, n::Integer,
    alpha::Number, beta::Number, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, t::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldt = max(1,stride(t,2))
    info = Ref{BlasInt}()

    ccall((:mb01wd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong, Clong), dico, uplo, trans, hess, n,
            alpha, beta, r, ldr, a, lda, t, ldt, info, 1, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01xd!(uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01xd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), uplo,
            n, a, lda, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01xy!(uplo::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb01xy_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), uplo,
            n, a, lda, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01yd!(uplo::AbstractChar, trans::AbstractChar, n::Integer,
    k::Integer, l::Integer, alpha::Number, beta::Number,
    a::AbstractMatrix{Float64}, c::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()

    ccall((:mb01yd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), uplo, trans,
            n, k, l, alpha, beta, a, lda, c, ldc, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb01zd!(side::AbstractChar, uplo::AbstractChar,
    transt::AbstractChar, diag::AbstractChar, m::Integer, n::Integer,
    l::Integer, alpha::Number, t::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64})

    ldt = max(1,stride(t,2))
    ldh = max(1,stride(h,2))
    info = Ref{BlasInt}()

    ccall((:mb01zd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong), side, uplo, transt, diag, m, n, l, alpha,
            t, ldt, h, ldh, info, 1, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02cd!(job::AbstractChar, typet::AbstractChar, k::Integer,
    n::Integer, t::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64}, r::AbstractMatrix{Float64},
    l::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    lcs::Integer, ldwork::Integer)

    ldt = max(1,stride(t,2))
    ldg = max(1,stride(g,2))
    ldr = max(1,stride(r,2))
    ldl = max(1,stride(l,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02cd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            job, typet, k, n, t, ldt, g, ldg, r, ldr, l, ldl, cs,
            lcs, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rnk,  info)
"""
function mb02cu!(typeg::AbstractChar, k::Integer, p::Integer,
    q::Integer, nb::Integer, a1::AbstractMatrix{Float64},
    a2::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    ipvt::AbstractVector{BlasInt}, cs::AbstractVector{Float64},
    tol::Number, ldwork::Integer)

    lda1 = max(1,stride(a1,2))
    lda2 = max(1,stride(a2,2))
    ldb = max(1,stride(b,2))
    rnk = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02cu_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), typeg, k, p, q, nb, a1, lda1, a2, lda2, b, ldb,
            rnk, ipvt, cs, tol, dwork, ldwork, info, 1)
    chkargsok(info[])

    return rnk[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02cv!(typeg::AbstractChar, strucg::AbstractChar,
    k::Integer, n::Integer, p::Integer, q::Integer, nb::Integer,
    rnk::Integer, a1::AbstractMatrix{Float64},
    a2::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    f1::AbstractMatrix{Float64}, f2::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    ldwork::Integer)

    lda1 = max(1,stride(a1,2))
    lda2 = max(1,stride(a2,2))
    ldb = max(1,stride(b,2))
    ldf1 = max(1,stride(f1,2))
    ldf2 = max(1,stride(f2,2))
    ldg = max(1,stride(g,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02cv_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), typeg,
            strucg, k, n, p, q, nb, rnk, a1, lda1, a2, lda2, b, ldb,
            f1, ldf1, f2, ldf2, g, ldg, cs, dwork, ldwork, info, 1,
            1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02cx!(typet::AbstractChar, p::Integer, q::Integer,
    k::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    lcs::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02cx_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), typet,
            p, q, k, a, lda, b, ldb, cs, lcs, dwork, ldwork, info,
            1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02cy!(typet::AbstractChar, strucg::AbstractChar,
    p::Integer, q::Integer, n::Integer, k::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    lcs::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldh = max(1,stride(h,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02cy_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            typet, strucg, p, q, n, k, a, lda, b, ldb, h, ldh, cs,
            lcs, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02dd!(job::AbstractChar, typet::AbstractChar, k::Integer,
    m::Integer, n::Integer, ta::AbstractMatrix{Float64},
    t::AbstractMatrix{Float64}, g::AbstractMatrix{Float64},
    r::AbstractMatrix{Float64}, l::AbstractMatrix{Float64},
    cs::AbstractVector{Float64}, ldwork::Integer)

    ldta = max(1,stride(ta,2))
    ldt = max(1,stride(t,2))
    ldg = max(1,stride(g,2))
    ldr = max(1,stride(r,2))
    ldl = max(1,stride(l,2))
    lcs = length(cs)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02dd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), job, typet,
            k, m, n, ta, ldta, t, ldt, g, ldg, r, ldr, l, ldl, cs,
            lcs, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02ed!(typet::AbstractChar, k::Integer, n::Integer,
    nrhs::Integer, t::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, ldwork::Integer)

    ldt = max(1,stride(t,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02ed_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), typet, k, n, nrhs, t, ldt, b, ldb,
            dwork, ldwork, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02fd!(typet::AbstractChar, k::Integer, n::Integer,
    p::Integer, s::Integer, t::AbstractMatrix{Float64},
    r::AbstractMatrix{Float64}, ldwork::Integer)

    ldt = max(1,stride(t,2))
    ldr = max(1,stride(r,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02fd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong), typet, k, n, p, s,
            t, ldt, r, ldr, dwork, ldwork, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02gd!(typet::AbstractChar, triu::AbstractChar, k::Integer,
    n::Integer, nl::Integer, p::Integer, s::Integer,
    t::AbstractMatrix{Float64}, rb::AbstractMatrix{Float64})

    ldt = max(1,stride(t,2))
    ldrb = max(1,stride(rb,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02gd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), typet, triu, k, n, nl, p, s, t, ldt, rb,
            ldrb, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02hd!(triu::AbstractChar, k::Integer, l::Integer,
    m::Integer, ml::Integer, n::Integer, nu::Integer, p::Integer,
    s::Integer, tc::AbstractMatrix{Float64},
    tr::AbstractMatrix{Float64}, rb::AbstractMatrix{Float64})

    ldtc = max(1,stride(tc,2))
    ldtr = max(1,stride(tr,2))
    ldrb = max(1,stride(rb,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02hd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), triu, k, l, m, ml, n, nu, p, s, tc, ldtc, tr,
            ldtr, rb, ldrb, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02id!(job::AbstractChar, k::Integer, l::Integer,
    m::Integer, n::Integer, rb::Integer, rc::Integer,
    tc::AbstractMatrix{Float64}, tr::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, c::AbstractMatrix{Float64})

    ldtc = max(1,stride(tc,2))
    ldtr = max(1,stride(tr,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02id_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), job, k, l, m, n, rb, rc, tc, ldtc, tr, ldtr, b,
            ldb, c, ldc, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02jd!(job::AbstractChar, k::Integer, l::Integer,
    m::Integer, n::Integer, p::Integer, s::Integer,
    tc::AbstractMatrix{Float64}, tr::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, r::AbstractMatrix{Float64})

    ldtc = max(1,stride(tc,2))
    ldtr = max(1,stride(tr,2))
    ldq = max(1,stride(q,2))
    ldr = max(1,stride(r,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02jd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), job, k, l, m, n, p, s, tc, ldtc, tr, ldtr, q,
            ldq, r, ldr, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rnk,  info)
"""
function mb02jx!(job::AbstractChar, k::Integer, l::Integer,
    m::Integer, n::Integer, tc::AbstractMatrix{Float64},
    tr::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    r::AbstractMatrix{Float64}, jpvt::AbstractVector{BlasInt},
    tol1::Number, tol2::Number, ldwork::Integer)

    ldtc = max(1,stride(tc,2))
    ldtr = max(1,stride(tr,2))
    ldq = max(1,stride(q,2))
    ldr = max(1,stride(r,2))
    rnk = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02jx_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong), job, k, l, m, n, tc,
            ldtc, tr, ldtr, rnk, q, ldq, r, ldr, jpvt, tol1, tol2,
            dwork, ldwork, info, 1)
    chkargsok(info[])

    return rnk[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02kd!(ldblk::AbstractChar, trans::AbstractChar, k::Integer,
    l::Integer, m::Integer, n::Integer, r::Integer, alpha::Number,
    beta::Number, tc::AbstractMatrix{Float64},
    tr::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64})

    ldtc = max(1,stride(tc,2))
    ldtr = max(1,stride(tr,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02kd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), ldblk, trans,
            k, l, m, n, r, alpha, beta, tc, ldtc, tr, ldtr, b, ldb,
            c, ldc, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank, info,  iwarn)
"""
function mb02md!(job::AbstractChar, m::Integer, n::Integer,
    l::Integer, ini_rank::Integer, c::AbstractMatrix{Float64},
    s::AbstractVector{Float64}, x::AbstractMatrix{Float64},
    tol::Number)

    ldc = max(1,stride(c,2))
    ldx = max(1,stride(x,2))
    rank = Ref{BlasInt}(ini_rank)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, l)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02md_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong), job, m, n, l, rank,
            c, ldc, s, x, ldx, tol, iwork, dwork, ldwork, iwarn,
            info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return rank[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (rank, theta, info,  iwarn)
"""
function mb02nd!(m::Integer, n::Integer, l::Integer,
    ini_rank::Integer, ini_theta::Number, c::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, q::AbstractVector{Float64},
    inul::AbstractVector{BlasBool}, tol::Number, reltol::Number)

    ldc = max(1,stride(c,2))
    ldx = max(1,stride(x,2))
    rank = Ref{BlasInt}(ini_rank)
    theta = Ref{Float64}(ini_theta)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n+2*l)
    bwork = Vector{BlasBool}(undef, n+l)
    ldwork = BlasInt(-1)
    # some of this is used even in the workspace query
    dwork = Vector{Float64}(undef, 64)

    local jlres
    for iwq in 1:2
        ccall((:mb02nd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{BlasBool}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt},
            Ptr{BlasInt}), m, n, l, rank, theta, c, ldc, x, ldx, q,
            inul, tol, reltol, iwork, dwork, ldwork, bwork, iwarn,
            info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return rank[], theta[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb02ny!(updatu::Bool, updatv::Bool, m::Integer, n::Integer,
    i::Integer, k::Integer, q::AbstractVector{Float64},
    e::AbstractVector{Float64}, u::AbstractMatrix{Float64},
    v::AbstractMatrix{Float64})

    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    dwork = Vector{Float64}(undef, max(1,ldwork))

    ccall((:mb02ny_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}), updatu,
            updatv, m, n, i, k, q, e, u, ldu, v, ldv, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (rcond,  info)
"""
function mb02od!(side::AbstractChar, uplo::AbstractChar,
    trans::AbstractChar, diag::AbstractChar, norm::AbstractChar,
    m::Integer, n::Integer, alpha::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, tol::Number)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    rcond = Ref{Float64}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, k)
    dwork = Vector{Float64}(undef, 3*k)

    ccall((:mb02od_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong, Clong), side, uplo, trans, diag, norm, m,
            n, alpha, a, lda, b, ldb, rcond, tol, iwork, dwork,
            info, 1, 1, 1, 1, 1)
    chkargsok(info[])

    return rcond[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rcond,  info)
"""
function mb02pd!(fact::AbstractChar, trans::AbstractChar, n::Integer,
    nrhs::Integer, a::AbstractMatrix{Float64},
    af::AbstractMatrix{Float64}, ipiv::AbstractVector{BlasInt},
    equed::AbstractChar, r::AbstractVector{Float64},
    c::AbstractVector{Float64}, b::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, ferr::AbstractVector{Float64},
    berr::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldaf = max(1,stride(af,2))
    ldb = max(1,stride(b,2))
    ldx = max(1,stride(x,2))
    rcond = Ref{Float64}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)
    dwork = Vector{Float64}(undef, 4*n)

    ccall((:mb02pd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af,
            ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr,
            berr, iwork, dwork, info, 1, 1, 1)
    chkargsok(info[])

    return rcond[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb02qd!(job::AbstractChar, iniper::AbstractChar, m::Integer,
    n::Integer, nrhs::Integer, rcond::Number, svlmax::Number,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64}, y::AbstractVector{Float64},
    jpvt::AbstractVector{BlasInt}, sval::AbstractVector{Float64},
    ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02qd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), job, iniper, m, n, nrhs, rcond, svlmax,
            a, lda, b, ldb, y, jpvt, rank, sval, dwork, ldwork,
            info, 1, 1)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02qy!(m::Integer, n::Integer, nrhs::Integer, rank::Integer,
    a::AbstractMatrix{Float64},
    jpvt::AbstractVector{BlasInt}, b::AbstractMatrix{Float64},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(a,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  ldwork )
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02qy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), m, n, nrhs,
            rank, a, lda, jpvt, b, ldb, tau, dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02rd!(trans::AbstractChar, n::Integer, nrhs::Integer,
    h::AbstractMatrix{Float64}, ipiv::AbstractVector{BlasInt},
    b::AbstractMatrix{Float64})

    ldh = max(1,stride(h,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb02rd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), trans,
            n, nrhs, h, ldh, ipiv, b, ldb, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02rz!(trans::AbstractChar, n::Integer, nrhs::Integer,
    h::AbstractMatrix{ComplexF64}, ipiv::AbstractVector{BlasInt},
    b::AbstractMatrix{ComplexF64})

    ldh = max(1,stride(h,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb02rz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), trans, n, nrhs, h, ldh, ipiv, b,
            ldb, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02sd!(n::Integer, h::AbstractMatrix{Float64},
    ipiv::AbstractVector{BlasInt})

    ldh = max(1,stride(h,2))
    info = Ref{BlasInt}()

    ccall((:mb02sd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}), n, h, ldh,
            ipiv, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02sz!(n::Integer, h::AbstractMatrix{ComplexF64},
    ipiv::AbstractVector{BlasInt})

    ldh = max(1,stride(h,2))
    info = Ref{BlasInt}()

    ccall((:mb02sz_, libslicot), Cvoid, (Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}), n, h, ldh, ipiv, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rcond,  info)
"""
function mb02td!(norm::AbstractChar, n::Integer, hnorm::Number,
    h::AbstractMatrix{Float64}, ipiv::AbstractVector{BlasInt})

    ldh = max(1,stride(h,2))
    rcond = Ref{Float64}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)
    dwork = Vector{Float64}(undef, 3*n)

    ccall((:mb02td_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong), norm, n, hnorm, h, ldh, ipiv, rcond, iwork,
            dwork, info, 1)
    chkargsok(info[])

    return rcond[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rcond,  info)
"""
function mb02tz!(norm::AbstractChar, n::Integer, hnorm::Number,
    h::AbstractMatrix{ComplexF64}, ipiv::AbstractVector{BlasInt})

    ldh = max(1,stride(h,2))
    rcond = Ref{Float64}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)
    zwork = Vector{ComplexF64}(undef, 2*n)

    ccall((:mb02tz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{ComplexF64}, Ptr{BlasInt}, Clong), norm, n, hnorm,
            h, ldh, ipiv, rcond, dwork, zwork, info, 1)
    chkargsok(info[])

    return rcond[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02ud!(fact::AbstractChar, side::AbstractChar,
    trans::AbstractChar, jobp::AbstractChar, m::Integer, n::Integer,
    alpha::Number, rcond::Number, rank::Integer,
    r::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    sv::AbstractVector{Float64}, b::AbstractMatrix{Float64},
    rp::AbstractMatrix{Float64})

    ldr = max(1,stride(r,2))
    ldq = max(1,stride(q,2))
    ldb = max(1,stride(b,2))
    ldrp = max(1,stride(rp,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb02ud_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong), fact, side, trans, jobp, m, n, alpha,
            rcond, rank, r, ldr, q, ldq, sv, b, ldb, rp, ldrp,
            dwork, ldwork, info, 1, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns scale
"""
function mb02uu!(n::Integer, a::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64}, ipiv::AbstractVector{BlasInt},
    jpiv::AbstractVector{BlasInt})

    lda = max(1,stride(a,2))
    scale = Ref{Float64}()

    ccall((:mb02uu_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}), n, a, lda, rhs, ipiv, jpiv, scale)

    return scale[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02uv!(n::Integer, a::AbstractMatrix{Float64},
    ipiv::AbstractVector{BlasInt}, jpiv::AbstractVector{BlasInt})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb02uv_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
            n, a, lda, ipiv, jpiv, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (scale,  iwarn)
"""
function mb02uw!(ltrans::Bool, n::Integer, m::Integer,
    par::AbstractVector{Float64}, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    scale = Ref{Float64}()
    iwarn = Ref{BlasInt}()

    ccall((:mb02uw_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}),
            ltrans, n, m, par, a, lda, b, ldb, scale, iwarn)

    return scale[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02vd!(trans::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64}, ipiv::AbstractVector{BlasInt},
    b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()

    ccall((:mb02vd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), trans,
            m, n, a, lda, ipiv, b, ldb, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb02yd!(cond::AbstractChar, n::Integer,
    r::AbstractMatrix{Float64}, ipvt::AbstractVector{BlasInt},
    diag::AbstractVector{Float64}, qtb::AbstractVector{Float64},
    rank::Integer, x::AbstractVector{Float64}, tol::Number,
    ldwork::Integer)

    ldr = max(1,stride(r,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb02yd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong), cond,
            n, r, ldr, ipvt, diag, qtb, rank, x, tol, dwork, ldwork,
            info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ab!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3},
    w1::Number, w2::Number)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()

    ccall((:mb03ab_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Clong), shft, k, n, amap, s, sinv, a,
            lda1, lda2, w1, w2, c1, s1, c2, s2, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ad!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()

    ccall((:mb03ad_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Clong), shft,
            k, n, amap, s, sinv, a, lda1, lda2, c1, s1, c2, s2, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ae!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()

    ccall((:mb03ae_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Clong), shft,
            k, n, amap, s, sinv, a, lda1, lda2, c1, s1, c2, s2, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03af!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()

    ccall((:mb03af_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Clong), shft,
            k, n, amap, s, sinv, a, lda1, lda2, c1, s1, c2, s2, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ag!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()
    iwork = Vector{BlasInt}(undef, 2*n)
    dwork = Vector{Float64}(undef, 2*n*n)

    ccall((:mb03ag_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Clong), shft, k, n, amap, s, sinv, a,
            lda1, lda2, c1, s1, c2, s2, iwork, dwork, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ah!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()

    ccall((:mb03ah_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Clong), shft,
            k, n, amap, s, sinv, a, lda1, lda2, c1, s1, c2, s2, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns (c1, s1, c2,  s2)
"""
function mb03ai!(shft::AbstractChar, k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    c1 = Ref{Float64}()
    s1 = Ref{Float64}()
    c2 = Ref{Float64}()
    s2 = Ref{Float64}()
    dwork = Vector{Float64}(undef, n*(n+2))

    ccall((:mb03ai_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Clong), shft, k, n, amap, s, sinv, a, lda1, lda2, c1,
            s1, c2, s2, dwork, 1)

    return c1[], s1[], c2[], s2[]
end


"""
$(TYPEDSIGNATURES)
returns smult
"""
function mb03ba!(k::Integer, h::Integer, s::AbstractVector{BlasInt},
    amap::AbstractVector{BlasInt}, qmap::AbstractVector{BlasInt})

    smult = Ref{BlasInt}()

    ccall((:mb03ba_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
            k, h, s, smult, amap, qmap)

    return smult[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03bb!(base::Number, lgbas::Number, ulp::Number, k::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, scal::AbstractVector{BlasInt})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, 8*k)

    ccall((:mb03bb_, libslicot), Cvoid, (Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}), base, lgbas, ulp, k, amap,
            s, sinv, a, lda1, lda2, alphar, alphai, beta, scal,
            dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb03bc!(k::Integer, amap::AbstractVector{BlasInt},
    s::AbstractVector{BlasInt}, sinv::Integer, a::Array{Float64,3},
    macpar::AbstractVector{Float64},
    cv::AbstractVector{Float64}, sv::AbstractVector{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    dwork = Vector{Float64}(undef, 3*(k-1))

    ccall((:mb03bc_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}), k, amap, s, sinv, a, lda1, lda2, macpar,
            cv, sv, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (info,  iwarn)
"""
function mb03bd!(job::AbstractChar, defl::AbstractChar,
    compq::AbstractChar, qind::AbstractVector{BlasInt}, k::Integer,
    n::Integer, h::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{Float64,3},
    q::Array{Float64,3},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, scal::AbstractVector{BlasInt},
    liwork::Integer, ldwork::Integer)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03bd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), job, defl, compq, qind, k, n, h,
            ilo, ihi, s, a, lda1, lda2, q, ldq1, ldq2, alphar,
            alphai, beta, scal, iwork, liwork, dwork, ldwork, iwarn,
            info, 1, 1, 1)
    chkargsok(info[])

    return info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb03be!(k::Integer, amap::AbstractVector{BlasInt},
    s::AbstractVector{BlasInt}, sinv::Integer, a::Array{Float64,3})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)

    ccall((:mb03be_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}), k, amap, s, sinv, a, lda1, lda2)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb03bf!(k::Integer, amap::AbstractVector{BlasInt},
    s::AbstractVector{BlasInt}, sinv::Integer, a::Array{Float64,3},
    ulp::Number)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)

    ccall((:mb03bf_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}), k, amap, s, sinv, a, lda1,
            lda2, ulp)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03bg!(k::Integer, n::Integer,
    amap::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    sinv::Integer, a::Array{Float64,3},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    info = Ref{BlasInt}()

    ccall((:mb03bg_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}),
            k, n, amap, s, sinv, a, lda1, lda2, wr, wi)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03bz!(job::AbstractChar, compq::AbstractChar, k::Integer,
    n::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{ComplexF64,3}, q::Array{ComplexF64,3},
    alpha::AbstractVector{ComplexF64},
    beta::AbstractVector{ComplexF64}, scal::AbstractVector{BlasInt},
    ldwork::Integer, lzwork::Integer)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)
    zwork = Vector{ComplexF64}(undef, lzwork)

    ccall((:mb03bz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Clong,
            Clong), job, compq, k, n, ilo, ihi, s, a, lda1, lda2, q,
            ldq1, ldq2, alpha, beta, scal, dwork, ldwork, zwork,
            lzwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (n1, n2,  info)
"""
function mb03cd!(uplo::AbstractChar, ini_n1::Integer, ini_n2::Integer,
    prec::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    q1::AbstractMatrix{Float64}, q2::AbstractMatrix{Float64},
    q3::AbstractMatrix{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldd = max(1,stride(d,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldq3 = max(1,stride(q3,2))
    n1 = Ref{BlasInt}(ini_n1)
    n2 = Ref{BlasInt}(ini_n2)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03cd_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ptr{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), uplo, n1, n2, prec, a, lda, b,
            ldb, d, ldd, q1, ldq1, q2, ldq2, q3, ldq3, dwork,
            ldwork, info, 1)
    chkargsok(info[])

    return n1[], n2[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (co1, si1, co2, si2, co3,  si3)
"""
function mb03cz!(a::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, d::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldd = max(1,stride(d,2))
    co1 = Ref{Float64}()
    si1 = Ref{ComplexF64}()
    co2 = Ref{Float64}()
    si2 = Ref{ComplexF64}()
    co3 = Ref{Float64}()
    si3 = Ref{ComplexF64}()

    ccall((:mb03cz_, libslicot), Cvoid, (Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64},
            Ptr{Float64}, Ptr{ComplexF64}), a, lda, b, ldb, d, ldd,
            co1, si1, co2, si2, co3, si3)

    return co1[], si1[], co2[], si2[], co3[], si3[]
end


"""
$(TYPEDSIGNATURES)
returns (n1, n2,  info)
"""
function mb03dd!(uplo::AbstractChar, ini_n1::Integer, ini_n2::Integer,
    prec::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    n1 = Ref{BlasInt}(ini_n1)
    n2 = Ref{BlasInt}(ini_n2)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03dd_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ptr{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), uplo, n1, n2, prec, a, lda, b,
            ldb, q1, ldq1, q2, ldq2, dwork, ldwork, info, 1)
    chkargsok(info[])

    return n1[], n2[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (co1, si1, co2,  si2)
"""
function mb03dz!(a::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    co1 = Ref{Float64}()
    si1 = Ref{ComplexF64}()
    co2 = Ref{Float64}()
    si2 = Ref{ComplexF64}()

    ccall((:mb03dz_, libslicot), Cvoid, (Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
            Ptr{ComplexF64}), a, lda, b, ldb, co1, si1, co2, si2)

    return co1[], si1[], co2[], si2[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03ed!(n::Integer, prec::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    q1::AbstractMatrix{Float64}, q2::AbstractMatrix{Float64},
    q3::AbstractMatrix{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldd = max(1,stride(d,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldq3 = max(1,stride(q3,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03ed_, libslicot), Cvoid, (Ref{BlasInt}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), n, prec, a,
            lda, b, ldb, d, ldd, q1, ldq1, q2, ldq2, q3, ldq3,
            dwork, ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03fd!(n::Integer, prec::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03fd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), n, prec, a,
            lda, b, ldb, q1, ldq1, q2, ldq2, dwork, ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03fz!(compq::AbstractChar, compu::AbstractChar,
    orth::AbstractChar, n::Integer, z::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, fg::AbstractMatrix{ComplexF64},
    d::AbstractMatrix{ComplexF64}, c::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, u::AbstractMatrix{ComplexF64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer,
    bwork::AbstractVector{BlasBool})

    ldz = max(1,stride(z,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldd = max(1,stride(d,2))
    ldc = max(1,stride(c,2))
    ldq = max(1,stride(q,2))
    ldu = max(1,stride(u,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03fz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt},
            Clong, Clong, Clong), compq, compu, orth, n, z, ldz, b,
            ldb, fg, ldfg, neig, d, ldd, c, ldc, q, ldq, u, ldu,
            alphar, alphai, beta, iwork, liwork, dwork, ldwork,
            zwork, lzwork, bwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03gd!(n::Integer, b::AbstractMatrix{Float64},
    d::AbstractMatrix{Float64}, macpar::AbstractVector{Float64},
    q::AbstractMatrix{Float64}, u::AbstractMatrix{Float64},
    ldwork::Integer)

    ldb = max(1,stride(b,2))
    ldd = max(1,stride(d,2))
    ldq = max(1,stride(q,2))
    ldu = max(1,stride(u,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03gd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), n, b, ldb, d,
            ldd, macpar, q, ldq, u, ldu, dwork, ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (co1, si1, co2,  si2)
"""
function mb03gz!(z11::Complex, z12::Complex, z22::Complex,
    h11::Complex, h12::Complex)

    co1 = Ref{Float64}()
    si1 = Ref{ComplexF64}()
    co2 = Ref{Float64}()
    si2 = Ref{ComplexF64}()

    ccall((:mb03gz_, libslicot), Cvoid, (Ref{ComplexF64},
            Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
            Ref{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64},
            Ptr{Float64}, Ptr{ComplexF64}), z11, z12, z22, h11, h12,
            co1, si1, co2, si2)

    return co1[], si1[], co2[], si2[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03hd!(n::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, macpar::AbstractVector{Float64},
    q::AbstractMatrix{Float64}, dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()

    ccall((:mb03hd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}),
            n, a, lda, b, ldb, macpar, q, ldq, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (co,  si)
"""
function mb03hz!(s11::Complex, s12::Complex, h11::Complex,
    h12::Complex)

    co = Ref{Float64}()
    si = Ref{ComplexF64}()

    ccall((:mb03hz_, libslicot), Cvoid, (Ref{ComplexF64},
            Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
            Ptr{Float64}, Ptr{ComplexF64}), s11, s12, h11, h12, co,
            si)

    return co[], si[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03id!(compq::AbstractChar, compu::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, c::AbstractMatrix{Float64},
    d::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    f::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    u1::AbstractMatrix{Float64}, u2::AbstractMatrix{Float64},
    liwork::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldc = max(1,stride(c,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03id_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            compq, compu, n, a, lda, c, ldc, d, ldd, b, ldb, f, ldf,
            q, ldq, u1, ldu1, u2, ldu2, neig, iwork, liwork, dwork,
            ldwork, info, 1, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03iz!(compq::AbstractChar, compu::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, c::AbstractMatrix{ComplexF64},
    d::AbstractMatrix{ComplexF64}, b::AbstractMatrix{ComplexF64},
    f::AbstractMatrix{ComplexF64}, q::AbstractMatrix{ComplexF64},
    u1::AbstractMatrix{ComplexF64}, u2::AbstractMatrix{ComplexF64},
    tol::Number)

    lda = max(1,stride(a,2))
    ldc = max(1,stride(c,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb03iz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
            Ref{Float64}, Ptr{BlasInt}, Clong, Clong), compq, compu,
            n, a, lda, c, ldc, d, ldd, b, ldb, f, ldf, q, ldq, u1,
            ldu1, u2, ldu2, neig, tol, info, 1, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03jd!(compq::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, f::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, liwork::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03jd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), compq, n, a, lda, d, ldd, b, ldb, f, ldf, q,
            ldq, neig, iwork, liwork, dwork, ldwork, info, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03jp!(compq::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, f::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, liwork::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03jp_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), compq, n, a, lda, d, ldd, b, ldb, f, ldf, q,
            ldq, neig, iwork, liwork, dwork, ldwork, info, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03jz!(compq::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, d::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, f::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, tol::Number)

    lda = max(1,stride(a,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb03jz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ptr{BlasInt},
            Clong), compq, n, a, lda, d, ldd, b, ldb, f, ldf, q,
            ldq, neig, tol, info, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (ifst, ilst,  info)
"""
function mb03ka!(compq::AbstractChar, whichq::AbstractVector{BlasInt},
    ws::Bool, k::Integer, nc::Integer, kschur::Integer,
    ini_ifst::Integer, ini_ilst::Integer, n::AbstractVector{BlasInt},
    ni::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    t::AbstractVector{Float64}, ldt::AbstractVector{BlasInt},
    ixt::AbstractVector{BlasInt}, q::AbstractVector{Float64},
    ldq::AbstractVector{BlasInt}, ixq::AbstractVector{BlasInt},
    tol::AbstractVector{Float64})
    # NOTE: intentionally strange signature: `t,q` are vector holding 3d array

    ifst = Ref{BlasInt}(ini_ifst)
    ilst = Ref{BlasInt}(ini_ilst)
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, 4*k)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03ka_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ref{BlasBool}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), compq, whichq, ws, k, nc, kschur, ifst, ilst, n,
            ni, s, t, ldt, ixt, q, ldq, ixq, tol, iwork, dwork,
            ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return ifst[], ilst[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03kb!(compq::AbstractChar, whichq::AbstractVector{BlasInt},
    ws::Bool, k::Integer, nc::Integer, kschur::Integer, j1::Integer,
    n1::Integer, n2::Integer, n::AbstractVector{BlasInt},
    ni::AbstractVector{BlasInt}, s::AbstractVector{BlasInt},
    t::AbstractVector{Float64}, ldt::AbstractVector{BlasInt},
    ixt::AbstractVector{BlasInt}, q::AbstractVector{Float64},
    ldq::AbstractVector{BlasInt}, ixq::AbstractVector{BlasInt},
    tol::AbstractVector{Float64})
    # NOTE: intentionally strange signature: `t,q` are vector holding 3d array

    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, 4*k)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03kb_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ref{BlasBool}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), compq, whichq, ws, k, nc, kschur,
            j1, n1, n2, n, ni, s, t, ldt, ixt, q, ldq, ixq, tol,
            iwork, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb03kc!(k::Integer, khess::Integer, n::Integer, r::Integer,
    s::AbstractVector{BlasInt}, a::AbstractVector{Float64},
    lda::Integer, v::AbstractVector{Float64},
    tau::AbstractVector{Float64})
    # NOTE: intentionally strange signature: `a` is vector holding 3d array

    ccall((:mb03kc_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}), k, khess, n,
            r, s, a, lda, v, tau)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (m,  info)
"""
function mb03kd!(compq::AbstractChar, whichq::AbstractVector{BlasInt},
    strong::AbstractChar, k::Integer, nc::Integer, kschur::Integer,
    n::AbstractVector{BlasInt}, ni::AbstractVector{BlasInt},
    s::AbstractVector{BlasInt}, select::AbstractVector{BlasBool},
    t::AbstractVector{Float64}, ldt::AbstractVector{BlasInt},
    ixt::AbstractVector{BlasInt}, q::AbstractVector{Float64},
    ldq::AbstractVector{BlasInt}, ixq::AbstractVector{BlasInt},
    tol::Number)

    # NOTE: intentionally strange signature: `t,q` are vector holding 3d array
    m = Ref{BlasInt}()
    info = Ref{BlasInt}(0)
    iwork = Vector{BlasInt}(undef, 4*k)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03kd_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasBool},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), compq, whichq, strong, k, nc, kschur, n,
            ni, s, select, t, ldt, ixt, q, ldq, ixq, m, tol, iwork,
            dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return m[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (scale,  info)
"""
function mb03ke!(trana::Bool, tranb::Bool, isgn::Integer, k::Integer,
    m::Integer, n::Integer, prec::Number, smin::Number,
    s::AbstractVector{BlasInt}, a::AbstractVector{Float64},
    b::AbstractVector{Float64}, c::AbstractVector{Float64})

    scale = Ref{Float64}()
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03ke_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), trana, tranb, isgn, k, m,
            n, prec, smin, s, a, b, c, scale, dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return scale[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03ld!(compq::AbstractChar, orth::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    liwork::Integer)

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    bwork = Vector{BlasBool}(undef, n÷2)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03ld_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool},
            Ptr{BlasInt}, Clong, Clong), compq, orth, n, a, lda, de,
            ldde, b, ldb, fg, ldfg, neig, q, ldq, alphar, alphai,
            beta, iwork, liwork, dwork, ldwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig, info,  iwarn)
"""
function mb03lf!(compq::AbstractChar, compu::AbstractChar,
    orth::AbstractChar, n::Integer, z::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, u::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer)

    ldz = max(1,stride(z,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    ldu = max(1,stride(u,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    bwork = Vector{BlasBool}(undef, n÷2)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03lf_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasBool}, Ptr{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong), compq, compu, orth, n, z, ldz, b, ldb, fg, ldfg,
            neig, q, ldq, u, ldu, alphar, alphai, beta, iwork,
            liwork, dwork, ldwork, bwork, iwarn, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return neig[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03lp!(compq::AbstractChar, orth::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    liwork::Integer)

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    bwork = Vector{BlasBool}(undef, n÷2)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03lp_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool},
            Ptr{BlasInt}, Clong, Clong), compq, orth, n, a, lda, de,
            ldde, b, ldb, fg, ldfg, neig, q, ldq, alphar, alphai,
            beta, iwork, liwork, dwork, ldwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb03lz!(compq::AbstractChar, orth::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, de::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, fg::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    bwork::AbstractVector{BlasBool})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n+1)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03lz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt}, Clong, Clong),
            compq, orth, n, a, lda, de, ldde, b, ldb, fg, ldfg,
            neig, q, ldq, alphar, alphai, beta, iwork, dwork,
            ldwork, zwork, lzwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (l, theta, info,  iwarn)
"""
function mb03md!(n::Integer, ini_l::Integer, ini_theta::Number,
    q::AbstractVector{Float64}, e::AbstractVector{Float64},
    q2::AbstractVector{Float64}, e2::AbstractVector{Float64},
    pivmin::Number, tol::Number, reltol::Number)

    l = Ref{BlasInt}(ini_l)
    theta = Ref{Float64}(ini_theta)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mb03md_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}), n, l, theta, q, e, q2, e2,
            pivmin, tol, reltol, iwarn, info)
    chkargsok(info[])

    return l[], theta[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb03my!(nx::Integer, x::AbstractVector{Float64},
    incx::Integer)


    jlres = ccall((:mb03my_, libslicot), Float64, (Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}), nx, x, incx)

    return jlres
end


"""
$(TYPEDSIGNATURES)
returns (result,  info)
"""
function mb03nd!(n::Integer, theta::Number,
    q2::AbstractVector{Float64}, e2::AbstractVector{Float64},
    pivmin::Number)

    info = Ref{BlasInt}()

    jlres = ccall((:mb03nd_, libslicot), BlasInt, (Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
            Ptr{BlasInt}), n, theta, q2, e2, pivmin, info)
    chkargsok(info[])

    return jlres, info[]
end


"""
$(TYPEDSIGNATURES)
returns (result,  info)
"""
function mb03ny!(n::Integer, omega::Number,
    a::AbstractMatrix{Float64},
    s::AbstractVector{Float64}, ldwork::Integer, lcwork::Integer)

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)
    cwork = Vector{ComplexF64}(undef, lcwork)

    jlres = ccall((:mb03ny_, libslicot), Float64, (Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasInt}), n, omega, a, lda, s, dwork,
            ldwork, cwork, lcwork, info)
    chkargsok(info[])

    return jlres, info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb03od!(jobqr::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64},
    jpvt::AbstractVector{BlasInt}, rcond::Number, svlmax::Number,
    tau::AbstractVector{Float64}, sval::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03od_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), jobqr, m, n, a, lda, jpvt, rcond, svlmax, tau,
            rank, sval, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb03oy!(m::Integer, n::Integer, a::AbstractMatrix{Float64},
    rcond::Number, svlmax::Number,
    sval::AbstractVector{Float64}, jpvt::AbstractVector{BlasInt},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  3*n-1 )

    ccall((:mb03oy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), m, n, a, lda, rcond,
            svlmax, rank, sval, jpvt, tau, dwork, info)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb03pd!(jobrq::AbstractChar, m::Integer, n::Integer,
    a::AbstractMatrix{Float64},
    jpvt::AbstractVector{BlasInt}, rcond::Number, svlmax::Number,
    tau::AbstractVector{Float64}, sval::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  ldwork )

    ccall((:mb03pd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Clong), jobrq,
            m, n, a, lda, jpvt, rcond, svlmax, tau, rank, sval,
            dwork, info, 1)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb03py!(m::Integer, n::Integer, a::AbstractMatrix{Float64},
    rcond::Number, svlmax::Number,
    sval::AbstractVector{Float64}, jpvt::AbstractVector{BlasInt},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  3*m-1 )

    ccall((:mb03py_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), m, n, a, lda, rcond,
            svlmax, rank, sval, jpvt, tau, dwork, info)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (ndim,  info)
"""
function mb03qd!(dico::AbstractChar, stdom::AbstractChar,
    jobu::AbstractChar, n::Integer, nlow::Integer, nsup::Integer,
    alpha::Number, a::AbstractMatrix{Float64},
    u::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldu = max(1,stride(u,2))
    ndim = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03qd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong, Clong, Clong), dico, stdom, jobu, n, nlow, nsup,
            alpha, a, lda, u, ldu, ndim, dwork, info, 1, 1, 1)
    chkargsok(info[])

    return ndim[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (ndim,  info)
"""
function mb03qg!(dico::AbstractChar, stdom::AbstractChar,
    jobu::AbstractChar, jobv::AbstractChar, n::Integer, nlow::Integer,
    nsup::Integer, alpha::Number, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, u::AbstractMatrix{Float64},
    v::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    ndim = Ref{BlasInt}()
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03qg_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong, Clong),
            dico, stdom, jobu, jobv, n, nlow, nsup, alpha, a, lda,
            e, lde, u, ldu, v, ldv, ndim, dwork, ldwork, info, 1, 1,
            1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return ndim[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03qv!(n::Integer, s::AbstractMatrix{Float64}, lds::Integer,
    t::AbstractMatrix{Float64}, ldt::Integer,
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64})

    lds = max(1,stride(s,2))
    ldt = max(1,stride(t,2))
    info = Ref{BlasInt}()

    ccall((:mb03qv_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}), n, s, lds, t,
            ldt, alphar, alphai, beta, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03qw!(n::Integer, l::Integer, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, u::AbstractMatrix{Float64},
    v::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    info = Ref{BlasInt}()

    ccall((:mb03qw_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            n, l, a, lda, e, lde, u, ldu, v, ldv, alphar, alphai,
            beta, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03qx!(n::Integer, t::AbstractMatrix{Float64},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64})

    ldt = max(1,stride(t,2))
    info = Ref{BlasInt}()

    ccall((:mb03qx_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            n, t, ldt, wr, wi, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03qy!(n::Integer, l::Integer, a::AbstractMatrix{Float64},
    u::AbstractMatrix{Float64}, e1::Number, e2::Number)

    lda = max(1,stride(a,2))
    ldu = max(1,stride(u,2))
    info = Ref{BlasInt}()

    ccall((:mb03qy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ptr{BlasInt}), n, l, a, lda,
            u, ldu, e1, e2, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (nblcks,  info)
"""
function mb03rd!(jobx::AbstractChar, sort::AbstractChar, n::Integer,
    pmax::Number, a::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, blsize::AbstractVector{BlasInt},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    tol::Number)

    lda = max(1,stride(a,2))
    ldx = max(1,stride(x,2))
    nblcks = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03rd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
            Ptr{BlasInt}, Clong, Clong), jobx, sort, n, pmax, a,
            lda, x, ldx, nblcks, blsize, wr, wi, tol, dwork, info,
            1, 1)
    chkargsok(info[])

    return nblcks[], info[]
end


"""
$(TYPEDSIGNATURES)
returns ku
"""
function mb03rx!(jobv::AbstractChar, n::Integer, kl::Integer,
    ini_ku::Integer, a::AbstractMatrix{Float64},
    x::AbstractMatrix{Float64}, wr::AbstractVector{Float64},
    wi::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldx = max(1,stride(x,2))
    ku = Ref{BlasInt}(ini_ku)
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03rx_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Clong), jobv, n, kl, ku, a, lda, x, ldx,
            wr, wi, dwork, 1)

    return ku[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03ry!(m::Integer, n::Integer, pmax::Number,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()

    ccall((:mb03ry_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            m, n, pmax, a, lda, b, ldb, c, ldc, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03sd!(jobscl::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    ldwork::Integer)

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03sd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), jobscl, n, a, lda, qg, ldqg, wr,
            wi, dwork, ldwork, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (m,  info)
"""
function mb03td!(typ::AbstractChar, compu::AbstractChar,
    select::AbstractVector{BlasBool}, lower::AbstractVector{BlasBool},
    n::Integer, a::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64}, u1::AbstractMatrix{Float64},
    u2::AbstractMatrix{Float64}, wr::AbstractVector{Float64},
    wi::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldg = max(1,stride(g,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    m = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03td_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ptr{BlasBool}, Ptr{BlasBool}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), typ, compu, select, lower,
            n, a, lda, g, ldg, u1, ldu1, u2, ldu2, wr, wi, m, dwork,
            ldwork, info, 1, 1)
    chkargsok(info[])

    return m[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03ts!(isham::Bool, wantu::Bool, n::Integer,
    a::AbstractMatrix{Float64}, g::AbstractMatrix{Float64},
    u1::AbstractMatrix{Float64}, u2::AbstractMatrix{Float64},
    j1::Integer, n1::Integer, n2::Integer)

    lda = max(1,stride(a,2))
    ldg = max(1,stride(g,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03ts_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}), isham, wantu, n, a, lda, g,
            ldg, u1, ldu1, u2, ldu2, j1, n1, n2, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03ud!(jobq::AbstractChar, jobp::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    sv::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03ud_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), jobq, jobp, n, a, lda, q,
            ldq, sv, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03vd!(n::Integer, p::Integer, ilo::Integer, ihi::Integer,
    a::Array{Float64,3}, tau::AbstractMatrix{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldtau = max(1,stride(tau,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03vd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}), n, p, ilo, ihi, a, lda1, lda2, tau,
            ldtau, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03vy!(n::Integer, p::Integer, ilo::Integer, ihi::Integer,
    a::Array{Float64,3}, tau::AbstractMatrix{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldtau = max(1,stride(tau,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03vy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), n, p, ilo, ihi, a, lda1,
            lda2, tau, ldtau, dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03wa!(wantq::Bool, wantz::Bool, n1::Integer, n2::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, z::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()

    ccall((:mb03wa_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), wantq, wantz,
            n1, n2, a, lda, b, ldb, q, ldq, z, ldz, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03wd!(job::AbstractChar, compz::AbstractChar, n::Integer,
    p::Integer, ilo::Integer, ihi::Integer, iloz::Integer,
    ihiz::Integer, h::Array{Float64,3}, z::Array{Float64,3},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    ldwork::Integer)

    ldh1 = max(1,stride(h,2))
    ldh2 = max(1,stride(h,3)÷ldh1)
    ldz1 = max(1,stride(z,2))
    ldz2 = max(1,stride(z,3)÷ldz1)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03wd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), job, compz, n, p, ilo, ihi,
            iloz, ihiz, h, ldh1, ldh2, z, ldz1, ldz2, wr, wi, dwork,
            ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03wx!(n::Integer, p::Integer, t::Array{Float64,3},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64})

    ldt1 = max(1,stride(t,2))
    ldt2 = max(1,stride(t,3)÷ldt1)
    info = Ref{BlasInt}()

    ccall((:mb03wx_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), n, p, t, ldt1, ldt2, wr,
            wi, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo,  info)
"""
function mb03xd!(balanc::AbstractChar, job::AbstractChar,
    jobu::AbstractChar, jobv::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    t::AbstractMatrix{Float64}, u1::AbstractMatrix{Float64},
    u2::AbstractMatrix{Float64}, v1::AbstractMatrix{Float64},
    v2::AbstractMatrix{Float64}, wr::AbstractVector{Float64},
    wi::AbstractVector{Float64}, scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ldt = max(1,stride(t,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03xd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong, Clong), balanc, job, jobu, jobv, n,
            a, lda, qg, ldqg, t, ldt, u1, ldu1, u2, ldu2, v1, ldv1,
            v2, ldv2, wr, wi, ilo, scale, dwork, ldwork, info, 1, 1,
            1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return ilo[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03xp!(job::AbstractChar, compq::AbstractChar,
    compz::AbstractChar, n::Integer, ilo::Integer, ihi::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, z::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03xp_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong), job,
            compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z,
            ldz, alphar, alphai, beta, dwork, ldwork, info, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03xs!(jobu::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    u1::AbstractMatrix{Float64}, u2::AbstractMatrix{Float64},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03xs_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), jobu, n, a, lda, qg, ldqg, u1,
            ldu1, u2, ldu2, wr, wi, dwork, ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb03xu!(ltra::Bool, ltrb::Bool, n::Integer, k::Integer,
    nb::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, g::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, xa::AbstractMatrix{Float64},
    xb::AbstractMatrix{Float64}, xg::AbstractMatrix{Float64},
    xq::AbstractMatrix{Float64}, ya::AbstractMatrix{Float64},
    yb::AbstractMatrix{Float64}, yg::AbstractMatrix{Float64},
    yq::AbstractMatrix{Float64}, csl::AbstractVector{Float64},
    csr::AbstractVector{Float64}, taul::AbstractVector{Float64},
    taur::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldg = max(1,stride(g,2))
    ldq = max(1,stride(q,2))
    ldxa = max(1,stride(xa,2))
    ldxb = max(1,stride(xb,2))
    ldxg = max(1,stride(xg,2))
    ldxq = max(1,stride(xq,2))
    ldya = max(1,stride(ya,2))
    ldyb = max(1,stride(yb,2))
    ldyg = max(1,stride(yg,2))
    ldyq = max(1,stride(yq,2))
    dwork = Vector{Float64}(undef, 5*nb)

    ccall((:mb03xu_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
            ltra, ltrb, n, k, nb, a, lda, b, ldb, g, ldg, q, ldq,
            xa, ldxa, xb, ldxb, xg, ldxg, xq, ldxq, ya, ldya, yb,
            ldyb, yg, ldyg, yq, ldyq, csl, csr, taul, taur, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (ilo,  info)
"""
function mb03xz!(balanc::AbstractChar, job::AbstractChar,
    jobu::AbstractChar, n::Integer, a::AbstractMatrix{ComplexF64},
    qg::AbstractMatrix{ComplexF64}, u1::AbstractMatrix{ComplexF64},
    u2::AbstractMatrix{ComplexF64}, wr::AbstractVector{Float64},
    wi::AbstractVector{Float64}, scale::AbstractVector{Float64},
    bwork::AbstractVector{BlasBool})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03xz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt}, Clong, Clong,
            Clong), balanc, job, jobu, n, a, lda, qg, ldqg, u1,
            ldu1, u2, ldu2, wr, wi, ilo, scale, dwork, ldwork,
            zwork, lzwork, bwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return ilo[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03ya!(wantt::Bool, wantq::Bool, wantz::Bool, n::Integer,
    ilo::Integer, ihi::Integer, iloq::Integer, ihiq::Integer,
    pos::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()

    ccall((:mb03ya_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasBool}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            wantt, wantq, wantz, n, ilo, ihi, iloq, ihiq, pos, a,
            lda, b, ldb, q, ldq, z, ldz, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb03yd!(wantt::Bool, wantq::Bool, wantz::Bool, n::Integer,
    ilo::Integer, ihi::Integer, iloq::Integer, ihiq::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, z::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03yd_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasBool}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            wantt, wantq, wantz, n, ilo, ihi, iloq, ihiq, a, lda, b,
            ldb, q, ldq, z, ldz, alphar, alphai, beta, dwork,
            ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (csl, snl, csr,  snr)
"""
function mb03yt!(a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    csl = Ref{Float64}()
    snl = Ref{Float64}()
    csr = Ref{Float64}()
    snr = Ref{Float64}()

    ccall((:mb03yt_, libslicot), Cvoid, (Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}), a, lda, b, ldb, alphar, alphai, beta,
            csl, snl, csr, snr)

    return csl[], snl[], csr[], snr[]
end


"""
$(TYPEDSIGNATURES)
returns (m,  info)
"""
function mb03za!(compc::AbstractChar, compu::AbstractChar,
    compv::AbstractChar, compw::AbstractChar, which::AbstractChar,
    select::AbstractVector{BlasBool}, n::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, u1::AbstractMatrix{Float64},
    u2::AbstractMatrix{Float64}, v1::AbstractMatrix{Float64},
    v2::AbstractMatrix{Float64}, w::AbstractMatrix{Float64},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    ldw = max(1,stride(w,2))
    m = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03za_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasBool},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong, Clong), compc, compu, compv, compw, which,
            select, n, a, lda, b, ldb, c, ldc, u1, ldu1, u2, ldu2,
            v1, ldv1, v2, ldv2, w, ldw, wr, wi, m, dwork, ldwork,
            info, 1, 1, 1, 1, 1)
    chkargsok(info[])

    return m[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (m,  info)
"""
function mb03zd!(which::AbstractChar, meth::AbstractChar,
    stab::AbstractChar, balanc::AbstractChar, ortbal::AbstractChar,
    select::AbstractVector{BlasBool}, n::Integer, mm::Integer,
    ilo::Integer, scale::AbstractVector{Float64},
    s::AbstractMatrix{Float64}, t::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64}, u1::AbstractMatrix{Float64},
    u2::AbstractMatrix{Float64}, v1::AbstractMatrix{Float64},
    v2::AbstractMatrix{Float64},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    us::AbstractMatrix{Float64}, uu::AbstractMatrix{Float64},
    iwork::AbstractVector{BlasInt})

    lds = max(1,stride(s,2))
    ldt = max(1,stride(t,2))
    ldg = max(1,stride(g,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    ldus = max(1,stride(us,2))
    lduu = max(1,stride(uu,2))
    m = Ref{BlasInt}()
    info = Ref{BlasInt}()
    lwork = Vector{BlasBool}(undef, 2*n)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03zd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong, Clong,
            Clong), which, meth, stab, balanc, ortbal, select, n,
            mm, ilo, scale, s, lds, t, ldt, g, ldg, u1, ldu1, u2,
            ldu2, v1, ldv1, v2, ldv2, m, wr, wi, us, ldus, uu, lduu,
            lwork, iwork, dwork, ldwork, info, 1, 1, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return m[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04ad!(job::AbstractChar, compq1::AbstractChar,
    compq2::AbstractChar, compu1::AbstractChar, compu2::AbstractChar,
    n::Integer, z::AbstractMatrix{Float64},
    h::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, u11::AbstractMatrix{Float64},
    u12::AbstractMatrix{Float64}, u21::AbstractMatrix{Float64},
    u22::AbstractMatrix{Float64}, t::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer)

    ldz = max(1,stride(z,2))
    ldh = max(1,stride(h,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldu11 = max(1,stride(u11,2))
    ldu12 = max(1,stride(u12,2))
    ldu21 = max(1,stride(u21,2))
    ldu22 = max(1,stride(u22,2))
    ldt = max(1,stride(t,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04ad_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong, Clong,
            Clong), job, compq1, compq2, compu1, compu2, n, z, ldz,
            h, ldh, q1, ldq1, q2, ldq2, u11, ldu11, u12, ldu12, u21,
            ldu21, u22, ldu22, t, ldt, alphar, alphai, beta, iwork,
            liwork, dwork, ldwork, info, 1, 1, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04az!(job::AbstractChar, compq::AbstractChar,
    compu::AbstractChar, n::Integer, z::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, fg::AbstractMatrix{ComplexF64},
    d::AbstractMatrix{ComplexF64}, c::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, u::AbstractMatrix{ComplexF64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer,
    bwork::AbstractVector{BlasBool})

    ldz = max(1,stride(z,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldd = max(1,stride(d,2))
    ldc = max(1,stride(c,2))
    ldq = max(1,stride(q,2))
    ldu = max(1,stride(u,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04az_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt},
            Clong, Clong, Clong), job, compq, compu, n, z, ldz, b,
            ldb, fg, ldfg, d, ldd, c, ldc, q, ldq, u, ldu, alphar,
            alphai, beta, iwork, liwork, dwork, ldwork, zwork,
            lzwork, bwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04bd!(job::AbstractChar, compq1::AbstractChar,
    compq2::AbstractChar, n::Integer, a::AbstractMatrix{Float64},
    de::AbstractMatrix{Float64}, c1::AbstractMatrix{Float64},
    vw::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    f::AbstractMatrix{Float64}, c2::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldc1 = max(1,stride(c1,2))
    ldvw = max(1,stride(vw,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldc2 = max(1,stride(c2,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04bd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), job, compq1, compq2, n, a, lda,
            de, ldde, c1, ldc1, vw, ldvw, q1, ldq1, q2, ldq2, b,
            ldb, f, ldf, c2, ldc2, alphar, alphai, beta, iwork,
            liwork, dwork, ldwork, info, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (info,  iwarn)
"""
function mb04bp!(job::AbstractChar, compq1::AbstractChar,
    compq2::AbstractChar, n::Integer, a::AbstractMatrix{Float64},
    de::AbstractMatrix{Float64}, c1::AbstractMatrix{Float64},
    vw::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    f::AbstractMatrix{Float64}, c2::AbstractMatrix{Float64},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, liwork::Integer, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldc1 = max(1,stride(c1,2))
    ldvw = max(1,stride(vw,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldc2 = max(1,stride(c2,2))
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04bp_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), job, compq1, compq2, n, a, lda,
            de, ldde, c1, ldc1, vw, ldvw, q1, ldq1, q2, ldq2, b,
            ldb, f, ldf, c2, ldc2, alphar, alphai, beta, iwork,
            liwork, dwork, ldwork, info, 1, 1, 1)
    chkargsok(info[])

    return info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04bz!(job::AbstractChar, compq::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, de::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, fg::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    bwork::AbstractVector{BlasBool})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, 2*n+4)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04bz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasBool},
            Ptr{BlasInt}, Clong, Clong), job, compq, n, a, lda, de,
            ldde, b, ldb, fg, ldfg, q, ldq, alphar, alphai, beta,
            iwork, dwork, ldwork, zwork, lzwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04cd!(compq1::AbstractChar, compq2::AbstractChar,
    compq3::AbstractChar, n::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    q1::AbstractMatrix{Float64}, q2::AbstractMatrix{Float64},
    q3::AbstractMatrix{Float64}, liwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldd = max(1,stride(d,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    ldq3 = max(1,stride(q3,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    bwork = Vector{BlasBool}(undef, n÷2)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04cd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt},
            Clong, Clong, Clong), compq1, compq2, compq3, n, a, lda,
            b, ldb, d, ldd, q1, ldq1, q2, ldq2, q3, ldq3, iwork,
            liwork, dwork, ldwork, bwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04db!(job::AbstractChar, sgn::AbstractChar, n::Integer,
    ilo::Integer, lscale::AbstractVector{Float64},
    rscale::AbstractVector{Float64}, m::Integer,
    v1::AbstractMatrix{Float64}, v2::AbstractMatrix{Float64})

    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    info = Ref{BlasInt}()

    ccall((:mb04db_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), job, sgn, n,
            ilo, lscale, rscale, m, v1, ldv1, v2, ldv2, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo,  info)
"""
function mb04dd!(job::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb04dd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Clong), job,
            n, a, lda, qg, ldqg, ilo, scale, info, 1)
    chkargsok(info[])

    return ilo[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04di!(job::AbstractChar, sgn::AbstractChar, n::Integer,
    ilo::Integer, scale::AbstractVector{Float64}, m::Integer,
    v1::AbstractMatrix{Float64}, v2::AbstractMatrix{Float64})

    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    info = Ref{BlasInt}()

    ccall((:mb04di_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), job, sgn, n, ilo, scale, m,
            v1, ldv1, v2, ldv2, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo, ihi, info,  iwarn)
"""
function mb04dl!(job::AbstractChar, n::Integer, thresh::Number,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    lscale::AbstractVector{Float64}, rscale::AbstractVector{Float64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ilo = Ref{BlasInt}()
    ihi = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mb04dl_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Clong), job, n, thresh, a, lda, b, ldb, ilo, ihi,
            lscale, rscale, dwork, iwarn, info, 1)
    chkargsok(info[])

    return ilo[], ihi[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo, info,  iwarn)
"""
function mb04dp!(job::AbstractChar, n::Integer, thresh::Number,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, vw::AbstractMatrix{Float64},
    lscale::AbstractVector{Float64}, rscale::AbstractVector{Float64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldc = max(1,stride(c,2))
    ldvw = max(1,stride(vw,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mb04dp_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Clong), job,
            n, thresh, a, lda, de, ldde, c, ldc, vw, ldvw, ilo,
            lscale, rscale, dwork, iwarn, info, 1)
    chkargsok(info[])

    return ilo[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo,  info)
"""
function mb04ds!(job::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb04ds_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Clong), job,
            n, a, lda, qg, ldqg, ilo, scale, info, 1)
    chkargsok(info[])

    return ilo[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04dy!(jobscl::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    d::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb04dy_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Clong),
            jobscl, n, a, lda, qg, ldqg, d, dwork, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo,  info)
"""
function mb04dz!(job::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, qg::AbstractMatrix{ComplexF64},
    scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()

    ccall((:mb04dz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong), job, n, a, lda, qg, ldqg, ilo, scale, info, 1)
    chkargsok(info[])

    return ilo[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04ed!(job::AbstractChar, compq::AbstractChar,
    compu::AbstractChar, n::Integer, z::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, u1::AbstractMatrix{Float64},
    u2::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    liwork::Integer)

    ldz = max(1,stride(z,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04ed_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong), job,
            compq, compu, n, z, ldz, b, ldb, fg, ldfg, q, ldq, u1,
            ldu1, u2, ldu2, alphar, alphai, beta, iwork, liwork,
            dwork, ldwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04fd!(job::AbstractChar, compq::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n÷2+1)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04fd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), job, compq,
            n, a, lda, de, ldde, b, ldb, fg, ldfg, q, ldq, alphar,
            alphai, beta, iwork, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04fp!(job::AbstractChar, compq::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, de::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, fg::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n÷2+1)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04fp_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), job, compq,
            n, a, lda, de, ldde, b, ldb, fg, ldfg, q, ldq, alphar,
            alphai, beta, iwork, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04gd!(m::Integer, n::Integer, a::AbstractMatrix{Float64},
    jpvt::AbstractVector{BlasInt}, tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, 3*m)

    ccall((:mb04gd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), m, n, a, lda, jpvt, tau,
            dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04hd!(compq1::AbstractChar, compq2::AbstractChar,
    n::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, liwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    bwork = Vector{BlasBool}(undef, n÷2)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04hd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt}, Clong, Clong),
            compq1, compq2, n, a, lda, b, ldb, q1, ldq1, q2, ldq2,
            iwork, liwork, dwork, ldwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04id!(n::Integer, m::Integer, p::Integer, l::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04id_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), n, m, p, l, a, lda, b, ldb,
            tau, dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04iy!(side::AbstractChar, trans::AbstractChar, n::Integer,
    m::Integer, k::Integer, p::Integer, a::AbstractMatrix{Float64},
    tau::AbstractVector{Float64}, c::AbstractMatrix{Float64},
    ldwork::Integer)

    lda = max(1,stride(a,2))
    ldc = max(1,stride(c,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04iy_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), side, trans, n, m, k, p, a, lda, tau, c,
            ldc, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04iz!(n::Integer, m::Integer, p::Integer, l::Integer,
    a::AbstractMatrix{ComplexF64}, b::AbstractMatrix{ComplexF64},
    tau::AbstractVector{ComplexF64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04iz_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}), n, m, p, l, a, lda, b, ldb, tau, zwork,
            lzwork, info)
        chkargsok(info[])
        if iwq == 1
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04jd!(n::Integer, m::Integer, p::Integer, l::Integer,
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04jd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), n, m, p, l, a, lda, b, ldb,
            tau, dwork, ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04kd!(uplo::AbstractChar, n::Integer, m::Integer,
    p::Integer, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, tau::AbstractVector{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    dwork = Vector{Float64}(undef, n)

    ccall((:mb04kd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Clong), uplo, n, m, p, r, ldr, a, lda, b, ldb, c, ldc,
            tau, dwork, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04ld!(uplo::AbstractChar, n::Integer, m::Integer,
    p::Integer, l::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, tau::AbstractVector{Float64})

    ldl = max(1,stride(l,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    dwork = Vector{Float64}(undef, n)

    ccall((:mb04ld_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Clong), uplo, n, m, p, l, ldl, a, lda, b, ldb, c, ldc,
            tau, dwork, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns (maxred,  info)
"""
function mb04md!(n::Integer, ini_maxred::Number,
    a::AbstractMatrix{Float64}, scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    maxred = Ref{Float64}(ini_maxred)
    info = Ref{BlasInt}()

    ccall((:mb04md_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}),
            n, maxred, a, lda, scale, info)
    chkargsok(info[])

    return maxred[], info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04nd!(uplo::AbstractChar, n::Integer, m::Integer,
    p::Integer, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, tau::AbstractVector{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    dwork = Vector{Float64}(undef, max(n-1,m))

    ccall((:mb04nd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Clong), uplo, n, m, p, r, ldr, a, lda, b, ldb, c, ldc,
            tau, dwork, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04ny!(m::Integer, n::Integer, v::AbstractVector{Float64},
    incv::Integer, tau::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    dwork = Vector{Float64}(undef, m)

    ccall((:mb04ny_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}),
            m, n, v, incv, tau, a, lda, b, ldb, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04od!(uplo::AbstractChar, n::Integer, m::Integer,
    p::Integer, r::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, tau::AbstractVector{Float64})

    ldr = max(1,stride(r,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))
    dwork = Vector{Float64}(undef, max(n-1,m))

    ccall((:mb04od_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Clong), uplo, n, m, p, r, ldr, a, lda, b, ldb, c, ldc,
            tau, dwork, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04ow!(m::Integer, n::Integer, p::Integer,
    a::AbstractMatrix{Float64}, t::AbstractMatrix{Float64},
    x::AbstractVector{Float64}, incx::Integer,
    b::AbstractMatrix{Float64}, c::AbstractMatrix{Float64},
    d::AbstractVector{Float64}, incd::Integer)

    lda = max(1,stride(a,2))
    ldt = max(1,stride(t,2))
    ldb = max(1,stride(b,2))
    ldc = max(1,stride(c,2))

    ccall((:mb04ow_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}), m, n, p, a, lda, t, ldt, x, incx, b, ldb,
            c, ldc, d, incd)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04ox!(n::Integer, a::AbstractMatrix{Float64},
    x::AbstractVector{Float64}, incx::Integer)

    lda = max(1,stride(a,2))

    ccall((:mb04ox_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, a, lda, x,
            incx)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04oy!(m::Integer, n::Integer, v::AbstractVector{Float64},
    tau::Number, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    dwork = Vector{Float64}(undef, n)

    ccall((:mb04oy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}), m, n, v, tau,
            a, lda, b, ldb, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04pa!(lham::Bool, n::Integer, k::Integer, nb::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    xa::AbstractMatrix{Float64}, xg::AbstractMatrix{Float64},
    xq::AbstractMatrix{Float64}, ya::AbstractMatrix{Float64},
    cs::AbstractVector{Float64}, tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ldxa = max(1,stride(xa,2))
    ldxg = max(1,stride(xg,2))
    ldxq = max(1,stride(xq,2))
    ldya = max(1,stride(ya,2))
    dwork = Vector{Float64}(undef, 3*nb)

    ccall((:mb04pa_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}), lham, n, k, nb, a, lda, qg, ldqg, xa,
            ldxa, xg, ldxg, xq, ldxq, ya, ldya, cs, tau, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04pb!(n::Integer, ilo::Integer, a::AbstractMatrix{Float64},
    qg::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04pb_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), n, ilo, a, lda, qg, ldqg, cs, tau, dwork,
            ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04pu!(n::Integer, ilo::Integer, a::AbstractMatrix{Float64},
    qg::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04pu_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), n, ilo, a, lda, qg, ldqg, cs, tau, dwork,
            ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04py!(side::AbstractChar, m::Integer, n::Integer,
    v::AbstractVector{Float64}, tau::Number,
    c::AbstractMatrix{Float64}, dwork::AbstractVector{Float64})

    ldc = max(1,stride(c,2))

    ccall((:mb04py_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n, v, tau,
            c, ldc, dwork, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04qb!(tranc::AbstractChar, trand::AbstractChar,
    tranq::AbstractChar, storev::AbstractChar, storew::AbstractChar,
    m::Integer, n::Integer, k::Integer, v::AbstractMatrix{Float64},
    w::AbstractMatrix{Float64}, c::AbstractMatrix{Float64},
    d::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64})

    ldv = max(1,stride(v,2))
    ldw = max(1,stride(w,2))
    ldc = max(1,stride(c,2))
    ldd = max(1,stride(d,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04qb_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong, Clong), tranc, trand, tranq, storev,
            storew, m, n, k, v, ldv, w, ldw, c, ldc, d, ldd, cs,
            tau, dwork, ldwork, info, 1, 1, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04qc!(strab::AbstractChar, trana::AbstractChar,
    tranb::AbstractChar, tranq::AbstractChar, direct::AbstractChar,
    storev::AbstractChar, storew::AbstractChar, m::Integer,
    n::Integer, k::Integer, v::AbstractMatrix{Float64},
    w::AbstractMatrix{Float64}, rs::AbstractMatrix{Float64},
    ldrs::Integer, t::AbstractMatrix{Float64},
    a::AbstractMatrix{Float64}, b::AbstractMatrix{Float64},
    dwork::AbstractVector{Float64})

    ldv = max(1,stride(v,2))
    ldw = max(1,stride(w,2))
    ldt = max(1,stride(t,2))
    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))

    ccall((:mb04qc_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Clong, Clong, Clong, Clong, Clong, Clong,
            Clong), strab, trana, tranb, tranq, direct, storev,
            storew, m, n, k, v, ldv, w, ldw, rs, ldrs, t, ldt, a,
            lda, b, ldb, dwork, 1, 1, 1, 1, 1, 1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04qf!(direct::AbstractChar, storev::AbstractChar,
    storew::AbstractChar, n::Integer, k::Integer,
    v::AbstractMatrix{Float64}, w::AbstractMatrix{Float64},
    cs::AbstractVector{Float64}, tau::AbstractVector{Float64},
    rs::AbstractMatrix{Float64}, t::AbstractMatrix{Float64})

    ldv = max(1,stride(v,2))
    ldw = max(1,stride(w,2))
    ldrs = max(1,stride(rs,2))
    ldt = max(1,stride(t,2))
    dwork = Vector{Float64}(undef, 3*k)

    ccall((:mb04qf_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Clong, Clong, Clong),
            direct, storev, storew, n, k, v, ldv, w, ldw, cs, tau,
            rs, ldrs, t, ldt, dwork, 1, 1, 1)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04qs!(tranc::AbstractChar, trand::AbstractChar,
    tranu::AbstractChar, m::Integer, n::Integer, ilo::Integer,
    v::AbstractMatrix{Float64}, w::AbstractMatrix{Float64},
    c::AbstractMatrix{Float64}, d::AbstractMatrix{Float64},
    cs::AbstractVector{Float64}, tau::AbstractVector{Float64})

    ldv = max(1,stride(v,2))
    ldw = max(1,stride(w,2))
    ldc = max(1,stride(c,2))
    ldd = max(1,stride(d,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04qs_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong, Clong), tranc, trand, tranu,
            m, n, ilo, v, ldv, w, ldw, c, ldc, d, ldd, cs, tau,
            dwork, ldwork, info, 1, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04qu!(tranc::AbstractChar, trand::AbstractChar,
    tranq::AbstractChar, storev::AbstractChar, storew::AbstractChar,
    m::Integer, n::Integer, k::Integer, v::AbstractMatrix{Float64},
    w::AbstractMatrix{Float64}, c::AbstractMatrix{Float64},
    d::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    ldv = max(1,stride(v,2))
    ldw = max(1,stride(w,2))
    ldc = max(1,stride(c,2))
    ldd = max(1,stride(d,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04qu_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong,
            Clong, Clong, Clong), tranc, trand, tranq, storev,
            storew, m, n, k, v, ldv, w, ldw, c, ldc, d, ldd, cs,
            tau, dwork, ldwork, info, 1, 1, 1, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04rb!(n::Integer, ilo::Integer, a::AbstractMatrix{Float64},
    qg::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04rb_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), n, ilo, a, lda, qg, ldqg, cs, tau, dwork,
            ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04ru!(n::Integer, ilo::Integer, a::AbstractMatrix{Float64},
    qg::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04ru_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), n, ilo, a, lda, qg, ldqg, cs, tau, dwork,
            ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04su!(m::Integer, n::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04su_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), m, n, a, lda, b, ldb, cs, tau, dwork,
            ldwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04tb!(trana::AbstractChar, tranb::AbstractChar, n::Integer,
    ilo::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, g::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, csl::AbstractVector{Float64},
    csr::AbstractVector{Float64}, taul::AbstractVector{Float64},
    taur::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldg = max(1,stride(g,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04tb_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), trana, tranb, n, ilo, a,
            lda, b, ldb, g, ldg, q, ldq, csl, csr, taul, taur,
            dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04ts!(trana::AbstractChar, tranb::AbstractChar, n::Integer,
    ilo::Integer, a::AbstractMatrix{Float64},
    b::AbstractMatrix{Float64}, g::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, csl::AbstractVector{Float64},
    csr::AbstractVector{Float64}, taul::AbstractVector{Float64},
    taur::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ldg = max(1,stride(g,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04ts_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), trana, tranb, n, ilo, a,
            lda, b, ldb, g, ldg, q, ldq, csl, csr, taul, taur,
            dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns rank
"""
function mb04tt!(updatq::Bool, updatz::Bool, m::Integer, n::Integer,
    ifira::Integer, ifica::Integer, nca::Integer,
    a::AbstractMatrix{Float64}, e::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, z::AbstractMatrix{Float64},
    istair::AbstractVector{BlasInt}, tol::Number)

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    rank = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)

    ccall((:mb04tt_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{Float64},
            Ptr{BlasInt}), updatq, updatz, m, n, ifira, ifica, nca,
            a, lda, e, lde, q, ldq, z, ldz, istair, rank, tol,
            iwork)

    return rank[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04tu!(n::Integer, x::AbstractVector{Float64},
    incx::Integer, y::AbstractVector{Float64}, incy::Integer,
    c::Number, s::Number)


    ccall((:mb04tu_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}), n, x, incx, y, incy, c, s)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04tv!(updatz::Bool, n::Integer, nra::Integer, nca::Integer,
    ifira::Integer, ifica::Integer, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, z::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldz = max(1,stride(z,2))

    ccall((:mb04tv_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}), updatz, n, nra, nca, ifira,
            ifica, a, lda, e, lde, z, ldz)

    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function mb04tw!(updatq::Bool, m::Integer, n::Integer, nre::Integer,
    nce::Integer, ifire::Integer, ifice::Integer, ifica::Integer,
    a::AbstractMatrix{Float64}, e::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))

    ccall((:mb04tw_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
            updatq, m, n, nre, nce, ifire, ifice, ifica, a, lda, e,
            lde, q, ldq)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns nblcks
"""
function mb04tx!(updatq::Bool, updatz::Bool, m::Integer, n::Integer,
    ini_nblcks::Integer, inuk::AbstractVector{BlasInt},
    imuk::AbstractVector{BlasInt}, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64}, mnei::AbstractVector{BlasInt})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    nblcks = Ref{BlasInt}(ini_nblcks)

    ccall((:mb04tx_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), updatq, updatz, m, n,
            nblcks, inuk, imuk, a, lda, e, lde, q, ldq, z, ldz,
            mnei)

    return nblcks[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04ty!(updatq::Bool, updatz::Bool, m::Integer, n::Integer,
    nblcks::Integer, inuk::AbstractVector{BlasInt},
    imuk::AbstractVector{BlasInt}, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    info = Ref{BlasInt}()

    ccall((:mb04ty_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), updatq, updatz, m, n,
            nblcks, inuk, imuk, a, lda, e, lde, q, ldq, z, ldz,
            info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ranke,  info)
"""
function mb04ud!(jobq::AbstractChar, jobz::AbstractChar, m::Integer,
    n::Integer, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64}, istair::AbstractVector{BlasInt},
    tol::Number)

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    ranke = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, max(m,n))

    ccall((:mb04ud_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Clong, Clong),
            jobq, jobz, m, n, a, lda, e, lde, q, ldq, z, ldz, ranke,
            istair, tol, dwork, info, 1, 1)
    chkargsok(info[])

    return ranke[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (nblcks, nblcki,  info)
"""
function mb04vd!(mode::AbstractChar, jobq::AbstractChar,
    jobz::AbstractChar, m::Integer, n::Integer, ranke::Integer,
    a::AbstractMatrix{Float64}, e::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64}, z::AbstractMatrix{Float64},
    istair::AbstractVector{BlasInt}, imuk::AbstractVector{BlasInt},
    inuk::AbstractVector{BlasInt}, imuk0::AbstractVector{BlasInt},
    mnei::AbstractVector{BlasInt}, tol::Number)

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))
    nblcks = Ref{BlasInt}()
    nblcki = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)

    ccall((:mb04vd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong), mode,
            jobq, jobz, m, n, ranke, a, lda, e, lde, q, ldq, z, ldz,
            istair, nblcks, nblcki, imuk, inuk, imuk0, mnei, tol,
            iwork, info, 1, 1, 1)
    chkargsok(info[])

    return nblcks[], nblcki[], info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04vx!(updatq::Bool, updatz::Bool, m::Integer, n::Integer,
    nblcks::Integer, inuk::AbstractVector{BlasInt},
    imuk::AbstractVector{BlasInt}, a::AbstractMatrix{Float64},
    e::AbstractMatrix{Float64}, q::AbstractMatrix{Float64},
    z::AbstractMatrix{Float64}, mnei::AbstractVector{BlasInt})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldq = max(1,stride(q,2))
    ldz = max(1,stride(z,2))

    ccall((:mb04vx_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), updatq, updatz, m, n,
            nblcks, inuk, imuk, a, lda, e, lde, q, ldq, z, ldz,
            mnei)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04wd!(tranq1::AbstractChar, tranq2::AbstractChar,
    m::Integer, n::Integer, k::Integer, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64})

    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04wd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), tranq1, tranq2, m, n, k, q1, ldq1, q2,
            ldq2, cs, tau, dwork, ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04wp!(n::Integer, ilo::Integer,
    u1::AbstractMatrix{Float64}, u2::AbstractMatrix{Float64},
    cs::AbstractVector{Float64}, tau::AbstractVector{Float64})

    ldu1 = max(1,stride(u1,2))
    ldu2 = max(1,stride(u2,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04wp_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}), n, ilo, u1, ldu1, u2, ldu2, cs, tau,
            dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04wr!(job::AbstractChar, trans::AbstractChar, n::Integer,
    ilo::Integer, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64})

    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04wr_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong),
            job, trans, n, ilo, q1, ldq1, q2, ldq2, cs, tau, dwork,
            ldwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04wu!(tranq1::AbstractChar, tranq2::AbstractChar,
    m::Integer, n::Integer, k::Integer, q1::AbstractMatrix{Float64},
    q2::AbstractMatrix{Float64}, cs::AbstractVector{Float64},
    tau::AbstractVector{Float64}, ldwork::Integer)

    ldq1 = max(1,stride(q1,2))
    ldq2 = max(1,stride(q2,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04wu_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong, Clong), tranq1, tranq2, m, n, k, q1, ldq1, q2,
            ldq2, cs, tau, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank, theta, info,  iwarn)
"""
function mb04xd!(jobu::AbstractChar, jobv::AbstractChar, m::Integer,
    n::Integer, ini_rank::Integer, ini_theta::Number,
    a::AbstractMatrix{Float64}, u::AbstractMatrix{Float64},
    v::AbstractMatrix{Float64}, q::AbstractVector{Float64},
    inul::AbstractVector{BlasBool}, tol::Number, reltol::Number)

    lda = max(1,stride(a,2))
    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    rank = Ref{BlasInt}(ini_rank)
    theta = Ref{Float64}(ini_theta)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb04xd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasBool},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong, Clong), jobu, jobv,
            m, n, rank, theta, a, lda, u, ldu, v, ldv, q, inul, tol,
            reltol, dwork, ldwork, iwarn, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return rank[], theta[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04xy!(jobu::AbstractChar, jobv::AbstractChar, m::Integer,
    n::Integer, x::AbstractMatrix{Float64},
    taup::AbstractVector{Float64}, tauq::AbstractVector{Float64},
    u::AbstractMatrix{Float64}, v::AbstractMatrix{Float64},
    inul::AbstractVector{BlasBool})

    ldx = max(1,stride(x,2))
    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    info = Ref{BlasInt}()

    ccall((:mb04xy_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt},
            Clong, Clong), jobu, jobv, m, n, x, ldx, taup, tauq, u,
            ldu, v, ldv, inul, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank, theta, info,  iwarn)
"""
function mb04yd!(jobu::AbstractChar, jobv::AbstractChar, m::Integer,
    n::Integer, ini_rank::Integer, ini_theta::Number,
    q::AbstractVector{Float64}, e::AbstractVector{Float64},
    u::AbstractMatrix{Float64}, v::AbstractMatrix{Float64},
    inul::AbstractVector{BlasBool}, tol::Number, reltol::Number,
    ldwork::Integer)

    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    rank = Ref{BlasInt}(ini_rank)
    theta = Ref{Float64}(ini_theta)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb04yd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasBool}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Clong, Clong), jobu, jobv, m, n, rank,
            theta, q, e, u, ldu, v, ldv, inul, tol, reltol, dwork,
            ldwork, iwarn, info, 1, 1)
    chkargsok(info[])

    return rank[], theta[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
"""
function mb04yw!(qrit::Bool, updatu::Bool, updatv::Bool, m::Integer,
    n::Integer, l::Integer, k::Integer, shift::Number,
    d::AbstractVector{Float64}, e::AbstractVector{Float64},
    u::AbstractMatrix{Float64}, v::AbstractMatrix{Float64})

    ldu = max(1,stride(u,2))
    ldv = max(1,stride(v,2))
    dwork = Vector{Float64}(undef, max(1,ldwork))

    ccall((:mb04yw_, libslicot), Cvoid, (Ref{BlasBool}, Ref{BlasBool},
            Ref{BlasBool}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}), qrit, updatu, updatv, m, n, l, k, shift,
            d, e, u, ldu, v, ldv, dwork)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb04zd!(compu::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, qg::AbstractMatrix{Float64},
    u::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldqg = max(1,stride(qg,2))
    ldu = max(1,stride(u,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, 2*n)

    ccall((:mb04zd_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Clong), compu, n, a, lda, qg, ldqg, u, ldu, dwork, info,
            1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb05md!(balanc::AbstractChar, n::Integer, delta::Number,
    a::AbstractMatrix{Float64}, v::AbstractMatrix{Float64},
    y::AbstractMatrix{Float64}, valr::AbstractVector{Float64},
    vali::AbstractVector{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    ldv = max(1,stride(v,2))
    ldy = max(1,stride(y,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb05md_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong), balanc, n, delta, a, lda, v, ldv,
            y, ldy, valr, vali, iwork, dwork, ldwork, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb05my!(balanc::AbstractChar, n::Integer,
    a::AbstractMatrix{Float64}, wr::AbstractVector{Float64},
    wi::AbstractVector{Float64}, r::AbstractMatrix{Float64},
    q::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    ldr = max(1,stride(r,2))
    ldq = max(1,stride(q,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb05my_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong),
            balanc, n, a, lda, wr, wi, r, ldr, q, ldq, dwork,
            ldwork, info, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb05nd!(n::Integer, delta::Number,
    a::AbstractMatrix{Float64}, ex::AbstractMatrix{Float64},
    exint::AbstractMatrix{Float64}, tol::Number,
    ldwork::Integer)

    lda = max(1,stride(a,2))
    ldex = max(1,stride(ex,2))
    ldexin = max(1,stride(exint,2))
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb05nd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{Float64},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), n, delta, a,
            lda, ex, ldex, exint, ldexin, tol, iwork, dwork, ldwork,
            info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (mdig, idig, info,  iwarn)
"""
function mb05od!(balanc::AbstractChar, n::Integer, ndiag::Integer,
    delta::Number, a::AbstractMatrix{Float64}, ldwork::Integer)

    lda = max(1,stride(a,2))
    mdig = Ref{BlasInt}()
    idig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb05od_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Clong),
            balanc, n, ndiag, delta, a, lda, mdig, idig, iwork,
            dwork, ldwork, iwarn, info, 1)
    chkargsok(info[])

    return mdig[], idig[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb05oy!(job::AbstractChar, n::Integer, low::Integer,
    igh::Integer, a::AbstractMatrix{Float64},
    scale::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    info = Ref{BlasInt}()

    ccall((:mb05oy_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}, Clong), job, n, low, igh, a,
            lda, scale, info, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb3jzp!(compq::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, d::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, f::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, tol::Number)

    lda = max(1,stride(a,2))
    ldd = max(1,stride(d,2))
    ldb = max(1,stride(b,2))
    ldf = max(1,stride(f,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n÷2)
    zwork = Vector{ComplexF64}(undef, n÷2)

    ccall((:mb3jzp_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ptr{Float64},
            Ptr{ComplexF64}, Ptr{BlasInt}, Clong), compq, n, a, lda,
            d, ldd, b, ldb, f, ldf, q, ldq, neig, tol, dwork, zwork,
            info, 1)
    chkargsok(info[])

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (neig,  info)
"""
function mb3lzp!(compq::AbstractChar, orth::AbstractChar, n::Integer,
    a::AbstractMatrix{ComplexF64}, de::AbstractMatrix{ComplexF64},
    b::AbstractMatrix{ComplexF64}, fg::AbstractMatrix{ComplexF64},
    q::AbstractMatrix{ComplexF64}, alphar::AbstractVector{Float64},
    alphai::AbstractVector{Float64}, beta::AbstractVector{Float64},
    bwork::AbstractVector{BlasBool})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldb = max(1,stride(b,2))
    ldfg = max(1,stride(fg,2))
    ldq = max(1,stride(q,2))
    neig = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, n+1)
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)
    lzwork = BlasInt(-1)
    zwork = Vector{ComplexF64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb3lzp_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{BlasBool}, Ptr{BlasInt}, Clong, Clong),
            compq, orth, n, a, lda, de, ldde, b, ldb, fg, ldfg,
            neig, q, ldq, alphar, alphai, beta, iwork, dwork,
            ldwork, zwork, lzwork, bwork, info, 1, 1)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
            lzwork = BlasInt(real(zwork[1]))
            resize!(zwork, lzwork)
        end
    end

    return neig[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb3oyz!(m::Integer, n::Integer,
    a::AbstractMatrix{ComplexF64}, rcond::Number,
    svlmax::Number, sval::AbstractVector{Float64},
    jpvt::AbstractVector{BlasInt}, tau::AbstractVector{ComplexF64})

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  2*n )
    zwork = Vector{ComplexF64}(undef,  3*n-1 )

    ccall((:mb3oyz_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64},
            Ptr{BlasInt}), m, n, a, lda, rcond, svlmax, rank, sval,
            jpvt, tau, dwork, zwork, info)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (rank,  info)
"""
function mb3pyz!(m::Integer, n::Integer,
    a::AbstractMatrix{ComplexF64}, rcond::Number,
    svlmax::Number, sval::AbstractVector{Float64},
    jpvt::AbstractVector{BlasInt}, tau::AbstractVector{ComplexF64})

    lda = max(1,stride(a,2))
    rank = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef,  2*m )
    zwork = Vector{ComplexF64}(undef,  3*m-1 )

    ccall((:mb3pyz_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64},
            Ptr{BlasInt}), m, n, a, lda, rcond, svlmax, rank, sval,
            jpvt, tau, dwork, zwork, info)
    chkargsok(info[])

    return rank[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mb4dbz!(job::AbstractChar, sgn::AbstractChar, n::Integer,
    ilo::Integer, lscale::AbstractVector{Float64},
    rscale::AbstractVector{Float64}, m::Integer,
    v1::AbstractMatrix{ComplexF64}, v2::AbstractMatrix{ComplexF64})

    ldv1 = max(1,stride(v1,2))
    ldv2 = max(1,stride(v2,2))
    info = Ref{BlasInt}()

    ccall((:mb4dbz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Clong,
            Clong), job, sgn, n, ilo, lscale, rscale, m, v1, ldv1,
            v2, ldv2, info, 1, 1)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo, ihi, info,  iwarn)
"""
function mb4dlz!(job::AbstractChar, n::Integer, thresh::Number,
    a::AbstractMatrix{ComplexF64}, b::AbstractMatrix{ComplexF64},
    lscale::AbstractVector{Float64}, rscale::AbstractVector{Float64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldb = max(1,stride(b,2))
    ilo = Ref{BlasInt}()
    ihi = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mb4dlz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong), job, n, thresh, a,
            lda, b, ldb, ilo, ihi, lscale, rscale, dwork, iwarn,
            info, 1)
    chkargsok(info[])

    return ilo[], ihi[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (ilo, info,  iwarn)
"""
function mb4dpz!(job::AbstractChar, n::Integer, thresh::Number,
    a::AbstractMatrix{ComplexF64}, de::AbstractMatrix{ComplexF64},
    c::AbstractMatrix{ComplexF64}, vw::AbstractMatrix{ComplexF64},
    lscale::AbstractVector{Float64}, rscale::AbstractVector{Float64},
    dwork::AbstractVector{Float64})

    lda = max(1,stride(a,2))
    ldde = max(1,stride(de,2))
    ldc = max(1,stride(c,2))
    ldvw = max(1,stride(vw,2))
    ilo = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mb4dpz_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong), job, n, thresh, a,
            lda, de, ldde, c, ldc, vw, ldvw, ilo, lscale, rscale,
            dwork, iwarn, info, 1)
    chkargsok(info[])

    return ilo[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01md!(dp::Integer, alpha::Number, k::Integer,
    p::AbstractVector{Float64}, q::AbstractVector{Float64})

    info = Ref{BlasInt}()

    ccall((:mc01md_, libslicot), Cvoid, (Ref{BlasInt}, Ref{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            dp, alpha, k, p, q, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (vr, vi,  info)
"""
function mc01nd!(dp::Integer, xr::Number, xi::Number,
    p::AbstractVector{Float64})

    vr = Ref{Float64}()
    vi = Ref{Float64}()
    info = Ref{BlasInt}()

    ccall((:mc01nd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}), dp, xr, xi, p, vr, vi, info)
    chkargsok(info[])

    return vr[], vi[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01od!(k::Integer, rez::AbstractVector{Float64},
    imz::AbstractVector{Float64}, rep::AbstractVector{Float64},
    imp::AbstractVector{Float64})

    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, 2*k+2)

    ccall((:mc01od_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{BlasInt}), k, rez, imz, rep, imp, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01pd!(k::Integer, rez::AbstractVector{Float64},
    imz::AbstractVector{Float64}, p::AbstractVector{Float64})

    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, k+1)

    ccall((:mc01pd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            k, rez, imz, p, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01py!(k::Integer, rez::AbstractVector{Float64},
    imz::AbstractVector{Float64}, p::AbstractVector{Float64})

    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, k)

    ccall((:mc01py_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            k, rez, imz, p, dwork, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (db, info,  iwarn)
"""
function mc01qd!(da::Integer, ini_db::Integer,
    a::AbstractVector{Float64}, b::AbstractVector{Float64},
    rq::AbstractVector{Float64})

    db = Ref{BlasInt}(ini_db)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()

    ccall((:mc01qd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{BlasInt}), da, db, a, b, rq, iwarn, info)
    chkargsok(info[])

    return db[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (dp3,  info)
"""
function mc01rd!(dp1::Integer, dp2::Integer, ini_dp3::Integer,
    alpha::Number, p1::AbstractVector{Float64},
    p2::AbstractVector{Float64}, p3::AbstractVector{Float64})

    dp3 = Ref{BlasInt}(ini_dp3)
    info = Ref{BlasInt}()

    ccall((:mc01rd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), dp1, dp2, dp3, alpha, p1,
            p2, p3, info)
    chkargsok(info[])

    return dp3[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (s, t,  info)
"""
function mc01sd!(dp::Integer, p::AbstractVector{Float64},
    mant::AbstractVector{Float64}, e::AbstractVector{BlasInt})

    s = Ref{BlasInt}()
    t = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, dp+1)

    ccall((:mc01sd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}), dp, p, s, t, mant, e,
            iwork, info)
    chkargsok(info[])

    return s[], t[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (m,  e)
"""
function mc01sw!(a::Number, b::Integer)

    m = Ref{Float64}()
    e = Ref{BlasInt}()

    ccall((:mc01sw_, libslicot), Cvoid, (Ref{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{BlasInt}), a, b, m, e)

    return m[], e[]
end


"""
$(TYPEDSIGNATURES)
"""
function mc01sx!(lb::Integer, ub::Integer, e::AbstractVector{BlasInt},
    mant::AbstractVector{Float64})


    jlres = ccall((:mc01sx_, libslicot), BlasInt, (Ref{BlasInt},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}), lb, ub, e,
            mant)

    return jlres
end


"""
$(TYPEDSIGNATURES)
returns (a,  ovflow)
"""
function mc01sy!(m::Number, e::Integer, b::Integer)

    a = Ref{Float64}()
    ovflow = Ref{BlasBool}()

    ccall((:mc01sy_, libslicot), Cvoid, (Ref{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasBool}), m, e, b, a,
            ovflow)

    return a[], ovflow[]
end


"""
$(TYPEDSIGNATURES)
returns (dp, stable, nz, info,  iwarn)
"""
function mc01td!(dico::AbstractChar, ini_dp::Integer,
    p::AbstractVector{Float64})

    dp = Ref{BlasInt}(ini_dp)
    stable = Ref{BlasBool}()
    nz = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, 2*ini_dp+2)

    ccall((:mc01td_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt},
            Ptr{Float64}, Ptr{BlasBool}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}, Ptr{BlasInt}, Clong), dico, dp, p, stable,
            nz, dwork, iwarn, info, 1)
    chkargsok(info[])

    return dp[], stable[], nz[], info[], iwarn[]
end


"""
$(TYPEDSIGNATURES)
returns (z1re, z1im, z2re, z2im,  info)
"""
function mc01vd!(a::Number, b::Number, c::Number)

    z1re = Ref{Float64}()
    z1im = Ref{Float64}()
    z2re = Ref{Float64}()
    z2im = Ref{Float64}()
    info = Ref{BlasInt}()

    ccall((:mc01vd_, libslicot), Cvoid, (Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}), a, b, c, z1re, z1im, z2re,
            z2im, info)
    chkargsok(info[])

    return z1re[], z1im[], z2re[], z2im[], info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01wd!(dp::Integer, p::AbstractVector{Float64}, u1::Number,
    u2::Number, q::AbstractVector{Float64})

    info = Ref{BlasInt}()

    ccall((:mc01wd_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{Float64},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}),
            dp, p, u1, u2, q, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc01xd!(alpha::Number, beta::Number, gamma::Number,
    delta::Number, evr::AbstractVector{Float64},
    evi::AbstractVector{Float64}, evq::AbstractVector{Float64})

    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    # some of this is used even in the workspace query
    dwork = Vector{Float64}(undef, 64)

    local jlres
    for iwq in 1:2
        ccall((:mc01xd_, libslicot), Cvoid, (Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            alpha, beta, gamma, delta, evr, evi, evq, dwork, ldwork,
            info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (dp3,  info)
"""
function mc03md!(rp1::Integer, cp1::Integer, cp2::Integer,
    dp1::Integer, dp2::Integer, ini_dp3::Integer, alpha::Number,
    p1::Array{Float64,3}, p2::Array{Float64,3}, p3::Array{Float64,3})

    ldp11 = max(1,stride(p1,2))
    ldp12 = max(1,stride(p1,3)÷ldp11)
    ldp21 = max(1,stride(p2,2))
    ldp22 = max(1,stride(p2,3)÷ldp21)
    ldp31 = max(1,stride(p3,2))
    ldp32 = max(1,stride(p3,3)÷ldp31)
    dp3 = Ref{BlasInt}(ini_dp3)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, cp1)

    ccall((:mc03md_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}),
            rp1, cp1, cp2, dp1, dp2, dp3, alpha, p1, ldp11, ldp12,
            p2, ldp21, ldp22, p3, ldp31, ldp32, dwork, info)
    chkargsok(info[])

    return dp3[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (dk,  info)
"""
function mc03nd!(mp::Integer, np::Integer, dp::Integer,
    p::Array{Float64,3},
    gam::AbstractVector{BlasInt}, nullsp::AbstractMatrix{Float64},
    ker::Array{Float64,3}, tol::Number, iwork::AbstractVector{BlasInt},
    ldwork::Integer)

    ldp1 = max(1,stride(p,2))
    ldp2 = max(1,stride(p,3)÷ldp1)
    ldnull = max(1,stride(nullsp,2))
    ldker1 = max(1,stride(ker,2))
    ldker2 = max(1,stride(ker,3)÷ldker1)
    dk = Ref{BlasInt}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mc03nd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
            mp, np, dp, p, ldp1, ldp2, dk, gam, nullsp, ldnull, ker,
            ldker1, ldker2, tol, iwork, dwork, ldwork, info)
    chkargsok(info[])

    return dk[], info[]
end


"""
$(TYPEDSIGNATURES)
"""
function mc03nx!(mp::Integer, np::Integer, dp::Integer,
    p::Array{Float64,3},
    a::AbstractMatrix{Float64}, e::AbstractMatrix{Float64})

    ldp1 = max(1,stride(p,2))
    ldp2 = max(1,stride(p,3)÷ldp1)
    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))

    ccall((:mc03nx_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
            mp, np, dp, p, ldp1, ldp2, a, lda, e, lde)

    return nothing
end


"""
$(TYPEDSIGNATURES)
returns info
"""
function mc03ny!(nblcks::Integer, nra::Integer, nca::Integer,
    a::AbstractMatrix{Float64}, e::AbstractMatrix{Float64},
    imuk::AbstractVector{BlasInt}, inuk::AbstractVector{BlasInt},
    veps::AbstractMatrix{Float64})

    lda = max(1,stride(a,2))
    lde = max(1,stride(e,2))
    ldveps = max(1,stride(veps,2))
    info = Ref{BlasInt}()

    ccall((:mc03ny_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), nblcks, nra, nca, a, lda,
            e, lde, imuk, inuk, veps, ldveps, info)
    chkargsok(info[])

    return info[]
end


"""
$(TYPEDSIGNATURES)
returns (gnorm,  info)
"""
function md03ba!(n::Integer, ipar::AbstractVector{BlasInt},
    lipar::Integer, fnorm::Number, j::AbstractVector{Float64},
    e::AbstractVector{Float64},
    jnorms::AbstractVector{Float64}, ipvt::AbstractVector{BlasInt},
    ldwork::Integer)

    ldj = max(1,stride(j,2))
    gnorm = Ref{Float64}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:md03ba_, libslicot), Cvoid, (Ref{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}), n, ipar,
            lipar, fnorm, j, ldj, e, jnorms, gnorm, ipvt, dwork,
            ldwork, info)
    chkargsok(info[])

    return gnorm[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (par,  info)
"""
function md03bb!(cond::AbstractChar, n::Integer,
    ipar::AbstractVector{BlasInt}, lipar::Integer,
    r::AbstractMatrix{Float64}, ipvt::AbstractVector{BlasInt},
    diag::AbstractVector{Float64}, qtb::AbstractVector{Float64},
    delta::Number, ini_par::Number, ranks::AbstractVector{BlasInt},
    x::AbstractVector{Float64}, rx::AbstractVector{Float64},
    tol::Number, ldwork::Integer)

    ldr = max(1,stride(r,2))
    par = Ref{Float64}(ini_par)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:md03bb_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
            Clong), cond, n, ipar, lipar, r, ldr, ipvt, diag, qtb,
            delta, par, ranks, x, rx, tol, dwork, ldwork, info, 1)
    chkargsok(info[])

    return par[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (gnorm,  info)
"""
function md03bx!(m::Integer, n::Integer, fnorm::Number,
    j::AbstractVector{Float64},
    e::AbstractVector{Float64}, jnorms::AbstractVector{Float64},
    ipvt::AbstractVector{BlasInt}, ldwork::Integer)

    ldj = max(1,stride(j,2))
    gnorm = Ref{Float64}()
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:md03bx_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), m, n, fnorm, j, ldj, e,
            jnorms, gnorm, ipvt, dwork, ldwork, info)
    chkargsok(info[])

    return gnorm[], info[]
end


"""
$(TYPEDSIGNATURES)
returns (par,  info)
"""
function md03by!(cond::AbstractChar, n::Integer,
    r::AbstractMatrix{Float64}, ipvt::AbstractVector{BlasInt},
    diag::AbstractVector{Float64}, qtb::AbstractVector{Float64},
    delta::Number, ini_par::Number, rank::Integer,
    x::AbstractVector{Float64}, rx::AbstractVector{Float64},
    tol::Number, ldwork::Integer)

    ldr = max(1,stride(r,2))
    par = Ref{Float64}(ini_par)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:md03by_, libslicot), Cvoid, (Ref{UInt8}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}, Clong), cond, n, r, ldr,
            ipvt, diag, qtb, delta, par, rank, x, rx, tol, dwork,
            ldwork, info, 1)
    chkargsok(info[])

    return par[], info[]
end
