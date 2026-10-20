const EPS = 1e-9
using Printf

function escolhe_entrante(z::Vector{Float64}, nb_idxs::Vector{Int})
    best = 0
    best_idx = 0
    for j in nb_idxs
        if z[j] > EPS
            if best == 0 || j < best
                best = j
                best_idx = j
            end
        end
    end
    return best_idx
end

function escolhe_sainte(col::Vector{Float64}, b::Vector{Float64}, bas_idxs::Vector{Int})
    best_row = 0
    best_val = 0.0
    for (r, br) in enumerate(bas_idxs)
        a = col[r]
        if a > EPS
            ratio = b[r] / a
            if best_row == 0 || ratio < best_val - EPS || (abs(ratio - best_val) ≤ EPS && br < bas_idxs[best_row])
                best_row = r
                best_val = ratio
            end
        end
    end
    return best_row
end

function pivot!(T::Array{Float64,2}, r::Int, c::Int)
    m, n = size(T)
    piv = T[r,c]
    for j in 1:n
        T[r,j] /= piv
    end
    
    for i in 1:m
        i == r && continue
        f = T[i,c]
        if abs(f) > EPS
            for j in 1:n
                T[i,j] -= f * T[r,j]
            end
        end
    end
end

function simplex_fase2!(T::Array{Float64,2}, bas::Vector{Int}, nb::Vector{Int})
    m, n = size(T)
    while true
        z = T[m, :] 
        j_entra = escolhe_entrante(z, nb)
        if j_entra == 0
            return (:otima, T, bas, nb)
        end
        
        col = T[1:m-1, j_entra]
        b = T[1:m-1, n]  
        r = escolhe_sainte(col, b, bas)
        if r == 0
            return (:ilimitado, T, bas, nb)
        end
        pivot!(T, r, j_entra)

        antigo = bas[r]
        bas[r] = j_entra
        
        for k in 1:length(nb)
            if nb[k] == j_entra
                nb[k] = antigo
                break
            end
        end
    end
end


function monta_auxiliar(A::Array{Float64,2}, b::Vector{Float64}, c::Vector{Float64})
    n, m = size(A)
    totalvars = m + n
    T = zeros(Float64, n+1, totalvars + 1)

    for i in 1:n
        for j in 1:m
            T[i,j] = A[i,j]
        end
        T[i, m+i] = 1.0
        T[i, end] = b[i]
    end

    for i in 1:n
        T[end, m+i] = -1.0
    end

    for i in 1:n
        for j in 1:(totalvars+1)
            T[end,j] += T[i,j]
        end
    end

    bas = [m+i for i in 1:n]
    nb  = [j for j in 1:m]

    return T, bas, nb
end

function constroi_fase2!(T1::Array{Float64,2}, bas::Vector{Int}, nb::Vector{Int}, c::Vector{Float64}, m::Int)
    nlin, ncol = size(T1)
    n = nlin - 1
    T2 = zeros(Float64, n+1, m+1)

    for i in 1:n
        for j in 1:m
            T2[i,j] = T1[i,j]
        end
        T2[i,end] = T1[i,end]
    end

    bas2 = copy(bas)
    nb2  = Int[]
    for j in 1:m
        push!(nb2, j)
    end

    for b in bas2
        if b ≤ m
            filter!(x -> x != b, nb2)
        end
    end

    for i in 1:n
        if bas2[i] > m
            mudou = false
            for j in nb2
                if abs(T2[i,j]) > EPS
                    pivot!(T2, i, j)
                    antigo = bas2[i]
                    bas2[i] = j
                    filter!(x -> x != j, nb2)
                    mudou = true
                    break
                end
            end
        end
    end

    for j in 1:m
        T2[end,j] = c[j]
    end
    T2[end,end] = 0.0

    for i in 1:n
        bvar = bas2[i]
        if 1 ≤ bvar ≤ m
            cb = c[bvar]
            if abs(cb) > EPS
                for j in 1:(m+1)
                    T2[end,j] -= cb * T2[i,j]
                end
            end
        end
    end

    nb2 = Int[]
    for j in 1:m
        if !(j in bas2)
            push!(nb2, j)
        end
    end

    return T2, bas2, nb2
end

function resolve_pl(A::Array{Float64,2}, b::Vector{Float64}, c::Vector{Float64})
    n, m = size(A)

    T1, bas1, nb1 = monta_auxiliar(A, b, c)
    status1, T1, bas1, nb1 = simplex_fase2!(T1, bas1, nb1)

    if abs(T1[end,end]) > 1e-6
        return :inviavel, nothing, nothing
    end

    T2, bas2, nb2 = constroi_fase2!(T1, bas1, nb1, c, m)
    status2, T2, bas2, nb2 = simplex_fase2!(T2, bas2, nb2)

    if status2 == :ilimitado
        return :ilimitado, nothing, nothing
    else
        valor = T2[end,end]
        return :otima, valor, T2
    end
end

function main()
    fin = ARGS[1]
    fout = ARGS[2]

    data = open(fin, "r") do io
        read(io, String)
    end
    nums = split(data)
    pos = 1
    n = parse(Int, nums[pos]); pos+=1
    m = parse(Int, nums[pos]); pos+=1

    c = zeros(Float64, m)
    for j in 1:m
        c[j] = parse(Float64, nums[pos]); pos+=1
    end

    A = zeros(Float64, n, m)
    b = zeros(Float64, n)
    for i in 1:n
        for j in 1:m
            A[i,j] = parse(Float64, nums[pos]); pos+=1
        end
        b[i] = parse(Float64, nums[pos]); pos+=1
	    if b[i] < 0
            for j in 1:m
                A[i, j] = -A[i, j]
            end
            b[i] = -b[i]
        end
    end

    status, valor, _ = resolve_pl(A,b,c)

    open(fout, "w") do io
    	if status == :inviavel
            println(io, "inviavel")
        elseif status == :ilimitado
            println(io, "ilimitada")
        else
            @printf(io, "otima\n%.3f\n", -valor)
        end
    end
end

main()
