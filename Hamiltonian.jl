using AlgebraicPetri
using LinearAlgebra

function periodicHamiltonian(pn::AbstractLabelledReactionNet) 
    tm = TransitionMatrices(pn) 
    scc = tm.output - tm.input
    r = rates(pn)
    c = concentrations(pn)

    function H(k::Vector{Float64}) 
        n = 2^rank(scc)
        H = zeros(Complex{Float64}, n, n)

        for i in 1:nt(pn)
            for j in 1:n
                H[j,j] -= r[i] * prod(c[k]^tm.input[i,k] for k in 1:ns(pn))
            end
        end

        for i in 1:rank(scc)
            for j in findall(x->(1<= x%(2^i) <= 2^(i-1)), collect(1:n))
                H[j,j+2^(i-1)] += r[i] * prod(c[k]^tm.input[i,k] for k in 1:ns(pn)) 
                H[j+2^(i-1),j] += r[i] * prod(c[k]^tm.input[i,k] for k in 1:ns(pn)) * exp(im*k[i]) 
            end

            rev_i = i + rank(scc)

            for j in findall(x->(1<= x%(2^i) <= 2^(i-1)), collect(1:n))
                H[j+2^(i-1), j] += r[rev_i] * prod(c[k]^tm.input[rev_i,k] for k in 1:ns(pn)) 
                H[j, j+2^(i-1)] += r[rev_i] * prod(c[k]^tm.input[rev_i,k] for k in 1:ns(pn)) * exp(-im*k[i]) 
            end
        end

        return H
    end
    return H 
end

function bandplot(H, n::Int64) 
    k = collect(-π:0.01:π)

    f(k) = begin
        arg = zeros(n)
        arg[1] += k
        arg
    end

    p = plot(k, k->eigvals(H(f(k)))[1], xlims=[-π, π])

    for i in 2:2^n
        plot!(p, k, k->eigvals(H(f(k)))[i], xlims=[-π, π])
    end
    p

    # TODO: complex eigenvalue plotting
end
