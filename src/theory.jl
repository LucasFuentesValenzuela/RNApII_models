"""
Theoretical model from Lakatos, Chou, 2003.
"""


"""
    J(α, β, γ, L)

Theoretical model for the current density. 
    `α`: initiation rate
    `β`: elongation rate
    `γ`: termination rate
    `L`: particle footprint

Note: Eq (22) in Lakatos, Chou expresses quantities relative to the
elongation rate
"""
function J(α, β, γ, L)

    γ = check_γ(γ)

    α_ = α/β
    γ_ = γ/β

    regime = get_regime(α, β, γ, L)

	if regime == "entry-limited"
		return β*(α_ * (1-α_)/(1+(L-1)*α_))
    elseif regime == "exit-limited"
		return β * (γ_*(1-γ_)/(1+ (L-1)*γ_))
	elseif regime == "maxJ"
		return β*1/(L^(1/2) + 1)^2
	end
end

"""
    get_regime(α, β, γ, L)

Extract regime the system is in, according the Lakatos, Chou 2003.
    `α`: initiation rate
    `β`: elongation rate
    `γ`: termination rate
    `L`: particle footprint
"""
function get_regime(α, β, γ, L)

    γ = check_γ(γ)

    # critical values
	αc = 1/(1+L^(1/2))
    γc = αc

    α_ = α/β
    γ_ = γ/β

	if (α_ < αc) & (γ_ > α_) # entry-limited
        regime = "entry-limited"
    elseif (γ_ < γc) & (α_ > γ_) # exit-limited
        regime = "exit-limited"
	else # (α, γ) >= (αc, γc) : maximum current regime
        regime = "maxJ"
	end

    return regime

end

"""
    ρ(α, β, γ, L)

Theoretical model for the particle density. (not profiles, limit densities)
    `α`: initiation rate
    `β`: elongation rate
    `γ`: termination rate
    `L`: particle footprint

Note: Eq (25) in Lakatos, Chou expresses quantities relative to the
elongation rate
"""
function ρ(α_, β_, γ_, L)

    γ_ = check_γ(γ_)

    J_crt = J(α_, β_, γ_, L)/β_

    # critical values
	αc = 1/(1+L^(1/2))
    γc = αc

    α = α_/β_
    γ = γ_/β_


	if (α < αc) & (γ > α) # entry-limited

        ρL = α/(1+(L-1)*α)
        ρN = (1+(L-1) * J_crt  - sqrt((1+(L-1)*J_crt)^2 - 4L*J_crt))/2L
        ρR = α*(1-α)/γ/(1+(L-1)*α)
        
    elseif (γ < γc) & (α > γ) # exit-limiteimage.pngd

        ρL = 1/L - γ*(1-γ)/(α*L)/(1+(L-1)*γ)
        ρN = (1+(L-1) * J_crt  - sqrt((1+(L-1)*J_crt)^2 - 4L*J_crt))/2L
		ρR = (1-γ)/(1+(L-1)*γ)

	else # (α, γ) >= (αc, γc) : maximum current regime

        ρL = ((L^(1/2)+1)^2 - 1) / (α*L*(L^(1/2)+1)^2)
		ρN = 1/(L^(1/2) * (L^(1/2) + 1))
        ρR = 1/γ/(L^(1/2) + 1)^2 

	end

    return ρL, ρN, ρR
end

"""
Effective initiation rate with a promoter
"""
effective_α(kon, koff, α) = kon * α / (kon + koff + α)

"""
Assign value to γ if it is equal to `nothing`
"""
check_γ(γ) = γ === nothing ? LARGE_γ : γ

"""
Promoter occupancy
"""
ρp(kon, koff, α) = kon / (kon + koff + α)


"""
    ρp_MC(T)

Promoter occupancy based on a MC approximation for a transition matrix `T`.
"""
function ρp_MC(T)
    n = size(T, 1)
    sum(mean(T^1e3; dims=2)[Int(n/2)+1:end])
end

"""
    ρp_MC(α, β, kon, koff, Δt; order=4

Promoter occupancy based on a MC approximation with initiation rate `\alpha`, 
elongation rate `β`, `kon` and `koff` and timestep `Δt`. The order refers to the number of sites
that are modeled in the MC approximation
"""
ρp_MC(α, β, kon, koff, Δt; order=4) = ρp_MC(TransitionMatrix(α, β, kon, koff, Δt; order=order))

ρp_MC(α, β, kon, koff; order=4) = ρp_MC(α, β, kon, koff, set_Δt(α, β, β, kon, koff, nothing); order=order)

"""
    TransitionMatrix(α, β, kon, koff, Δt; order)

Transition Matrix for an MC approximation with initiation rate `\alpha`, 
elongation rate `β`, `kon` and `koff` and timestep `Δt`. The order refers to the number of sites
that are modeled in the MC approximation
"""
function TransitionMatrix(α, β, kon, koff, Δt; order)
	n = 2^order
	T = zeros(n, n)

	# 0-sites
	T[1, 1] = 1-kon * Δt
	T[Int(n/2)+1, 1] = kon * Δt
	for k in 2:Int(n/2)
		T[k-1, k] = β * Δt
		T[k, k] = 1 - (β + kon) * Δt
		T[Int(n/2) + k, k] = kon * Δt
	end

	# 1-sites
	T[1, Int(n/2)+1] = koff * Δt
	T[Int(n/4) + 1, Int(n/2) + 1] = α * Δt
	T[Int(n/2) + 1, Int(n/2) + 1] = 1 - (α + koff) * Δt
	for k in Int(n/2) + 2: Int(3*n/4)
		T[Int(n/4) + k - Int(n/2), k] = α * Δt
		T[Int(n/4) + k - Int(n/2) + 1, k] = β * Δt
		T[Int(n/4) + k - Int(n/2) + 2, k] = 1 - (α+β) * Δt
	end
	for k in Int(3n/4)+1:n
		T[k-1, k] = β*Δt
		T[k, k] = 1-β*Δt
	end
	return T
end
