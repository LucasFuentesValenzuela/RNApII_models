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
function effective_α(kon, koff, α)
    return kon * α / (kon + koff + α)
end

check_γ(γ) = γ === nothing ? LARGE_γ : γ