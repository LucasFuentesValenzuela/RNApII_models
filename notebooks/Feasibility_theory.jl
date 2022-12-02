### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 6fa6e566-eb31-4803-924f-8e4e1a10b9f6
begin
	using Revise
	import Pkg
	Pkg.activate()
	using RNApIIModels
end

# ╔═╡ 798003cc-abd9-42f2-8f28-a52713ac72f7
md"""# Description

We want to replicate the same type of analysis we do with the simulations purely from a theoretical basis. This allows faster iterations, results verification, and more fine-grained analysis

"""

# ╔═╡ 77bf5141-2715-4016-bd3d-e7aeaeee3a51
md"""# Load"""

# ╔═╡ 3cb614e0-cbfa-40d6-b2fd-5e2a925a3216
begin
	using Plots
	using PlutoUI
end

# ╔═╡ Cell order:
# ╠═798003cc-abd9-42f2-8f28-a52713ac72f7
# ╠═77bf5141-2715-4016-bd3d-e7aeaeee3a51
# ╠═6fa6e566-eb31-4803-924f-8e4e1a10b9f6
# ╠═3cb614e0-cbfa-40d6-b2fd-5e2a925a3216
