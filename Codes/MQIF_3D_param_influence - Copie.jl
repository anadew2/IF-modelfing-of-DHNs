### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ a389152e-ca91-4e87-b032-416b45376de5
begin
	using DifferentialEquations
	using LinearAlgebra
	using Plots
	using Statistics
end

# ╔═╡ 651f4fe0-a054-11eb-102f-5596e5a831f7
md" #### Packages"

# ╔═╡ 05814d8c-92a6-4a28-865c-d5b55d8a10d6
md" #### Problem parameters"

# ╔═╡ 0fffa3a2-76a1-4470-a5d2-4b78789cf6df
p=(4.7,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ abe05fb9-a447-41ce-8981-7e49dc35acb9
p_regen=(4.7,-40.0,-38.4,-30.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ acef61f0-7ad2-4eca-8049-537b1f283b61
function change_I(I,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==1
			p2[i]=I
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 98821bce-ccb4-403c-8af9-23021322e059
p_I_l = change_I(0.5,p)

# ╔═╡ 6496c76f-f3e3-4b2b-a1b9-6edf9d6d330c
p_I_h = change_I(11,p)

# ╔═╡ 089d1752-141c-46cb-8b64-23f0ad02c336
p_regen_I_l=change_I(0.5,p_regen)

# ╔═╡ 9d452ac9-c876-494a-bd55-92ba30ae71b6
p_regen_I_h=change_I(10,p_regen)

# ╔═╡ b55387e2-a88e-48d1-8c5e-fc4e0dddc14e
function change_gs(gs,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==7
			p2[i]=gs
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 21f17a1f-ba58-42ea-b3da-19e7f83ea34f
p_gs_l = change_gs(0.1,p)

# ╔═╡ a1fc8af4-2cd4-457e-a25c-ebcd78fe5636
p_gs_h = change_gs(0.9,p)

# ╔═╡ 2595be40-b236-46ef-9e88-f2c4225b3197
p_regen_gs_l = change_gs(0.1,p_regen)

# ╔═╡ 96f8c76b-1156-4b5c-b0c4-fd8e820662c6
p_regen_gs_h = change_gs(0.9,p_regen)

# ╔═╡ 352523e4-026b-40b8-a850-8610e2646a67
function change_gus(gus,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==8
			p2[i]=gus
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 9355dcb7-c3a1-4787-a026-5858ef83c52d
p_gus_h = change_gus(0.03,p)

# ╔═╡ 1884ccd5-b52e-439b-b0aa-fffad18c4cac
p_gus_l = change_gus(0.005,p)

# ╔═╡ 8d4e1961-ab4c-40fb-97a6-784c4a6057af
p_regen_gus_h = change_gus(0.018,p)

# ╔═╡ 841be4cd-afd5-4fa3-a560-3607049bb23f
p_regen_gus_l = change_gus(0.012,p)

# ╔═╡ 8bf784c8-3b88-46b4-8b62-dbb1a581402b
function change_vs0(vs0,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==3
			p2[i]=vs0
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ eb855792-b818-46b3-a02e-f4b97ccebfc4
p_vs0_l = change_vs0(-41,p)

# ╔═╡ e1726250-f07d-45b6-b82b-74f543f01215
p_vs0_h = change_vs0(-38.2,p)

# ╔═╡ 8c75d7e7-db77-43a4-b221-3adea64d2ef7
p_regen_vs0_l = change_vs0(-41,p_regen)

# ╔═╡ e72b85df-5951-4275-81cc-1f068d8946a0
p_regen_vs0_h = change_vs0(-38,p_regen)

# ╔═╡ ed1d0cdc-a5bc-48f5-b81c-d9e9d12bdf36
md" #### Nullclines functions"

# ╔═╡ 7e8b53bc-4800-46ba-bcff-55aca599d1cb
function Vnullcline1(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I) >= 0
		v1 = v0 + sqrt((gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I)/gf)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ 3cf9b4c3-bfe1-489e-a0b3-a0713114dab6
function Vnullcline2(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I) >= 0
		v2 = v0 - sqrt((gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I)/gf)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ a0726536-711c-4d21-a2c4-26767e00eeb7
function Vsnullcline(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	v = vs
	return v
end

# ╔═╡ 0cad5c79-4f97-4331-8781-e9544a08523e
function Vusnullcline(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	v = vus
	return v
end

# ╔═╡ 80235e6a-ef03-4068-84ec-02c618d60e4b
function V_Vs_nullcline1(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	US = gus*(vus-vus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US)

	if delta >= 0 #there is 1 (2) fixed point
		v1 = ((gf*v0)-(gs*vs0)+sqrt(delta))/(gf-gs)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ 86f8d8f0-2f6b-4ec1-900f-1d4e1f743324
function V_Vs_nullcline2(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	US = gus*(vus-vus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US)

	if delta >= 0 #there is 1 (2) fixed point
		v2 = ((gf*v0)-(gs*vs0)-sqrt(delta))/(gf-gs)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ ae5c617e-aad0-4321-b76e-5ff38b554646
function V_Vus_nullcline1(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	S = gs*(vs-vs0)^2
	delta = gf*gus*(v0-vus0)^2 - (gf-gus)*(I-S)

	if delta >= 0 #there is 1 (2) fixed point
		v1 = ((gf*v0)-(gus*vus0)+sqrt(delta))/(gf-gus)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ 9ee2c719-b34a-41be-b951-212303f2e58d
function V_Vus_nullcline2(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	S = gs*(vs-vs0)^2
	delta = gf*gus*(v0-vus0)^2 - (gf-gus)*(I-S)

	if delta >= 0 #there is 1 (2) fixed point
		v2 = ((gf*v0)-(gus*vus0)-sqrt(delta))/(gf-gus)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ 0da094a5-2899-4dc0-a03d-be4997b1e507
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# ╔═╡ 4af9fe12-05fa-4dfa-84af-9fe5b739410c
md" #### Bifurcation functions defined based on analytical expression"

# ╔═╡ 6e504105-ac3a-419d-9ee1-89f15ae2445c
function bifI3D(p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	Ibif = (gf*gs*(v0-vs0)^2+gf*gus*(v0-vus0)^2-gs*gus*(vs0-vus0)^2)/(gf-gs-gus)
	return Ibif
end

# ╔═╡ 4218be3f-9b82-476f-b8c1-21a90dad00bb
function jacobian(v,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	J = zeros(3,3)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[1,3] = -2*gus*(v-vus0)
	J[2,1] = 1
	J[2,2] = -1
	J[2,3] = 0
	J[3,1] = 1
	J[3,2] = 0
	J[3,3] = -1

	return J
end

# ╔═╡ b34deb54-5b25-4050-98ba-8147a1f27d8c
function fp_stability_3D(fp,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	stability = ones(1,2)*66 #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable
	for i=1:length(fp)
		J = jacobian(fp[i],p)
		lambda1,lambda2,lambda3 = eigvals(J)
		if real(lambda1)>0 && real(lambda2)>0 && real(lambda3)>0
			stability[i]=1.0
		end
		if (real(lambda1)>0 && (real(lambda2)<0 || real(lambda3)<0)) || (real(lambda1)<0 && (real(lambda2)>0 || real(lambda3)>0))
			stability[i]=0.0
		end
		if real(lambda1)<0 && real(lambda2)<0 && real(lambda3)<0
			stability[i]=-1.0
		end
	end
	return stability
end

# ╔═╡ 86755dcd-d1a5-423d-ac23-8d2615afddda
function global_I_bifurcation(p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	common_term = gf*v0 - gs*vs0 - gus*vus0
	delta = gf*gs*(v0-vs0)^2 + gf*gus*(v0-vus0)^2 - gs*gus*(vs0-vus0)^2 -(gf-gs-gus)*I

	if delta >= 0
		numerator1 = common_term + sqrt(delta)
		numerator2 = common_term - sqrt(delta)
		denominator = gf-gs-gus
		V_fp_1 = numerator1/denominator
		V_fp_2 = numerator2/denominator
		V_fp = [V_fp_1,V_fp_2]

		stability = fp_stability_3D(V_fp,p)
	else
		V_fp_1 = NaN
		V_fp_2 = NaN
		V_fp = [V_fp_1,V_fp_2]

		stability = [NaN,NaN]
	end
	return V_fp,stability
end

# ╔═╡ bead8814-21a5-4a7a-8488-d5493e14b391
md" #### ODE problem functions"

# ╔═╡ 53e4acac-b36e-48a0-b32a-4ff6bb68f67b
function MQIF_3D!(du,u,p,t)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 - gus*(u[3]-vus0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
	du[3] = (u[1]-u[3])/tus
end

# ╔═╡ f2107a48-c7ee-46db-9acd-44e16488ca0a
begin
	Vmax = 30 #30
	Vr = -40 #-40
	Vsr = -35 #-20
	DVusr = 3

	md"""Voltages used for reset"""
end

# ╔═╡ 8d01dc21-6666-42eb-9b7f-ac05be91b1d8
function spike(x)  # spikes when spike(x) goes from negative to positive
    (x[1] - Vmax)
end

# ╔═╡ 6d652924-4231-4dc4-81d1-f6a438aea675
# event when event_f(u,t) == 0
function condition(x,t,integrator) #
    spike(x)
end

# ╔═╡ c43b5adc-0253-4b9d-835b-3ec44362ed03
function reset!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
	x[3] = x[3] + DVusr
end

# ╔═╡ 24065204-66b6-4dbb-9656-0c252b014cce
# when condition == 0 and upcrossing (from negative to positive)
function affect!(integrator)
    reset!(integrator.u)
end

# ╔═╡ 8035aaab-6b4f-49ef-874a-06989d53cc41
md" #### Simulations"

# ╔═╡ b342899c-0cd3-4b35-b261-0807216de37e
begin
	v_min = -100
	v_max = 20
	vs_min = -100
	vs_max = 0
	vus_min = -100
	vus_max = 0
	vus_min_ = -80
	vus_max_ = -20

	lim_hole = 1
	vs_inter_l = collect(range(vs_min,stop=p[3]-4,length=150))
	vs_inter_h = collect(range(p[3]+4,stop=vs_max,length=150))
	vs_hole = collect(range(p[3]-4,stop=p[3]+4,length=600))

	lim_vus_hole = -50 + sqrt(lim_hole+p[1]/p[8]) + sqrt(lim_hole+p[1]/p[8])/10
	vus_inter_l = collect(range(vus_min,stop=-50,length=150))
	vus_inter_h = collect(range(lim_vus_hole,stop=vus_max,length=150))
	vus_hole = collect(range(-50,stop=lim_vus_hole,length=600))

	vs = zeros(length(vs_inter_l)+length(vs_inter_h)+length(vs_hole))
	vus = zeros(length(vus_inter_l)+length(vus_inter_h)+length(vus_hole))

	for i=1:(length(vs))
		if i <= length(vs_inter_l)
			vs[i] = vs_inter_l[i]
		else
			if i<= length(vs_inter_l)+length(vs_hole)
				vs[i] = vs_hole[i-length(vs_inter_l)]
			else
				vs[i] = vs_inter_h[i-(length(vs_inter_l)+length(vs_hole))]
			end
		end
	end

	for i=1:(length(vus))
		if i <= length(vus_inter_l)
			vus[i] = vus_inter_l[i]
		else
			if i<= length(vus_inter_l)+length(vus_hole)
				vus[i] = vus_hole[i-length(vus_inter_l)]
			else
				vus[i] = vus_inter_h[i-(length(vus_inter_l)+length(vus_hole))]
			end
		end
	end
	md"""Nullclines calculation"""
end

# ╔═╡ 29bb5f71-a408-479e-83ed-69a62e6bba29
md"###### Current variation for restorarive ultraslow feedback and regenerative slow feedback"

# ╔═╡ 01252cf7-d74c-4f02-977a-fcbf743dc056
begin
	u0=[-40,-40,-40]
	tspan = (0.0,8000.0)
	cb   = ContinuousCallback(condition,affect!,nothing)
	prob = ODEProblem(MQIF_3D!,u0,tspan,p,callback=cb)
	sol = solve(prob,dense=false)
end

# ╔═╡ 82eb0286-c992-4459-9d65-529e8ba30362
gr()

# ╔═╡ 2eb57d95-e17c-4f71-8f9f-34a6df5b0392
begin
	plt_time = plot(sol.t,sol[1,:],linecolor="blue",label="V(t)")
	plot!(sol.t,sol[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol.t,sol[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol.t,p[3]ones(size(sol.t)),linecolor="red",label="Vs0")
	plot!(sol.t,p[4]ones(size(sol.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("I=$(p[1])")
end

# ╔═╡ 02927149-ba9f-4da7-83b1-390a34d6f62f
md"""Due to the reset, vus is sometimes brought a little bit too high making the current more negative. Thus, the trajectory is sucked on the left part of the v nullcline for It<0 making v to decrease a little more. To contrast with the previous case, the gradient of Vus is higher here due to the high I. Therefore, the current I_t varies faster. The trajectory continues to evolve slowly with v increases but stay lower than V0 (left part of the phase plane). At some point, due to the highly varying current, the trajectory is repelled by the saddle and the system spikes.  """

# ╔═╡ 17518825-5527-4478-a4c3-52f7ff864a87
md"""Changing I moves the system on a SN bifurcation"""

# ╔═╡ 70af372c-963d-4e53-b162-f3b2ca733b48
begin
	V1_null_(a,b) = Vnullcline1(b,a,p)
	V2_null_(a,b) = Vnullcline2(b,a,p)
	Vs_null_(a,b) = Vsnullcline(b,a,p)
	Vus_null_(a,b) = Vusnullcline(b,a,p)
	
	inter_V_Vs_null1 = zeros(length(vus))
	inter_V_Vs_null2 = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1[i] = V_Vs_nullcline1(NaN,vus[i],p)
		inter_V_Vs_null2[i] = V_Vs_nullcline2(NaN,vus[i],p)
	end

	inter_V_Vus_null1 = zeros(size(vs))
	inter_V_Vus_null2 = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1[i] = V_Vus_nullcline1(vs[i],NaN,p)
		inter_V_Vus_null2[i] = V_Vus_nullcline2(vs[i],NaN,p)
	end


	plt_v_null = plot(vus,vs,V1_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20)) #yellow-orange
	plot!(vus,vs,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1],[vs],[inter_V_Vus_null1],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol[3,1:end],sol[2,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p[1])")

end

# ╔═╡ 1cfe6b28-6d83-4bfa-920a-b4eb387590cb
md"""The small jump in the trajectory before being sucked into the hole seems to be due to the fact that the total current I_t becomes negative and make the system converge to the Vnullcline as at this moment, the gradient of Vs is still high enough. It will then decreases as the trajectory is brought close to the Vs nullcline. The system reaches a plateau as the place where the gradient of Vs becomes close to zero is the same place where the gradient of Vus becomes close to zero, the system is therefore "blocked" for a while""" 

# ╔═╡ cfa819c9-a58a-44b3-bffd-3405fabaf6cc
md"""Changing I increases the size of the hole where the system is excitable and therefore moves a lot the intersection with the Vs nullcline"""

# ╔═╡ 49b1f686-bcf7-49ed-93a3-aa74d0081110
gr()

# ╔═╡ 6ec8f2d4-8de5-480b-b8e4-1a3b49fa4a4a
md"###### Current variation for regenerative ultraslow feedback and regenerative slow feedback"

# ╔═╡ f1aba699-0d2c-4835-ab43-42bea662a141
begin
	prob_regen = ODEProblem(MQIF_3D!,u0,tspan,p_regen,callback=cb)
	sol_regen = solve(prob_regen,dense=false)
end

# ╔═╡ b990cc17-1ea4-468a-8a2e-39f897cade99
begin
	plt_time_regen = plot(sol_regen.t,sol_regen[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen.t,sol_regen[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen.t,sol_regen[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen.t,p_regen[3]ones(size(sol_regen.t)),linecolor="red",label="Vs0")
	plot!(sol_regen.t,p_regen[4]ones(size(sol_regen.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("I=$(p_regen[1])")
end

# ╔═╡ c67e6b50-ee12-4361-ac57-1e602641fad6
md"""When I increases, the limit cycle frequency increases"""

# ╔═╡ 6a379d9f-05c3-4a3c-935b-b2b8bc96d43f
begin
	V1_null_regen_(a,b) = Vnullcline1(b,a,p_regen)
	V2_null_regen_(a,b) = Vnullcline2(b,a,p_regen)
	Vs_null_regen_(a,b) = Vsnullcline(b,a,p_regen)
	Vus_null_regen_(a,b) = Vusnullcline(b,a,p_regen)
	
	inter_V_Vs_null1_regen = zeros(length(vus))
	inter_V_Vs_null2_regen = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen[i] = V_Vs_nullcline1(NaN,vus[i],p_regen)
		inter_V_Vs_null2_regen[i] = V_Vs_nullcline2(NaN,vus[i],p_regen)
	end

	inter_V_Vus_null1_regen = zeros(size(vs))
	inter_V_Vus_null2_regen = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_regen[i] = V_Vus_nullcline1(vs[i],NaN,p_regen)
		inter_V_Vus_null2_regen[i] = V_Vus_nullcline2(vs[i],NaN,p_regen)
	end


	plt_v_null_regen = plot(vus,vs,V1_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen[3,1:end],sol_regen[2,1:end],sol_regen[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p_regen[1])")

end

# ╔═╡ 7371f342-f293-497d-b5c4-bbd24f02967c
gr()

# ╔═╡ b7a260fd-9f58-4fd8-a5c2-45997bf62bc8
md"""Changing I increases the size of the hole where the system is excitable and therefore moves a lot the intersection with the Vs nullcline"""

# ╔═╡ 0365925b-e208-4f13-8908-747f238a2a4f
md"###### GS variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ e05e543d-a6ba-4a26-a61c-9c05daba95c5
md"""When gs increases, the V nullcline is brought closer to the plane Vs=Vs0

When gs decreases, it means that the V nullcline becomes more sensitive to V than to Vs or Vus which is seen as the fact that the nullcline contracts towards the plane V=V0"""

# ╔═╡ c7875827-d1a8-4058-be16-e4a1437cbc6a
md"###### GS variation for regenerative ultraslow feedback and regenerative slow feedback"

# ╔═╡ 9d900972-fc6d-49b0-b1e3-3a980deb2b97
md"""gs seems to have an impact mostly on the burst frequency before converging towards the limit cycle"""

# ╔═╡ ee174df3-a596-441e-a3d0-808013360e31
gr()

# ╔═╡ 298c4acd-c3a5-4f88-89e6-8d5b5fb249b0
md"###### gus variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ 509fa73f-5b82-42ab-a417-82ab648c1eed
begin
	prob_gus_l = ODEProblem(MQIF_3D!,u0,tspan,p_gus_l,callback=cb)
	sol_gus_l = solve(prob_gus_l,dense=false)
end

# ╔═╡ a7648315-fbe4-4bd2-a3e2-0ad038fefbbe
begin
	prob_gus_h = ODEProblem(MQIF_3D!,u0,tspan,p_gus_h,callback=cb)
	sol_gus_h = solve(prob_gus_h,dense=false)
end

# ╔═╡ 1768e5d7-82de-46cf-8206-a61bb392c781
begin
	plot(sol_gus_l.t,sol_gus_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gus_l.t,sol_gus_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gus_l.t,sol_gus_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gus_l.t,p[3]ones(size(sol_gus_l.t)),linecolor="red",label="Vs0")
	plot!(sol_gus_l.t,p[4]ones(size(sol_gus_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_gus_l[8])")
end

# ╔═╡ 5f82aa48-4287-41bb-914e-ad9f62a2a270
begin
	plot(plt_time)
	title!("gus=$(p[8])")
end

# ╔═╡ dadccbbb-a331-424f-9c97-a7753b8cf92a
begin
	plot(sol_gus_h.t,sol_gus_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gus_h.t,sol_gus_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gus_h.t,sol_gus_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gus_h.t,p[3]ones(size(sol_gus_h.t)),linecolor="red",label="Vs0")
	plot!(sol_gus_h.t,p[4]ones(size(sol_gus_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_gus_h[8])")
end

# ╔═╡ 3dc77dcb-b68f-4006-bc0d-9dd39cbf2f38
md"""The overall effect of modifying gus is that it changes the range of value for the excitability hole. Gus is inversely porpotional to its size. If it decreases, vus will tend to increase more which accelerates the ocsillations. On the contrary, if it increases too much, the hole becomes too small and the system have a stable node"""

# ╔═╡ b5d1876f-c572-4e77-b694-f8270162f03f
begin
	V1_null_gus_l_(a,b) = Vnullcline1(b,a,p_gus_l)
	V2_null_gus_l_(a,b) = Vnullcline2(b,a,p_gus_l)
	Vs_null_gus_l_(a,b) = Vsnullcline(b,a,p_gus_l)
	Vus_null_gus_l_(a,b) = Vusnullcline(b,a,p_gus_l)
	
	inter_V_Vs_null1_gus_l = zeros(length(vus))
	inter_V_Vs_null2_gus_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_gus_l[i] = V_Vs_nullcline1(NaN,vus[i],p_gus_l)
		inter_V_Vs_null2_gus_l[i] = V_Vs_nullcline2(NaN,vus[i],p_gus_l)
	end

	inter_V_Vus_null1_gus_l = zeros(size(vs))
	inter_V_Vus_null2_gus_l = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_gus_l[i] = V_Vus_nullcline1(vs[i],NaN,p_gus_l)
		inter_V_Vus_null2_gus_l[i] = V_Vus_nullcline2(vs[i],NaN,p_gus_l)
	end


	plot(vus,vs,V1_null_gus_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_gus_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_gus_l,inter_V_Vs_null2_gus_l],[inter_V_Vs_null1_gus_l,inter_V_Vs_null2_gus_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_gus_l,inter_V_Vus_null2_gus_l],[vs,vs],[inter_V_Vus_null1_gus_l,inter_V_Vus_null2_gus_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_gus_l[3,1:end],sol_gus_l[2,1:end],sol_gus_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gus=$(p_gus_l[8])")

end

# ╔═╡ ee04cc5e-9fdc-4468-a589-7645343b7682
begin
	plot(plt_v_null)
	title!("V nullcline when gus=$(p[8])")
end

# ╔═╡ 672a0a95-82e3-4d42-ab3f-f2cf26510652
begin
	V1_null_gus_h_(a,b) = Vnullcline1(b,a,p_gus_h)
	V2_null_gus_h_(a,b) = Vnullcline2(b,a,p_gus_h)
	Vs_null_gus_h_(a,b) = Vsnullcline(b,a,p_gus_h)
	Vus_null_gus_h_(a,b) = Vusnullcline(b,a,p_gus_h)
	
	inter_V_Vs_null1_gus_h = zeros(length(vus))
	inter_V_Vs_null2_gus_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_gus_h[i] = V_Vs_nullcline1(NaN,vus[i],p_gus_h)
		inter_V_Vs_null2_gus_h[i] = V_Vs_nullcline2(NaN,vus[i],p_gus_h)
	end

	inter_V_Vus_null1_gus_h = zeros(size(vs))
	inter_V_Vus_null2_gus_h = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_gus_h[i] = V_Vus_nullcline1(vs[i],NaN,p_gus_h)
		inter_V_Vus_null2_gus_h[i] = V_Vus_nullcline2(vs[i],NaN,p_gus_h)
	end


	plot(vus,vs,V1_null_gus_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_gus_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_gus_h,inter_V_Vs_null2_gus_h],[inter_V_Vs_null1_gus_h,inter_V_Vs_null2_gus_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_gus_h,inter_V_Vus_null2_gus_h],[vs,vs],[inter_V_Vus_null1_gus_h,inter_V_Vus_null2_gus_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_gus_h[3,1:end],sol_gus_h[2,1:end],sol_gus_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gus=$(p_gus_h[8])")

end

# ╔═╡ 0ed3d289-8c77-47e5-8972-a5fceb9a5fe2
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_gus_h,inter_V_Vus_null2_gus_h],[inter_V_Vus_null1_gus_h,inter_V_Vus_null2_gus_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; gus=$(p_gus_h[8])",legend=:outertopright) #pink
	plot!([inter_V_Vus_null1_gus_l,inter_V_Vus_null2_gus_l],[inter_V_Vus_null1_gus_l,inter_V_Vus_null2_gus_l],[vs,vs],linecolor=RGB(0.0,0.5,1),linewidth=3,label="dV/dt =dVus/dt = 0; gus=$(p_gus_l[8])",legend=:outertopright)
	title!("Regenerative slow feedback")
end

# ╔═╡ 5e07a6ae-a8f5-480e-8d54-369b94282cee
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_gus_h,inter_V_Vs_null2_gus_h],[inter_V_Vs_null1_gus_h,inter_V_Vs_null2_gus_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; gus=$(p_gus_h[8])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_gus_l,inter_V_Vs_null2_gus_l],[inter_V_Vs_null1_gus_l,inter_V_Vs_null2_gus_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; gus=$(p_gus_l[8])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ 86e7ab9d-4544-400d-b7b0-6c4fea65e264
begin
	prob_regen_gus_l = ODEProblem(MQIF_3D!,u0,tspan,p_regen_gus_l,callback=cb)
	sol_regen_gus_l = solve(prob_regen_gus_l,dense=false)
end

# ╔═╡ 879293d4-1d92-44da-923b-3e2cc9e5bc60
begin
	prob_regen_gus_h = ODEProblem(MQIF_3D!,u0,tspan,p_regen_gus_h,callback=cb)
	sol_regen_gus_h = solve(prob_regen_gus_h,dense=false)
end

# ╔═╡ a0bfcfa1-5437-468a-ae71-d5d36c7e21b7
begin
	plot(sol_regen_gus_l.t,sol_regen_gus_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_gus_l.t,sol_regen_gus_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_gus_l.t,sol_regen_gus_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_gus_l.t,p[3]ones(size(sol_regen_gus_l.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_gus_l.t,p[4]ones(size(sol_regen_gus_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_regen_gus_l[8])")
end

# ╔═╡ 403bbf12-9c0d-459e-89d0-54b9bf6a524d
begin
	plot(plt_time_regen)
	title!("gus=$(p_regen[8])")
end

# ╔═╡ e15ecbf6-59ad-4675-87ea-f9803d5883f1
begin
	plot(sol_regen_gus_h.t,sol_regen_gus_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_gus_h.t,sol_regen_gus_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_gus_h.t,sol_regen_gus_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_gus_h.t,p[3]ones(size(sol_regen_gus_h.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_gus_h.t,p[4]ones(size(sol_regen_gus_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_regen_gus_h[8])")
end

# ╔═╡ e42297bc-471f-4abf-9939-28a31f67f5a2
md"""The overall effect of modifying gus is that it changes the range of value for the excitability hole. Gus is inversely porpotional to its size. If it decreases, vus will tend to increase more which accelerates the ocsillations. On the contrary, if it increases too much, the hole becomes too small and the system have a stable node"""

# ╔═╡ f241d2d3-38a9-441d-916b-51289a4a1e61
begin
	V1_null_regen_gus_l_(a,b) = Vnullcline1(b,a,p_regen_gus_l)
	V2_null_regen_gus_l_(a,b) = Vnullcline2(b,a,p_regen_gus_l)
	Vs_null_regen_gus_l_(a,b) = Vsnullcline(b,a,p_regen_gus_l)
	Vus_null_regen_gus_l_(a,b) = Vusnullcline(b,a,p_regen_gus_l)
	
	inter_V_Vs_null1_regen_gus_l = zeros(length(vus))
	inter_V_Vs_null2_regen_gus_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_gus_l[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_gus_l)
		inter_V_Vs_null2_regen_gus_l[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_gus_l)
	end

	inter_V_Vus_null1_regen_gus_l = zeros(size(vs))
	inter_V_Vus_null2_regen_gus_l = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_regen_gus_l[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_gus_l)
		inter_V_Vus_null2_regen_gus_l[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_gus_l)
	end


	plot(vus,vs,V1_null_regen_gus_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_gus_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_gus_l,inter_V_Vs_null2_regen_gus_l],[inter_V_Vs_null1_regen_gus_l,inter_V_Vs_null2_regen_gus_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_gus_l,inter_V_Vus_null2_regen_gus_l],[vs,vs],[inter_V_Vus_null1_regen_gus_l,inter_V_Vus_null2_regen_gus_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_gus_l[3,1:end],sol_regen_gus_l[2,1:end],sol_regen_gus_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gus=$(p_regen_gus_l[8])")

end

# ╔═╡ d5c6139f-7125-468f-a21b-576970c7b478
begin
	plot(plt_v_null_regen)
	title!("V nullcline when gus=$(p_regen[8])")
end

# ╔═╡ 1f77716b-61f6-48ca-ae45-1616771cd5db
begin
	V1_null_regen_gus_h_(a,b) = Vnullcline1(b,a,p_regen_gus_h)
	V2_null_regen_gus_h_(a,b) = Vnullcline2(b,a,p_regen_gus_h)
	Vs_null_regen_gus_h_(a,b) = Vsnullcline(b,a,p_regen_gus_h)
	Vus_null_regen_gus_h_(a,b) = Vusnullcline(b,a,p_regen_gus_h)
	
	inter_V_Vs_null1_regen_gus_h = zeros(length(vus))
	inter_V_Vs_null2_regen_gus_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_gus_h[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_gus_h)
		inter_V_Vs_null2_regen_gus_h[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_gus_h)
	end

	inter_V_Vus_null1_regen_gus_h = zeros(size(vs))
	inter_V_Vus_null2_regen_gus_h = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_regen_gus_h[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_gus_h)
		inter_V_Vus_null2_regen_gus_h[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_gus_h)
	end


	plot(vus,vs,V1_null_regen_gus_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_gus_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_gus_h,inter_V_Vs_null2_regen_gus_h],[inter_V_Vs_null1_regen_gus_h,inter_V_Vs_null2_regen_gus_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_gus_h,inter_V_Vus_null2_regen_gus_h],[vs,vs],[inter_V_Vus_null1_regen_gus_h,inter_V_Vus_null2_regen_gus_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_gus_h[3,1:end],sol_regen_gus_h[2,1:end],sol_regen_gus_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gus=$(p_regen_gus_h[8])")

end

# ╔═╡ 55bcf2d5-0b82-45ce-9689-fe0165513784
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_regen_gus_h,inter_V_Vus_null2_regen_gus_h],[inter_V_Vus_null1_regen_gus_h,inter_V_Vus_null2_regen_gus_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; gus=$(p_regen_gus_h[8])",legend=:outertopright) #pink
	plot!([inter_V_Vus_null1_regen_gus_l,inter_V_Vus_null2_regen_gus_l],[inter_V_Vus_null1_regen_gus_l,inter_V_Vus_null2_regen_gus_l],[vs,vs],linecolor=RGB(0.0,0.5,1),linewidth=3,label="dV/dt =dVus/dt = 0; gus=$(p_regen_gus_l[8])",legend=:outertopright)
	title!("Regenerative slow feedback")
end

# ╔═╡ 002db113-6b87-47ec-9a3d-009a6143ec52
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_regen_gus_h,inter_V_Vs_null2_regen_gus_h],[inter_V_Vs_null1_regen_gus_h,inter_V_Vs_null2_regen_gus_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; gus=$(p_regen_gus_h[8])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_regen_gus_l,inter_V_Vs_null2_regen_gus_l],[inter_V_Vs_null1_regen_gus_l,inter_V_Vs_null2_regen_gus_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; gus=$(p_regen_gus_l[8])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ 23bbe085-afec-42e4-99ee-b863ce68d62f
md"###### Vs0 variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ cb445d6e-ac17-4137-985b-f73ac538dfbf
begin
	prob_vs0_l = ODEProblem(MQIF_3D!,u0,tspan,p_vs0_l,callback=cb)
	sol_vs0_l = solve(prob_vs0_l,dense=false)
end

# ╔═╡ e75132c6-17ef-4bb2-ab33-9550bc88091e
begin
	prob_vs0_h = ODEProblem(MQIF_3D!,u0,tspan,p_vs0_h,callback=cb)
	sol_vs0_h = solve(prob_vs0_h,dense=false)
end

# ╔═╡ 04500988-7e91-4a99-9262-7e46d02a732e
begin
	plot(sol_vs0_l.t,sol_vs0_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_vs0_l.t,sol_vs0_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_vs0_l.t,sol_vs0_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_vs0_l.t,p[3]ones(size(sol_vs0_l.t)),linecolor="red",label="Vs0")
	plot!(sol_vs0_l.t,p[4]ones(size(sol_vs0_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("vs0=$(p_vs0_l[3])")
end

# ╔═╡ 13ed704b-4ed1-47d0-a384-c66a8de1d0a5
begin
	plot(plt_time)
	title!("vs0=$(p[3])")
end

# ╔═╡ 0da78619-4c28-4c41-a5ef-2fd0fc349e35
begin
	plot(sol_vs0_h.t,sol_vs0_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_vs0_h.t,sol_vs0_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_vs0_h.t,sol_vs0_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_vs0_h.t,p[3]ones(size(sol_vs0_h.t)),linecolor="red",label="Vs0")
	plot!(sol_vs0_h.t,p[4]ones(size(sol_vs0_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("vs0=$(p_vs0_h[3])")
end

# ╔═╡ 129badd2-31e9-4dbd-aa02-f6b0df5b10d6
md"""Changing vs0 moves the whole plane to the right if vs0 increases and to the left otherwise on the vs axis. As the other plane where dVs/dt =0 or dVus/dt=0 do not move, the projection of the nullclines inersections in the plane V;Vs and V;Vus move. 

Having a slow feedback too regenerative seems to slow the system as it can make the 3 nullclines to intersect. If the slow feedback becomes more and more resorative, the system accelerates due to the instable node created seen in the view V;Vs and then converge towards a stable state due to the stable node created""" 

# ╔═╡ 0df08dc3-3dc9-4255-895e-b20f21243bf5
begin
	V1_null_vs0_l_(a,b) = Vnullcline1(b,a,p_vs0_l)
	V2_null_vs0_l_(a,b) = Vnullcline2(b,a,p_vs0_l)
	Vs_null_vs0_l_(a,b) = Vsnullcline(b,a,p_vs0_l)
	Vus_null_vs0_l_(a,b) = Vusnullcline(b,a,p_vs0_l)
	
	inter_V_Vs_null1_vs0_l = zeros(length(vus))
	inter_V_Vs_null2_vs0_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_vs0_l[i] = V_Vs_nullcline1(NaN,vus[i],p_vs0_l)
		inter_V_Vs_null2_vs0_l[i] = V_Vs_nullcline2(NaN,vus[i],p_vs0_l)
	end

	inter_V_Vus_null1_vs0_l = zeros(size(vs))
	inter_V_Vus_null2_vs0_l = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_vs0_l[i] = V_Vus_nullcline1(vs[i],NaN,p_vs0_l)
		inter_V_Vus_null2_vs0_l[i] = V_Vus_nullcline2(vs[i],NaN,p_vs0_l)
	end


	plot(vus,vs,V1_null_vs0_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_vs0_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_vs0_l,inter_V_Vs_null2_vs0_l],[inter_V_Vs_null1_vs0_l,inter_V_Vs_null2_vs0_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_vs0_l,inter_V_Vus_null2_vs0_l],[vs,vs],[inter_V_Vus_null1_vs0_l,inter_V_Vus_null2_vs0_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_vs0_l[3,1:end],sol_vs0_l[2,1:end],sol_vs0_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when vs0=$(p_vs0_l[3])")

end

# ╔═╡ 8bde95e1-f76c-4943-8d1e-5b43ac2defdf
begin
	plot(plt_v_null)
	title!("V nullcline when vs0=$(p[3])")
end

# ╔═╡ 9a6033ce-60df-40a3-a823-49ff3b3c5261
begin
	V1_null_vs0_h_(a,b) = Vnullcline1(b,a,p_vs0_h)
	V2_null_vs0_h_(a,b) = Vnullcline2(b,a,p_vs0_h)
	Vs_null_vs0_h_(a,b) = Vsnullcline(b,a,p_vs0_h)
	Vus_null_vs0_h_(a,b) = Vusnullcline(b,a,p_vs0_h)
	
	inter_V_Vs_null1_vs0_h = zeros(length(vus))
	inter_V_Vs_null2_vs0_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_vs0_h[i] = V_Vs_nullcline1(NaN,vus[i],p_vs0_h)
		inter_V_Vs_null2_vs0_h[i] = V_Vs_nullcline2(NaN,vus[i],p_vs0_h)
	end

	inter_V_Vus_null1_vs0_h = zeros(size(vs))
	inter_V_Vus_null2_vs0_h = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_vs0_h[i] = V_Vus_nullcline1(vs[i],NaN,p_vs0_h)
		inter_V_Vus_null2_vs0_h[i] = V_Vus_nullcline2(vs[i],NaN,p_vs0_h)
	end


	plot(vus,vs,V1_null_vs0_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_vs0_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_vs0_h,inter_V_Vs_null2_vs0_h],[inter_V_Vs_null1_vs0_h,inter_V_Vs_null2_vs0_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_vs0_h,inter_V_Vus_null2_vs0_h],[vs,vs],[inter_V_Vus_null1_vs0_h,inter_V_Vus_null2_vs0_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_vs0_h[3,1:end],sol_vs0_h[2,1:end],sol_vs0_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when vs0=$(p_vs0_h[3])")

end

# ╔═╡ 97bac848-795e-4770-8ea0-3cb021823bf1
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_vs0_h,inter_V_Vus_null2_vs0_h],[inter_V_Vus_null1_vs0_h,inter_V_Vus_null2_vs0_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; vs0=$(p_vs0_h[3])",legend=:outertopright) #pink
	plot!([inter_V_Vus_null1_vs0_l,inter_V_Vus_null2_vs0_l],[inter_V_Vus_null1_vs0_l,inter_V_Vus_null2_vs0_l],[vs,vs],linecolor=RGB(0.0,0.5,1),linewidth=3,label="dV/dt =dVus/dt = 0; vs0=$(p_vs0_l[3])",legend=:outertopright)
	title!("Regenerative slow feedback")
end

# ╔═╡ f380fa87-c7a9-4fac-9be0-7987d8ca8960
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_vs0_h,inter_V_Vs_null2_vs0_h],[inter_V_Vs_null1_vs0_h,inter_V_Vs_null2_vs0_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; vs0=$(p_vs0_h[3])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_vs0_l,inter_V_Vs_null2_vs0_l],[inter_V_Vs_null1_vs0_l,inter_V_Vs_null2_vs0_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; vs0=$(p_vs0_l[3])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ 88c7243f-d325-459d-9297-abc999bbb269
begin
	prob_regen_vs0_l = ODEProblem(MQIF_3D!,u0,tspan,p_regen_vs0_l,callback=cb)
	sol_regen_vs0_l = solve(prob_regen_vs0_l,dense=false)
end

# ╔═╡ 6a48b3b3-6204-4127-b431-c7125321adcf
begin
	prob_regen_vs0_h = ODEProblem(MQIF_3D!,u0,tspan,p_regen_vs0_h,callback=cb)
	sol_regen_vs0_h = solve(prob_regen_vs0_h,dense=false)
end

# ╔═╡ 5a474cef-3084-45ed-8b83-bd65adfbfe6b
begin
	plot(sol_regen_vs0_l.t,sol_regen_vs0_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_vs0_l.t,sol_regen_vs0_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_vs0_l.t,sol_regen_vs0_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_vs0_l.t,p[3]ones(size(sol_regen_vs0_l.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_vs0_l.t,p[4]ones(size(sol_regen_vs0_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("vs0=$(p_regen_vs0_l[3])")
end

# ╔═╡ 1090cbfa-667c-474d-b00f-e2db4e37cd90
begin
	plot(plt_time_regen)
	title!("vs0=$(p_regen[3])")
end

# ╔═╡ 0a89f507-e3ae-4d4d-8eff-f6bbdb822a92
begin
	plot(sol_regen_vs0_h.t,sol_regen_vs0_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_vs0_h.t,sol_regen_vs0_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_vs0_h.t,sol_regen_vs0_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_vs0_h.t,p[3]ones(size(sol_regen_vs0_h.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_vs0_h.t,p[4]ones(size(sol_regen_vs0_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("vs0=$(p_regen_vs0_h[3])")
end

# ╔═╡ 32449dbb-0486-458e-812a-bf343f23bfc4
md"""Changing vs0 moves the whole plane to the right if vs0 increases and to the left otherwise on the vs axis. As the other plane where dVs/dt =0 or dVus/dt=0 do not move, the projection of the nullclines inersections in the plane V;Vs and V;Vus move. 

Having a slow feedback too regenerative seems to slow the system as it can make the 3 nullclines to intersect. If the slow feedback becomes more and more resorative, the system accelerates due to the instable node created seen in the view V;Vs and then converge towards a stable state due to the stable node created""" 

# ╔═╡ f56a539c-72d5-4856-a5e9-73b20f13b0a3
begin
	V1_null_regen_vs0_l_(a,b) = Vnullcline1(b,a,p_regen_vs0_l)
	V2_null_regen_vs0_l_(a,b) = Vnullcline2(b,a,p_regen_vs0_l)
	Vs_null_regen_vs0_l_(a,b) = Vsnullcline(b,a,p_regen_vs0_l)
	Vus_null_regen_vs0_l_(a,b) = Vusnullcline(b,a,p_regen_vs0_l)
	
	inter_V_Vs_null1_regen_vs0_l = zeros(length(vus))
	inter_V_Vs_null2_regen_vs0_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_vs0_l[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_vs0_l)
		inter_V_Vs_null2_regen_vs0_l[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_vs0_l)
	end

	inter_V_Vus_null1_regen_vs0_l = zeros(size(vs))
	inter_V_Vus_null2_regen_vs0_l = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_regen_vs0_l[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_vs0_l)
		inter_V_Vus_null2_regen_vs0_l[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_vs0_l)
	end


	plot(vus,vs,V1_null_regen_vs0_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_vs0_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_vs0_l,inter_V_Vs_null2_regen_vs0_l],[inter_V_Vs_null1_regen_vs0_l,inter_V_Vs_null2_regen_vs0_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_vs0_l,inter_V_Vus_null2_regen_vs0_l],[vs,vs],[inter_V_Vus_null1_regen_vs0_l,inter_V_Vus_null2_regen_vs0_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_vs0_l[3,1:end],sol_regen_vs0_l[2,1:end],sol_regen_vs0_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when vs0=$(p_regen_vs0_l[3])")

end

# ╔═╡ 2a17b374-9b13-4c37-ac06-702835a4087a
begin
	plot(plt_v_null_regen)
	title!("V nullcline when vs0=$(p_regen[3])")
end

# ╔═╡ c7a96ad0-f6fa-49c2-8054-6ce8c08d5529
begin
	V1_null_regen_vs0_h_(a,b) = Vnullcline1(b,a,p_regen_vs0_h)
	V2_null_regen_vs0_h_(a,b) = Vnullcline2(b,a,p_regen_vs0_h)
	Vs_null_regen_vs0_h_(a,b) = Vsnullcline(b,a,p_regen_vs0_h)
	Vus_null_regen_vs0_h_(a,b) = Vusnullcline(b,a,p_regen_vs0_h)
	
	inter_V_Vs_null1_regen_vs0_h = zeros(length(vus))
	inter_V_Vs_null2_regen_vs0_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_vs0_h[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_vs0_h)
		inter_V_Vs_null2_regen_vs0_h[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_vs0_h)
	end

	inter_V_Vus_null1_regen_vs0_h = zeros(size(vs))
	inter_V_Vus_null2_regen_vs0_h = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_regen_vs0_h[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_vs0_h)
		inter_V_Vus_null2_regen_vs0_h[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_vs0_h)
	end


	plot(vus,vs,V1_null_regen_vs0_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_vs0_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_vs0_h,inter_V_Vs_null2_regen_vs0_h],[inter_V_Vs_null1_regen_vs0_h,inter_V_Vs_null2_regen_vs0_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_vs0_h,inter_V_Vus_null2_regen_vs0_h],[vs,vs],[inter_V_Vus_null1_regen_vs0_h,inter_V_Vus_null2_regen_vs0_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_vs0_h[3,1:end],sol_regen_vs0_h[2,1:end],sol_regen_vs0_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when vs0=$(p_regen_vs0_h[3])")

end

# ╔═╡ 0256fbec-2fa6-4c23-bbbf-e7da53dbc52f
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_regen_vs0_h,inter_V_Vus_null2_regen_vs0_h],[inter_V_Vus_null1_regen_vs0_h,inter_V_Vus_null2_regen_vs0_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; vs0=$(p_regen_vs0_h[3])",legend=:outertopright) #pink
	plot!([inter_V_Vus_null1_regen_vs0_l,inter_V_Vus_null2_regen_vs0_l],[inter_V_Vus_null1_regen_vs0_l,inter_V_Vus_null2_regen_vs0_l],[vs,vs],linecolor=RGB(0.0,0.5,1),linewidth=3,label="dV/dt =dVus/dt = 0; vs0=$(p_regen_vs0_l[3])",legend=:outertopright)
	title!("Regenerative slow feedback")
end

# ╔═╡ 2dc59388-02bc-492c-92c0-cfd055e65eec
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_regen_vs0_h,inter_V_Vs_null2_regen_vs0_h],[inter_V_Vs_null1_regen_vs0_h,inter_V_Vs_null2_regen_vs0_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; vs0=$(p_regen_vs0_h[3])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_regen_vs0_l,inter_V_Vs_null2_regen_vs0_l],[inter_V_Vs_null1_regen_vs0_l,inter_V_Vs_null2_regen_vs0_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; vs0=$(p_regen_vs0_l[3])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ 1570f53a-efc9-4c0e-83b0-ffb172a69f9a
md"#### Summary"

# ╔═╡ 81c2ab95-5fb4-4cf2-bf9f-e209643a4fcf
md"###### Change in I "

# ╔═╡ 591a8483-5065-4113-9bfe-b1fb48e2b22d
md"""The overall effect of I in the 3D phase plane is that the hole size is modifified. 
If I increases, then the hole size increases. 
The contrary effect is true. The projections of the nullclines intersection in V;Vs and V;Vus show that I both change the extrama posiions on only one axis at the time depending on the sign of I.

A change in I will move the system around a saddle-node bifurcation. For a low value of I, the hole is too small and the system cannot burst. 
When I increases, the limit cycle frequency increases and the plateau disappear. """


# ╔═╡ 1210ff7f-73ee-4b52-9be8-6515ffe0d2c0
md"###### Change in gs "

# ╔═╡ fd609ef0-66e2-4599-971c-f06307185f24
md""" For a restorative ultraslow feedback, the convergence of the system is not the same
if we change gs. If gs is decreased, the system converges towards a limit cycle without 
plateau potentials and burst between. If gs is increases, the system is slowed down and converges towards a stable state.
When we increase gs, the nullclines in the plane (V;Vs) are less sharp and the extrema are closer. 
Therefore, the gradient is reduced in between (which slow down the system) and the nullclines might intersect, leading to a resting state.
When gs is decreased, the gradient is higher in the bottleneck whre the stable node takes place. 
Indeed, in the plane V;Vus, we see that the distance between the extrema increases as gs decreases. 

For a regenerative ultraslow feedback, only the speed at which 
the system converges towards the limit cycle seems to be impacted by a change in gs.
If gs decreases, the system converge fastly towards the stable limit cycle. 
This may be due to the fact that the hole width increases so the gradient of Vs is higher
as the trajectory is on the v nullcline. 

In the 3D phase plane, the effect of gs is on the hole size and the nullcline slope.
Indeed, the width of the hole is inversely proportional to gs.
This is consistent with the fact that the V nullcline plane tends to contract on Vs=Vs0* when gs increases.
"""

# ╔═╡ d5f310fb-5b1d-41a2-a897-aa4f0234899d
md"###### Change in gus "

# ╔═╡ 32e320d3-e31b-45ac-8760-ac938374fcd2
md"""The overall effect of modifying gus is that it changes the range of value for the excitability hole. Gus is inversely porpotional to its size. If it decreases, vus will tend to increase more which accelerates the ocsillations. On the contrary, if it increases too much, the hole becomes too small and the system have a stable node

Consequently, an increase in gus slow down the system (and make it to converge towards a resting state) while a decrease in gus accerelates the system and increases the equilibrium limit cycle frequency as the corresponding vus is higher, making |dVus/dt| bigger. For a regenerative US feedback, the effect is not that straightforward. It seems that frequency increases until a given value and decreases after. 

In the phase plane projections in V;Vus, we see that the the nullclines dVs/dt=dV/dt=0 contract on the symmetry axis V=V0* and the distance between its extrema increases as gus decreases 

In the phase plane projections in V;Vs, the effect is less clear but it seems that if the hole size is too low, the nullcine dVus/dt=dV/dt=0 is similar to what we had when the current is negative in the 2D model, making the system not excitable. """

# ╔═╡ eec55970-8d99-4a94-8968-7a20c5467f10
md"""		It=0 <-> Vus-Vus0 = +- sqrt(I/gus)"""

# ╔═╡ f503b8d5-a973-44ee-a4e7-b15935b88b2b
md"###### Change in Vs0 "

# ╔═╡ 3bcfea38-7819-41bf-891d-571f1dfc224f
md"""Changing vs0 moves the whole V nullcline plane to the right if vs0 increases and to the left otherwise on the vs axis. As the other planes where dVs/dt =0 or dVus/dt=0 do not move, the projection of the nullclines inersections in the plane V;Vs and V;Vus move. Consistently, the projection of the nullclines dVus/dt=0=dV/dt in the plane (V;Vs) move vertically according to vs0

Having a slow feedback too regenerative seems to slow the system as it can make the 3 nullclines to intersect. If the slow feedback becomes more and more resorative, the system accelerates due to the instable node created seen in the view V;Vs and then converge towards a stable state due to the stable node created. If vs0 is further decreased, then this point becomes stable and the system evolves toxards a rest state.""" 

# ╔═╡ dcb0fca8-a9f0-4d3b-ac61-19261afadb6a


# ╔═╡ 77b66127-095b-4191-96c0-f07235741fd7
@bind vs0_level html"<input type=range min=-40 max=-37 step=0.2>"

# ╔═╡ cff1b1ff-aa54-43ad-ae24-4e43f7e1837d
begin
	p_vs0_level = change_vs0(vs0_level,p)
	prob_vs0_level = ODEProblem(MQIF_3D!,u0,tspan,p_vs0_level,callback=cb)
	sol_vs0_level = solve(prob_vs0_level,dense=false)

	sub_vs0_lev = plot(sol_vs0_level.t,sol_vs0_level[1,:],linecolor="blue",label="V(t)")
	plot!(sol_vs0_level.t,sol_vs0_level[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_vs0_level.t,sol_vs0_level[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_vs0_level.t,p[3]ones(size(sol_vs0_level.t)),linecolor="red",label="Vs0")
	plot!(sol_vs0_level.t,p[4]ones(size(sol_vs0_level.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("Vs0=$(p_vs0_level[3])")
	
	
	sub_vs0 = plot(plt_time)
	title!("Vs0=$(p[3])")
	
	plot(sub_vs0_lev,sub_vs0,layout=(2,1))
end

# ╔═╡ 1dd6b220-b4be-49b7-930c-9f6f6a096dad
plotly()

# ╔═╡ a06e2104-bf28-4dac-aead-2a36d76eb36e
sub_vs0_lev

# ╔═╡ b8044619-7fcf-4e7a-936e-3a84932633c2
@bind vs0_level_r html"<input type=range min=-43 max=-37 step=0.5>"

# ╔═╡ ed25e322-b669-4a5e-9e9b-66a919368c36
begin
	p_vs0_level_r = change_vs0(vs0_level_r,p_regen)
	prob_vs0_level_r = ODEProblem(MQIF_3D!,u0,tspan,p_vs0_level_r,callback=cb)
	sol_vs0_level_r = solve(prob_vs0_level_r,dense=false)

	sub_vs0_lev_r = plot(sol_vs0_level_r.t,sol_vs0_level_r[1,:],linecolor="blue",label="V(t)")
	plot!(sol_vs0_level_r.t,sol_vs0_level_r[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_vs0_level_r.t,sol_vs0_level_r[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_vs0_level_r.t,p[3]ones(size(sol_vs0_level_r.t)),linecolor="red",label="Vs0")
	plot!(sol_vs0_level_r.t,p[4]ones(size(sol_vs0_level_r.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("Vs0=$(p_vs0_level_r[3])")
	
	
	sub_vs0_r = plot(plt_time_regen)
	title!("Vs0=$(p_regen[3])")
	
	plot(sub_vs0_lev_r,sub_vs0_r,layout=(2,1))
end

# ╔═╡ 7bd8401f-6cd3-407c-9621-dd69f3576764
p_regen

# ╔═╡ c9a43de3-0bb7-4e28-a7de-431b9953e08d
@bind gus_level html"<input type=range min=0.005 max=0.025 step=0.0025>"

# ╔═╡ 51c38624-17e6-4a67-afbf-af7e1f558e9d
begin
	p_gus_level = change_gus(gus_level,p)
	prob_gus_level = ODEProblem(MQIF_3D!,u0,tspan,p_gus_level,callback=cb)
	sol_gus_level = solve(prob_gus_level,dense=false)

	sub_gus_lev = plot(sol_gus_level.t,sol_gus_level[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gus_level.t,sol_gus_level[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gus_level.t,sol_gus_level[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gus_level.t,p[3]ones(size(sol_gus_level.t)),linecolor="red",label="Vs0")
	plot!(sol_gus_level.t,p[4]ones(size(sol_gus_level.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_gus_level[8])")
	
	
	sub_gus = plot(plt_time)
	title!("gus=$(p[8])")
	
	plot(sub_gus_lev,sub_gus,layout=(2,1))
end

# ╔═╡ aaf73312-3425-4c75-a00a-9cd3bd62ca79
@bind gus_level_r html"<input type=range min=0.005 max=0.035 step=0.0025>"

# ╔═╡ 927d6784-ac9e-447c-a27a-aa2fd258d749
begin
	p_gus_level_r = change_gus(gus_level_r,p_regen)
	prob_gus_level_r = ODEProblem(MQIF_3D!,u0,tspan,p_gus_level_r,callback=cb)
	sol_gus_level_r = solve(prob_gus_level_r,dense=false)

	sub_gus_lev_r = plot(sol_gus_level_r.t,sol_gus_level_r[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gus_level_r.t,sol_gus_level_r[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gus_level_r.t,sol_gus_level_r[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gus_level_r.t,p[3]ones(size(sol_gus_level_r.t)),linecolor="red",label="Vs0")
	plot!(sol_gus_level_r.t,p[4]ones(size(sol_gus_level_r.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gus=$(p_gus_level_r[8])")
	
	
	sub_gus_r = plot(plt_time_regen)
	title!("gus=$(p_regen[8])")
	
	plot(sub_gus_lev_r,sub_gus_r,layout=(2,1))
end

# ╔═╡ Cell order:
# ╟─651f4fe0-a054-11eb-102f-5596e5a831f7
# ╠═a389152e-ca91-4e87-b032-416b45376de5
# ╟─05814d8c-92a6-4a28-865c-d5b55d8a10d6
# ╠═0fffa3a2-76a1-4470-a5d2-4b78789cf6df
# ╠═98821bce-ccb4-403c-8af9-23021322e059
# ╠═6496c76f-f3e3-4b2b-a1b9-6edf9d6d330c
# ╠═21f17a1f-ba58-42ea-b3da-19e7f83ea34f
# ╠═a1fc8af4-2cd4-457e-a25c-ebcd78fe5636
# ╠═9355dcb7-c3a1-4787-a026-5858ef83c52d
# ╠═1884ccd5-b52e-439b-b0aa-fffad18c4cac
# ╠═eb855792-b818-46b3-a02e-f4b97ccebfc4
# ╠═e1726250-f07d-45b6-b82b-74f543f01215
# ╠═abe05fb9-a447-41ce-8981-7e49dc35acb9
# ╠═089d1752-141c-46cb-8b64-23f0ad02c336
# ╠═9d452ac9-c876-494a-bd55-92ba30ae71b6
# ╠═2595be40-b236-46ef-9e88-f2c4225b3197
# ╠═96f8c76b-1156-4b5c-b0c4-fd8e820662c6
# ╠═8d4e1961-ab4c-40fb-97a6-784c4a6057af
# ╠═841be4cd-afd5-4fa3-a560-3607049bb23f
# ╠═8c75d7e7-db77-43a4-b221-3adea64d2ef7
# ╠═e72b85df-5951-4275-81cc-1f068d8946a0
# ╟─acef61f0-7ad2-4eca-8049-537b1f283b61
# ╟─b55387e2-a88e-48d1-8c5e-fc4e0dddc14e
# ╟─352523e4-026b-40b8-a850-8610e2646a67
# ╟─8bf784c8-3b88-46b4-8b62-dbb1a581402b
# ╟─ed1d0cdc-a5bc-48f5-b81c-d9e9d12bdf36
# ╟─7e8b53bc-4800-46ba-bcff-55aca599d1cb
# ╟─3cf9b4c3-bfe1-489e-a0b3-a0713114dab6
# ╟─a0726536-711c-4d21-a2c4-26767e00eeb7
# ╟─0cad5c79-4f97-4331-8781-e9544a08523e
# ╟─80235e6a-ef03-4068-84ec-02c618d60e4b
# ╟─86f8d8f0-2f6b-4ec1-900f-1d4e1f743324
# ╟─ae5c617e-aad0-4321-b76e-5ff38b554646
# ╟─9ee2c719-b34a-41be-b951-212303f2e58d
# ╟─0da094a5-2899-4dc0-a03d-be4997b1e507
# ╟─4af9fe12-05fa-4dfa-84af-9fe5b739410c
# ╟─6e504105-ac3a-419d-9ee1-89f15ae2445c
# ╟─86755dcd-d1a5-423d-ac23-8d2615afddda
# ╟─4218be3f-9b82-476f-b8c1-21a90dad00bb
# ╟─b34deb54-5b25-4050-98ba-8147a1f27d8c
# ╟─bead8814-21a5-4a7a-8488-d5493e14b391
# ╟─53e4acac-b36e-48a0-b32a-4ff6bb68f67b
# ╟─8d01dc21-6666-42eb-9b7f-ac05be91b1d8
# ╟─c43b5adc-0253-4b9d-835b-3ec44362ed03
# ╟─6d652924-4231-4dc4-81d1-f6a438aea675
# ╟─24065204-66b6-4dbb-9656-0c252b014cce
# ╠═f2107a48-c7ee-46db-9acd-44e16488ca0a
# ╟─8035aaab-6b4f-49ef-874a-06989d53cc41
# ╟─b342899c-0cd3-4b35-b261-0807216de37e
# ╟─29bb5f71-a408-479e-83ed-69a62e6bba29
# ╠═01252cf7-d74c-4f02-977a-fcbf743dc056
# ╠═82eb0286-c992-4459-9d65-529e8ba30362
# ╟─2eb57d95-e17c-4f71-8f9f-34a6df5b0392
# ╟─02927149-ba9f-4da7-83b1-390a34d6f62f
# ╟─17518825-5527-4478-a4c3-52f7ff864a87
# ╟─70af372c-963d-4e53-b162-f3b2ca733b48
# ╟─1cfe6b28-6d83-4bfa-920a-b4eb387590cb
# ╟─cfa819c9-a58a-44b3-bffd-3405fabaf6cc
# ╠═49b1f686-bcf7-49ed-93a3-aa74d0081110
# ╟─6ec8f2d4-8de5-480b-b8e4-1a3b49fa4a4a
# ╟─f1aba699-0d2c-4835-ab43-42bea662a141
# ╟─b990cc17-1ea4-468a-8a2e-39f897cade99
# ╟─c67e6b50-ee12-4361-ac57-1e602641fad6
# ╟─6a379d9f-05c3-4a3c-935b-b2b8bc96d43f
# ╠═7371f342-f293-497d-b5c4-bbd24f02967c
# ╟─b7a260fd-9f58-4fd8-a5c2-45997bf62bc8
# ╟─0365925b-e208-4f13-8908-747f238a2a4f
# ╟─e05e543d-a6ba-4a26-a61c-9c05daba95c5
# ╟─c7875827-d1a8-4058-be16-e4a1437cbc6a
# ╟─9d900972-fc6d-49b0-b1e3-3a980deb2b97
# ╠═ee174df3-a596-441e-a3d0-808013360e31
# ╟─298c4acd-c3a5-4f88-89e6-8d5b5fb249b0
# ╟─509fa73f-5b82-42ab-a417-82ab648c1eed
# ╟─a7648315-fbe4-4bd2-a3e2-0ad038fefbbe
# ╟─1768e5d7-82de-46cf-8206-a61bb392c781
# ╟─5f82aa48-4287-41bb-914e-ad9f62a2a270
# ╟─dadccbbb-a331-424f-9c97-a7753b8cf92a
# ╟─3dc77dcb-b68f-4006-bc0d-9dd39cbf2f38
# ╟─b5d1876f-c572-4e77-b694-f8270162f03f
# ╟─ee04cc5e-9fdc-4468-a589-7645343b7682
# ╟─672a0a95-82e3-4d42-ab3f-f2cf26510652
# ╟─0ed3d289-8c77-47e5-8972-a5fceb9a5fe2
# ╟─5e07a6ae-a8f5-480e-8d54-369b94282cee
# ╟─86e7ab9d-4544-400d-b7b0-6c4fea65e264
# ╟─879293d4-1d92-44da-923b-3e2cc9e5bc60
# ╟─a0bfcfa1-5437-468a-ae71-d5d36c7e21b7
# ╟─403bbf12-9c0d-459e-89d0-54b9bf6a524d
# ╟─e15ecbf6-59ad-4675-87ea-f9803d5883f1
# ╟─e42297bc-471f-4abf-9939-28a31f67f5a2
# ╟─f241d2d3-38a9-441d-916b-51289a4a1e61
# ╟─d5c6139f-7125-468f-a21b-576970c7b478
# ╟─1f77716b-61f6-48ca-ae45-1616771cd5db
# ╟─55bcf2d5-0b82-45ce-9689-fe0165513784
# ╟─002db113-6b87-47ec-9a3d-009a6143ec52
# ╟─23bbe085-afec-42e4-99ee-b863ce68d62f
# ╟─cb445d6e-ac17-4137-985b-f73ac538dfbf
# ╟─e75132c6-17ef-4bb2-ab33-9550bc88091e
# ╟─04500988-7e91-4a99-9262-7e46d02a732e
# ╟─13ed704b-4ed1-47d0-a384-c66a8de1d0a5
# ╟─0da78619-4c28-4c41-a5ef-2fd0fc349e35
# ╟─129badd2-31e9-4dbd-aa02-f6b0df5b10d6
# ╟─0df08dc3-3dc9-4255-895e-b20f21243bf5
# ╟─8bde95e1-f76c-4943-8d1e-5b43ac2defdf
# ╟─9a6033ce-60df-40a3-a823-49ff3b3c5261
# ╟─97bac848-795e-4770-8ea0-3cb021823bf1
# ╟─f380fa87-c7a9-4fac-9be0-7987d8ca8960
# ╟─88c7243f-d325-459d-9297-abc999bbb269
# ╟─6a48b3b3-6204-4127-b431-c7125321adcf
# ╟─5a474cef-3084-45ed-8b83-bd65adfbfe6b
# ╟─1090cbfa-667c-474d-b00f-e2db4e37cd90
# ╟─0a89f507-e3ae-4d4d-8eff-f6bbdb822a92
# ╟─32449dbb-0486-458e-812a-bf343f23bfc4
# ╟─f56a539c-72d5-4856-a5e9-73b20f13b0a3
# ╟─2a17b374-9b13-4c37-ac06-702835a4087a
# ╟─c7a96ad0-f6fa-49c2-8054-6ce8c08d5529
# ╟─0256fbec-2fa6-4c23-bbbf-e7da53dbc52f
# ╟─2dc59388-02bc-492c-92c0-cfd055e65eec
# ╟─1570f53a-efc9-4c0e-83b0-ffb172a69f9a
# ╟─81c2ab95-5fb4-4cf2-bf9f-e209643a4fcf
# ╟─591a8483-5065-4113-9bfe-b1fb48e2b22d
# ╟─1210ff7f-73ee-4b52-9be8-6515ffe0d2c0
# ╟─fd609ef0-66e2-4599-971c-f06307185f24
# ╟─d5f310fb-5b1d-41a2-a897-aa4f0234899d
# ╟─32e320d3-e31b-45ac-8760-ac938374fcd2
# ╟─eec55970-8d99-4a94-8968-7a20c5467f10
# ╟─f503b8d5-a973-44ee-a4e7-b15935b88b2b
# ╟─3bcfea38-7819-41bf-891d-571f1dfc224f
# ╟─dcb0fca8-a9f0-4d3b-ac61-19261afadb6a
# ╠═77b66127-095b-4191-96c0-f07235741fd7
# ╟─cff1b1ff-aa54-43ad-ae24-4e43f7e1837d
# ╠═1dd6b220-b4be-49b7-930c-9f6f6a096dad
# ╠═a06e2104-bf28-4dac-aead-2a36d76eb36e
# ╟─b8044619-7fcf-4e7a-936e-3a84932633c2
# ╟─ed25e322-b669-4a5e-9e9b-66a919368c36
# ╠═7bd8401f-6cd3-407c-9621-dd69f3576764
# ╟─c9a43de3-0bb7-4e28-a7de-431b9953e08d
# ╠═51c38624-17e6-4a67-afbf-af7e1f558e9d
# ╟─aaf73312-3425-4c75-a00a-9cd3bd62ca79
# ╟─927d6784-ac9e-447c-a27a-aa2fd258d749
