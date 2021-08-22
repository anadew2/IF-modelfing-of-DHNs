### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

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
p_I_l = change_I(-1,p)

# ╔═╡ 6496c76f-f3e3-4b2b-a1b9-6edf9d6d330c
p_I_h = change_I(10,p)

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
p_gus_h = change_gus(0.15,p)

# ╔═╡ 1884ccd5-b52e-439b-b0aa-fffad18c4cac
p_gus_l = change_gus(0.0015,p)

# ╔═╡ 8d4e1961-ab4c-40fb-97a6-784c4a6057af
p_regen_gus_h = change_gus(0.15,p)

# ╔═╡ 841be4cd-afd5-4fa3-a560-3607049bb23f
p_regen_gus_l = change_gus(0.0015,p)

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
p_vs0_l = change_vs0(-42,p)

# ╔═╡ e1726250-f07d-45b6-b82b-74f543f01215
p_vs0_h = change_vs0(-36,p)

# ╔═╡ 8c75d7e7-db77-43a4-b221-3adea64d2ef7
p_regen_vs0_l = change_vs0(-42,p_regen)

# ╔═╡ e72b85df-5951-4275-81cc-1f068d8946a0
p_regen_vs0_h = change_vs0(-28,p_regen)

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

# ╔═╡ 972118cb-6c8b-4618-9e66-dfdae5b343f6
begin
	prob_I_l = ODEProblem(MQIF_3D!,u0,tspan,p_I_l,callback=cb)
	sol_I_l = solve(prob_I_l,dense=false)
end

# ╔═╡ aefc61d1-50b0-4f86-9204-2935495de66d
begin
	prob_I_h = ODEProblem(MQIF_3D!,u0,tspan,p_I_h,callback=cb)
	sol_I_h = solve(prob_I_h,dense=false)
end

# ╔═╡ 82eb0286-c992-4459-9d65-529e8ba30362
gr()

# ╔═╡ 89fb57df-ffb6-4a97-acb5-d2c293b3f1d3
begin
	plot(sol_I_l.t,sol_I_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_I_l.t,sol_I_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_I_l.t,sol_I_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_I_l.t,p[3]ones(size(sol_I_l.t)),linecolor="red",label="Vs0")
	plot!(sol_I_l.t,p[4]ones(size(sol_I_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("I=$(p_I_l[1])")
end

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

# ╔═╡ 8698b02a-a34e-4e4f-9ecf-a709e57136a2
begin
	plot(sol_I_h.t,sol_I_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_I_h.t,sol_I_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_I_h.t,sol_I_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_I_h.t,p[3]ones(size(sol_I_h.t)),linecolor="red",label="Vs0")
	plot!(sol_I_h.t,p[4]ones(size(sol_I_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("I=$(p_I_h[1])")
end

# ╔═╡ 02927149-ba9f-4da7-83b1-390a34d6f62f
md"""Due to the reset, vus is sometimes brought a little bit too high making the current more negative. 
Thus, the trajectory is sucked on the left part of the v nullcline for It<0 making v to decrease a little more. 
To contrast with the previous case, the gradient of Vus is higher here due to the high I. 
Therefore, the current I_t varies faster. 
The trajectory continues to evolve slowly with v increases but stay lower than V0 (left part of the phase plane). 
At some point, due to the highly varying current, the trajectory is repelled by the saddle and the system spikes.  """

# ╔═╡ 17518825-5527-4478-a4c3-52f7ff864a87
md"""Changing I moves the system on a SN bifurcation"""

# ╔═╡ 371dc1f8-570b-4d77-a235-03b14efca425


# ╔═╡ 7ee5fec7-9299-4fb4-bc9a-8af9e94fec34
begin
	V1_null_I_l_(a,b) = Vnullcline1(b,a,p_I_l)
	V2_null_I_l_(a,b) = Vnullcline2(b,a,p_I_l)
	Vs_null_I_l_(a,b) = Vsnullcline(b,a,p_I_l)
	Vus_null_I_l_(a,b) = Vusnullcline(b,a,p_I_l)
	
	inter_V_Vs_null1_I_l = zeros(length(vus))
	inter_V_Vs_null2_I_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_I_l[i] = V_Vs_nullcline1(NaN,vus[i],p_I_l)
		inter_V_Vs_null2_I_l[i] = V_Vs_nullcline2(NaN,vus[i],p_I_l)
	end

	inter_V_Vus_null1_I_l = zeros(size(vs))
	inter_V_Vus_null2_I_l = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_I_l[i] = V_Vus_nullcline1(vs[i],NaN,p_I_l)
		inter_V_Vus_null2_I_l[i] = V_Vus_nullcline2(vs[i],NaN,p_I_l)
	end


	plot(vus,vs,V1_null_I_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_I_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_I_l,inter_V_Vs_null2_I_l],[inter_V_Vs_null1_I_l,inter_V_Vs_null2_I_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_I_l,inter_V_Vus_null2_I_l],[vs,vs],[inter_V_Vus_null1_I_l,inter_V_Vus_null2_I_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	#plot!(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p_I_l[1])")

end

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

# ╔═╡ 2a7ef951-82a6-4433-8911-8936704b9424
begin
	V1_null_I_h_(a,b) = Vnullcline1(b,a,p_I_h)
	V2_null_I_h_(a,b) = Vnullcline2(b,a,p_I_h)
	Vs_null_I_h_(a,b) = Vsnullcline(b,a,p_I_h)
	Vus_null_I_h_(a,b) = Vusnullcline(b,a,p_I_h)
	
	inter_V_Vs_null1_I_h = zeros(length(vus))
	inter_V_Vs_null2_I_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_I_h[i] = V_Vs_nullcline1(NaN,vus[i],p_I_h)
		inter_V_Vs_null2_I_h[i] = V_Vs_nullcline2(NaN,vus[i],p_I_h)
	end

	inter_V_Vus_null1_I_h = zeros(size(vs))
	inter_V_Vus_null2_I_h = zeros(size(vs))
	for i=1:length(vs)
		inter_V_Vus_null1_I_h[i] = V_Vus_nullcline1(vs[i],NaN,p_I_h)
		inter_V_Vus_null2_I_h[i] = V_Vus_nullcline2(vs[i],NaN,p_I_h)
	end


	plot(vus,vs,V1_null_I_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_I_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],[inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],[vs,vs],[inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_I_h[3,1:end],sol_I_h[2,1:end],sol_I_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p_I_h[1])")

end

# ╔═╡ cfa819c9-a58a-44b3-bffd-3405fabaf6cc
md"""Changing I increases the size of the hole where the system is excitable and therefore moves a lot the intersection with the Vs nullcline"""

# ╔═╡ 49b1f686-bcf7-49ed-93a3-aa74d0081110
gr()

# ╔═╡ f494a04b-0a96-411e-bc22-1c7e8da6341c
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],[inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; I=$(p_I_h[1])",legend=:outertopright) #pink
	title!("Regenerative slow feedback")
end

# ╔═╡ 871c279c-7579-458b-b110-286389f710ad
gr()

# ╔═╡ 18ac0378-77c7-4570-a0d4-f687ec76d420
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],[inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; I=$(p_I_h[1])",legend=:outertopright) #pink
	title!("Restorative ultraslow feedback")
end

# ╔═╡ e7d5ccae-c62d-4aae-910d-0cb3c172d7b8
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0",legend=:outertopright) #pink
	xaxis!((-80,vs_max),"V")
	yaxis!("Vus")
end

# ╔═╡ 6ec8f2d4-8de5-480b-b8e4-1a3b49fa4a4a
md"###### Current variation for regenerative ultraslow feedback and regenerative slow feedback"

# ╔═╡ eedb2f8e-f389-4a80-9d7f-6dd3b5833110
begin
	prob_regen_I_l = ODEProblem(MQIF_3D!,u0,tspan,p_regen_I_l,callback=cb)
	sol_regen_I_l = solve(prob_regen_I_l,dense=false)
end

# ╔═╡ f1aba699-0d2c-4835-ab43-42bea662a141
begin
	prob_regen = ODEProblem(MQIF_3D!,u0,tspan,p_regen,callback=cb)
	sol_regen = solve(prob_regen,dense=false)
end

# ╔═╡ 243fe1ba-9b96-4028-a8cb-4928153a7f19
begin
	prob_regen_I_h = ODEProblem(MQIF_3D!,u0,tspan,p_regen_I_h,callback=cb)
	sol_regen_I_h = solve(prob_regen_I_h,dense=false)
end

# ╔═╡ 8ea57601-f46f-4488-8dd5-ec0189c5b751
begin
	plot(sol_regen_I_l.t,sol_regen_I_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_I_l.t,sol_regen_I_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_I_l.t,sol_regen_I_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_I_l.t,p_regen[3]ones(size(sol_regen_I_l.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_I_l.t,p_regen[4]ones(size(sol_regen_I_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("I=$(p_regen_I_l[1])")
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

# ╔═╡ 410eee64-5670-4b57-be26-d45b9b4e5899
begin
	plot(sol_regen_I_h.t,sol_regen_I_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_I_h.t,sol_regen_I_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_I_h.t,sol_regen_I_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_I_h.t,p_regen[3]ones(size(sol_regen_I_h.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_I_h.t,p_regen[4]ones(size(sol_regen_I_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("I=$(p_regen_I_h[1])")
end

# ╔═╡ 5464ca79-59ea-4b55-adf9-8c17ec5afa7b
begin
	plot(sol_regen[1,Int(round(end/2)):end],sol_regen[2,Int(round(end/2)):end],label="I=$(p_regen[1])")
	plot!(sol_regen_I_h[1,Int(round(end/2)):end],sol_regen_I_h[2,Int(round(end/2)):end],label="I=$(p_regen_I_h[1])")
	plot!(size=(200,200))
	
	xaxis!("V")
	yaxis!("Vs")
	title!("Limit cycle")
end

# ╔═╡ c67e6b50-ee12-4361-ac57-1e602641fad6
md"""When I increases, the limit cycle frequency increases"""

# ╔═╡ 7011c21d-da15-4b7e-b9ac-f8d06bc08511


# ╔═╡ 7fe6cfa0-1f7d-498d-a8db-2afc98fdd838
begin
	V1_null_regen_I_l_(a,b) = Vnullcline1(b,a,p_regen_I_l)
	V2_null_regen_I_l_(a,b) = Vnullcline2(b,a,p_regen_I_l)
	Vs_null_regen_I_l_(a,b) = Vsnullcline(b,a,p_regen_I_l)
	Vus_null_regen_I_l_(a,b) = Vusnullcline(b,a,p_regen_I_l)
	
	inter_V_Vs_null1_regen_I_l = zeros(length(vus))
	inter_V_Vs_null2_regen_I_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_I_l[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_I_l)
		inter_V_Vs_null2_regen_I_l[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_I_l)
	end

	inter_V_Vus_null1_regen_I_l = zeros(size(vs))
	inter_V_Vus_null2_regen_I_l = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_regen_I_l[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_I_l)
		inter_V_Vus_null2_regen_I_l[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_I_l)
	end


	plot(vus,vs,V1_null_regen_I_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_I_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_I_l,inter_V_Vs_null2_regen_I_l],[inter_V_Vs_null1_regen_I_l,inter_V_Vs_null2_regen_I_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_I_l,inter_V_Vus_null2_regen_I_l],[vs,vs],[inter_V_Vus_null1_regen_I_l,inter_V_Vus_null2_regen_I_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	#plot!(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p_regen_I_l[1])")

end

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

# ╔═╡ 91bb6de8-df2f-4598-a231-1f855591c65e
begin
	V1_null_regen_I_h_(a,b) = Vnullcline1(b,a,p_regen_I_h)
	V2_null_regen_I_h_(a,b) = Vnullcline2(b,a,p_regen_I_h)
	Vs_null_regen_I_h_(a,b) = Vsnullcline(b,a,p_regen_I_h)
	Vus_null_regen_I_h_(a,b) = Vusnullcline(b,a,p_regen_I_h)
	
	inter_V_Vs_null1_regen_I_h = zeros(length(vus))
	inter_V_Vs_null2_regen_I_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_I_h[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_I_h)
		inter_V_Vs_null2_regen_I_h[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_I_h)
	end

	inter_V_Vus_null1_regen_I_h = zeros(size(vs))
	inter_V_Vus_null2_regen_I_h = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_regen_I_h[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_I_h)
		inter_V_Vus_null2_regen_I_h[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_I_h)
	end


	plot(vus,vs,V1_null_regen_I_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_I_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_I_h,inter_V_Vs_null2_regen_I_h],[inter_V_Vs_null1_regen_I_h,inter_V_Vs_null2_regen_I_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_I_h,inter_V_Vus_null2_regen_I_h],[vs,vs],[inter_V_Vus_null1_regen_I_h,inter_V_Vus_null2_regen_I_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	#plot!(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when I=$(p_regen_I_h[1])")

end

# ╔═╡ 7371f342-f293-497d-b5c4-bbd24f02967c
gr()

# ╔═╡ 8f0158b1-9664-4c27-b778-912c480e3ae8
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_regen_I_h,inter_V_Vus_null2_regen_I_h],[inter_V_Vus_null1_regen_I_h,inter_V_Vus_null2_regen_I_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; I=$(p_I_h[1])",legend=:outertopright) #pink
	title!("Regenerative slow feedback")
end

# ╔═╡ b1c9b459-7d75-4082-9876-a3ab480714c7
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_regen_I_h,inter_V_Vs_null2_regen_I_h],[inter_V_Vs_null1_regen_I_h,inter_V_Vs_null2_regen_I_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; I=$(p_I_h[1])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_regen_I_l,inter_V_Vs_null2_regen_I_l],[inter_V_Vs_null1_regen_I_l,inter_V_Vs_null2_regen_I_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; I=$(p_I_l[1])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ b7a260fd-9f58-4fd8-a5c2-45997bf62bc8
md"""Changing I increases the size of the hole where the system is excitable and therefore moves a lot the intersection with the Vs nullcline"""

# ╔═╡ 0365925b-e208-4f13-8908-747f238a2a4f
md"###### GS variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ 59d2716d-3001-4499-b213-f3b3848bdb0f
begin
	prob_gs_l = ODEProblem(MQIF_3D!,u0,tspan,p_gs_l,callback=cb)
	sol_gs_l = solve(prob_gs_l,dense=false)
end

# ╔═╡ d0985722-90cc-4c25-8826-9a12a9eef8c2
begin
	prob_gs_h = ODEProblem(MQIF_3D!,u0,tspan,p_gs_h,callback=cb)
	sol_gs_h = solve(prob_gs_h,dense=false)
end

# ╔═╡ c69a5267-9b2f-4531-9b5f-27d8f56940e2
begin
	plot(sol_gs_l.t,sol_gs_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gs_l.t,sol_gs_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gs_l.t,sol_gs_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gs_l.t,p[3]ones(size(sol_gs_l.t)),linecolor="red",label="Vs0")
	plot!(sol_gs_l.t,p[4]ones(size(sol_gs_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gs=$(p_gs_l[7])")
end

# ╔═╡ 97492156-9228-42f4-bab0-3175948559d0
begin
	plot(plt_time)
	title!("gs=$(p[7])")
end

# ╔═╡ bd571ae8-74e1-480e-ac7f-be65d20ef89d
begin
	plot(sol_gs_h.t,sol_gs_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_gs_h.t,sol_gs_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_gs_h.t,sol_gs_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_gs_h.t,p[3]ones(size(sol_gs_h.t)),linecolor="red",label="Vs0")
	plot!(sol_gs_h.t,p[4]ones(size(sol_gs_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage ")
	title!("gs=$(p_gs_h[7])")
end

# ╔═╡ 7520f4b2-8948-474b-988b-23573b9311e6
gr()

# ╔═╡ eeb780aa-3943-4bc9-8108-c5c4dfe6b773


# ╔═╡ e156e461-1bac-427d-b6d1-e8ea5ce05433
begin
	V1_null_gs_l_(a,b) = Vnullcline1(b,a,p_gs_l)
	V2_null_gs_l_(a,b) = Vnullcline2(b,a,p_gs_l)
	Vs_null_gs_l_(a,b) = Vsnullcline(b,a,p_gs_l)
	Vus_null_gs_l_(a,b) = Vusnullcline(b,a,p_gs_l)
	
	inter_V_Vs_null1_gs_l = zeros(length(vus))
	inter_V_Vs_null2_gs_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_gs_l[i] = V_Vs_nullcline1(NaN,vus[i],p_gs_l)
		inter_V_Vs_null2_gs_l[i] = V_Vs_nullcline2(NaN,vus[i],p_gs_l)
	end

	inter_V_Vus_null1_gs_l = zeros(size(vs))
	inter_V_Vus_null2_gs_l = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_gs_l[i] = V_Vus_nullcline1(vs[i],NaN,p_gs_l)
		inter_V_Vus_null2_gs_l[i] = V_Vus_nullcline2(vs[i],NaN,p_gs_l)
	end


	plot(vus,vs,V1_null_gs_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_gs_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_gs_l,inter_V_Vs_null2_gs_l],[inter_V_Vs_null1_gs_l,inter_V_Vs_null2_gs_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_gs_l,inter_V_Vus_null2_gs_l],[vs,vs],[inter_V_Vus_null1_gs_l,inter_V_Vus_null2_gs_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_gs_l[3,1:end],sol_gs_l[2,1:end],sol_gs_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gs=$(p_gs_l[7])")

end

# ╔═╡ 4f1874ab-569b-4145-93c6-1688c1d3aadf
begin
	plot(plt_v_null)
	title!("V nullcline when gs=$(p[7])")
end

# ╔═╡ 3a35b2f1-ba94-4597-bb1b-fba90a72930b
begin
	V1_null_gs_h_(a,b) = Vnullcline1(b,a,p_gs_h)
	V2_null_gs_h_(a,b) = Vnullcline2(b,a,p_gs_h)
	Vs_null_gs_h_(a,b) = Vsnullcline(b,a,p_gs_h)
	Vus_null_gs_h_(a,b) = Vusnullcline(b,a,p_gs_h)
	
	inter_V_Vs_null1_gs_h = zeros(length(vus))
	inter_V_Vs_null2_gs_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_gs_h[i] = V_Vs_nullcline1(NaN,vus[i],p_gs_h)
		inter_V_Vs_null2_gs_h[i] = V_Vs_nullcline2(NaN,vus[i],p_gs_h)
	end

	inter_V_Vus_null1_gs_h = zeros(size(vs))
	inter_V_Vus_null2_gs_h = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_gs_h[i] = V_Vus_nullcline1(vs[i],NaN,p_gs_h)
		inter_V_Vus_null2_gs_h[i] = V_Vus_nullcline2(vs[i],NaN,p_gs_h)
	end


	plot(vus,vs,V1_null_gs_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_gs_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_gs_h,inter_V_Vs_null2_gs_h],[inter_V_Vs_null1_gs_h,inter_V_Vs_null2_gs_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_gs_h,inter_V_Vus_null2_gs_h],[vs,vs],[inter_V_Vus_null1_gs_h,inter_V_Vus_null2_gs_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	#plot!(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gs=$(p_gs_h[7])")

end

# ╔═╡ 772cc897-4364-4e80-97cc-65b064547884
plotly()

# ╔═╡ 48f59077-f18e-469a-877d-9f7a4c5c4d10
begin
	
	plot(vus,vs,V1_null_gs_l_,st=:surface,c=:Purples_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(-40,-30),ylims=(-45,-35))  
	#plot!(vus,vs,V2_null_gs_l_,st=:surface,c=:Purples_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(-40,-30),ylims=(-45,-35))  
	plot!(vus,vs,V1_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(-40,-30),ylims=(-45,-35))  
	#plot!(vus,vs,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(-40,-30),ylims=(-45,-35))  
	plot!(vus,vs,V1_null_gs_h_,st=:surface,c=:Oranges_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(-40,-30),ylims=(-45,-35),zlims=(-45,-35),camera=(-90,0))
	plot!(vus,vs,V2_null_gs_l_,st=:surface,c=:Oranges_3,colorbar_entry=false) 
end

# ╔═╡ e05e543d-a6ba-4a26-a61c-9c05daba95c5
md"""When gs increases, the V nullcline is brought closer to the plane Vs=Vs0

When gs decreases, it means that the V nullcline becomes more sensitive to V than to Vs or Vus which is seen 
as the fact that the nullcline contracts towards the plane V=V0"""

# ╔═╡ f4ae01bc-4092-414e-a70d-7abe0c38cd6c
gr()

# ╔═╡ 53d6b75b-41e1-44bb-a3ca-95863677b92b
begin
	plot([inter_V_Vs_null1,inter_V_Vs_null2],[inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1,inter_V_Vus_null2],[vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],[inter_V_Vs_null1_I_h,inter_V_Vs_null2_I_h],[vus,vus],linecolor=RGB(1,0.3,0),linewidth=3,label="dV/dt =dVs/dt = 0; I=$(p_I_h[1])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_gs_h,inter_V_Vs_null2_gs_h],[inter_V_Vs_null1_gs_h,inter_V_Vs_null2_gs_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; gs=$(p_gs_h[7])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_gs_l,inter_V_Vs_null2_gs_l],[inter_V_Vs_null1_gs_l,inter_V_Vs_null2_gs_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; gs=$(p_gs_l[7])",legend=:outertopright) #pink
	title!("Restorative ultraslow feedback")
end

# ╔═╡ 20b8ef9a-18e8-4aed-8a7e-ae927938fd7b
p_I_h

# ╔═╡ c7875827-d1a8-4058-be16-e4a1437cbc6a
md"###### GS variation for regenerative ultraslow feedback and regenerative slow feedback"

# ╔═╡ fe6c6934-d5f6-4d05-88dc-845b9759122d
begin
	prob_regen_gs_l = ODEProblem(MQIF_3D!,u0,tspan,p_regen_gs_l,callback=cb)
	sol_regen_gs_l = solve(prob_regen_gs_l,dense=false)
end

# ╔═╡ dd575c93-296e-445c-bc2b-b27ab3f89b45
begin
	prob_regen_gs_h = ODEProblem(MQIF_3D!,u0,tspan,p_regen_gs_h,callback=cb)
	sol_regen_gs_h = solve(prob_regen_gs_h,dense=false)
end

# ╔═╡ 3878b997-3ec2-42f7-9d71-195028344152
begin
	plot(sol_regen_gs_l.t,sol_regen_gs_l[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_gs_l.t,sol_regen_gs_l[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_gs_l.t,sol_regen_gs_l[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_gs_l.t,p_regen[3]ones(size(sol_regen_gs_l.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_gs_l.t,p_regen[4]ones(size(sol_regen_gs_l.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("gs=$(p_regen_gs_l[7])")
end

# ╔═╡ 628d70e4-4d6d-4d65-aa0c-f028ef6317d4
begin
	plot(plt_time_regen)
	title!("gs=$(p_regen[7])")
end

# ╔═╡ d80fdc73-c420-49a0-9e6d-c8cfebfd7d60
begin
	plot(sol_regen_gs_h.t,sol_regen_gs_h[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen_gs_h.t,sol_regen_gs_h[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen_gs_h.t,sol_regen_gs_h[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen_gs_h.t,p_regen[3]ones(size(sol_regen_gs_h.t)),linecolor="red",label="Vs0")
	plot!(sol_regen_gs_h.t,p_regen[4]ones(size(sol_regen_gs_h.t)),linecolor="magenta",label="Vus0")
	plot!(size=(300,300))
	
	xaxis!("Time (ms)")
	yaxis!("Voltage")
	title!("gs=$(p_regen_gs_h[7])")
end

# ╔═╡ 9d900972-fc6d-49b0-b1e3-3a980deb2b97
md"""gs seems to have an impact mostly on the burst frequency before converging towards the limit cycle"""

# ╔═╡ 081428d2-c3ca-421c-8f9e-e45d03b24407
begin
	plot(sol_regen_gs_h[1,Int(round(end/2)):end],sol_regen_gs_h[2,Int(round(end/2)):end],label="gs=$(p_regen_gs_h[7])")
	plot!(sol_regen[1,Int(round(end/2)):end],sol_regen[2,Int(round(end/2)):end],label="gs=$(p_regen[7])")
	plot!(sol_regen_gs_l[1,Int(round(end/2)):end],sol_regen_gs_l[2,Int(round(end/2)):end],label="gs=$(p_regen_gs_l[7])")
	plot!(size=(200,200))
	
	xaxis!("V")
	yaxis!("Vs")
	title!("Limit cycle")
end

# ╔═╡ 07a1ef09-5a86-4b05-9e8e-5355af702b1c
begin
	V1_null_regen_gs_l_(a,b) = Vnullcline1(b,a,p_regen_gs_l)
	V2_null_regen_gs_l_(a,b) = Vnullcline2(b,a,p_regen_gs_l)
	Vs_null_regen_gs_l_(a,b) = Vsnullcline(b,a,p_regen_gs_l)
	Vus_null_regen_gs_l_(a,b) = Vusnullcline(b,a,p_regen_gs_l)
	
	inter_V_Vs_null1_regen_gs_l = zeros(length(vus))
	inter_V_Vs_null2_regen_gs_l = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_gs_l[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_gs_l)
		inter_V_Vs_null2_regen_gs_l[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_gs_l)
	end

	inter_V_Vus_null1_regen_gs_l = zeros(size(vs))
	inter_V_Vus_null2_regen_gs_l = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_regen_gs_l[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_gs_l)
		inter_V_Vus_null2_regen_gs_l[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_gs_l)
	end


	plot(vus,vs,V1_null_regen_gs_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_gs_l_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_gs_l,inter_V_Vs_null2_regen_gs_l],[inter_V_Vs_null1_regen_gs_l,inter_V_Vs_null2_regen_gs_l],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_gs_l,inter_V_Vus_null2_regen_gs_l],[vs,vs],[inter_V_Vus_null1_regen_gs_l,inter_V_Vus_null2_regen_gs_l],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_gs_l[3,1:end],sol_regen_gs_l[2,1:end],sol_regen_gs_l[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gs=$(p_gs_l[7])")

end

# ╔═╡ 742dd2c3-040b-425d-b8a2-85e208af7b1e
begin
	plot(plt_v_null_regen)
	title!("V nullcline when gs=$(p_regen[7])")
end

# ╔═╡ 97d8c973-4841-422c-be66-34fb55115c5a
begin
	V1_null_regen_gs_h_(a,b) = Vnullcline1(b,a,p_regen_gs_h)
	V2_null_regen_gs_h_(a,b) = Vnullcline2(b,a,p_regen_gs_h)
	Vs_null_regen_gs_h_(a,b) = Vsnullcline(b,a,p_regen_gs_h)
	Vus_null_regen_gs_h_(a,b) = Vusnullcline(b,a,p_regen_gs_h)
	
	inter_V_Vs_null1_regen_gs_h = zeros(length(vus))
	inter_V_Vs_null2_regen_gs_h = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_regen_gs_h[i] = V_Vs_nullcline1(NaN,vus[i],p_regen_gs_h)
		inter_V_Vs_null2_regen_gs_h[i] = V_Vs_nullcline2(NaN,vus[i],p_regen_gs_h)
	end

	inter_V_Vus_null1_regen_gs_h = zeros(size(vs))
	inter_V_Vus_null2_regen_gs_h = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1_regen_gs_h[i] = V_Vus_nullcline1(vs[i],NaN,p_regen_gs_h)
		inter_V_Vus_null2_regen_gs_h[i] = V_Vus_nullcline2(vs[i],NaN,p_regen_gs_h)
	end


	plot(vus,vs,V1_null_regen_gs_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vus", ylabel="Vs",zlabel="V",xlims=(vus_min,vus_max),ylims=(vs_min,vs_max),zlims=(-55,-20))  #yellow-orange
	plot!(vus,vs,V2_null_regen_gs_h_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([vus,vus],[inter_V_Vs_null1_regen_gs_h,inter_V_Vs_null2_regen_gs_h],[inter_V_Vs_null1_regen_gs_h,inter_V_Vs_null2_regen_gs_h],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen_gs_h,inter_V_Vus_null2_regen_gs_h],[vs,vs],[inter_V_Vus_null1_regen_gs_h,inter_V_Vus_null2_regen_gs_h],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0") #blue
	plot!(sol_regen_gs_h[3,1:end],sol_regen_gs_h[2,1:end],sol_regen_gs_h[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
	title!("V nullcline when gs=$(p_gs_l[7])")

end

# ╔═╡ c43f1aac-7cc6-47d5-9fae-bac5a20c70a4
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=1,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],[inter_V_Vus_null1_I_h,inter_V_Vus_null2_I_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; I=$(p_I_h[1])",legend=:outertopright)
	plot!([inter_V_Vus_null1_regen_gs_h,inter_V_Vus_null2_gs_h],[inter_V_Vus_null1_gs_h,inter_V_Vus_null2_gs_h],[vs,vs],linecolor=RGB(0.2,0.2,0.9),linewidth=1,label="dV/dt =dVus/dt = 0; gs=$(p_gs_h[7])",legend=:outertopright)
	plot!([inter_V_Vus_null1_gs_l,inter_V_Vus_null2_gs_l],[inter_V_Vus_null1_gs_l,inter_V_Vus_null2_gs_l],[vs,vs],linecolor=RGB(0.0,0.0,0.5),linewidth=3,label="dV/dt =dVus/dt = 0; gs=$(p_gs_l[7])",legend=:outertopright,xlims=(-50,-25),zlims=(-50,-25))
	plot!((-39.8).*ones(length(vus)),vus,vus,linecolor="green")
	title!("Regenerative slow feedback")
end

# ╔═╡ ee174df3-a596-441e-a3d0-808013360e31
gr()

# ╔═╡ c7c54efe-a200-4818-adfa-2f6178837112
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vus",zlabel="Vs",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vus_null1_regen_gs_h,inter_V_Vus_null2_regen_gs_h],[inter_V_Vus_null1_regen_gs_h,inter_V_Vus_null2_regen_gs_h],[vs,vs],linecolor=RGB(0.2,0.8,1),linewidth=3,label="dV/dt =dVus/dt = 0; gs=$(p_gs_h[7])",legend=:outertopright) #pink
	plot!([inter_V_Vus_null1_regen_gs_l,inter_V_Vus_null2_regen_gs_l],[inter_V_Vus_null1_regen_gs_l,inter_V_Vus_null2_regen_gs_l],[vs,vs],linecolor=RGB(0.0,0.5,1),linewidth=3,label="dV/dt =dVus/dt = 0; gs=$(p_gs_l[7])",legend=:outertopright)
	title!("Regenerative slow feedback")
end

# ╔═╡ f15bc695-c9c5-4c5e-9867-0150235a1565
begin
	plot([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	
	plot!([inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[vs,vs],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V", ylabel="Vs",zlabel="Vus",label="dV/dt = dVus/dt = 0",camera=(0,0)) #blue
	plot!([inter_V_Vs_null1_regen_gs_h,inter_V_Vs_null2_regen_gs_h],[inter_V_Vs_null1_regen_gs_h,inter_V_Vs_null2_regen_gs_h],[vus,vus],linecolor="pink",linewidth=3,label="dV/dt =dVs/dt = 0; gs=$(p_regen_gs_h[7])",legend=:outertopright) #pink
	plot!([inter_V_Vs_null1_regen_gs_l,inter_V_Vs_null2_regen_gs_l],[inter_V_Vs_null1_regen_gs_l,inter_V_Vs_null2_regen_gs_l],[vus,vus],linecolor=RGB(1,0.7,1),linewidth=3,label="dV/dt =dVs/dt = 0; gs=$(p_regen_gs_l[7])",legend=:outertopright) #pink
	title!("Regenerative ultraslow feedback")
end

# ╔═╡ 298c4acd-c3a5-4f88-89e6-8d5b5fb249b0
md"###### gus variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ 23bbe085-afec-42e4-99ee-b863ce68d62f
md"###### Vs0 variation for restorative ultraslow feedback and regenerative slow feedback"

# ╔═╡ c7a96ad0-f6fa-49c2-8054-6ce8c08d5529


# ╔═╡ 0256fbec-2fa6-4c23-bbbf-e7da53dbc52f


# ╔═╡ 2dc59388-02bc-492c-92c0-cfd055e65eec


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
# ╟─01252cf7-d74c-4f02-977a-fcbf743dc056
# ╟─972118cb-6c8b-4618-9e66-dfdae5b343f6
# ╟─aefc61d1-50b0-4f86-9204-2935495de66d
# ╠═82eb0286-c992-4459-9d65-529e8ba30362
# ╟─89fb57df-ffb6-4a97-acb5-d2c293b3f1d3
# ╟─2eb57d95-e17c-4f71-8f9f-34a6df5b0392
# ╟─8698b02a-a34e-4e4f-9ecf-a709e57136a2
# ╟─02927149-ba9f-4da7-83b1-390a34d6f62f
# ╟─17518825-5527-4478-a4c3-52f7ff864a87
# ╠═371dc1f8-570b-4d77-a235-03b14efca425
# ╟─7ee5fec7-9299-4fb4-bc9a-8af9e94fec34
# ╟─70af372c-963d-4e53-b162-f3b2ca733b48
# ╟─1cfe6b28-6d83-4bfa-920a-b4eb387590cb
# ╟─2a7ef951-82a6-4433-8911-8936704b9424
# ╟─cfa819c9-a58a-44b3-bffd-3405fabaf6cc
# ╠═49b1f686-bcf7-49ed-93a3-aa74d0081110
# ╟─f494a04b-0a96-411e-bc22-1c7e8da6341c
# ╠═871c279c-7579-458b-b110-286389f710ad
# ╠═18ac0378-77c7-4570-a0d4-f687ec76d420
# ╠═e7d5ccae-c62d-4aae-910d-0cb3c172d7b8
# ╟─6ec8f2d4-8de5-480b-b8e4-1a3b49fa4a4a
# ╟─eedb2f8e-f389-4a80-9d7f-6dd3b5833110
# ╟─f1aba699-0d2c-4835-ab43-42bea662a141
# ╟─243fe1ba-9b96-4028-a8cb-4928153a7f19
# ╟─8ea57601-f46f-4488-8dd5-ec0189c5b751
# ╟─b990cc17-1ea4-468a-8a2e-39f897cade99
# ╟─410eee64-5670-4b57-be26-d45b9b4e5899
# ╟─5464ca79-59ea-4b55-adf9-8c17ec5afa7b
# ╟─c67e6b50-ee12-4361-ac57-1e602641fad6
# ╠═7011c21d-da15-4b7e-b9ac-f8d06bc08511
# ╟─7fe6cfa0-1f7d-498d-a8db-2afc98fdd838
# ╟─6a379d9f-05c3-4a3c-935b-b2b8bc96d43f
# ╟─91bb6de8-df2f-4598-a231-1f855591c65e
# ╠═7371f342-f293-497d-b5c4-bbd24f02967c
# ╟─8f0158b1-9664-4c27-b778-912c480e3ae8
# ╟─b1c9b459-7d75-4082-9876-a3ab480714c7
# ╟─b7a260fd-9f58-4fd8-a5c2-45997bf62bc8
# ╟─0365925b-e208-4f13-8908-747f238a2a4f
# ╟─59d2716d-3001-4499-b213-f3b3848bdb0f
# ╟─d0985722-90cc-4c25-8826-9a12a9eef8c2
# ╟─c69a5267-9b2f-4531-9b5f-27d8f56940e2
# ╟─97492156-9228-42f4-bab0-3175948559d0
# ╟─bd571ae8-74e1-480e-ac7f-be65d20ef89d
# ╠═7520f4b2-8948-474b-988b-23573b9311e6
# ╠═eeb780aa-3943-4bc9-8108-c5c4dfe6b773
# ╟─e156e461-1bac-427d-b6d1-e8ea5ce05433
# ╟─4f1874ab-569b-4145-93c6-1688c1d3aadf
# ╟─3a35b2f1-ba94-4597-bb1b-fba90a72930b
# ╠═772cc897-4364-4e80-97cc-65b064547884
# ╠═48f59077-f18e-469a-877d-9f7a4c5c4d10
# ╟─e05e543d-a6ba-4a26-a61c-9c05daba95c5
# ╠═f4ae01bc-4092-414e-a70d-7abe0c38cd6c
# ╟─c43f1aac-7cc6-47d5-9fae-bac5a20c70a4
# ╟─53d6b75b-41e1-44bb-a3ca-95863677b92b
# ╠═20b8ef9a-18e8-4aed-8a7e-ae927938fd7b
# ╟─c7875827-d1a8-4058-be16-e4a1437cbc6a
# ╟─fe6c6934-d5f6-4d05-88dc-845b9759122d
# ╟─dd575c93-296e-445c-bc2b-b27ab3f89b45
# ╟─3878b997-3ec2-42f7-9d71-195028344152
# ╟─628d70e4-4d6d-4d65-aa0c-f028ef6317d4
# ╟─d80fdc73-c420-49a0-9e6d-c8cfebfd7d60
# ╟─9d900972-fc6d-49b0-b1e3-3a980deb2b97
# ╟─081428d2-c3ca-421c-8f9e-e45d03b24407
# ╟─07a1ef09-5a86-4b05-9e8e-5355af702b1c
# ╠═742dd2c3-040b-425d-b8a2-85e208af7b1e
# ╟─97d8c973-4841-422c-be66-34fb55115c5a
# ╠═ee174df3-a596-441e-a3d0-808013360e31
# ╟─c7c54efe-a200-4818-adfa-2f6178837112
# ╟─f15bc695-c9c5-4c5e-9867-0150235a1565
# ╟─298c4acd-c3a5-4f88-89e6-8d5b5fb249b0
# ╟─23bbe085-afec-42e4-99ee-b863ce68d62f
# ╠═c7a96ad0-f6fa-49c2-8054-6ce8c08d5529
# ╠═0256fbec-2fa6-4c23-bbbf-e7da53dbc52f
# ╠═2dc59388-02bc-492c-92c0-cfd055e65eec
