### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ bd35d5de-7544-11eb-2023-9345c9da0b0f
begin
	using DifferentialEquations
	using LinearAlgebra
	using Plots
	using Statistics
end

# ╔═╡ 792f7770-95e7-11eb-23a5-45e952a738da
md" #### Packages"

# ╔═╡ 05913e12-7545-11eb-2eec-875820f00015
gr()

# ╔═╡ 8da532ce-95e7-11eb-1d64-a37503684718
md" #### Problem parameters"

# ╔═╡ 009dbc70-7546-11eb-1ae1-8dbac638d9e9
p=(4,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ 79a51975-67f2-4992-8b7c-d12176075249
p_regen=(4.5,-40.0,-38.4,-30.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ 9967a94d-ef13-45a0-b936-48bca8a64526
p_t=(5.5,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ 262a2b7f-7f42-4002-b91b-4bea4dea074c
p_regen_t=(5.5,-40.0,-38.4,-30.0,1.0,1.0,0.5,0.015,10.0,1000.0)

# ╔═╡ da7703c2-95ee-11eb-2ce4-399430fe71e5
function change_I(I)
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

# ╔═╡ 2b0581dc-f37c-4a44-bbe1-a16944a5e1ab
function change_I_(I,p)
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

# ╔═╡ 756351fa-5648-4cfb-a804-3fe5bda48cb9
function p3Dto2D(p3D)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p3D
	p2D=(I,v0,vs0,C,gf,gs,ts)
	return p2D
end

# ╔═╡ 616560f0-95e7-11eb-0ae3-ddcbdbe0be8a
md" #### Nullclines functions"

# ╔═╡ 63a575ae-7555-11eb-0cb8-43a1d04dabde
function Vnullcline1(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I) >= 0
		v1 = v0 + sqrt((gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I)/gf)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ 4fd67ae0-79e5-11eb-0278-838411a5c092
function Vnullcline2(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I) >= 0
		v2 = v0 - sqrt((gs*(vs-vs0)^2 +gus*(vus-vus0)^2 - I)/gf)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ 20f2eb20-79e3-11eb-268e-bbd0237d7af2
function Vsnullcline(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	v = vs
	return v
end

# ╔═╡ a9a5bb00-79e3-11eb-0ecf-2397fdc5fa21
function Vusnullcline(vs,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	v = vus
	return v
end

# ╔═╡ 939c4bf2-7ac2-11eb-1bc3-1bd1119508d6
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

# ╔═╡ d60571d0-7b60-11eb-0e78-93ffbb86bc10
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

# ╔═╡ e6f8e440-7b60-11eb-07b0-cba945d97277
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

# ╔═╡ 310f1bd2-7b61-11eb-3e3e-65e1334b589f
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

# ╔═╡ 6d492f43-e36b-4725-82da-905d02e83d73
md"###### Updated for x=V ; y=Vus ; z=Vs"

# ╔═╡ e4594047-eb27-44ec-aea2-c2df1defba77
function Vnullcline1_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gf*(v-v0)^2 -gus*(vus-vus0)^2 + I) >= 0
		vs1 = vs0 + sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 + I)/gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 3444011c-122f-44ce-b53a-ddf04e46f1dc
function Vnullcline2_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if (gf*(v-v0)^2 -gus*(vus-vus0)^2 + I) >= 0
		vs2 = vs0 - sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 + I)/gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 841b5bf9-8b23-469d-bf4d-31764ef5d900
function V_Vs_nullcline1_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	US = gus*(vus-vus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US)

	if delta >= 0 #there is 1 (2) fixed point
		vs1 = ((gf*v0)-(gs*vs0)+sqrt(delta))/(gf-gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ e2e956e1-e100-4315-83f1-557a2634c3f9
function V_Vs_nullcline2_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	#compute US component
	US = gus*(vus-vus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US)

	if delta >= 0 #there is 1 (2) fixed point
		vs2 = ((gf*v0)-(gs*vs0)-sqrt(delta))/(gf-gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ b2850211-edf0-4f48-aeff-942b23948ad8
function V_Vus_nullcline1_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if v == vus
		if gf*(v-v0)^2 - gus*(v-vus0)^2 +I >=0
			vs1 = vs0 + sqrt((gf*(v-v0)^2 - gus*(v-vus0)^2 +I)/gs)
		else
			vs1 = NaN
		end
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 1abb4876-62be-4df3-bd98-555ef0c3fd2a
function V_Vus_nullcline2_vs(v,vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	if v == vus
		if gf*(v-v0)^2 - gus*(v-vus0)^2 +I >=0
			vs2 = vs0 - sqrt((gf*(v-v0)^2 - gus*(v-vus0)^2 +I)/gs)
		else
			vs2 = NaN
		end
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 5d5293ed-0060-4615-b068-3acf16b2c1b5
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# ╔═╡ 6ffd59cc-13e2-41c7-8ac6-4031438f8780
vs_meshvsus, vus_meshvsus = meshgrid(-100.0:1.0:0.0, -100.0:1.0:0.0)

# ╔═╡ d5115900-95e7-11eb-34e4-ff02328e7981
md" #### Bifurcation functions defined based on analytical expression"

# ╔═╡ 5702cbc0-948e-11eb-2b59-1d53d6a15304
function bifI3D(p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	Ibif = (gf*gs*(v0-vs0)^2+gf*gus*(v0-vus0)^2-gs*gus*(vs0-vus0)^2)/(gf-gs-gus)
	return Ibif
end

# ╔═╡ 97c1b30e-95e7-11eb-2047-27179c477642
begin
	Ibif = bifI3D(p)
	md""" The current where we have a bifurcation of the 3D system (from a stable resting state to a limit cycle behavior) is $(round(Ibif*100)/100) for the parameters chosen for the restorative feedback

	For the parameters chosen for the regeneratove feedback, the bifurcation current is $(round(bifI3D(p_regen)*100)/100)"""
end

# ╔═╡ 5f26b7f2-9b35-4390-aa0e-29926da44b82
md" #### MQIF 2D functions"

# ╔═╡ e5b70336-c1c2-4ec4-858f-33a8d481599b
function Vnullcline1_2D(v,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gf*(v-v0)^2 + I) >= 0
		vs1 = vs0 + sqrt((gf*(v-v0)^2 + I)/gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 6bdb76b4-7b9e-4657-8b00-213e8453c8e4
function Vnullcline2_2D(v,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gf*(v-v0)^2 + I) >= 0
		vs2 = vs0 - sqrt((gf*(v-v0)^2 + I)/gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 14e62c79-d9b3-4f2d-9a8c-f16e2f0fd6b1
function Vnullcline1_2D_(vs,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gs*(vs-vs0)^2 - I) >= 0
		v1 = v0 + sqrt((gs*(vs-vs0)^2 - I)/gf)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ e1037f19-64a6-42f8-bef5-5a470b8e438d
function Vnullcline2_2D_(vs,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gs*(vs-vs0)^2 - I) >= 0
		v2 = v0 - sqrt((gs*(vs-vs0)^2 - I)/gf)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ 3a8758c2-84e3-47dc-bdac-52fe7bec5800
function Vsnullcline_2D(v,p_2D)
	return v
end

# ╔═╡ af7f680e-5209-46a3-a6b3-5323d1b66729
function bifI2D(p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D
	Ibif = (gf*gs*(v0-vs0)^2)/(gf-gs)
	return Ibif
end

# ╔═╡ df44c90b-868a-4da5-833d-2a508a3b8e91
function max_2D_low_V_null(I)
	gs = 0.5
	vs0 = -38.4

	if I>=0
		Vs_max = vs0 - sqrt(I/gs)
	else
		Vs_max = NaN
	end

	return Vs_max
end

# ╔═╡ d60d0341-703b-4908-9017-996b4d222c37
function fixedpoints_only_2D(p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D
	#with the vsnullcline, we know that vs = v, we then use this equality in the V-
	#nullcline and solve the system using the discriminant method

	a = (gf-gs)
	b = 2*(gs*vs0 - gf*v0)
	c = gf*v0^2 - gs*vs0^2 + I

	delta = b^2 - 4*a*c

	if delta <0
		#no real solution
		vfp1 = NaN
		vfp2 = NaN
	else
		vfp1 = (-b+sqrt(delta))/(2*a)
		vfp2 = (-b-sqrt(delta))/(2*a)
		stability = zeros(1,2) #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable
		fp = zeros(1,2)
		fp[1] = vfp1
		fp[2] = vfp2
	end

	fp1=ones(1,2).*vfp1
	fp2=ones(1,2).*vfp2

	return fp1,fp2
end

# ╔═╡ a3edd5c2-00bc-4925-a18d-e4b51b54f96e
function change_I_2D(I,p_2D)
	p=zeros(length(p_2D))
	for i=1:length(p)
		if i==1
			p[i]=I
		else
			p[i]=p_2D[i]
		end
	end
	return p
end

# ╔═╡ 2131d4e4-5c7a-4b6e-b6be-05d5e78ce8f9
function MQIF_2D!(du,u,p,t)
	I,v0,vs0,C,gf,gs,ts = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
end

# ╔═╡ 56d892b0-95e7-11eb-1b32-6db8dc815056
md" #### ODE problem functions"

# ╔═╡ 78f557b0-7545-11eb-2810-9bd1d2f85b33
function MQIF_3D!(du,u,p,t)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 - gus*(u[3]-vus0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
	du[3] = (u[1]-u[3])/tus
end

# ╔═╡ ddc72af0-7546-11eb-3d4e-f9ff68884993
begin
	Vmax = 30 #30
	Vr = -40 #-40
	Vsr = -35 #-20
	DVusr = 3

	md"""Voltages used for reset"""
end

# ╔═╡ 4fcc4093-0b2b-4645-82ba-b6bee0213695
function reset_2D!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
end

# ╔═╡ 8aa1c774-8026-4cac-bd58-f9ac89c5390c
# when condition == 0 and upcrossing (from negative to positive)
function affect_2D!(integrator)
    reset_2D!(integrator.u)
end

# ╔═╡ 8a1d2b20-7546-11eb-2d7a-494ae827de36
function spike(x)  # spikes when spike(x) goes from negative to positive
    (x[1] - Vmax)
end

# ╔═╡ bb367180-7546-11eb-361d-53889265d7be
# event when event_f(u,t) == 0
function condition(x,t,integrator) #
    spike(x)
end

# ╔═╡ fb9b0595-c11c-4be4-89df-6c07e0f9a02d
cb_2D = ContinuousCallback(condition,affect_2D!,nothing)

# ╔═╡ 9173ad8e-7546-11eb-1c3e-51967fc5a4ce
function reset!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
	x[3] = x[3] + DVusr
end

# ╔═╡ c0b88210-7546-11eb-1343-31a62c6e60b6
# when condition == 0 and upcrossing (from negative to positive)
function affect!(integrator)
    reset!(integrator.u)
end

# ╔═╡ e83ae550-95e7-11eb-11d3-0bb50817f70e
md" #### Simulations"

# ╔═╡ 93032350-9618-11eb-0a3f-95b734f703bb
plotly()

# ╔═╡ dc5c0f90-79e3-11eb-34f6-bfda108d916f
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
	vs_inter_l = collect(range(vs_min,stop=p[3]-4,length=30))
	vs_inter_h = collect(range(p[3]+4,stop=vs_max,length=30))
	vs_hole = collect(range(p[3]-4,stop=p[3]+4,length=600))

	lim_vus_hole = -50 + sqrt(lim_hole+p[1]/p[8]) + sqrt(lim_hole+p[1]/p[8])/10
	vus_inter_l = collect(range(vus_min,stop=-50,length=30))
	vus_inter_h = collect(range(lim_vus_hole,stop=vus_max,length=30))
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

# ╔═╡ e07ed7c0-7b63-11eb-3b80-8d22fcd902e3
gr()

# ╔═╡ a7ae4c37-3547-44ee-a9cc-8d3d3f5df221
md"""
Intuitition about 3D phase plane

V nullcline (yellow): (bottom part attractive, top part repulsive (excitable))

	Above the plane / under the plane / in the hole(basin) joining the 2 : dV/dt >0
	Between the 2 planes : dV/dt <0

Vs nullcline (pink): (attractive)

	Above the plane : dVs/dt >0
	Under the plane : dVs/dt <0

Vus nullcline (blue): (attractive)

	Above the plane : dVus/dt >0
	Under the plane : dVus/dt <0
"""

# ╔═╡ c28bddf0-79e5-11eb-1c82-014cf00fadd9
begin
	V1_null_(a,b) = Vnullcline1(a,b,p)
	V2_null_(a,b) = Vnullcline2(a,b,p)
	Vs_null_(a,b) = Vsnullcline(a,b,p)
	Vus_null_(a,b) = Vusnullcline(a,b,p)

	plot(vs,vus,V1_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	V_null_plt = plot!(vs,vus,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",zlims=(-55,-20)) #yellow-orange
	title!("dV/dt =0")

	Vs_null_plt = plot(vs,vus,Vs_null_,st=:surface,c=:RdPu_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",zlims=(-55,-20)) #pink
	title!("dVs/dt =0")

	Vus_null_plt = plot(vs,vus,Vus_null_,st=:surface,c=:PuBu_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",zlims=(-55,-20)) #blue
	title!("dVus/dt =0")

	null_plt_split = plot(V_null_plt,Vs_null_plt,Vus_null_plt,layout=(1,3),size=(2000,500),xlabel="Vs")

end

# ╔═╡ 4cb043a9-20aa-488a-8020-d9383d83b774
gr()

# ╔═╡ 55b08c31-a827-4ea1-bea6-e5e9ecf7d639
v=collect(range(-55,stop=-25,step=0.01))

# ╔═╡ b598cd61-be3c-4e09-85b2-0518a1b85054
v_=collect(range(-55,stop=-25,step=0.001))

# ╔═╡ 25c6a491-28a0-4e19-829f-7b790fc97d5f
vus_=collect(range(-100,stop=0,step=0.25))

# ╔═╡ ae0b4f6d-7269-44e5-9acc-4c618875b084
#plotly()

# ╔═╡ f1111e66-63d0-4a3e-bb59-49b729d0281c
begin 
	p2D = p3Dto2D(p)
	
	vus_1 = -80
	It_1 = p[1] - p[8]*(vus_1-p[4])^2
	
	vus_2 = -50
	It_2 = p[1] - p[8]*(vus_2-p[4])^2
	
	vus_3 = -20
	It_3 = p[1] - p[8]*(vus_3-p[4])^2
	
	Vnull1_2D_1 = zeros(size(v_))
	Vnull2_2D_1 = zeros(size(v_))
	Vnull1_2D_2 = zeros(size(v_))
	Vnull2_2D_2 = zeros(size(v_))
	Vnull1_2D_3 = zeros(size(v_))
	Vnull2_2D_3 = zeros(size(v_))
	for i=1:length(v_)
		Vnull1_2D_1[i] =Vnullcline1_2D(v_[i],change_I_(It_1,p2D))
		Vnull2_2D_1[i] =Vnullcline2_2D(v_[i],change_I_(It_1,p2D))
		Vnull1_2D_2[i] =Vnullcline1_2D(v_[i],change_I_(It_2,p2D))
		Vnull2_2D_2[i] =Vnullcline2_2D(v_[i],change_I_(It_2,p2D))
		Vnull1_2D_3[i] =Vnullcline1_2D(v_[i],change_I_(It_3,p2D))
		Vnull2_2D_3[i] =Vnullcline2_2D(v_[i],change_I_(It_3,p2D))
	end
	plot2D_1 = plot(v_,Vnull1_2D_1,c=:orange,colorbar_entry=false,linewidth=3)
	plot!(v_,Vnull2_2D_1,c=:orange,colorbar_entry=false,linewidth=3,legend=:none)
	xaxis!("V")
	yaxis!("Vs")
	title!("Vus = $vus_1")
	
	plot2D_2 = plot(v_,Vnull1_2D_2,c=:orange,colorbar_entry=false,linewidth=3)
	plot!(v_,Vnull2_2D_2,c=:orange,colorbar_entry=false,linewidth=3,legend=:none)
	xaxis!("V")
	yaxis!("Vs")
	title!("Vus = $vus_2")
	
	plot2D_3 = plot(v_,Vnull1_2D_3,c=:orange,colorbar_entry=false,linewidth=3)
	plot!(v_,Vnull2_2D_3,c=:orange,colorbar_entry=false,linewidth=3,legend=:none)
	xaxis!("V")
	yaxis!("Vs")
	title!("Vus = $vus_3")
	
	plot(plot2D_1,plot2D_2,plot2D_3,layout=(1,3),size=(1000,750))
	#savefig("projphaseplane3Dfor3vus.pdf")
end 

# ╔═╡ 31f42cfb-1e08-4120-b48d-d195788eadcb
change_I_(It,p2D)

# ╔═╡ 5f8af83d-577f-4a6b-acb4-0b67dbbdf113
begin
	V1_null_vs_(a,b) = Vnullcline1_vs(a,b,p)
	V2_null_vs_(a,b) = Vnullcline2_vs(a,b,p)
	
	inter_V_Vs_null1_vs = zeros(length(vus_))
	inter_V_Vs_null2_vs = zeros(length(vus_))
	for i=1:length(vus_)
		inter_V_Vs_null1_vs[i] = V_Vs_nullcline1_vs(NaN,vus_[i],p)
		inter_V_Vs_null2_vs[i] = V_Vs_nullcline2_vs(NaN,vus_[i],p)
	end

	inter_V_Vus_null1_vs = zeros(size(v))
	inter_V_Vus_null2_vs = zeros(size(v))
	for i=1:length(v)
		inter_V_Vus_null1_vs[i] = V_Vus_nullcline1_vs(v[i],v[i],p)
		inter_V_Vus_null2_vs[i] = V_Vus_nullcline2_vs(v[i],v[i],p)
	end

	plot(v,vus_,V1_null_vs_,st=:surface,c=:orange,colorbar_entry=false,xlabel="V", ylabel="Vus",zlabel="Vs",camera=(60,20)) #yellow-orange
	plot!(v,vus_,V2_null_vs_,st=:surface,c=:orange,colorbar_entry=false,label="dV/dt = 0") #yellow-orange
	plot!([inter_V_Vs_null1_vs,inter_V_Vs_null2_vs],[vus_,vus_],[inter_V_Vs_null1_vs,inter_V_Vs_null2_vs],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink
	plot!([v,v],[v,v],[inter_V_Vus_null1_vs,inter_V_Vus_null2_vs],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0",xlims=(-55,-20),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	#title!("V nullcline and its intersections")

end

# ╔═╡ 362b5e4f-f9ba-4e41-9164-e22dcb1a6731
gr()

# ╔═╡ 81513370-79eb-11eb-2257-f39b0b5789af
begin
	inter_V_Vs_null1 = zeros(length(vus))
	inter_V_Vs_null2 = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1[i] = V_Vs_nullcline1(NaN,vus[i],p)
		inter_V_Vs_null2[i] = V_Vs_nullcline2(NaN,vus[i],p)
	end



	inter_V_Vus_null1 = zeros(size(vs))
	inter_V_Vus_null2 = zeros(size(vs))
	for i=1:length(vus)
		inter_V_Vus_null1[i] = V_Vus_nullcline1(vs[i],NaN,p)
		inter_V_Vus_null2[i] = V_Vus_nullcline2(vs[i],NaN,p)
	end


	plot(vs,vus,V1_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	plot!(vs,vus,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #pink
	plot!([vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	#plot!(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")


end

# ╔═╡ 8e4240b0-95e8-11eb-1d42-cdae28c3935c
md"""Observations :

When the current is below the bifurcation value, all nullclines intersect which creates a stable node making the system converge towards a stable resting state.

When the current is above this specific value, the system evolves towards a stable limit cycle. When the current stays low but above the bifurcation value, the intsntataneous frequency of spikes decreases and the trajectory becomes attracted in the bottleneck region where the stable node is when the current is lower.

When the current is further increased, the bottleneck region is extended and therefore we do not observe long time intervals where the system does not oscillate as we did when the current is above but close to the bifurcation value. """

# ╔═╡ 12d22342-f5ec-4a88-ab9e-5c4d76f0f92a
gr()

# ╔═╡ 7ba8807a-3f9f-4919-a267-70539fa955a4
plotly()

# ╔═╡ 3d22c39f-d723-4c99-87b6-6b60f861d816
p_regen

# ╔═╡ 381d77d7-0723-47e9-b724-fed87d2f93a9
begin
	V1_null_regen_(a,b) = Vnullcline1(a,b,p_regen)
	V2_null_regen_(a,b) = Vnullcline2(a,b,p_regen)
	Vs_null_regen_(a,b) = Vsnullcline(a,b,p_regen)
	Vus_null_regen_(a,b) = Vusnullcline(a,b,p_regen)

	plot(vs,vus,V1_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	V_null_plt_regen = plot!(vs,vus,V2_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",zlims=(-55,-20)) #yellow-orange
	title!("dV/dt =0")

	Vs_null_plt_regen = plot(vs,vus,Vs_null_regen_,st=:surface,c=:RdPu_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",zlims=(-55,-20)) #pink
	title!("dVs/dt =0")

	Vus_null_plt_regen = plot(vs,vus,Vus_null_regen_,st=:surface,c=:PuBu_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",zlims=(-55,-20)) #blue
	title!("dVus/dt =0")

	null_plt_split_regen = plot(V_null_plt_regen,Vs_null_plt_regen,Vus_null_plt_regen,layout=(1,3),size=(2000,500),xlabel="Vs")
end

# ╔═╡ 20fa3ea1-7a0e-43ea-947b-ce78a8929e83
begin 
	V1_null_t_(a,b) = Vnullcline1(a,b,p_t)
	plot(vs,vus,V1_null_t_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	plot!(vs,vus,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright,size=(200,200)) #yellow-orange
end

# ╔═╡ 13feea03-5093-4c99-bae0-d88da15a2925
begin
	V1_null_regen_vs_(a,b) = Vnullcline1_vs(a,b,p_regen)
	V2_null_regen_vs_(a,b) = Vnullcline2_vs(a,b,p_regen)
	
	inter_V_Vs_null1_regen_vs = zeros(length(vus_))
	inter_V_Vs_null2_regen_vs = zeros(length(vus_))
	for i=1:length(vus_)
		inter_V_Vs_null1_regen_vs[i] = V_Vs_nullcline1_vs(NaN,vus_[i],p_regen)
		inter_V_Vs_null2_regen_vs[i] = V_Vs_nullcline2_vs(NaN,vus_[i],p_regen)
	end

	inter_V_Vus_null1_regen_vs = zeros(size(v))
	inter_V_Vus_null2_regen_vs = zeros(size(v))
	for i=1:length(v)
		inter_V_Vus_null1_regen_vs[i] = V_Vus_nullcline1_vs(v[i],v[i],p_regen)
		inter_V_Vus_null2_regen_vs[i] = V_Vus_nullcline2_vs(v[i],v[i],p_regen)
	end

	plot(v,vus_,V1_null_regen_vs_,st=:surface,c=:orange,colorbar_entry=false,xlabel="V", ylabel="Vus",zlabel="Vs",camera=(60,20)) #yellow-orange
	plot!(v,vus_,V2_null_regen_vs_,st=:surface,c=:orange,colorbar_entry=false,label="dV/dt = 0") #yellow-orange
	plot!([inter_V_Vs_null1_regen_vs,inter_V_Vs_null2_regen_vs],[vus_,vus_],[inter_V_Vs_null1_regen_vs,inter_V_Vs_null2_regen_vs],linecolor=RGB(1,0.4,0.7),linewidth=5,label="dVs/dt = 0") #pink 
	plot!([v,v],[v,v],[inter_V_Vus_null1_regen_vs,inter_V_Vus_null2_regen_vs],linecolor=RGB(0.31,0.66,1),linewidth=5,label="dVus/dt = 0",xlims=(-55,-20),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	title!("V nullcline and its intersections")

end

# ╔═╡ ec0e5f62-2114-43b3-931f-324cb45046a5
md"""
Observations :

If we change the ultraslow feedback to be regenerative, we see that the fixed points appear in the right part (Vus low) of the nullcline intersection dV/dt = 0 & dVs/dt =0. For a restorative ultrasolw feedback, the fixed points appeared in the left part (Vus high) of this nullcline. Concretely, changing the value of Vus0 moves the place where the basin between the 2 planes of the V nullcline sits.

In fact, when the feedback is regenerative, the nullcline intersection dV/dt = 0 & dVus/dt =0(blue) is much far away that the trajectory (is at the other side of the basin) making the gradient of Vus much bigger making it able to increase enough to join back the basin and to produce a spike again. Thus, there is no period of rest. When the feedback is restorative, the nullcline intersection is on the good side (since the trajectory follows closely the V nullcline)(also, the system tends to spike and thus increase Vus is not attracted towards a stable resting state)

Based on this, it makes also sense that the system automatically converges towards a stable resting state as soon as there is a stable node in the system for a restorative ultraslow feedback."""

# ╔═╡ e2726064-c076-4b44-98d7-f80ce2d572aa
begin
	vs_regen = zeros(length(vs_inter_l)+length(vs_inter_h)+length(vs_hole))
	vus_regen = zeros(length(vus_inter_l)+length(vus_inter_h)+length(vus_hole))

	for i=1:(length(vs_regen))
		if i <= length(vs_inter_l)
			vs_regen[i] = vs_inter_l[i]
		else
			if i<= length(vs_inter_l)+length(vs_hole)
				vs_regen[i] = vs_hole[i-length(vs_inter_l)]
			else
				vs_regen[i] = vs_inter_h[i-(length(vs_inter_l)+length(vs_hole))]
			end
		end
	end

	for i=1:(length(vus))
		if i <= length(vus_inter_l)
			vus_regen[i] = vus_inter_l[i]
		else
			if i<= length(vus_inter_l)+length(vus_hole)
				vus_regen[i] = vus_hole[i-length(vus_inter_l)]
			else
				vus_regen[i] = vus_inter_h[i-(length(vus_inter_l)+length(vus_hole))]
			end
		end
	end


	md"""Nullclines calculation"""
end

# ╔═╡ b19620bd-f7eb-4e53-b4e7-4c7b42675d31
begin 
	V1_null_regen_t_(a,b) = Vnullcline1(a,b,p_regen_t)
	plot(vs,vus,V1_null_regen_t_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	plot!(vs_regen,vus_regen,V2_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright,size=(200,200)) #yellow-orange
end

# ╔═╡ 7ebed900-95e8-11eb-0e9d-d7bd5f6748b1
gr()

# ╔═╡ 82582c60-95e8-11eb-1a32-3998de2baa18
md" #### Bifurcation with I "

# ╔═╡ f512e270-95ef-11eb-0ba6-91e5e91857c2
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

# ╔═╡ 357243f0-95f1-11eb-37be-b58f6942268f
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

# ╔═╡ 2497dbf0-95ea-11eb-2f36-dfcad08b86da
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

# ╔═╡ 3cce0240-960c-11eb-205e-bb2765422b75
begin
	V_fp_simu,stab_fp_simu = global_I_bifurcation(p)
	md"""Fixed points calculation"""
end

# ╔═╡ ded847d4-16e7-4e77-ab39-8f2ab65d66f1
begin
	V_fp_simu_regen,stab_fp_simu_regen = global_I_bifurcation(p_regen)
	md"""Fixed points calculation"""
end

# ╔═╡ 243d6310-95ee-11eb-3a7e-c19eb0280760
md"""###### Global bifurcation"""

# ╔═╡ 7938db90-960e-11eb-1815-b7ddd5a1853c
md""" It is assumed that there is no limit cycle when a stable node exists in the system. Indeed, we see that for a current higher than the bifurcation current, the trajectory has to go in the bottleneck between the nullclines of Vs and Vus to restart the limit cycle"""

# ╔═╡ 342d6f3e-95ee-11eb-1e02-f15eeda830f5
begin
	I_global_bif = collect(range(0.025,stop=50,step=0.025))
	V_fp_global_bif = zeros(length(I_global_bif),2)
	stab_fp_global_bif = zeros(length(I_global_bif),2)
	
	V_fp_global_bif_r = zeros(length(I_global_bif),2)
	stab_fp_global_bif_r = zeros(length(I_global_bif),2)

	for i=1:length(I_global_bif)
		p_global_bif = change_I(I_global_bif[i])
		p_global_bif_r = change_I_(I_global_bif[i],p_regen)
		V_fp_global_bif[i,:],stab_fp_global_bif[i,:]=global_I_bifurcation(p_global_bif)
		V_fp_global_bif_r[i,:],stab_fp_global_bif_r[i,:]=global_I_bifurcation(p_global_bif_r)
	end

	ind_bif = findfirst(x-> x>Ibif,I_global_bif)-1
	Vbif = mean(V_fp_global_bif[ind_bif,:])
	
	Ibif_r = bifI3D(p_regen)
	ind_bif_r = findfirst(x-> x>Ibif_r,I_global_bif)-1
	Vbif_r = mean(V_fp_global_bif_r[ind_bif_r,:])
	md"""Fixed point computation for a range of current values"""
end

# ╔═╡ 0dc4c81e-7547-11eb-3439-3d4a5750c2cb
begin
	u0=[Vbif,Vbif,Vbif]
	tspan = (0.0,5000.0)
	cb   = ContinuousCallback(condition,affect!,nothing)
	prob = ODEProblem(MQIF_3D!,u0,tspan,p,callback=cb)
	sol = solve(prob,dense=false)
end

# ╔═╡ 56bf7d90-7547-11eb-054a-6ff0654fd396
begin
	plot(sol.t,sol[1,:],linecolor="blue",label="V(t)")
	plot!(sol.t,sol[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol.t,sol[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol.t,p[3]ones(size(sol.t)),linecolor="red",label="Vs0")
	plot!(sol.t,p[4]ones(size(sol.t)),linecolor="magenta",label="Vus0")
end

# ╔═╡ eb8f7590-9489-11eb-353d-dff43f1b970c
begin
plot(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory",xlabel="Vs", ylabel="Vus",zlabel="V")
	md"""Sytem trajectory for restorative feedabck in the 3D phase plane"""
end

# ╔═╡ d8db632e-9490-11eb-27d0-73a81a9eeab6
begin
	plot(sol[2,1:end],sol[3,1:end],sol[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory",xlabel="Vs", ylabel="Vus",zlabel="V")
	plot!([inter_V_Vs_null1,inter_V_Vs_null2],[vus,vus],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #pink
	plot!([vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	for i=1:length(stab_fp_simu)
		if stab_fp_simu[i] > 0
			plot!([V_fp_simu[i]],[V_fp_simu[i]], [V_fp_simu[i]],
     seriestype=:scatter,markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
		end
		if stab_fp_simu[i] < 0
			plot!([V_fp_simu[i]],[V_fp_simu[i]], [V_fp_simu[i]],
     seriestype=:scatter,markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
		end
		if stab_fp_simu[i] == 0
			plot!([V_fp_simu[i]],[V_fp_simu[i]], [V_fp_simu[i]],
     seriestype=:scatter,markershape = :xcross,markersize = 2,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
		end
	end
	plot!()
	md"""System trajectory and nullclines intersection in 3D phase plane for param p"""
end

# ╔═╡ 0a06f401-42db-430d-acce-6cb609d654d9
begin
	u0_regen =[-10.0,-20.0,-50.0]
	tspan_regen=(0.0,20000.0)
	prob_regen = ODEProblem(MQIF_3D!,u0_regen,tspan_regen,p_regen,callback=cb)
	sol_regen = solve(prob_regen,dense=false)
end

# ╔═╡ 37e4ba5e-33c3-4083-a483-f20b111c217d
begin
	plot(sol_regen.t,sol_regen[1,:],linecolor="blue",label="V(t)")
	plot!(sol_regen.t,sol_regen[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_regen.t,sol_regen[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_regen.t,p_regen[3]ones(size(sol_regen.t)),linecolor="red",label="Vs0")
	plot!(sol_regen.t,p_regen[4]ones(size(sol_regen.t)),linecolor="magenta",label="Vus0")
end

# ╔═╡ 8aae255d-a738-4236-a4a9-15ad137eff0a
begin
	inter_V_Vs_null1_regen = zeros(length(vus_regen))
	inter_V_Vs_null2_regen = zeros(length(vus_regen))
	for i=1:length(vus_regen)
		inter_V_Vs_null1_regen[i] = V_Vs_nullcline1(NaN,vus_regen[i],p_regen)
		inter_V_Vs_null2_regen[i] = V_Vs_nullcline2(NaN,vus_regen[i],p_regen)
	end



	inter_V_Vus_null1_regen = zeros(size(vs_regen))
	inter_V_Vus_null2_regen = zeros(size(vs_regen))
	for i=1:length(vus_regen)
		inter_V_Vus_null1_regen[i] = V_Vus_nullcline1(vs_regen[i],NaN,p_regen)
		inter_V_Vus_null2_regen[i] = V_Vus_nullcline2(vs_regen[i],NaN,p_regen)
	end


	plot(vs_regen,vus_regen,V1_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	plot!(vs_regen,vus_regen,V2_null_regen_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus_regen,vus_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #pink
	plot!([vs_regen,vs_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	plot!(sol_regen[2,1:end],sol_regen[3,1:end],sol_regen[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
end

# ╔═╡ 781d46f6-0a91-408c-a2bb-4c78cbb91150
plot(sol_regen[2,1:end],sol_regen[3,1:end],sol_regen[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory",xlabel="Vs", ylabel="Vus",zlabel="V")

# ╔═╡ 2f225060-6402-4fae-bced-192328ecf4f7
begin
	plot(sol_regen[2,1:end],sol_regen[3,1:end],sol_regen[1,1:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory",xlabel="Vs", ylabel="Vus",zlabel="V")
	plot!([inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],[vus_regen,vus_regen],[inter_V_Vs_null1_regen,inter_V_Vs_null2_regen],linecolor=RGB(1,0.4,0.7),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #pink
	plot!([vs_regen,vs_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],[inter_V_Vus_null1_regen,inter_V_Vus_null2_regen],linecolor=RGB(0.31,0.66,1),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	for i=1:length(stab_fp_simu_regen)
		if stab_fp_simu_regen[i] > 0
			plot!([V_fp_simu_regen[i]],[V_fp_simu_regen[i]], [V_fp_simu_regen[i]],
     seriestype=:scatter,markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
		end
		if stab_fp_simu_regen[i] < 0
			plot!([V_fp_simu_regen[i]],[V_fp_simu_regen[i]], [V_fp_simu_regen[i]],
     seriestype=:scatter,markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
		end
		if stab_fp_simu_regen[i] == 0
			plot!([V_fp_simu_regen[i]],[V_fp_simu_regen[i]], [V_fp_simu_regen[i]],
     seriestype=:scatter,markershape = :xcross,markersize = 2,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
		end
	end
	plot!()
end

# ╔═╡ f30c88e3-1bc0-4eaa-96a0-37aa7ce9a25c
begin 	plot(vs,vus,V1_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	plot!(vs,vus,V2_null_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",legend=:outertopright) #yellow-orange
	plot!([inter_V_Vs_null1,inter_V_Vs_null2],[vus_regen,vus_regen],[inter_V_Vs_null1,inter_V_Vs_null2],linecolor=RGB(1,0.4,0.7),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVs/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #pink
	plot!([vs,vs],[inter_V_Vus_null1,inter_V_Vus_null2],[inter_V_Vus_null1,inter_V_Vus_null2],linecolor=RGB(0.31,0.66,1),linewidth=5,xlabel="Vs", ylabel="Vus",zlabel="V",label="dVus/dt = 0",xlims=(vs_min,vs_max),ylims=(vus_min,vus_max),zlims=(-55,-20)) #blue
	plot!(sol[2,end-100:end],sol[3,end-100:end],sol[1,end-100:end],linecolor=RGB(0.1,0.8,0.5),linewidth=2,label="Trajectory")
end

# ╔═╡ 63230109-92cf-44bc-8b95-4489053cb92a
function find_Vspikes(solt,solv,tspan_max)
	i_tmin = findfirst(x -> x>=50*tspan_max/100, solt)
	solv_=solv[i_tmin:end]
	i_spike_ = findall(x -> x >= (Vmax-1), solv_)
	i_spike=[]
	for i=1:length(i_spike_)
		if solv_[i_spike_[i]+1]<solv_[i_spike_[i]]
			append!(i_spike,i_spike_[i])
		end
	end
	return i_spike.+ (i_tmin-1)
end

# ╔═╡ ff4c958c-1a63-4a3f-a3b8-5e8aad4e5992
function find_Vusspikes(solt,solvus,tspan_max)
	i_tmin = findfirst(x -> x>=50*tspan_max/100, solt)
	solvus_=solvus[i_tmin:end]
	i_spike_ = findall(x -> x <= minimum(solvus_)+1, solvus_)
	i_spike=[]
	for i=1:length(i_spike_)
		if i_spike_[i]<length(solvus_)
			if i_spike_[i]>1
				if solvus_[i_spike_[i]+1]>solvus_[i_spike_[i]] && solvus_[i_spike_[i]-1]>solvus_[i_spike_[i]]
					append!(i_spike,i_spike_[i])
				end
			end
		end
	end
	return i_spike.+ (i_tmin-1)
end

# ╔═╡ a6cfb516-8f6c-43a2-a25d-672994a754eb
function burst_frequency(cycle,i,sol,tspan)
	if cycle==true 
			ind_spikeV = find_Vspikes(sol.t,sol[1,:],maximum(tspan)) 
			ind_spikeVus = find_Vusspikes(sol.t,sol[3,:],maximum(tspan))
			
			delta_t_Vus = sol.t[ind_spikeVus[2:end]] - sol.t[ind_spikeVus[1:end-1]]
			delta_i_Vus = ind_spikeVus[end]-ind_spikeVus[end-1]
			t_period_inter= mean(delta_t_Vus)
			
			
			delta_t_V = sol.t[ind_spikeV[2:end]] - sol.t[ind_spikeV[1:end-1]]
			i_p_n_1 =findfirst(x-> x>=ind_spikeVus[end-1],ind_spikeV)
			i_p_n_2 =length(ind_spikeV) - (findfirst(x-> x<=ind_spikeVus[end],reverse(ind_spikeV))-1)
			
			n_spike= i_p_n_2-i_p_n_1+1
			
			if n_spike >1
				t_period_intra1=maximum(delta_t_V[i_p_n_1:i_p_n_2-1])
				t_period_intra2=minimum(delta_t_V[i_p_n_1:i_p_n_2-1])
			
				ar = delta_t_V[i_p_n_1:i_p_n_2-1]
				t_period_intra_av = mean(ar) 
			else
				t_period_intra1=NaN
				t_period_intra2=NaN
				t_period_intra_av = NaN
			end

	else
		t_period_inter = NaN	
		t_period_intra1=NaN
		t_period_intra2=NaN
		t_period_intra_av = NaN
		n_spike= NaN
	end	
	return t_period_inter,t_period_intra1,t_period_intra2,t_period_intra_av, n_spike
end

# ╔═╡ ce60c8d0-960e-11eb-17d1-a1855ad6a7c4
begin
	I_global_bif_cycle = collect(range(3.75,stop=50,step=0.025))
	V_cycle_global_bif = zeros(length(I_global_bif_cycle),2) #[max,min]
	t_period_inter = zeros(length(I_global_bif_cycle)) #[max,min]
	
	t_period_intra = zeros(length(I_global_bif_cycle),2) #[max,min]
	t_period_intra_av = zeros(length(I_global_bif_cycle))
	n_spike = zeros(length(I_global_bif_cycle))
	for i=1:length(I_global_bif_cycle)
		p_global_bif_cycle=change_I(I_global_bif_cycle[i])
		tspan = (0.0,12000.0)
		u0=[-40.0,-40.0,-40.0]
		prob = ODEProblem(MQIF_3D!,u0,tspan,p_global_bif_cycle,callback=cb)
		sol = solve(prob,reltol=1e-5,abstol=1e-5)
		ind_tmin = findfirst(x -> x>=50*maximum(tspan)/100, sol.t)
		v_bob = sol[1,ind_tmin:end]
		if maximum(v_bob)-minimum(v_bob) >40
			V_cycle_global_bif[i,1]=maximum(v_bob)
			V_cycle_global_bif[i,2]=minimum(v_bob)
			
			cycle=true
		else
			V_cycle_global_bif[i,1]=NaN
			V_cycle_global_bif[i,2]=NaN
			
			cycle=false
		end
		t_period_inter[i],t_period_intra[i,1],t_period_intra[i,2],t_period_intra_av[i], n_spike[i] = burst_frequency(cycle,i,sol,tspan)
	end
	md"""Try to calculate max and min of the cycle in the 3D model after the SN bif but not smooth"""
end

# ╔═╡ 37ebe428-4359-4203-a2e9-93cdcc821c58
function burst_frequency_(cycle,i,sol,tspan)
	if cycle==true 
			ind_spikeV = find_Vspikes(sol.t,sol[1,:],maximum(tspan)) 
			ind_spikeVus = find_Vusspikes(sol.t,sol[3,:],maximum(tspan))
			
			delta_t_Vus = sol.t[ind_spikeVus[2:end]] - sol.t[ind_spikeVus[1:end-1]]
			delta_i_Vus = ind_spikeVus[end]-ind_spikeVus[end-1]
			t_period_inter= mean(delta_t_Vus)
			
			
			delta_t_V = sol.t[ind_spikeV[2:end]] - sol.t[ind_spikeV[1:end-1]]
			i_p_n_1 =findfirst(x-> x>=ind_spikeVus[end-1],ind_spikeV)
			i_p_n_2 =length(ind_spikeV) - (findfirst(x-> x<=ind_spikeVus[end],reverse(ind_spikeV))-1)
			
			n_spike= i_p_n_2-i_p_n_1+1
			t_period_intra_full =[]
			
			for j=i_p_n_1:i_p_n_2-1
				append!(t_period_intra_full,delta_t_V[j])
			end
		
			if n_spike >1
				t_period_intra1=maximum(delta_t_V[i_p_n_1:i_p_n_2-1])
				t_period_intra2=minimum(delta_t_V[i_p_n_1:i_p_n_2-1])
			
				ar = delta_t_V[i_p_n_1:i_p_n_2-1]
				t_period_intra_av = mean(ar) 
			else
				t_period_intra1=NaN
				t_period_intra2=NaN
				t_period_intra_av = NaN
			end

	else
		t_period_inter = NaN	
		t_period_intra1=NaN
		t_period_intra2=NaN
		t_period_intra_av = NaN
		n_spike= NaN
		t_period_intra_full = [NaN]	
	end	
	return t_period_inter,t_period_intra1,t_period_intra2,t_period_intra_av, n_spike,t_period_intra_full
end

# ╔═╡ 03e6717a-4959-4c54-bd2b-50c78e7282c7
function frequency_from_period(T)
	F=zeros(length(T))
	for i=1:length(T)
		F[i]=1/T[i]
	end
	return F
end

# ╔═╡ 84475afd-4453-4c85-83b3-951198d85fba
begin
	tspan_ref = (0.0,12000.0)
	prob_ref = ODEProblem(MQIF_3D!,u0,tspan_ref,p,callback=cb)
	sol_ref = solve(prob_ref,dense=false)
	
	t_period_inter_ref,t_period_intra1_ref,t_period_intra2_ref,t_period_intra_av_ref, n_spike_ref,t_period_intra_full_ref = burst_frequency_(true,NaN,sol_ref,tspan_ref)
	
	f_intra_full_ref=frequency_from_period(t_period_intra_full_ref[:]).*1000
	
	t_burst = ceil(sum(t_period_intra_full_ref))
	t_period = collect(range(0.0,stop=t_burst,step=0.1))
	f_intra = zeros(size(t_period))
	ind_t_period_s = [1]
	for j=1:length(t_period_intra_full_ref)
		ind_t_period_s_ = findfirst(x -> x>=sum(t_period_intra_full_ref[1:j]),t_period)
		f_intra[ind_t_period_s[end]:ind_t_period_s_].= f_intra_full_ref[j]
		append!(ind_t_period_s,ind_t_period_s_)
	end
end

# ╔═╡ ed9077c2-b76b-40e7-8d11-7c1a0b144df6
begin
	plot(t_period,f_intra,label="Burst instant frequency",size=(400,200))
	xaxis!("Time (ms)")
	yaxis!("Frequency [Hz]")
	#savefig("MQIF_3D_intraburst_f.pdf")
end

# ╔═╡ e0ff370f-d23e-44bb-b14a-9b6f1972421d
begin
	f_inter = frequency_from_period(t_period_inter).*1000
	f_intra_min = frequency_from_period(t_period_intra[:,1]).*1000
	f_intra_max = frequency_from_period(t_period_intra[:,2]).*1000
	f_intra_av = frequency_from_period(t_period_intra_av).*1000
end

# ╔═╡ e914d015-7dd2-4df0-af05-d298fbd262f5
begin
	plot(I_global_bif_cycle,f_intra_av,label="Intraburst average frequency")
	plot!(I_global_bif_cycle,f_intra_max,fill = (f_intra_min, 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Min max interburst frequency",size=(450,350),legend=:bottomright)
	#plot!(I_global_bif,f_intra_av_r,label="Intraburst average frequency")
	#plot!(I_global_bif,f_intra_max_r,fill = (f_intra_min_r, 0.5, :orange),linecolor=:orange,linealpha=0.1,label="Min max interburst frequency",size=(450,350),legend=:bottomright)
	yaxis!("Frequency [Hz]")
	xaxis!("Current",(0,40))
	#savefig("MQIF_3D_intraburst_I_f.pdf")
end

# ╔═╡ 0f28b311-cc76-40a8-8cea-cb357e450834
begin 
	plot_fr_inter = plot(I_global_bif_cycle,f_inter,linealpha=0.5,label="Interburst frequency      ")
	#plot!(I_global_bif,f_inter_r,linealpha=0.5,label="Interburst frequency regen  ")
	yaxis!("Frequency [Hz]")
	xaxis!("Current",(0,40))
	plot_n_fr = scatter(I_global_bif_cycle,n_spike,markershape = :circle,markersize = 1,markeralpha = 0.6,markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Intraburst spike number")
	#plot!(I_global_bif,n_spike_r,linealpha=0.5,label="Intraburst spike number reg")
	yaxis!("Number of spike")
	xaxis!("Current",(0,40))
	plot(plot_fr_inter,plot_n_fr,layout=(2,1),legend=:outertopright)
	#savefig("MQIF_3D_interburst_I_f_.pdf")
end

# ╔═╡ b1fb80f0-f9b2-470e-a628-c14a0b6a4dbe
begin
	V_cycle_global_bif_r = zeros(length(I_global_bif),2) #[max,min]
	V_period_r = zeros(length(I_global_bif),2) #[max,min]
	t_period_inter_r = zeros(length(I_global_bif)) #[max,min]
	
	t_period_intra_r = zeros(length(I_global_bif),2) #[max,min]
	t_period_intra_av_r = zeros(length(I_global_bif))
	n_spike_r = zeros(length(I_global_bif))
	for i=1:length(I_global_bif)
		p_global_bif_r=change_I_(I_global_bif[i],p_regen)
		tspan = (0.0,15000.0)
		u0=[p_global_bif_r[2],p_global_bif_r[3],p_global_bif_r[4]+10]
		prob = ODEProblem(MQIF_3D!,u0,tspan,p_global_bif_r,callback=cb)
		sol = solve(prob,reltol=1e-5,abstol=1e-5)
		ind_tmin = findfirst(x -> x>=50*maximum(tspan)/100, sol.t)
		v_bob = sol[1,ind_tmin:end]
		if maximum(v_bob)-minimum(v_bob) >40
			V_cycle_global_bif_r[i,1]=maximum(v_bob)
			V_cycle_global_bif_r[i,2]=minimum(v_bob)
			
			cycle=true
		else
			V_cycle_global_bif_r[i,1]=NaN
			V_cycle_global_bif_r[i,2]=NaN
			
			cycle=false
		end
		t_period_inter_r[i],t_period_intra_r[i,1],t_period_intra_r[i,2],t_period_intra_av_r[i], n_spike_r[i] = burst_frequency(cycle,i,sol,tspan)
	end
	md"""Try to calculate max and min of the cycle in the 3D model after the SN bif but not smooth"""
end

# ╔═╡ c550a58f-615e-442e-862a-baca5f741227
begin
	f_inter_r = frequency_from_period(t_period_inter_r).*1000
	f_intra_min_r = frequency_from_period(t_period_intra_r[:,1]).*1000
	f_intra_max_r = frequency_from_period(t_period_intra_r[:,2]).*1000
	f_intra_av_r = frequency_from_period(t_period_intra_av_r).*1000
end

# ╔═╡ cb6382ce-9f1d-40b6-99b9-e9f3f5133eed
begin
	p_global_bif_cycledebug=change_I_(4,p)
	p_global_bif_cycledebug2=change_I_(5,p)
	p_global_bif_cycledebug3=change_I_(5,p_regen)
	tspandebug = (0.0,4000.0)
	u0debug=[-40.0,-40.0,-40.0]
	probdebug = ODEProblem(MQIF_3D!,u0debug,tspandebug,p_global_bif_cycledebug,callback=cb)
	soldebug = solve(probdebug,reltol=1e-5,abstol=1e-5)
	probdebug2 = ODEProblem(MQIF_3D!,u0debug,tspandebug,p_global_bif_cycledebug2,callback=cb)
	soldebug2 = solve(probdebug2,reltol=1e-6,abstol=1e-6)
	probdebug3 = ODEProblem(MQIF_3D!,u0debug,tspandebug,p_global_bif_cycledebug3,callback=cb)
	soldebug3 = solve(probdebug3,reltol=1e-5,abstol=1e-5)
	ind_tmindebug3 = findfirst(x -> x>=50*maximum(tspandebug)/100, soldebug3.t)
	ind_tmindebug2 = findfirst(x -> x>=50*maximum(tspandebug)/100, soldebug2.t)
	ind_tmindebug = findfirst(x -> x>=50*maximum(tspandebug)/100, soldebug.t)
	v_bobdebug3 = soldebug3[1,ind_tmindebug3:end]
	v_bobdebug2 = soldebug2[1,ind_tmindebug2:end]
	v_bobdebug = soldebug[1,ind_tmindebug:end]
	minimum(v_bobdebug),minimum(v_bobdebug2)
	md"""Debug simulations"""
end

# ╔═╡ 7c3fb682-0460-4014-934a-be1d2c0f7acc
begin
	plot(soldebug.t,soldebug[1,:],label="I = $(p_global_bif_cycledebug[1])")
	plot!(soldebug2.t,soldebug2[1,:],label="I = $(p_global_bif_cycledebug2[1])")
	plot!(soldebug3.t,soldebug3[1,:],label="I = $(p_global_bif_cycledebug3[1])",legend=:outertopright,size=(400,300))
	yaxis!((-45,-37),"V")
	xaxis!((3000,4000))
end

# ╔═╡ 3b8c3790-aec6-45ca-90b9-5e48758c70e1
begin
	plot(soldebug.t,soldebug[3,:],label="I = $(p_global_bif_cycledebug[1])")
	plot!(soldebug2.t,soldebug2[3,:],label="I = $(p_global_bif_cycledebug2[1])")
	plot!(soldebug3.t,soldebug3[3,:],label="I = $(p_global_bif_cycledebug3[1])",legend=:outertopright,size=(400,300))
	yaxis!("Vus")

end

# ╔═╡ cbff4346-cef2-429e-86b9-d79e92f1389f
begin
	plot(soldebug[2,200:end],soldebug[3,200:end],soldebug[1,200:end])
	plot!(soldebug2[2,200:end],soldebug2[3,200:end],soldebug2[1,200:end],xlabel="Vs",ylabel="Vus",zlabel="V",size=(300,300))
end

# ╔═╡ f6388d7c-b102-4cc8-a7a9-fc9aa24e13bf
begin
	begin
	V1_null_d_(a,b) = Vnullcline1(a,b,p_global_bif_cycledebug2)
	V2_null_d_(a,b) = Vnullcline2(a,b,p_global_bif_cycledebug2)
	Vs_null_d_(a,b) = Vsnullcline(a,b,p_global_bif_cycledebug2)
	Vus_null_d_(a,b) = Vusnullcline(a,b,p_global_bif_cycledebug2)

	plot(vs,vus,V1_null_d_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs", ylabel="Vus",zlabel="V") #yellow-orange
	V_null_plt_d = plot!(vs,vus,V2_null_d_,st=:surface,c=:YlOrBr_3,colorbar_entry=false,label="dV/dt = 0",zlims=(-55,-20)) #yellow-orange
	title!("dV/dt =0")

	plot!(soldebug[2,200:end],soldebug[3,200:end],soldebug[1,200:end])
	plot!(soldebug2[2,200:end],soldebug2[3,200:end],soldebug2[1,200:end],xlabel="Vs",ylabel="Vus",zlabel="V",size=(300,300))
end
end

# ╔═╡ b0bf9322-dba5-4084-99df-f6b77e39273f
begin
	inter_V_Vs_null1_d1_vs = zeros(length(vus))
	inter_V_Vs_null2_d1_vs = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_d1_vs[i] = V_Vs_nullcline1_vs(NaN,vus[i],p_global_bif_cycledebug)
		inter_V_Vs_null2_d1_vs[i] = V_Vs_nullcline2_vs(NaN,vus[i],p_global_bif_cycledebug)
	end

	inter_V_Vus_null1_d1_vs = zeros(size(v))
	inter_V_Vus_null2_d1_vs = zeros(size(v))
	for i=1:length(v)
		inter_V_Vus_null1_d1_vs[i] = V_Vus_nullcline1_vs(v[i],v[i],p_global_bif_cycledebug)
		inter_V_Vus_null2_d1_vs[i] = V_Vus_nullcline2_vs(v[i],v[i],p_global_bif_cycledebug)
	end
	
	

	u0_d1 =[-40.0,-40.0,-40.0]
	tspan_d1=tspandebug
	prob_d1 = ODEProblem(MQIF_3D!,u0_d1,tspan_d1,p_global_bif_cycledebug,callback=cb)
	sol_d1 = solve(prob_d1,dense=false)
	
	proj_d1_v_vs = plot([inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],[inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[inter_V_Vus_null1_d1_vs,inter_V_Vus_null2_d1_vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug[1,200:end],soldebug[2,200:end],line_z=soldebug.t[200:end],linewidth=3,linecolor=:green,label="Trajectory")
	xaxis!((-50,-30))
	yaxis!((-50,-30))
	
	proj_d1_v_vus = plot([inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[v,v],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vus",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug[1,200:end],soldebug[3,200:end],line_z=soldebug.t[200:end],linewidth=3,label="Trajectory")
	xaxis!((-60,-20))
	
	plot(proj_d1_v_vs,proj_d1_v_vus,layout=(1,2),legend=:none)
end

# ╔═╡ 3a3c7729-d8c3-4b4a-92b9-020da890792e
begin
	inter_V_Vs_null1_d2_vs = zeros(length(vus))
	inter_V_Vs_null2_d2_vs = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_d2_vs[i] = V_Vs_nullcline1_vs(NaN,vus[i],p_global_bif_cycledebug2)
		inter_V_Vs_null2_d2_vs[i] = V_Vs_nullcline2_vs(NaN,vus[i],p_global_bif_cycledebug2)
	end

	inter_V_Vus_null1_d2_vs = zeros(size(v))
	inter_V_Vus_null2_d2_vs = zeros(size(v))
	for i=1:length(v)
		inter_V_Vus_null1_d2_vs[i] = V_Vus_nullcline1_vs(v[i],v[i],p_global_bif_cycledebug2)
		inter_V_Vus_null2_d2_vs[i] = V_Vus_nullcline2_vs(v[i],v[i],p_global_bif_cycledebug2)
	end
	
	

	u0_d2 =[-40.0,-40.0,-40.0]
	tspan_d2=tspandebug
	prob_d2 = ODEProblem(MQIF_3D!,u0_d2,tspan_d2,p_global_bif_cycledebug2,callback=cb)
	sol_d2 = solve(prob_d2,dense=false)
	
	proj_d2_v_vs = plot([inter_V_Vs_null1_d2_vs,inter_V_Vs_null2_d2_vs],[inter_V_Vs_null1_d2_vs,inter_V_Vs_null2_d2_vs],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[inter_V_Vus_null1_d2_vs,inter_V_Vus_null2_d1_vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug2[1,200:end],soldebug2[2,200:end],line_z=soldebug2.t[200:end],linewidth=3,linecolor=:green,label="Trajectory")
	xaxis!((-50,-30))
	yaxis!((-50,-30))
	
	proj_d2_v_vus = plot([inter_V_Vs_null1_d2_vs,inter_V_Vs_null2_d2_vs],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[v,v],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vus",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug2[1,200:end],soldebug2[3,200:end],line_z=soldebug2.t[200:end],linewidth=3,label="Trajectory")
	xaxis!((-60,-20))
	
	plot(proj_d2_v_vs,proj_d2_v_vus,layout=(1,2),legend=:none)
end

# ╔═╡ 8726355e-4969-4eb4-b435-511b507bd6bd
md"###### Actual bifurcation diagram"

# ╔═╡ 2c643ebe-95f7-11eb-3174-59c327320e07
begin
	list_V_stable_global_bif=[]
	list_I_stable_global_bif=[]
	list_V_unstable_global_bif=[]
	list_I_unstable_global_bif=[]
	list_V_saddle_global_bif=[]
	list_I_saddle_global_bif=[]

	for i=1:length(I_global_bif)
		if stab_fp_global_bif[i,1] <0
			append!(list_V_stable_global_bif,V_fp_global_bif[i,1])
			append!(list_I_stable_global_bif,I_global_bif[i])
		else
			if stab_fp_global_bif[i,1] >0
				append!(list_V_unstable_global_bif,V_fp_global_bif[i,1])
				append!(list_I_unstable_global_bif,I_global_bif[i])
			else
				append!(list_V_saddle_global_bif,V_fp_global_bif[i,1])
				append!(list_I_saddle_global_bif,I_global_bif[i])
			end
		end
		if stab_fp_global_bif[i,2] <0
			append!(list_V_stable_global_bif,V_fp_global_bif[i,2])
			append!(list_I_stable_global_bif,I_global_bif[i])
		else
			if stab_fp_global_bif[i,2] >0
				append!(list_V_unstable_global_bif,V_fp_global_bif[i,2])
				append!(list_I_unstable_global_bif,I_global_bif[i])
			else
				append!(list_V_saddle_global_bif,V_fp_global_bif[i,2])
				append!(list_I_saddle_global_bif,I_global_bif[i])
			end
		end
	end
	md"""Separation of fixed points value based on their stability"""
end

# ╔═╡ ee87bfb5-68b0-40e7-8b8a-db2c2a4130a0
begin
	list_V_stable_global_bif_r=[]
	list_I_stable_global_bif_r=[]
	list_V_unstable_global_bif_r=[]
	list_I_unstable_global_bif_r=[]
	list_V_saddle_global_bif_r=[]
	list_I_saddle_global_bif_r=[]

	for i=1:length(I_global_bif)
		if stab_fp_global_bif_r[i,1] <0
			append!(list_V_stable_global_bif_r,V_fp_global_bif_r[i,1])
			append!(list_I_stable_global_bif_r,I_global_bif[i])
		else
			if stab_fp_global_bif_r[i,1] >0
				append!(list_V_unstable_global_bif_r,V_fp_global_bif_r[i,1])
				append!(list_I_unstable_global_bif_r,I_global_bif_r[i])
			else
				append!(list_V_saddle_global_bif_r,V_fp_global_bif_r[i,1])
				append!(list_I_saddle_global_bif_r,I_global_bif[i])
			end
		end
		if stab_fp_global_bif_r[i,2] <0
			append!(list_V_stable_global_bif_r,V_fp_global_bif_r[i,2])
			append!(list_I_stable_global_bif_r,I_global_bif[i])
		else
			if stab_fp_global_bif_r[i,2] >0
				append!(list_V_unstable_global_bif_r,V_fp_global_bif_r[i,2])
				append!(list_I_unstable_global_bif_r,I_global_bif[i])
			else
				append!(list_V_saddle_global_bif_r,V_fp_global_bif_r[i,2])
				append!(list_I_saddle_global_bif_r,I_global_bif[i])
			end
		end
	end
	md"""Separation of fixed points value based on their stability"""
end

# ╔═╡ b3318e49-5cee-4b9e-82b2-d9daacf0f480
plotly()

# ╔═╡ 20175480-95f8-11eb-025f-ed13fb0faec1
begin
	bifI_resto = plot(list_I_stable_global_bif,list_V_stable_global_bif,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(list_I_saddle_global_bif,list_V_saddle_global_bif,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	#plot!(list_I_unstable_global_bif,list_V_unstable_global_bif,linecolor="red",linewidth=1.5,label="Unstable node")
	
	#plot!(I_global_bif_cycle,V_cycle_global_bif[:,1])
	plot!(I_global_bif_cycle,V_cycle_global_bif[:,2],label="Limit cycle minimum",linecolor=:pink)
	scatter!([Ibif],[Vbif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	#yaxis!((-50,-30))
	xaxis!("I")
	yaxis!("V")
	#savefig("bifdiag3D.pdf")
end

# ╔═╡ afd86ff3-3457-4e1d-964b-b5d20ab65da9
begin
	bifI_regen = plot(list_I_stable_global_bif_r,list_V_stable_global_bif_r,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(list_I_saddle_global_bif_r,list_V_saddle_global_bif_r,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	plot!(list_I_unstable_global_bif_r,list_V_unstable_global_bif_r,linecolor="red",linewidth=1.5,label="Unstable node")
	scatter!([Ibif_r],[Vbif_r],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	#plot!(I_global_bif_cycle,V_cycle_global_bif[:,1])
	plot!(I_global_bif,V_cycle_global_bif_r[:,2],label="Limit cycle minimum")
	#yaxis!((-50,-30))
end

# ╔═╡ 83f85190-9619-11eb-2b29-cb5032c8b81d
gr()

# ╔═╡ b3a6d3be-4a94-40fb-aedf-a0fa3b5a8bf3
begin
	title1_bifI = plot(title = "A.  Restorative feedback : Vs0 < V0", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=15,titlefontcolor=RGB(0,0.4,0.95))
	title2_bifI = plot(title = "B.   Regenerative feedback : Vs0 > V0", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=15,titlefontcolor=RGB(0,0.4,0.95))
	plot(title1_bifI, 
	bifI_resto,
	title2_bifI,	
	bifI_regen,
	layout = @layout([A{0.04h}; B; C{0.04h}; D]),size=(800,800))
	#savefig("MQIF_3D-bif_I_rest_rege.pdf")
	md"""Bifurcation with I subplot"""
end

# ╔═╡ 37486ae0-95ee-11eb-1912-b1e2ea244526
md"""###### Bifurcation from 2D model : taus = 100 ms, regeneartive and restorative feedbacks"""

# ╔═╡ 64cac8b3-1d84-4152-bef7-aa4f7ee3052f
function jacobian_2D(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	J = zeros(2,2)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[2,1] = 1
	J[2,2] = -1

	return J
end

# ╔═╡ 7736932f-d9be-49e3-be14-888b780f3cca
function fp_stability_2D(fp,p)
	I,v0,vs0,C,gf,gs,ts = p

	stability = ones(1,2)*66 #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable
	for i=1:length(fp)
		J = jacobian_2D(fp[i],p)
		lambda1,lambda2 = eigvals(J)
		if real(lambda1)>0 && real(lambda2)>0
			stability[i]=1.0
		end
		if (real(lambda1)>0 && real(lambda2)<0 ) || (real(lambda1)<0 && real(lambda2)>0)
			stability[i]=0.0
		end
		if real(lambda1)<0 && real(lambda2)<0
			stability[i]=-1.0
		end
	end
	return stability
end

# ╔═╡ 35b7b20d-7a32-4bdf-a682-9e2f72eba110
function fixedpoints_2D(p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D
	#with the vsnullcline, we know that vs = v, we then use this equality in the V-
	#nullcline and solve the system using the discriminant method

	a = (gf-gs)
	b = 2*(gs*vs0 - gf*v0)
	c = gf*v0^2 - gs*vs0^2 + I

	delta = b^2 - 4*a*c

	if delta <0
		#no real solution
		vfp1 = NaN
		vfp2 = NaN
		stability = 0.0
	else
		vfp1 = (-b+sqrt(delta))/(2*a)
		vfp2 = (-b-sqrt(delta))/(2*a)
		stability = zeros(1,2) #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable
		fp = zeros(1,2)
		fp[1] = vfp1
		fp[2] = vfp2
		stability = fp_stability_2D(fp,p_2D)
	end

	fp1=ones(1,2).*vfp1
	fp2=ones(1,2).*vfp2

	return fp1,fp2,stability
end

# ╔═╡ 19188423-c2e3-403e-9e16-7ccfeaca4f45
function bifI_matrices_2D(I,p_2D)
	I_,v0,vs0,C,gf,gs,ts = p_2D

	list_stable = zeros(length(I)) #column 1 : V
	list_cycleV = zeros(length(I),2) #column 1 : max, column 2 : min

	list_saddle = zeros(length(I)) #column 1 : V
	list_unstable = zeros(length(I)) #column 1 : V

	for i=1:length(I)
		pI=change_I_2D(I[i],p_2D)
		fp1,fp2,stability= fixedpoints_2D(pI)
		fp=[fp1,fp2]
		if length(stability)>1
			for j=1:length(stability)
				if stability[j]<0
					list_stable[i]=fp[j][1]
					list_unstable[i]=NaN
				else
					if stability[j]==0
						list_saddle[i]=fp[j][1]
					else
						list_unstable[i]=fp[j][1]
						list_stable[i]=NaN
					end
				end
			end
		else
			list_stable[i]=NaN
			list_unstable[i]=NaN
			list_saddle[i]=NaN
		end

		u0 =[v0+20,vs0]
		tspan = (0.0,200.0)
		prob = ODEProblem(MQIF_2D!,u0,tspan,pI,callback=cb_2D)
		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=3*maximum(tspan)/10, sol.t)
		list_cycleV[i,1] = maximum(sol[1,ind_tmin:length(sol)])
		list_cycleV[i,2] = minimum(sol[1,ind_tmin:length(sol)])

		if list_cycleV[i,1]-list_cycleV[i,2]<= 30
			list_cycleV[i,:]=[NaN,NaN]
		end
	end

	return list_stable,list_saddle,list_unstable,list_cycleV
end

# ╔═╡ 71b19713-3e2c-463d-a44a-fcd7c46a12a1
function p_3D_to_2D(p3D)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p3D
	p2D = zeros(1,7)
	p2D[1]=I
	p2D[2]=v0
	p2D[3]=vs0
	p2D[4]=C
	p2D[5]=gf
	p2D[6]=gs
	p2D[7]=ts
	return p2D
end

# ╔═╡ 05107ac0-69f7-493a-a78a-64f84167f547
begin
	inter_V_Vs_null1_d3_vs = zeros(length(vus))
	inter_V_Vs_null2_d3_vs = zeros(length(vus))
	for i=1:length(vus)
		inter_V_Vs_null1_d3_vs[i] = V_Vs_nullcline1_vs(NaN,vus[i],p_global_bif_cycledebug3)
		inter_V_Vs_null2_d3_vs[i] = V_Vs_nullcline2_vs(NaN,vus[i],p_global_bif_cycledebug3)
	end

	inter_V_Vus_null1_d3_vs = zeros(size(v))
	inter_V_Vus_null2_d3_vs = zeros(size(v))
	bob=[]
	for i=1:length(v)
		inter_V_Vus_null1_d3_vs[i] = V_Vus_nullcline1_vs(v[i],v[i],p_global_bif_cycledebug3)
		inter_V_Vus_null2_d3_vs[i] = V_Vus_nullcline2_vs(v[i],v[i],p_global_bif_cycledebug3)
		append!(bob,Vnullcline1_2D(v[i],p_3D_to_2D(p_global_bif_cycledebug3)))
	end
	
	
	

	u0_d3 =[-40.0,-40.0,-40.0]
	tspan_d3=tspandebug
	prob_d3 = ODEProblem(MQIF_3D!,u0_d3,tspan_d3,p_global_bif_cycledebug3,callback=cb)
	sol_d3 = solve(prob_d3,dense=false)
	
	proj_d3_v_vs = plot([inter_V_Vs_null1_d3_vs,inter_V_Vs_null2_d3_vs],[inter_V_Vs_null1_d3_vs,inter_V_Vs_null2_d3_vs],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[inter_V_Vus_null1_d3_vs,inter_V_Vus_null2_d3_vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug3[1,400:end],soldebug3[2,400:end],line_z=soldebug3.t[400:end],linewidth=3,linecolor=:green,label="Trajectory")
	xaxis!((-50,-30))
	yaxis!((-50,-30))
	
	proj_d3_v_vus = plot([inter_V_Vs_null1_d3_vs,inter_V_Vs_null2_d3_vs],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([v,v],[v,v],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vus",label="dV/dt = dVus/dt = 0") #blue
	plot!(soldebug3[1,400:end],soldebug3[3,400:end],line_z=soldebug3.t[400:end],linewidth=3,label="Trajectory")
	xaxis!((-60,-20))
	
	plot(proj_d3_v_vs,proj_d3_v_vus,layout=(1,2),legend=:none)
end

# ╔═╡ 5dfa49e0-b3ec-4563-ba63-dc945cdff46d
begin
	 proj_d_v_vs_ = plot([inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],[inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0")#pink
	plot!([v,v],[inter_V_Vus_null1_d1_vs,inter_V_Vus_null2_d1_vs],linecolor=RGB(0.31,0.8,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #light blue
	plot!([v,v],[inter_V_Vus_null1_d2_vs,inter_V_Vus_null2_d2_vs],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #blue
	plot!([v,v],[inter_V_Vus_null1_d3_vs,inter_V_Vus_null2_d3_vs],linecolor=RGB(0.31,0.2,1),linewidth=3,xlabel="V",ylabel="Vs",label="dV/dt = dVus/dt = 0") #dark blue
	xaxis!((p_global_bif_cycledebug[2]-5,p_global_bif_cycledebug[2]+5))
	yaxis!((p_global_bif_cycledebug[3]-5,p_global_bif_cycledebug[3]+5))
	
	proj_d_v_vus_ = plot([inter_V_Vs_null1_d1_vs,inter_V_Vs_null2_d1_vs],[vus,vus],linecolor=RGB(1,0.7,0.7),linewidth=3,label="dV/dt = dVs/dt = 0") #pink
	plot!([inter_V_Vs_null1_d2_vs,inter_V_Vs_null2_d2_vs],[vus,vus],linecolor=RGB(1,0.4,0.7),linewidth=3,label="dV/dt = dVs/dt = 0")
	plot!([inter_V_Vs_null1_d3_vs,inter_V_Vs_null2_d3_vs],[vus,vus],linecolor=RGB(1,0.1,0.5),linewidth=3,label="dV/dt = dVs/dt = 0")
	plot!([v,v],[v,v],linecolor=RGB(0.31,0.66,1),linewidth=3,xlabel="V",ylabel="Vus",label="dV/dt = dVus/dt = 0") #blue
	xaxis!((p_global_bif_cycledebug[2]-5,p_global_bif_cycledebug[2]+5))
	yaxis!((-75,-25))
	
	plot(proj_d_v_vs_,proj_d_v_vus_,layout=(1,2),legend=:none,size=(600,300))
end

# ╔═╡ 243e2330-95ea-11eb-10c8-7707ea0c7402
function local_I_bif_from2D_model(p,vus)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	Ius = -gus*(vus-vus0)^2
	It = I + Ius

	list_cycleV =[0.0,0.0]

	common_term = gf*v0 - gs*vs0
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*It
	if delta >=0
		list_cycleV =[NaN,NaN]

		numerator1 = common_term + sqrt(delta)
		numerator2 = common_term - sqrt(delta)
		denominator = gf-gs
		V_fp_1 = numerator1/denominator
		V_fp_2 = numerator2/denominator
		V_fp = [V_fp_1,V_fp_2]

		convergence=V_fp
		conv_type=0

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		stability = fp_stability_2D(V_fp,p_2D)
	else
		V_fp_1 = NaN
		V_fp_2 = NaN
		V_fp = [V_fp_1,V_fp_2]

		#must simulate 2D model to get the min and max cycle value
		tspan_2D=(0.0,150.0)
		u0_2D =[v0+20,vs0]

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		prob_2D = ODEProblem(MQIF_2D!,u0_2D,tspan_2D,p_2D,callback=cb_2D)
		sol_2D = solve(prob_2D,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=maximum(tspan_2D)/10, sol_2D.t)
		list_cycleV[1] = maximum(sol_2D[1,ind_tmin:length(sol_2D)])
		list_cycleV[2] = minimum(sol_2D[1,ind_tmin:length(sol_2D)])

		convergence=list_cycleV
		conv_type=1

		stability = [NaN,NaN]
	end
	return convergence,conv_type,stability,It
end

# ╔═╡ 87273901-f5b8-4db1-8810-28c7f5e8d60a
function local_I_bif_from2D_model_decr(p,vus,vs,v)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	Ius = -gus*(vus-vus0)^2
	It = I + Ius

	list_cycleV =[0.0,0.0]

	common_term = gf*v0 - gs*vs0
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*It
	if delta >=0
		list_cycleV =[NaN,NaN]

		numerator1 = common_term + sqrt(delta)
		numerator2 = common_term - sqrt(delta)
		denominator = gf-gs
		V_fp_1 = numerator1/denominator
		V_fp_2 = numerator2/denominator
		V_fp = [V_fp_1,V_fp_2]

		convergence=V_fp
		conv_type=0

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		stability = fp_stability_2D(V_fp,p_2D)
	else
		V_fp_1 = NaN
		V_fp_2 = NaN
		V_fp = [V_fp_1,V_fp_2]

		#must simulate 2D model to get the min and max cycle value
		tspan_2D=(0.0,200.0)
		u0_2D =[v0+20,vs0]

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		prob_2D = ODEProblem(MQIF_2D!,u0_2D,tspan_2D,p_2D,callback=cb_2D)
		sol_2D = solve(prob_2D,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=maximum(tspan_2D)/10, sol_2D.t)

		Vs_max_low_V_null_2D=max_2D_low_V_null(It)
		list_cycleV[1] = maximum(sol_2D[1,ind_tmin:length(sol_2D)])
		list_cycleV[2] = minimum(sol_2D[1,ind_tmin:length(sol_2D)])

		if vs>Vs_max_low_V_null_2D
			convergence=list_cycleV
			conv_type=1
		else
			spike_latency = [vs,list_cycleV[2]]
			convergence = spike_latency
			conv_type=2
		end


		stability = [NaN,NaN]
	end
	return convergence,conv_type,stability,It
end

# ╔═╡ d5a7342d-d27d-4fef-b2ff-3a768a4a5b7c
function local_I_bif_from2D_model_incr(p,vus,vs,v)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p

	Ius = -gus*(vus-vus0)^2
	It = I + Ius

	list_cycleV =[0.0,0.0]

	common_term = gf*v0 - gs*vs0
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*It
	if delta >=0
		list_cycleV =[NaN,NaN]

		numerator1 = common_term + sqrt(delta)
		numerator2 = common_term - sqrt(delta)
		denominator = gf-gs
		V_fp_1 = numerator1/denominator
		V_fp_2 = numerator2/denominator
		V_fp = [V_fp_1,V_fp_2]

		tspan_2D=(0.0,200.0)
		u0_2D =[v0,vs0+10] #above upper nullcline

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		prob_2D = ODEProblem(MQIF_2D!,u0_2D,tspan_2D,p_2D,callback=cb_2D)
		sol_2D = solve(prob_2D,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=2*maximum(tspan_2D)/10, sol_2D.t)

		list_cycleV[1] = maximum(sol_2D[1,ind_tmin:length(sol_2D)])
		list_cycleV[2] = minimum(sol_2D[1,ind_tmin:length(sol_2D)])

		if list_cycleV[1]-list_cycleV[2]>=50
			convergence = list_cycleV
			conv_type=1
			stability = [NaN,NaN]
		else
			convergence=V_fp
			conv_type=0
			stability = fp_stability_2D(V_fp,p_2D)
		end

	else
		V_fp_1 = NaN
		V_fp_2 = NaN
		V_fp = [V_fp_1,V_fp_2]

		#must simulate 2D model to get the min and max cycle value
		tspan_2D=(0.0,200.0)
		u0_2D =[v0,vs0+10] #above upper nullcline

		p_It = change_I(It)
		p_2D=p_3D_to_2D(p_It)

		prob_2D = ODEProblem(MQIF_2D!,u0_2D,tspan_2D,p_2D,callback=cb_2D)
		sol_2D = solve(prob_2D,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=maximum(tspan_2D)/10, sol_2D.t)

		Vs_max_low_V_null_2D=max_2D_low_V_null(It)
		list_cycleV[1] = maximum(sol_2D[1,ind_tmin:length(sol_2D)])
		list_cycleV[2] = minimum(sol_2D[1,ind_tmin:length(sol_2D)])

		if vs>Vs_max_low_V_null_2D
			convergence=list_cycleV
			conv_type=1
		else
			spike_latency = [vs,list_cycleV[2]]
			convergence = spike_latency
			conv_type=2
		end

		stability = [NaN,NaN]
	end
	return convergence,conv_type,stability,It
end

# ╔═╡ 09266a70-9c17-4f5f-a1e9-8917bbbcc0de
function separate_local_bif_conv(conv_local_bif,conv_type_local_bif,It_local_bif,stab_fp_local_bif)
	#Separation of convergence value based on their type and stability
	V_stable_local_bif=[]
	I_stable_local_bif=[]
	V_unstable_local_bif=[]
	I_unstable_local_bif=[]
	V_saddle_local_bif=[]
	I_saddle_local_bif=[]
	V_cyclemax_local_bif=[]
	V_cyclemin_local_bif=[]
	I_cycle_local_bif=[]

	for i=1:length(It_local_bif)
		if conv_type_local_bif[i]==0
			if stab_fp_local_bif[i,1] <0
				append!(V_stable_local_bif,conv_local_bif[i,1])
				append!(I_stable_local_bif,It_local_bif[i])
			else
				if stab_fp_local_bif[i,1] >0
					append!(V_unstable_local_bif,conv_local_bif[i,1])
					append!(I_unstable_local_bif,It_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,1])
					append!(I_saddle_local_bif,It_local_bif[i])
				end
			end
			if stab_fp_local_bif[i,2] <0
				append!(V_stable_local_bif,conv_local_bif[i,2])
				append!(I_stable_local_bif,It_local_bif[i])
			else
				if stab_fp_local_bif[i,2] >0
					append!(V_unstable_local_bif,conv_local_bif[i,2])
					append!(I_unstable_local_bif,It_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,2])
					append!(I_saddle_local_bif,It_local_bif[i])
				end
			end
		else
			append!(V_cyclemax_local_bif,conv_local_bif[i,1])
			append!(V_cyclemin_local_bif,conv_local_bif[i,2])
			append!(I_cycle_local_bif,It_local_bif[i])
		end
	end
	return V_stable_local_bif,I_stable_local_bif,V_unstable_local_bif,I_unstable_local_bif,V_saddle_local_bif,I_saddle_local_bif,V_cyclemax_local_bif,V_cyclemin_local_bif,I_cycle_local_bif
end

# ╔═╡ 835e25e2-a335-40e5-af56-d75f0db67c8a
function separate_local_bif_conv_full(conv_local_bif,conv_type_local_bif,It_local_bif,t_local_bif,stab_fp_local_bif)
	#Separation of convergence value based on their type and stability
	V_stable_local_bif=[]
	I_stable_local_bif=[]
	t_stable_local_bif=[]
	V_unstable_local_bif=[]
	I_unstable_local_bif=[]
	t_unstable_local_bif=[]
	V_saddle_local_bif=[]
	I_saddle_local_bif=[]
	t_saddle_local_bif=[]
	V_cyclemax_local_bif=[]
	V_cyclemin_local_bif=[]
	I_cycle_local_bif=[]
	t_cycle_local_bif=[]
	V_spikelat_local_bif=[]
	V_sl_cmin_local_bif=[]
	I_spikelat_local_bif=[]
	t_spikelat_local_bif=[]

	for i=1:length(It_local_bif)
		if conv_type_local_bif[i]==0
			if stab_fp_local_bif[i,1] <0
				append!(V_stable_local_bif,conv_local_bif[i,1])
				append!(I_stable_local_bif,It_local_bif[i])
				append!(t_stable_local_bif,t_local_bif[i])
			else
				if stab_fp_local_bif[i,1] >0
					append!(V_unstable_local_bif,conv_local_bif[i,1])
					append!(I_unstable_local_bif,It_local_bif[i])
					append!(t_unstable_local_bif,t_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,1])
					append!(I_saddle_local_bif,It_local_bif[i])
					append!(t_saddle_local_bif,t_local_bif[i])
				end
			end
			if stab_fp_local_bif[i,2] <0
				append!(V_stable_local_bif,conv_local_bif[i,2])
				append!(I_stable_local_bif,It_local_bif[i])
				append!(t_stable_local_bif,t_local_bif[i])
			else
				if stab_fp_local_bif[i,2] >0
					append!(V_unstable_local_bif,conv_local_bif[i,2])
					append!(I_unstable_local_bif,It_local_bif[i])
					append!(t_unstable_local_bif,t_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,2])
					append!(I_saddle_local_bif,It_local_bif[i])
					append!(t_saddle_local_bif,t_local_bif[i])
				end
			end
		else
			if conv_type_local_bif[i]==1
				append!(V_cyclemax_local_bif,conv_local_bif[i,1])
				append!(V_cyclemin_local_bif,conv_local_bif[i,2])
				append!(I_cycle_local_bif,It_local_bif[i])
				append!(t_cycle_local_bif,t_local_bif[i])
			else
				append!(V_spikelat_local_bif,conv_local_bif[i,1])
				append!(V_sl_cmin_local_bif,conv_local_bif[i,2])
				append!(I_spikelat_local_bif,It_local_bif[i])
				append!(t_spikelat_local_bif,t_local_bif[i])
			end
		end
	end
	return V_stable_local_bif,I_stable_local_bif,t_stable_local_bif,V_unstable_local_bif,I_unstable_local_bif,t_unstable_local_bif,V_saddle_local_bif,I_saddle_local_bif,t_saddle_local_bif,V_cyclemax_local_bif,V_cyclemin_local_bif,I_cycle_local_bif,t_cycle_local_bif,V_spikelat_local_bif,V_sl_cmin_local_bif,I_spikelat_local_bif,t_spikelat_local_bif
end

# ╔═╡ 7f469400-632a-40c6-96c7-2b674324ed60
function separate_local_bif_conv_time(conv_local_bif,conv_type_local_bif,t_local_bif,stab_fp_local_bif)
	#Separation of convergence value based on their type and stability
	V_stable_local_bif=[]
	t_stable_local_bif=[]
	V_unstable_local_bif=[]
	t_unstable_local_bif=[]
	V_saddle_local_bif=[]
	t_saddle_local_bif=[]
	V_cyclemax_local_bif=[]
	V_cyclemin_local_bif=[]
	t_cycle_local_bif=[]

	for i=1:length(t_local_bif)
		if conv_type_local_bif[i]==0
			if stab_fp_local_bif[i,1] <0
				append!(V_stable_local_bif,conv_local_bif[i,1])
				append!(t_stable_local_bif,t_local_bif[i])
			else
				if stab_fp_local_bif[i,1] >0
					append!(V_unstable_local_bif,conv_local_bif[i,1])
					append!(t_unstable_local_bif,t_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,1])
					append!(t_saddle_local_bif,t_local_bif[i])
				end
			end
			if stab_fp_local_bif[i,2] <0
				append!(V_stable_local_bif,conv_local_bif[i,2])
				append!(t_stable_local_bif,t_local_bif[i])
			else
				if stab_fp_local_bif[i,2] >0
					append!(V_unstable_local_bif,conv_local_bif[i,2])
					append!(t_unstable_local_bif,t_local_bif[i])
				else
					append!(V_saddle_local_bif,conv_local_bif[i,2])
					append!(t_saddle_local_bif,t_local_bif[i])
				end
			end
		else
			append!(V_cyclemax_local_bif,conv_local_bif[i,1])
			append!(V_cyclemin_local_bif,conv_local_bif[i,2])
			append!(t_cycle_local_bif,t_local_bif[i])
		end
	end
	return V_stable_local_bif,t_stable_local_bif,V_unstable_local_bif,t_unstable_local_bif,V_saddle_local_bif,t_saddle_local_bif,V_cyclemax_local_bif,V_cyclemin_local_bif,t_cycle_local_bif
end

# ╔═╡ 9f0566ed-09a1-4c4f-ac7c-dbc1b2452f68
function compute_I_from_vus(vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	
	Ius = gus*(vus-vus0)^2
	
	It = I - Ius
	
	return It
end

# ╔═╡ 40097d95-9548-4291-8101-a74c6f405440
begin
	Itdebug = zeros(size(soldebug[3,:]))
	Itdebug2 = zeros(size(soldebug2[3,:]))	
	Itdebug3 = zeros(size(soldebug3[3,:]))
				
	for i=1:length(Itdebug)
		Itdebug[i] = compute_I_from_vus(soldebug[3,i],p_global_bif_cycledebug)		
	end
	for i=1:length(Itdebug2)
		Itdebug2[i] = compute_I_from_vus(soldebug2[3,i],p_global_bif_cycledebug2)	
	end
	for i=1:length(Itdebug3)
		Itdebug3[i] = compute_I_from_vus(soldebug3[3,i],p_global_bif_cycledebug3)	
	end
	
	plot(soldebug.t,Itdebug,label="I = $(p_global_bif_cycledebug[1])")
	plot!(soldebug2.t,Itdebug2,label="I = $(p_global_bif_cycledebug2[1])")
	plot!(soldebug3.t,Itdebug3,label="I = $(p_global_bif_cycledebug3[1])",size=(400,200),legend=:outertopright)
	yaxis!((-7,7),"It")
	xaxis!((3000,4000))
end

# ╔═╡ 9c5553f3-751b-4b28-abd1-c06ca67a80aa
begin
	##Simulation of 3D model to get Vus(t) for these parameters
	p_local_bif=(4,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,100.0)

	u0_local_bif =[-40.0,-40.0,-40.0]
	tspan_local_bif=(0.0,1200.0)
	prob_local_bif = ODEProblem(MQIF_3D!,u0_local_bif,tspan_local_bif,p_local_bif,callback=cb)
	sol_local_bif = solve(prob_local_bif,dense=false,dtmax=0.02,reltol=1e-6,abstol=1e-6)
	##Must retain only one cycle of Vus, chosen to be between the 2 last spikes observed in the simulation
	ind_tmin = findfirst(x -> x>=20*maximum(tspan_local_bif)/100, sol_local_bif.t)
	ind_spike_ = findall(x -> x>=maximum(sol_local_bif[3,ind_tmin:end])-0.01,sol_local_bif[3,ind_tmin:end])
	ind_spike=[]
	for i=1:length(ind_spike_)
		if i>1
			if ind_spike_[i]-ind_spike_[i-1]>10
				append!(ind_spike,ind_spike_[i-1])
			end
		end
	end
	append!(ind_spike,ind_spike_[end]) #add last spike

	vus_local_bif = sol_local_bif[3,(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin]
	vs_local_bif = sol_local_bif[2,(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin]
	v_local_bif = sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin]
	t_local_bif = sol_local_bif.t[(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin]
	
	It_local_bif= zeros(length(t_local_bif))
	for i=1:length(t_local_bif)
		It_local_bif[i] = compute_I_from_vus(vus_local_bif[i],p_local_bif)
	end

	##Separation of the Vus cycle into 2 parts : Vus is decreasing and Vus is increasing
	ind_min = findfirst(x -> x == minimum(vus_local_bif),vus_local_bif)

	vus_local_bif_decr = vus_local_bif[1:ind_min]
	vs_local_bif_decr = vs_local_bif[1:ind_min]
	v_local_bif_decr = v_local_bif[1:ind_min]
	t_local_bif_decr = t_local_bif[1:ind_min]

	vus_local_bif_incr = vus_local_bif[ind_min:end]
	vs_local_bif_incr = vs_local_bif[ind_min:end]
	v_local_bif_incr = v_local_bif[ind_min:end]
	t_local_bif_incr = t_local_bif[ind_min:end]

	##Calculations to determine the 2D system convergence for It = I +Ius, depending on the value of Vus
	conv_local_bif_decr = zeros(length(vus_local_bif_decr),2)
	conv_local_bif_incr = zeros(length(vus_local_bif_incr),2)
	stab_fp_local_bif_decr = zeros(length(vus_local_bif_decr),2)
	stab_fp_local_bif_incr = zeros(length(vus_local_bif_incr),2)
	conv_type_local_bif_decr = zeros(length(vus_local_bif_decr))
	conv_type_local_bif_incr = zeros(length(vus_local_bif_incr))
	It_local_bif_decr = zeros(length(vus_local_bif_decr))
	It_local_bif_incr = zeros(length(vus_local_bif_incr))

	for i=1:length(vus_local_bif_decr)
		#conv_local_bif_decr[i,:],conv_type_local_bif_decr[i],stab_fp_local_bif_decr[i,:],It_local_bif_decr[i]=local_I_bif_from2D_model(p_local_bif,vus_local_bif_decr[i])
		conv_local_bif_decr[i,:],conv_type_local_bif_decr[i],stab_fp_local_bif_decr[i,:],It_local_bif_decr[i]=local_I_bif_from2D_model_decr(p_local_bif,vus_local_bif_decr[i],vs_local_bif_decr[i],v_local_bif_decr[i])
	end
	for i=1:length(vus_local_bif_incr)
		#conv_local_bif_incr[i,:],conv_type_local_bif_incr[i],stab_fp_local_bif_incr[i,:],It_local_bif_incr[i]=local_I_bif_from2D_model(p_local_bif,vus_local_bif_incr[i])

conv_local_bif_incr[i,:],conv_type_local_bif_incr[i],stab_fp_local_bif_incr[i,:],It_local_bif_incr[i]=local_I_bif_from2D_model_incr(p_local_bif,vus_local_bif_incr[i],vs_local_bif_incr[i],v_local_bif_incr[i])
		
		md"""Main simulatino, split for a period study, computation of the equivalent state in the 2D model for each value of the total current"""
	end
end

# ╔═╡ 17656593-6727-447e-97b0-d7f67faae4e2
begin
	##Simulation of 3D model to get Vus(t) for these parameters
	p_local_bif_hI=(5,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,100.0)

	prob_local_bif_hI = ODEProblem(MQIF_3D!,u0_local_bif,tspan_local_bif,p_local_bif_hI,callback=cb)
	sol_local_bif_hI = solve(prob_local_bif_hI,dense=false,dtmax=0.1,reltol=1e-6,abstol=1e-6)
	##Must retain only one cycle of Vus, chosen to be between the 2 last spikes observed in the simulation
	ind_tmin_hI = findfirst(x -> x>=20*maximum(tspan_local_bif)/100, sol_local_bif_hI.t)
	ind_spike_hI_ = findall(x -> x>=maximum(sol_local_bif_hI[3,ind_tmin_hI:end])-0.01,sol_local_bif_hI[3,ind_tmin_hI:end])
	ind_spike_hI=[]
	for i=1:length(ind_spike_hI_)
		if i>1
			if ind_spike_hI_[i]-ind_spike_hI_[i-1]>10
				append!(ind_spike_hI,ind_spike_hI_[i-1])
			end
		end
	end
	append!(ind_spike_hI,ind_spike_hI_[end]) #add last spike

	vus_local_bif_hI = sol_local_bif_hI[3,(ind_spike_hI[end-1]-1)+ind_tmin_hI:(ind_spike_hI[end]-1)+ind_tmin_hI]
	vs_local_bif_hI = sol_local_bif_hI[2,(ind_spike_hI[end-1]-1)+ind_tmin_hI:(ind_spike_hI[end]-1)+ind_tmin_hI]
	v_local_bif_hI = sol_local_bif_hI[1,(ind_spike_hI[end-1]-1)+ind_tmin_hI:(ind_spike_hI[end]-1)+ind_tmin_hI]
	t_local_bif_hI = sol_local_bif_hI.t[(ind_spike_hI[end-1]-1)+ind_tmin_hI:(ind_spike_hI[end]-1)+ind_tmin_hI]
	
	It_local_bif_hI= zeros(length(t_local_bif_hI))
	for i=1:length(t_local_bif_hI)
		It_local_bif_hI[i] = compute_I_from_vus(vus_local_bif_hI[i],p_local_bif_hI)
	end

	##Separation of the Vus cycle into 2 parts : Vus is decreasing and Vus is increasing
	ind_min_hI = findfirst(x -> x == minimum(vus_local_bif_hI),vus_local_bif_hI)

	vus_local_bif_decr_hI = vus_local_bif_hI[1:ind_min_hI]
	vs_local_bif_decr_hI = vs_local_bif_hI[1:ind_min_hI]
	v_local_bif_decr_hI = v_local_bif_hI[1:ind_min_hI]
	t_local_bif_decr_hI = t_local_bif_hI[1:ind_min_hI]
	vus_local_bif_incr_hI = vus_local_bif_hI[ind_min_hI:end]
	vs_local_bif_incr_hI = vs_local_bif_hI[ind_min_hI:end]
	v_local_bif_incr_hI = v_local_bif_hI[ind_min_hI:end]
	t_local_bif_incr_hI = t_local_bif_hI[ind_min_hI:end]

	##Calculations to determine the 2D system convergence for It = I +Ius, depending on the value of Vus
	conv_local_bif_decr_hI = zeros(length(vus_local_bif_decr_hI),2)
	conv_local_bif_incr_hI = zeros(length(vus_local_bif_incr_hI),2)
	stab_fp_local_bif_decr_hI = zeros(length(vus_local_bif_decr_hI),2)
	stab_fp_local_bif_incr_hI = zeros(length(vus_local_bif_incr_hI),2)
	conv_type_local_bif_decr_hI = zeros(length(vus_local_bif_decr_hI))
	conv_type_local_bif_incr_hI = zeros(length(vus_local_bif_incr_hI))
	It_local_bif_decr_hI = zeros(length(vus_local_bif_decr_hI))
	It_local_bif_incr_hI = zeros(length(vus_local_bif_incr_hI))

	for i=1:length(vus_local_bif_decr_hI)
		#conv_local_bif_decr_hI[i,:],conv_type_local_bif_decr_hI[i],stab_fp_local_bif_decr_hI[i,:],It_local_bif_decr_hI[i]=local_I_bif_from2D_model(p_local_bif_hI,vus_local_bif_decr_hI[i])
		conv_local_bif_decr_hI[i,:],conv_type_local_bif_decr_hI[i],stab_fp_local_bif_decr_hI[i,:],It_local_bif_decr_hI[i]=local_I_bif_from2D_model_decr(p_local_bif_hI,vus_local_bif_decr_hI[i],vs_local_bif_decr_hI[i],v_local_bif_decr_hI[i])

	end
	for i=1:length(vus_local_bif_incr_hI)
		#conv_local_bif_incr_hI[i,:],conv_type_local_bif_incr_hI[i],stab_fp_local_bif_incr_hI[i,:],It_local_bif_incr_hI[i]=local_I_bif_from2D_model(p_local_bif_hI,vus_local_bif_incr_hI[i])
		conv_local_bif_incr_hI[i,:],conv_type_local_bif_incr_hI[i],stab_fp_local_bif_incr_hI[i,:],It_local_bif_incr_hI[i]=local_I_bif_from2D_model_incr(p_local_bif_hI,vus_local_bif_incr_hI[i],vs_local_bif_incr_hI[i],v_local_bif_incr_hI[i])

	end
	md"""Main simulatino, split for a period study, computation of the equivalent state in the 2D model for each value of the total current, for a higher current"""
end

# ╔═╡ 86854a97-660b-4b10-a769-341a17ac71ae
gr()

# ╔═╡ 7d998b8e-6eae-43f6-8e67-684d6b6df2f8
begin
	Ibif_2D_local = bifI2D(p_3D_to_2D(p_local_bif))
	ind_Ibif_2D_local = findfirst(x -> x>= Ibif_2D_local,It_local_bif_decr)
	ind_nbif_2D_local = findfirst(x -> x<= 0,It_local_bif_incr)
	md"""Bifurcations indices in ROI"""
end

# ╔═╡ 378a31a7-4c78-4f56-b225-73fa47921fa4
begin
	It_local_bif_full = zeros(length(sol_local_bif[3,:]))
	for i=1:length(sol_local_bif[3,:])
		It_local_bif_full[i] = compute_I_from_vus(sol_local_bif[3,i],p_local_bif)
	end
	
	plot_It_time_full = plot(sol_local_bif.t,It_local_bif_full,label="It",linewidth=1.5,linecolor=RGB(0.65,0.29,0.65))	
	plot!(sol_local_bif.t,Ibif_2D_local.*ones(size(It_local_bif_full)),label="2D SN bifurcation",linewidth=1.5,linecolor=RGB(0.9,0.9,0),legend=:topleft)
	plot!(sol_local_bif.t,zeros(size(It_local_bif_full)),label="Bifurcation to rest",linewidth=1.5,linecolor=RGB(0,0.5,1))
	
	plot!(sol_local_bif.t,(4).*ones(size(sol_local_bif.t)),fill = ([Ibif_2D_local.*ones(size(sol_local_bif.t))], 0.1, RGB(0,0.5,1)),linecolor=:purple,linealpha=0.1,label="Limit cycle")
	plot!(sol_local_bif.t,Ibif_2D_local.*ones(size(sol_local_bif.t)),fill = ([zeros(size(It_local_bif_full))], 0.1, RGB(1,1,0)),linecolor=:purple,linealpha=0.1,label="Bistable")
	plot!(sol_local_bif.t,zeros(size(It_local_bif_full)),fill = ([(-1.5).*ones(size(It_local_bif_full))], 0.1, RGB(0,1,0)),linecolor=:purple,linealpha=0.1,label="Stable")	
	
	plot!(t_local_bif,(maximum(It_local_bif_full)+1).*ones(size(t_local_bif)),fill = ((-1.5).*ones(size(t_local_bif)), 0.6, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest",size=(250,250))
	
	
	yaxis!("Current",(-1.5,3.5))
	xaxis!("Time (ms)")
	#savefig("simu3D_local_bif_fullIt.pdf")
end

# ╔═╡ 3ba6624a-2219-498b-a435-a9d49fc6db14
begin
	It_local_bif_sl_test = zeros(length(It_local_bif_decr[ind_Ibif_2D_local:end])+length(It_local_bif_incr[1:ind_nbif_2D_local-1]))
	t_local_bif_sl_test = zeros(length(It_local_bif_sl_test))
	Vs_max_low_V_null_2D = zeros(length(It_local_bif_sl_test))

	for i=1:length(It_local_bif_sl_test)
		if i<=length(It_local_bif_decr)-(ind_Ibif_2D_local-1)
			It_local_bif_sl_test[i] = It_local_bif_decr[i+(ind_Ibif_2D_local-1)]
			t_local_bif_sl_test[i] = t_local_bif_decr[i+(ind_Ibif_2D_local-1)]
		else
			It_local_bif_sl_test[i] = It_local_bif_incr[i-(length(It_local_bif_decr)-(ind_Ibif_2D_local-1))]
			t_local_bif_sl_test[i] = t_local_bif_incr[i-(length(It_local_bif_decr)-(ind_Ibif_2D_local-1))]
		end
		Vs_max_low_V_null_2D[i]=max_2D_low_V_null(It_local_bif_sl_test[i])
	end
	md"""Computation of the V nullcline maximum to determine when spike latency ends"""
end

# ╔═╡ 465e64e3-ccc0-4951-b3c9-37710919e8a7
begin
	ind_sl_ends = [-1]
	for i=0:length(Vs_max_low_V_null_2D)-1
		if Vs_max_low_V_null_2D[i+1]<=sol_local_bif[2,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local+i] && ind_sl_ends[1] <0
			ind_sl_ends[1] = i+1
		end
	end

	plot_sl_local_bif_It = plot(t_local_bif_sl_test,It_local_bif_sl_test,linecolor=:pink,label="It",linewidth=1)
	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[Ibif_2D_local],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:bottomleft)
	scatter!([t_local_bif_sl_test[ind_sl_ends[1]]],[It_local_bif_sl_test[ind_sl_ends[1]]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends")
	xaxis!("Time (ms)")
	yaxis!("It")

	plot_sl_local_bif_v = plot(sol_local_bif.t[(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin],sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin],label="V",linewidth=1,linecolor =:blue)
	plot!(sol_local_bif.t[(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin],p[2].*ones(size(sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin])),label="V0",linewidth=1,linecolor=RGB(0,0.5,1))
	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:bottomleft)
	scatter!([t_local_bif_sl_test[ind_sl_ends]],[sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local+ind_sl_ends[1]-1]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)
	yaxis!((-45,-35),"V")
	xaxis!("Time (ms)")

	plot_sl_local_bif_vs  = plot(sol_local_bif.t[(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin],sol_local_bif[2,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local:(ind_spike[end]-1)+ind_tmin],label="Vs",linewidth=1,linecolor=:orange)
	plot!(t_local_bif_sl_test,Vs_max_low_V_null_2D,label="Vs of bottom V nullcline maximum",linecolor=RGB(1,0.4,0))
	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[sol_local_bif[2,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:bottomleft)
	scatter!([t_local_bif_sl_test[ind_sl_ends]],[Vs_max_low_V_null_2D[ind_sl_ends]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)
	yaxis!((-45,-35),"Vs")
	xaxis!("Time (ms)")

	md"""Individual spike latency plots"""
end

# ╔═╡ a9cb8c87-dfd9-46bf-b4b9-ab9c004c2667
begin
	plot_volt_time_full = plot(sol_local_bif.t,sol_local_bif[1,:],label="V",linewidth=1,linecolor=:blue)
	plot!(sol_local_bif.t,sol_local_bif[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif.t,sol_local_bif[3,:],label="Vus",linewidth=1,linecolor=RGB(1,0,0.5))

	plot!(t_local_bif,(maximum(vus_local_bif)+1).*ones(size(t_local_bif)),fill = ((minimum(vus_local_bif)-1.5).*ones(size(t_local_bif)), 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest")

	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[vus_local_bif_decr[ind_Ibif_2D_local]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	scatter!([t_local_bif_incr[ind_nbif_2D_local]],[vus_local_bif_incr[ind_nbif_2D_local]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	scatter!([t_local_bif_sl_test[ind_sl_ends]],[Vs_max_low_V_null_2D[ind_sl_ends]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft,size=(250,250))

	yaxis!("Voltage",(-45,-20))
	xaxis!("Time (ms)")
	#savefig("simu3D_local_bif.pdf")
end

# ╔═╡ c5a958b2-434b-4d66-b5a6-63fb293fa295
begin 
	plot(plot_volt_time_full,plot_It_time_full,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_full.pdf")
end

# ╔═╡ 1788b635-ae11-4617-958e-975a4d106c80
begin
	plot_volt_time_roi = plot(sol_local_bif.t,sol_local_bif[1,:],label="V",linewidth=1,linecolor=:blue,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white)
	plot!(sol_local_bif.t,sol_local_bif[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif.t,sol_local_bif[3,:],label="Vus",linewidth=1,linecolor=RGB(1,0,0.5))

	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[vus_local_bif_decr[ind_Ibif_2D_local]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	scatter!([t_local_bif_incr[ind_nbif_2D_local]],[vus_local_bif_incr[ind_nbif_2D_local]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	scatter!([t_local_bif_sl_test[ind_sl_ends]],[Vs_max_low_V_null_2D[ind_sl_ends]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft,size=(250,250))

	xaxis!("Time (ms)",(t_local_bif_decr[1]-2,t_local_bif_incr[end]+2))
	yaxis!("Voltage",((minimum(vus_local_bif)-1.5),maximum(vus_local_bif)+1))

	#savefig("simu3D_local_bif_It.pdf")
end

# ╔═╡ 997f5fe7-d3e6-4968-89f3-ab1348ffa1d4
begin
	plot_It_time_roi = plot(sol_local_bif.t,It_local_bif_full,linecolor=RGB(0.65,0.29,0.65),linewidth=2,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white,label="It")

	scatter!([t_local_bif_decr[ind_Ibif_2D_local]],[Ibif_2D_local],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation")
	scatter!([t_local_bif_incr[ind_nbif_2D_local]],[It_local_bif_incr[ind_nbif_2D_local]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:topleft)
	scatter!([t_local_bif_sl_test[ind_sl_ends[1]]],[It_local_bif_sl_test[ind_sl_ends[1]]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",size=(250,250))

	xaxis!("Time (ms)",(t_local_bif_decr[1]-2,t_local_bif_incr[end]+2))
	yaxis!("Current",(-1.5,3.5))
end

# ╔═╡ 9cd94904-cb24-48b6-9576-b73c2d3aa391
begin
	pltest = plot_It_time_roi

	#savefig("simu3D_local_bif_It.pdf")
	md"""pltest for display of It(t)"""
end

# ╔═╡ 6824efca-1041-4e55-a660-db38962f330b
begin 
	#plot_It_time_roi = pltest
	plot(plot_volt_time_roi,plot_It_time_roi,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_roi.pdf")
end

# ╔═╡ 0c2c6410-3d15-4ed9-9adb-568d986d0cad
begin
	plot(plot_sl_local_bif_It,plot_sl_local_bif_v,plot_sl_local_bif_vs,layout=(3,1),size=(600,400))
	#savefig("simu2D_for_spikelat_3D_time.pdf")
end

# ╔═╡ f46c788e-dd72-4dde-95bc-f5ff15ded9bc
begin
	plot(sol_local_bif.t[(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin],sol_local_bif[2,(ind_spike[end-1]-1)+ind_tmin:(ind_spike[end]-1)+ind_tmin],label="Vs",linewidth=1,linecolor=:orange,size=(200,200))
	plot!(t_local_bif_sl_test,Vs_max_low_V_null_2D,label="Vs of bottom V nullcline maximum",linecolor=RGB(1,0.4,0))
end

# ╔═╡ 579383f0-ef43-42e7-b771-9a0b2dbf1f66
begin
	v_potential = collect(range(-80,0,step=0.5))
	v_null_2D_1 = zeros(size(v_potential))
	v_null_2D_2 = zeros(size(v_potential))
	vs_null_2D = zeros(size(v_potential))
	p_local_bif_2D = change_I_2D(3,p_3D_to_2D(p_local_bif))
	for i=1:length(v_potential)
		v_null_2D_1[i] = Vnullcline1_2D(v_potential[i],p_local_bif_2D)
		v_null_2D_2[i] = Vnullcline2_2D(v_potential[i],p_local_bif_2D)
		vs_null_2D[i] = Vsnullcline_2D(v_potential[i],p_local_bif_2D)
	end


	u0_2D_1_ = fixedpoints_only_2D(change_I_2D(1,p_3D_to_2D(p_local_bif)))


	u0_2D=[sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local],Vnullcline2_2D(sol_local_bif[1,(ind_spike[end-1]-1)+ind_tmin+ind_Ibif_2D_local],p_local_bif_2D)]
	tspan_2D=(0.0,200)
	prob_2D = ODEProblem(MQIF_2D!,u0_2D,tspan_2D,p_local_bif_2D,callback=cb_2D)
	sol_2D = solve(prob_2D,DP5(),reltol=1e-6,abstol=1e-6)

	max_low_v_null_2D_local_bif = Vnullcline2_2D(p_local_bif_2D[2],p_local_bif_2D)
	ind_sl_ends_pp_2D = findfirst(x -> x>=max_low_v_null_2D_local_bif,sol_2D[2,:])

	md"""Nullclines and trajectory computation for the 2D model for the highest current It (3), initialized at the stable node of a low current It (1)"""
end

# ╔═╡ 933002ad-57f5-4573-8f18-fa6c77bef9d3
begin
	plot(v_potential,v_null_2D_1,linewidth=2,linecolor=RGB(0.6,0.4,0.7),label="dV/dt=0")
	plot!(v_potential,v_null_2D_2,linewidth=2,linecolor=RGB(0.6,0.4,0.7),label="dV/dt=0")
	plot!(v_potential,vs_null_2D,linewidth=2,linecolor=RGB(0,0.6,0.7),label="dVs/dt=0")
	plot!(sol_2D[1,:],sol_2D[2,:],linewidth=2,linecolor=RGB(0.2,0.2,0.2),label="Trajectory",linestyle=:dash)
	scatter!([sol_2D[1,ind_sl_ends_pp_2D]],[sol_2D[2,ind_sl_ends_pp_2D]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft,size=(200,200))
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Phase plane for I = $(p_local_bif_2D[1]), V0 = $(p_local_bif_2D[2]), Vs0 = $(p_local_bif_2D[3])")

	#savefig("simu2D_for_spikelat_3D_phase.pdf")
end

# ╔═╡ db7220d3-d245-44e2-86a4-b5e6f8c9ef90
md"""If we initialize the system at the stable equilibrium at the point (V;Vs) where the total current creates a bifurcation (-42,5;-42,5), we can see that the spike latency for the maximum current It (=3) is already of 100ms, supposing that the highest current is applied directly. However, in the 3D simulation, we can see that this current is applied progressively (from the bifurcation current and dring more than 100ms).

The time between the SN bifurcation in the 2D model (where the resting state disappears) is $(round(t_local_bif_sl_test[ind_sl_ends[1]]-t_local_bif_decr[ind_Ibif_2D_local])). which is """

# ╔═╡ 298033d0-6e82-4be3-9d87-70a10f734d2b
begin
	dt = round((sol_2D.t[ind_sl_ends_pp_2D]-sol_2D.t[1])*10)/10
	plot(sol_2D.t,sol_2D[1,:],linewidth=2,linecolor=RGB(0.2,0.2,0.2),label="Trajectory",legend=:topleft)
	scatter!([sol_2D.t[ind_sl_ends_pp_2D]],[sol_2D[1,ind_sl_ends_pp_2D]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",size=(200,200))
	xaxis!("Time (ms)")
	yaxis!("V")
	title!("Time to finish spike latency for I = $(p_local_bif_2D[1]) with I0 = 1")
	#savefig("simu2D_for_spikelat_3D_time_2D.pdf")

end

# ╔═╡ 5ea20fe5-3c23-4cb0-95bb-04dd8672a08c
begin
	#V_stable_local_bif_decr,I_stable_local_bif_decr,V_unstable_local_bif_decr,I_unstable_local_bif_decr,V_saddle_local_bif_decr,I_saddle_local_bif_decr,V_cyclemax_local_bif_decr,V_cyclemin_local_bif_decr,I_cycle_local_bif_decr = separate_local_bif_conv(conv_local_bif_decr,conv_type_local_bif_decr,It_local_bif_decr,stab_fp_local_bif_decr)
	#V_stable_local_bif_incr,I_stable_local_bif_incr,V_unstable_local_bif_incr,I_unstable_local_bif_incr,V_saddle_local_bif_incr,I_saddle_local_bif_incr,V_cyclemax_local_bif_incr,V_cyclemin_local_bif_incr,I_cycle_local_bif_incr = separate_local_bif_conv(conv_local_bif_incr,conv_type_local_bif_incr,It_local_bif_incr,stab_fp_local_bif_incr)

	V_stable_local_bif_decr,I_stable_local_bif_decr,t_stable_local_bif_decr,V_unstable_local_bif_decr,I_unstable_local_bif_decr,t_unstable_local_bif_decr,V_saddle_local_bif_decr,I_saddle_local_bif_decr,t_saddle_local_bif_decr,V_cyclemax_local_bif_decr,V_cyclemin_local_bif_decr,I_cycle_local_bif_decr,t_cycle_local_bif_decr,V_spikelat_local_bif_decr,V_sl_cmin_local_bif_decr,I_spikelat_local_bif_decr,t_spikelat_local_bif_decr = separate_local_bif_conv_full(conv_local_bif_decr,conv_type_local_bif_decr,It_local_bif_decr,t_local_bif_decr,stab_fp_local_bif_decr)

V_stable_local_bif_incr,I_stable_local_bif_incr,t_stable_local_bif_incr,V_unstable_local_bif_incr,I_unstable_local_bif_incr,t_unstable_local_bif_incr,V_saddle_local_bif_incr,I_saddle_local_bif_incr,t_saddle_local_bif_incr,V_cyclemax_local_bif_incr,V_cyclemin_local_bif_incr,I_cycle_local_bif_incr,t_cycle_local_bif_incr,V_spikelat_local_bif_incr,V_sl_cmin_local_bif_incr,I_spikelat_local_bif_incr,t_spikelat_local_bif_incr = separate_local_bif_conv_full(conv_local_bif_incr,conv_type_local_bif_incr,It_local_bif_incr,t_local_bif_incr,stab_fp_local_bif_incr)


#tV_stable_local_bif_decr,t_stable_local_bif_decr,tV_unstable_local_bif_decr,t_unstable_local_bif_decr,tV_saddle_local_bif_decr,t_saddle_local_bif_decr,tV_cyclemax_local_bif_decr,tV_cyclemin_local_bif_decr,t_cycle_local_bif_decr = separate_local_bif_conv_time(conv_local_bif_decr,conv_type_local_bif_decr,t_local_bif_decr,stab_fp_local_bif_decr)
	#tV_stable_local_bif_incr,t_stable_local_bif_incr,tV_unstable_local_bif_incr,t_unstable_local_bif_incr,tV_saddle_local_bif_incr,t_saddle_local_bif_incr,tV_cyclemax_local_bif_incr,tV_cyclemin_local_bif_incr,t_cycle_local_bif_incr = separate_local_bif_conv_time(conv_local_bif_incr,conv_type_local_bif_incr,t_local_bif_incr,stab_fp_local_bif_incr)

	V_sl_cmin_local_bif = zeros(length(V_sl_cmin_local_bif_decr)+length(V_sl_cmin_local_bif_incr))
	I_sl_cmin_local_bif = zeros(size(V_sl_cmin_local_bif))
	for i=1:length(V_sl_cmin_local_bif)
		if i <= length(V_sl_cmin_local_bif_decr)
			V_sl_cmin_local_bif[i]=V_sl_cmin_local_bif_decr[i]
			I_sl_cmin_local_bif[i]=I_spikelat_local_bif_decr[i]
		else
			V_sl_cmin_local_bif[i]=V_sl_cmin_local_bif_incr[i-length(V_sl_cmin_local_bif_decr)]
			I_sl_cmin_local_bif[i]=I_spikelat_local_bif_incr[i-length(V_sl_cmin_local_bif_decr)]
		end
	end
	md"""Separation of convergence value based on their type and stability"""
end

# ╔═╡ 3daae77a-ef7a-4766-b354-cb06598e2047
begin
	V_stable_local_bif_decr_hI,I_stable_local_bif_decr_hI,t_stable_local_bif_decr_hI,V_unstable_local_bif_decr_hI,I_unstable_local_bif_decr_hI,t_unstable_local_bif_decr_hI,V_saddle_local_bif_decr_hI,I_saddle_local_bif_decr_hI,t_saddle_local_bif_decr_hI,V_cyclemax_local_bif_decr_hI,V_cyclemin_local_bif_decr_hI,I_cycle_local_bif_decr_hI,t_cycle_local_bif_decr_hI,V_spikelat_local_bif_decr_hI,V_sl_cmin_local_bif_decr_hI,I_spikelat_local_bif_decr_hI,t_spikelat_local_bif_decr_hI = separate_local_bif_conv_full(conv_local_bif_decr_hI,conv_type_local_bif_decr_hI,It_local_bif_decr_hI,t_local_bif_decr_hI,stab_fp_local_bif_decr_hI)

V_stable_local_bif_incr_hI,I_stable_local_bif_incr_hI,t_stable_local_bif_incr_hI,V_unstable_local_bif_incr_hI,I_unstable_local_bif_incr_hI,t_unstable_local_bif_incr_hI,V_saddle_local_bif_incr_hI,I_saddle_local_bif_incr_hI,t_saddle_local_bif_incr_hI,V_cyclemax_local_bif_incr_hI,V_cyclemin_local_bif_incr_hI,I_cycle_local_bif_incr_hI,t_cycle_local_bif_incr_hI,V_spikelat_local_bif_incr_hI,V_sl_cmin_local_bif_incr_hI,I_spikelat_local_bif_incr_hI,t_spikelat_local_bif_incr_hI = separate_local_bif_conv_full(conv_local_bif_incr_hI,conv_type_local_bif_incr_hI,It_local_bif_incr_hI,t_local_bif_incr_hI,stab_fp_local_bif_incr_hI)

	begin
	V_sl_cmin_local_bif_hI= zeros(length(V_sl_cmin_local_bif_decr_hI)+length(V_sl_cmin_local_bif_incr_hI))
	I_sl_cmin_local_bif_hI = zeros(size(V_sl_cmin_local_bif_hI))
	for i=1:length(V_sl_cmin_local_bif_hI)
		if i <= length(V_sl_cmin_local_bif_decr_hI)
			V_sl_cmin_local_bif_hI[i]=V_sl_cmin_local_bif_decr_hI[i]
			I_sl_cmin_local_bif_hI[i]=I_spikelat_local_bif_decr_hI[i]
		else
			V_sl_cmin_local_bif_hI[i]=V_sl_cmin_local_bif_incr_hI[i-length(V_sl_cmin_local_bif_decr_hI)]
			I_sl_cmin_local_bif_hI[i]=I_spikelat_local_bif_incr_hI[i-length(V_sl_cmin_local_bif_decr_hI)]
		end
	end
end
	md"""Separation of convergence value based on their type and stability"""
end

# ╔═╡ a516dbb5-57b5-42ef-a10b-1358694b9daa
begin
	I_local_bif = collect(range(-1.5,stop=5,step=0.05))

	list_stable,list_saddle,list_unstable,list_cycleV = bifI_matrices_2D(I_local_bif,p_3D_to_2D(p_local_bif))

	md"""Matrices calculation to characterize fixed points and limit cycle for a bifurcation diagram with I in 2D"""
end

# ╔═╡ abdf6c7f-e761-46e7-8f5b-d0d2d7a2d285
begin
	ind_bif_2D_I = findfirst(x -> x>= Ibif_2D_local,I_local_bif)-1
	v_bif_2D_local_bif = (list_saddle[ind_bif_2D_I]+list_stable[ind_bif_2D_I])/2

	plot(I_local_bif, list_stable,label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)

	ind_real_unstable_ = findall(x -> x>=0 || x<0,list_unstable)
	if length(ind_real_unstable_)>0
		plot!(I_local_bif, list_unstable,label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	end
	plot!(I_local_bif, list_saddle,label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)

	#plot!(I_local_bif, list_cycleV[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	ind_real_ = findall(x -> x>=0 || x<0,list_cycleV[:,1])
	plot!(I_local_bif[ind_real_], (-35).*ones(size(list_cycleV[ind_real_,1])),label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5)

	plot!(I_local_bif, list_cycleV[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,legend=:outertopright)

	scatter!([Ibif_2D_local],[v_bif_2D_local_bif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	yaxis!("V")
	xaxis!("I")

	##Trajectory
	plot!(I_stable_local_bif_decr,V_stable_local_bif_decr,linecolor=RGB(0.0,0.6,0.2),linewidth=5,linealpha=0.5,label="Decreasing Vus, It <= Ibif")
	plot!(I_sl_cmin_local_bif,V_sl_cmin_local_bif,linecolor=RGB(0.1,0.6,0.6),linewidth=5,linealpha=0.5,label="Low Vus, It > Ibif")
	plot!(I_cycle_local_bif_incr,(V_cyclemax_local_bif_incr./maximum(V_cyclemax_local_bif_incr)).*(-35),linecolor=RGB(0.8,0.4,0.8),linewidth=5,linealpha=0.5,label="Increasing Vus, It > 0")

	##arrow SN bif
	plot!([I_stable_local_bif_decr[end],I_sl_cmin_local_bif[1]],[V_stable_local_bif_decr[end],V_sl_cmin_local_bif[1]],linecolor=RGB(1,0.8,0.0),linewidth=3,linealpha=0.5,linestyle=:dashdotdot,label="Low Vus, It -> Ibif")
	plot!([I_sl_cmin_local_bif[1]-0.15,I_sl_cmin_local_bif[1],I_sl_cmin_local_bif[1]+0.15],[V_sl_cmin_local_bif[1]-0.15,V_sl_cmin_local_bif[1],V_sl_cmin_local_bif[1]-0.15],linecolor=RGB(1,0.8,0.0),linewidth=2,linealpha=0.5,label="Low Vus, It -> Ibif")

	##arrow spike latency to cycle
	plot!([I_sl_cmin_local_bif[end],I_cycle_local_bif_incr[1]],[V_sl_cmin_local_bif[end],(V_cyclemax_local_bif_incr[1]./maximum(V_cyclemax_local_bif_incr)).*(-35)],linecolor=RGB(1,0.7,0.9),linewidth=2,linealpha=0.5,linestyle=:dashdotdot,label="Spike latency ends")
	plot!([I_cycle_local_bif_incr[1]-0.15,I_cycle_local_bif_incr[1],I_cycle_local_bif_incr[1]+0.15],[(-35)-0.15,(-35),(-35)-0.15],linecolor=RGB(1,0.7,0.9),linewidth=2,linealpha=0.5,label="Spike latency ends")

	ind_cycle_for_min_I = findfirst(x -> x==minimum(I_cycle_local_bif_incr),I_cycle_local_bif_incr)
	##arrow cycle to rest
	plot!([I_cycle_local_bif_incr[ind_cycle_for_min_I],I_stable_local_bif_decr[1]],[(-35),V_stable_local_bif_decr[1]],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,linestyle=:dashdotdot,label="High Vus, It <0")
	plot!([I_stable_local_bif_decr[1]-0.15,I_stable_local_bif_decr[1],I_stable_local_bif_decr[1]+0.15],[V_stable_local_bif_decr[1]+0.15,V_stable_local_bif_decr[1],V_stable_local_bif_decr[1]+0.15],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,label="High Vus, It <0")

	#savefig("bifdiag2D_traj_corresp_to_vus_3D.pdf")
end

# ╔═╡ 771927b2-814a-4422-8d78-c4020cc89cda
begin
	ind_v_local_bif_mod = findall(x->x>0,v_local_bif)
	v_local_bif_mod = v_local_bif
	for i=1:length(ind_v_local_bif_mod)
		v_local_bif_mod[ind_v_local_bif_mod[i]]=0
	end
end

# ╔═╡ d3db0388-7f11-41ef-887e-46187b3e4b53
gr()

# ╔═╡ 733c2729-3558-4660-a059-93d9d8220ad4
begin
	
	plot_bif_IV_traj_cmap = plot(I_local_bif, list_stable,label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)

	if length(ind_real_unstable_)>0
		plot!(I_local_bif, list_unstable,label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	end
	plot!(I_local_bif, list_saddle,label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)

	#plot!(I_local_bif, list_cycleV[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	plot!(I_local_bif[ind_real_], (0).*ones(size(list_cycleV[ind_real_,1])),label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5)

	plot!(I_local_bif, list_cycleV[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,legend=:outertopleft)

	scatter!([Ibif_2D_local],[v_bif_2D_local_bif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	yaxis!("V",(-45,0.5))
	xaxis!("I",(-1.5,5))

	##Trajectory
	plot!(I_stable_local_bif_decr,V_stable_local_bif_decr,linecolor=RGB(0.0,0.6,0.2),linewidth=5,linealpha=0.3,label="Decreasing Vus, It <= Ibif")
	plot!(I_sl_cmin_local_bif,V_sl_cmin_local_bif,linecolor=RGB(0.1,0.6,0.6),linewidth=5,linealpha=0.3,label="Low Vus, It > Ibif")
	plot!(I_cycle_local_bif_incr,(V_cyclemax_local_bif_incr./maximum(V_cyclemax_local_bif_incr)).*(0),linecolor=RGB(0.8,0.4,0.8),linewidth=5,linealpha=0.3,label="Increasing Vus, It > 0")


	##arrow cycle to rest
	plot!([I_cycle_local_bif_incr[ind_cycle_for_min_I],I_stable_local_bif_decr[1]],[(0),V_stable_local_bif_decr[1]],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,linestyle=:dashdotdot,label="High Vus, It <0")

	##Trajectory V(It)
	#cool color maps
	#isoluminant_cgo_70_c39_n256
	#linear_ternary_blue_0_44_c57_n256
	#linear_bmy_10_95_c78_n256
	#neon
	t_local_bif_normed=t_local_bif.-t_local_bif[1]
	plot!(It_local_bif,v_local_bif_mod,line_z=t_local_bif_normed,c=:matter,linewidth=2,linealpha=1,label="Trajectory V(It)",size=(800,500),colorbar=:bottom)

	#savefig("bif_diag_2D_traj_IV.pdf")
end

# ╔═╡ 053be29e-3c94-4ad2-b89f-7dafd86b06f9
begin
	plot_local_bif_decr = plot(I_stable_local_bif_decr,V_stable_local_bif_decr,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(I_saddle_local_bif_decr,V_saddle_local_bif_decr,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	if length(I_unstable_local_bif_decr)>0
		plot!(I_unstable_local_bif_decr,V_unstable_local_bif_decr,linecolor="red",linewidth=1.5,label="Unstable node")
	end
	plot!(I_spikelat_local_bif_decr,V_spikelat_local_bif_decr,linecolor="cyan",linewidth=1.5,label="Spike latency")

	plot!(I_spikelat_local_bif_decr,V_sl_cmin_local_bif_decr,linecolor="purple",linewidth=1.5,label="Limit cycle minimum",linestyle=:dot)

	scatter!([I_stable_local_bif_decr[end]],[(V_stable_local_bif_decr[end]+V_saddle_local_bif_decr[end])/2],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation",legend=:topright)

	#plot!(I_cycle_local_bif_decr,V_cyclemax_local_bif_decr,linewidth=1.5,label="Maximum value of the limit cycle")
	plot!(I_cycle_local_bif_decr,V_cyclemin_local_bif_decr,linewidth=1.5,label="Limit cycle minimum")
	yaxis!("V")
	xaxis!("It")
	title!("Increasing It with time")

	plot_vus_local_bif_decr = plot(t_local_bif_decr,vus_local_bif_decr,linecolor=RGB(1,0,0.5),linewidth=2,label="Vus(t)")
	plot!(t_local_bif_decr,p_local_bif[4].*ones(size(t_local_bif_decr)),linecolor=:blue,linewidth=2,label="Vus0")
	xaxis!("Time (ms)")
	yaxis!("Vus")
	title!("Corresponding Vus ")

	plot_It_local_bif_decr = plot(t_local_bif_decr,It_local_bif_decr,linecolor=:pink,linewidth=2,label="It(t)")
	xaxis!("Time (ms)")
	yaxis!("It")
	title!("Corresponding It ")

	plot(plot_local_bif_decr,plot_vus_local_bif_decr,plot_It_local_bif_decr,layout=(1,3),size=(700,400),legend=:topright)
	#savefig("increasing3D_local_bif.pdf")
	md"""Bifurcation with It, decreasing vus, subplot 3x1"""
end

# ╔═╡ b88a07e8-97ce-476b-8e0d-f5de53fda620
begin
	plot_local_bif_decr_t = plot(t_stable_local_bif_decr,V_stable_local_bif_decr,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(t_saddle_local_bif_decr,V_saddle_local_bif_decr,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	if length(I_unstable_local_bif_decr)>0
		plot!(t_unstable_local_bif_decr,V_unstable_local_bif_decr,linecolor="red",linewidth=1.5,label="Unstable node")
	end
	plot!(t_spikelat_local_bif_decr,V_spikelat_local_bif_decr,linecolor="cyan",linewidth=1.5,label="Limit cycle minimum")
	plot!(t_spikelat_local_bif_decr,V_sl_cmin_local_bif_decr,linecolor="purple",linewidth=1.5,label="Spike latency",linestyle=:dot)
	#plot!(t_local_bif_decr,v_local_bif_decr,linecolor="black",linewidth=1.5,label="Trajectory")

	scatter!([t_stable_local_bif_decr[end]],[(V_stable_local_bif_decr[end]+V_saddle_local_bif_decr[end])/2],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation",legend=:topright)

	#plot!(t_cycle_local_bif_decr,tV_cyclemax_local_bif_decr,linewidth=1.5,label="Maximum value of the limit cycle")
	if length(V_cyclemin_local_bif_decr)>0
		plot!(t_cycle_local_bif_decr,V_cyclemin_local_bif_decr,linewidth=1.5,label="Limit cycle minimum")
	end
	yaxis!("V")
	xaxis!("Time (ms)")
	title!("Increasing It with time")


	plot(plot_local_bif_decr_t,plot_vus_local_bif_decr,plot_It_local_bif_decr,layout=(1,3),size=(700,400))
	
	md"""Bifurcation with t, decreasing vus, subplot 3x1"""
end

# ╔═╡ 2517434c-1aa9-4cb5-8ed5-9b5943632538
begin
	plot_local_bif_incr = plot(t_stable_local_bif_incr,V_stable_local_bif_incr,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(t_saddle_local_bif_incr,V_saddle_local_bif_incr,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	if length(t_unstable_local_bif_incr)>0
		plot!(t_unstable_local_bif_incr,V_unstable_local_bif_incr,linecolor="red",linewidth=1.5,label="Unstable node")
	end
	plot!(t_spikelat_local_bif_incr,V_spikelat_local_bif_incr,linecolor="cyan",linewidth=1.5,label="Limit cycle minimum")
	plot!(t_spikelat_local_bif_incr,V_sl_cmin_local_bif_incr,linecolor="purple",linewidth=1.5,label="Spike latency",linestyle=:dot)
	#plot!(t_local_bif_incr,v_local_bif_incr,linecolor="black",linewidth=1.5,label="Trajectory")

	#plot!(t_cycle_local_bif_incr,V_cyclemax_local_bif_incr,linewidth=1.5,label="Maximum value of the limit cycle")
	plot!(t_cycle_local_bif_incr,V_cyclemin_local_bif_incr,linewidth=1.5,label="Limit cycle minimum",linecolor=:purple,legend=:bottomright)
	yaxis!("V")
	xaxis!("Time")
	title!("Decreasing It with time")

	plot_vus_local_bif_incr = plot(t_local_bif_incr,vus_local_bif_incr,linecolor=RGB(1,0,0.5),linewidth=2,label="Vus(t)")
	plot!(t_local_bif_incr,p_local_bif[4].*ones(size(t_local_bif_incr)),linecolor=:blue,linewidth=2,label="Vus0")
	xaxis!("Time (ms)")
	yaxis!("Vus")
	title!("Corresponding Vus ")

	plot_It_local_bif_incr = plot(t_local_bif_incr,It_local_bif_incr,linecolor=:pink,linewidth=2,label="It(t)")
	xaxis!("Time (ms)")
	yaxis!("It")
	title!("Corresponding It ")

	plot(plot_local_bif_incr,plot_vus_local_bif_incr,plot_It_local_bif_incr,layout=(1,3),size=(700,400))
	
	md"""Bifurcation with t, increasing vus, subplot 3x1"""
end

# ╔═╡ 9064707b-c621-44be-94f8-e4e15242e447
begin

	plot(t_stable_local_bif_incr,V_stable_local_bif_incr,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(t_stable_local_bif_decr.+(t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]),V_stable_local_bif_decr,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(t_saddle_local_bif_incr,V_saddle_local_bif_incr,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	plot!(t_saddle_local_bif_decr.+(t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]),V_saddle_local_bif_decr,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	plot!(t_spikelat_local_bif_incr,V_spikelat_local_bif_incr,linecolor="cyan",linewidth=1.5,label="Spike Latency")
	plot!(t_spikelat_local_bif_incr,V_sl_cmin_local_bif_incr,linecolor="purple",linewidth=1.5,label="Limit cycle minimum",linestyle=:dot)

	if length(t_unstable_local_bif_incr)>0
		plot!(t_unstable_local_bif_incr,V_unstable_local_bif_incr,linecolor="red",linewidth=1.5,label="Unstable node")
	end
	scatter!([t_stable_local_bif_decr[end].+(t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1])],[(V_stable_local_bif_decr[end]+V_saddle_local_bif_decr[end])/2],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation",legend=:topright)

	plot!(t_cycle_local_bif_incr,V_cyclemin_local_bif_incr,linewidth=1.5,label="Limit cycle minimum",linecolor=:purple,legend=:none)
	plot!(t_spikelat_local_bif_decr.+(t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]),V_spikelat_local_bif_decr,linecolor="cyan",linewidth=1.5,label="Spike latency")
	plot!(t_spikelat_local_bif_decr.+(t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]),V_sl_cmin_local_bif_decr,linecolor="purple",linewidth=1.5,label="Limit cycle minimum",linestyle=:dot)
	plot!(t_spikelat_local_bif_incr.+(t_spikelat_local_bif_decr[end]+t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]-t_spikelat_local_bif_incr[1]),V_spikelat_local_bif_incr,linecolor="cyan",linewidth=1.5,label="Spike latency")
	plot!(t_spikelat_local_bif_incr.+(t_spikelat_local_bif_decr[end]+t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]-t_spikelat_local_bif_incr[1]),V_sl_cmin_local_bif_incr,linecolor="purple",linewidth=1.5,label="Limit cycle minimum",linestyle=:dot)
	plot!(t_cycle_local_bif_incr.+(t_spikelat_local_bif_decr[end]+t_stable_local_bif_incr[end]-t_stable_local_bif_decr[1]-t_spikelat_local_bif_incr[1]),V_cyclemin_local_bif_incr,linewidth=1.5,label="Spike latency",linecolor=:purple,legend=:outertopright)

	plot!(sol_local_bif.t,sol_local_bif[1,:])

	#plot!(sol_local_bif.t,sol_local_bif[1,:],linealpha=0.5,linecolor=:blue)

	yaxis!((-45,-38),"V")
	xaxis!((945,1330),"Time(ms)")
	title!("Full Vus cycle impact on V")

	#savefig("bifdiag2D_time_traj_corresp_to_vus_3D.pdf")
	md"""Bifurcation with t, cycle of vus"""
end

# ╔═╡ ef6bdd3f-00bc-43ff-b8b6-0d266abdfd57
md"###### 2D phase planes +real-time 3D traj display for gif"

# ╔═╡ f61a1af4-5f02-4fe8-9c51-912d8ff2de11
begin
	v_local_bif,vs_local_bif,vus_local_bif,t_local_bif,It_local_bif
	md"""Vectors to consider"""
end

# ╔═╡ 7358098d-88a6-42a7-b401-a3088c98b8fd
function reconstruc_sample_I(sample_It_local_bif,sample_index_local_bif,t_local_bif,n_sample_local_bif_decr,n_sample_local_bif_incr)
	
	recon_It_local_bif = ones(size(t_local_bif))
	for i=2:n_sample_local_bif_decr
		ind1 = sample_index_local_bif[i-1]
		ind2 = sample_index_local_bif[i]
		recon_It_local_bif[ind1:ind2].=sample_It_local_bif[i-1]
	end
	
	indtrans1 = sample_index_local_bif[n_sample_local_bif_decr]
	indtrand2 = sample_index_local_bif[1+n_sample_local_bif_decr]
	recon_It_local_bif[indtrans1:indtrand2].=sample_It_local_bif[n_sample_local_bif_decr]
	
	for i=2:n_sample_local_bif_incr
		ind1 = sample_index_local_bif[i-1+n_sample_local_bif_decr]
		ind2 = sample_index_local_bif[i+n_sample_local_bif_decr]
		#recon_It_local_bif[ind1:ind2].=sample_It_local_bif[i+n_sample_local_bif_decr]
		recon_It_local_bif[ind1:ind2].=sample_It_local_bif[i-1+n_sample_local_bif_decr]
	end
	return recon_It_local_bif
end

# ╔═╡ 4802d83b-9051-46ef-94aa-b60b83a0fcde
begin
	##sampling
	n_sample_local_bif_decr = 30
	n_sample_local_bif_incr = 7
	
	ref_ind = ind_Ibif_2D_local+ind_sl_ends[1]-1
	
	step_sample_local_bif_decr = Int(round(length(t_local_bif[1:ref_ind])/n_sample_local_bif_decr))
	step_sample_local_bif_incr = Int(round(length(t_local_bif[ref_ind:end])/n_sample_local_bif_incr))
	
	#add offset to be sure to take the last sample of It into account 
	offset_sample_local_bif_incr = length(t_local_bif[ref_ind:end])-step_sample_local_bif_incr*n_sample_local_bif_incr
	
	sample_It_local_bif = []
	sample_t_local_bif = []
	sample_v_local_bif = []
	sample_vs_local_bif = []
	sample_index_local_bif = []
	
	for i=0:n_sample_local_bif_decr-1
		append!(sample_It_local_bif,It_local_bif[1+i*step_sample_local_bif_decr])
		append!(sample_t_local_bif,t_local_bif[1+i*step_sample_local_bif_decr])
		append!(sample_v_local_bif,v_local_bif[1+i*step_sample_local_bif_decr])
		append!(sample_vs_local_bif,vs_local_bif[1+i*step_sample_local_bif_decr])
		append!(sample_index_local_bif,1+i*step_sample_local_bif_decr)
	end
	
	for i=1:n_sample_local_bif_incr
		append!(sample_It_local_bif,It_local_bif[ref_ind-1+offset_sample_local_bif_incr+i*step_sample_local_bif_incr])
		append!(sample_t_local_bif,t_local_bif[ref_ind-1+offset_sample_local_bif_incr+i*step_sample_local_bif_incr])
		append!(sample_v_local_bif,v_local_bif[ref_ind-1+offset_sample_local_bif_incr+i*step_sample_local_bif_incr])
		append!(sample_vs_local_bif,vs_local_bif[ref_ind-1+offset_sample_local_bif_incr+i*step_sample_local_bif_incr])
		append!(sample_index_local_bif,ref_ind-1+offset_sample_local_bif_incr+i*step_sample_local_bif_incr)
	end
	
	##reconstruction
	recon_It_local_bif = reconstruc_sample_I(sample_It_local_bif,sample_index_local_bif,t_local_bif,n_sample_local_bif_decr,n_sample_local_bif_incr)
	recon_list_index = [1]
	
	
	
	md"""Samples of v,vs,vus,t and It to show the evolution of the system in a 2D phase plane"""
end

# ╔═╡ bd4faf19-cbee-46fd-9fbf-84d30da07cbb
begin
	plot(t_local_bif,It_local_bif,label="It")
	plot!(t_local_bif,recon_It_local_bif,legend=:outertopright,label="Sampled It for gif",size=(250,250))
	#plot!(t_local_bif_decr,It_local_bif_decr,line_z=t_local_bif,linewidth=2)
end

# ╔═╡ 9c063b76-e3ef-40b1-9691-aadcdd46b82b
function give_me_phase_planes(v_potential,p3D,It,v,vs,sample_index)
	
	test = [1,10]
	v_null_1_list = []
	v_null_2_list = []
	traj_v_list = []
	traj_vs_list = []
	for i=1:2	
		#for a given I	
		v_null_2D_1 = zeros(size(v_potential))
		v_null_2D_2 = zeros(size(v_potential))
		p_2D = change_I_2D(It[sample_index[test[i]]],p_3D_to_2D(p3D))
		for i=1:length(v_potential)
			v_null_2D_1[i] = Vnullcline1_2D(v_potential[i],p_2D)
			v_null_2D_2[i] = Vnullcline2_2D(v_potential[i],p_2D)
		end
		
		if i<2
			traj_v = v[sample_index[test[i]]:sample_index[test[i+1]]]
			traj_vs = vs[sample_index[test[i]]:sample_index[test[i+1]]]
		else
			traj_v = [NaN]
			traj_vs = [NaN]
		end
		
		append!(v_null_1_list,v_null_2D_1)
		append!(v_null_2_list,v_null_2D_2)
		append!(traj_v_list,traj_v)
		append!(traj_vs_list,traj_vs)
	end
	return v_null_1_list,v_null_2_list,traj_v_list,traj_vs_list
end

# ╔═╡ 676cd0a8-4f10-407d-88e5-d777acce0758
function give_me_traj_v(v,sample_index,i)
	if i<length(sample_index)
		traj_v = v[1:sample_index[i+1]]
	else
		traj_v = ones(length(v[1:sample_index[end]])+1)
		traj_v[1:end-1] .= v[1:sample_index[end]]
		traj_v[end] = v[sample_index[1]]
	end
	return traj_v
end

# ╔═╡ c03e6b55-6d09-4c71-987e-99227614b2d1
function give_me_traj_vs(vs,sample_index,i)
	if i<length(sample_index)
		traj_vs = vs[1:sample_index[i+1]]
	else
		traj_vs = ones(length(vs[1:sample_index[end]])+1)
		traj_vs[1:end-1] .= vs[1:sample_index[end]]
		traj_vs[end] = vs[sample_index[1]]
	end
	return traj_vs
end

# ╔═╡ dd3e2e91-9dbc-45ad-88db-2fee49fa0ece
function give_me_v_null_1(v_potential,p3D,It,sample_index,i)
	v_null_2D_1 = zeros(size(v_potential))
	p_2D = change_I_2D(It[sample_index[i]],p_3D_to_2D(p3D))
	for j=1:length(v_potential)
		if It[sample_index[i]] >=0
			v_null_2D_1[j] = Vnullcline1_2D(v_potential[j],p_2D)
		else
			v_null_2D_1[j] = Vnullcline1_2D_(v_potential[j],p_2D)
		end
	end
	return v_null_2D_1
end

# ╔═╡ 7b1ec15f-f926-49f5-bbd8-fedaec8432c9
function give_me_v_null_2(v_potential,p3D,It,sample_index,i)
	v_null_2D_2 = zeros(size(v_potential))
	p_2D = change_I_2D(It[sample_index[i]],p_3D_to_2D(p3D))
	for j=1:length(v_potential)
		if It[sample_index[i]] >=0
			v_null_2D_2[j] = Vnullcline2_2D(v_potential[j],p_2D)
		else
			v_null_2D_2[j] = Vnullcline2_2D_(v_potential[j],p_2D)
		end
	end
	return v_null_2D_2
end

# ╔═╡ 262f3030-053a-457f-bebc-32ddaa19b154
function give_me_fp1_2D(p3D,It,sample_index,i)
	p_2D = change_I_2D(It[sample_index[i]],p_3D_to_2D(p3D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return fp1	
end

# ╔═╡ e2fea39c-6a9f-4e68-9be2-ffae78c9a807
function give_me_fp2_2D(p3D,It,sample_index,i)
	p_2D = change_I_2D(It[sample_index[i]],p_3D_to_2D(p3D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return fp2	
end

# ╔═╡ 840c7361-a9b1-4fdb-ab67-8116428da226
function give_me_stab_2D(p3D,It,sample_index,i)
	p_2D = change_I_2D(It[sample_index[i]],p_3D_to_2D(p3D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return stability	
end

# ╔═╡ 88c95811-b6cc-4512-9756-c796ddcb6ea2
function give_me_phase_planes_l(v_potential,p3D,It,v,vs,sample_index)
	
	v_null_1_list = [give_me_v_null_1(v_potential,p3D,It,sample_index,i) for i=1:length(sample_index)]
	v_null_2_list = [give_me_v_null_2(v_potential,p3D,It,sample_index,i) for i=1:length(sample_index)]
	traj_v_list = [give_me_traj_v(v,sample_index,i) for i=1:length(sample_index)]
	traj_vs_list = [give_me_traj_vs(vs,sample_index,i) for i=1:length(sample_index)]
	fp1_list = [give_me_fp1_2D(p3D,It,sample_index,i) for i=1:length(sample_index)]
	fp2_list = [give_me_fp2_2D(p3D,It,sample_index,i) for i=1:length(sample_index)]
	stab_list = [give_me_stab_2D(p3D,It,sample_index,i) for i=1:length(sample_index)]
	
	return v_null_1_list,v_null_2_list,traj_v_list,traj_vs_list,fp1_list,fp2_list, stab_list
end

# ╔═╡ b8fd420a-da3d-485a-a98d-c2271b24b6fa
function make_my_phase_plane(v_potential,v_null_1_list,v_null_2_list,vs_null_2D, traj_v_list,traj_vs_list,fp1_list,fp2_list,stab_list,recon_It_local_bif,t_local_bif,sample_index_local_bif,i)
	
	plot_test = plot(v_potential,vs_null_2D,linewidth=2,linecolor=RGB(0,0.6,0.7),label="dVs/dt=0")
	
	if recon_It_local_bif[sample_index_local_bif[i]]>=0
		plot!(v_potential,v_null_1_list[i],linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
		plot!(v_potential,v_null_2_list[i],linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
	else
		plot!(v_null_1_list[i],v_potential,linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
		plot!(v_null_2_list[i],v_potential,linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
	end

	plot!(traj_v_list[i],traj_vs_list[i],linewidth=2,line_z=t_local_bif,c=:cyclic_mygbm_30_95_c78_n256,label="Trajectory",linealpha=0.9,legend=:topleft)
	
	if length(stab_list[i])>1
		##FP1
		if (stab_list[i][1])>0
			scatter!([fp1_list[i][1]],[fp1_list[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
		end
		if (stab_list[i][1])<0
			scatter!([fp1_list[i][1]],[fp1_list[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
		end
		if (stab_list[i][1])==0.0
			scatter!([fp1_list[i][1]],[fp1_list[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:topleft)#saddle
		end
		
		##FP2
		if (stab_list[i][2])>0
			scatter!([fp2_list[i][1]],[fp2_list[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
		end
		if (stab_list[i][2])<0
			scatter!([fp2_list[i][1]],[fp2_list[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
		end
		if (stab_list[i][2])==0.0
			scatter!([fp2_list[i][1]],[fp2_list[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:topleft)#saddle
		end
	end
	
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	#title!("t =$(round(t_local_bif[sample_index_local_bif[i]]*100)/100), I = $(round((recon_It_local_bif[sample_index_local_bif[i]])*100)/100), V0 = $(p_local_bif_2D[2]), Vs0 = $(p_local_bif_2D[3])")
	title!("  A. ",titlelocation=:left,titlefontsize=18)
	return plot_test
end

# ╔═╡ e30ec005-016a-45a6-abf8-3d6eda1de911
begin
	v_null_1_list,v_null_2_list,traj_v_list,traj_vs_list,fp1_list,fp2_list,stab_list = give_me_phase_planes_l(v_potential,p_local_bif,recon_It_local_bif,v_local_bif,vs_local_bif,sample_index_local_bif)
end

# ╔═╡ d13d81a7-c7e8-4096-833a-cda6f5f99c32
begin
	plot_list = [make_my_phase_plane(v_potential,v_null_1_list,v_null_2_list,vs_null_2D,traj_v_list,traj_vs_list,fp1_list,fp2_list,stab_list,recon_It_local_bif,t_local_bif,sample_index_local_bif,i) for i=1:length(traj_v_list)]
	

	#savefig("simu2D_for_spikelat_3D_phase.pdf")
end

# ╔═╡ f2ce1d7b-8392-4f20-b567-3593889b73c8
begin
	in = 37
	plt = plot_list[in]
	#savefig("gif_part_$(in).png")
end

# ╔═╡ 3d6b6513-397b-4369-895b-7eee2678ec37
begin
	plot(plt)
	#savefig("gif_part_$(in).png")
	md"""save gif plot part"""
end

# ╔═╡ 8b41be33-fad7-4f3f-94b2-dc4be766102c
md""" For a higher input current, Vus(t) obviously changes but the 2D SN bifurcation does not. """

# ╔═╡ 25e64d32-8603-44d4-8a9a-e3396c9ed17a
md"""###### Bifurcation from 2D model : taus = 100 ms,  restorative feedbacks"""

# ╔═╡ 82ee2ffb-5375-4271-ad7b-f15c5f520d4a
begin
	Ibif_2D_local_hI = bifI2D(p_3D_to_2D(p_local_bif_hI))
	ind_Ibif_2D_local_hI = findfirst(x -> x>= Ibif_2D_local,It_local_bif_decr_hI)
	ind_nbif_2D_local_hI = length(It_local_bif_incr_hI)-(findfirst(x -> x>= 0,reverse(It_local_bif_incr_hI))-2)
	md"""Bifurcation indices"""
end

# ╔═╡ 73cd72e6-c389-4c40-803e-c4513995d438
begin
	It_local_bif_sl_test_hI = zeros(length(It_local_bif_decr_hI[ind_Ibif_2D_local_hI:end])+length(It_local_bif_incr_hI[1:ind_nbif_2D_local_hI-1]))
	t_local_bif_sl_test_hI = zeros(length(It_local_bif_sl_test_hI))
	Vs_max_low_V_null_2D_hI = zeros(length(It_local_bif_sl_test_hI))

	for i=1:length(It_local_bif_sl_test_hI)
		if i<=length(It_local_bif_decr_hI)-(ind_Ibif_2D_local_hI-1)
			It_local_bif_sl_test_hI[i] = It_local_bif_decr_hI[i+(ind_Ibif_2D_local_hI-1)]
			t_local_bif_sl_test_hI[i] = t_local_bif_decr_hI[i+(ind_Ibif_2D_local_hI-1)]
		else
			It_local_bif_sl_test_hI[i] = It_local_bif_incr_hI[i-(length(It_local_bif_decr_hI)-(ind_Ibif_2D_local_hI-1))]
			t_local_bif_sl_test_hI[i] = t_local_bif_incr_hI[i-(length(It_local_bif_decr_hI)-(ind_Ibif_2D_local_hI-1))]
		end
		Vs_max_low_V_null_2D_hI[i]=max_2D_low_V_null(It_local_bif_sl_test_hI[i])
	end
	md"""Index for end of spike latency : V nullcline maximum computation"""
end

# ╔═╡ e9b8a4ca-05cb-4aad-b2fd-3b22c3fc23c2
begin
	ind_sl_ends_hI = [-1]
	for i=0:length(Vs_max_low_V_null_2D_hI)-1
		if Vs_max_low_V_null_2D_hI[i+1]<=sol_local_bif_hI[2,(ind_spike_hI[end-1]-1)+ind_tmin_hI+ind_Ibif_2D_local_hI+i] && ind_sl_ends_hI[1] <0
			ind_sl_ends_hI[1] = i+1
		end
	end
	md"""Index for end of spike latency"""
end

# ╔═╡ e7637146-0496-495b-892a-08668620e00b
begin
	plot_volt_time_full_hI = plot(sol_local_bif_hI.t,sol_local_bif_hI[1,:],label="V",linewidth=1,linecolor=:blue)
	plot!(sol_local_bif_hI.t,sol_local_bif_hI[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_hI.t,sol_local_bif_hI[3,:],label="Vus",linewidth=1,linecolor=RGB(1.0,0.0,0.5))
	plot!(t_local_bif_hI,(maximum(vus_local_bif_hI)+1).*ones(size(t_local_bif_hI)),fill = ((minimum(v_local_bif_hI)-1).*ones(size(t_local_bif_hI)), 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest",size=(250,250))

	#scatter!([t_local_bif_decr_hI[ind_Ibif_2D_local_hI]],[vus_local_bif_decr_hI[ind_Ibif_2D_local_hI]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	#scatter!([t_local_bif_incr_hI[ind_nbif_2D_local_hI]],[vus_local_bif_incr_hI[ind_nbif_2D_local_hI]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	#scatter!([t_local_bif_sl_test_hI[ind_sl_ends_hI]],[Vs_max_low_V_null_2D_hI[ind_sl_ends_hI]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)

	yaxis!("Voltage",(-45,-20))
	xaxis!("Time (ms)")
end

# ╔═╡ 07b5937f-f95a-49e4-8c66-8c7ae791b7d9
begin
	It_local_bif_full_hI = zeros(length(sol_local_bif_hI[3,:]))
	for i=1:length(sol_local_bif_hI[3,:])
		It_local_bif_full_hI[i] = compute_I_from_vus(sol_local_bif_hI[3,i],p_local_bif_hI)
	end
	
	plot_It_time_full_hI = plot(sol_local_bif_hI.t,It_local_bif_full_hI,label="It",linewidth=1.5,linecolor=RGB(0.65,0.29,0.65))	
	plot!(sol_local_bif_hI.t,Ibif_2D_local_hI.*ones(size(It_local_bif_full_hI)),label="2D SN bifurcation",linewidth=1.5,linecolor=RGB(0.9,0.9,0),legend=:topleft)
	plot!(sol_local_bif_hI.t,zeros(size(It_local_bif_full_hI)),label="Bifurcation to rest",linewidth=1.5,linecolor=RGB(0,0.5,1))
	
	plot!(sol_local_bif_hI.t,(4.5).*ones(size(sol_local_bif_hI.t)),fill = ([Ibif_2D_local_hI.*ones(size(sol_local_bif_hI.t))], 0.1, RGB(0,0.5,1)),linecolor=:purple,linealpha=0.1,label="Limit cycle")
	plot!(sol_local_bif_hI.t,Ibif_2D_local_hI.*ones(size(sol_local_bif_hI.t)),fill = ([zeros(size(It_local_bif_full_hI))], 0.1, RGB(1,1,0)),linecolor=:purple,linealpha=0.1,label="Bistable")
	plot!(sol_local_bif_hI.t,zeros(size(It_local_bif_full_hI)),fill = ([(-1.5).*ones(size(It_local_bif_full_hI))], 0.1, RGB(0,1,0)),linecolor=:purple,linealpha=0.1,label="Stable")	
	
	plot!(t_local_bif_hI,(maximum(It_local_bif_full_hI)+2).*ones(size(t_local_bif_hI)),fill = ((-1.5).*ones(size(t_local_bif_hI)), 0.6, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest",legend=:topleft,size=(250,250))
	
	
	yaxis!("Current",(-1.5,4))
	xaxis!("Time (ms)")
	#savefig("simu3D_local_bif_fullIt.pdf")
end

# ╔═╡ 5dcaee47-e015-4592-8144-979d2616804e
begin
	plot_volt_time_roi_hI = plot(sol_local_bif_hI.t,sol_local_bif_hI[1,:],label="V",linewidth=1,linecolor=:blue,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white)
	plot!(sol_local_bif_hI.t,sol_local_bif_hI[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_hI.t,sol_local_bif_hI[3,:],label="Vus",linewidth=1,linecolor=RGB(1,0,0.5),size=(250,250))

	#scatter!([t_local_bif_decr_hI[ind_Ibif_2D_local_hI]],[vus_local_bif_decr_hI[ind_Ibif_2D_local_hI]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	#scatter!([t_local_bif_incr_hI[ind_nbif_2D_local_hI]],[vus_local_bif_incr_hI[ind_nbif_2D_local_hI]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	#scatter!([t_local_bif_sl_test_hI[ind_sl_ends_hI]],[Vs_max_low_V_null_2D_hI[ind_sl_ends_hI]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)

	xaxis!("Time (ms)",(t_local_bif_decr_hI[1]-2,t_local_bif_incr_hI[end]+2))
	yaxis!("Voltage",((minimum(v_local_bif_hI)-1.5),maximum(vus_local_bif_hI)+1))

	#savefig("simu3D_local_bif_It.pdf")
end

# ╔═╡ 0f0ae675-3b24-4b60-aabc-6d8dff897a2e
begin
	plot_It_time_roi_hI = plot(sol_local_bif_hI.t,It_local_bif_full_hI,linecolor=RGB(0.65,0.29,0.65),linewidth=2,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white,label="It",size=(250,250))

	#scatter!([t_local_bif_decr_hI[ind_Ibif_2D_local_hI]],[Ibif_2D_local_hI],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation")
	#scatter!([t_local_bif_incr_hI[ind_nbif_2D_local_hI]],[It_local_bif_incr_hI[ind_nbif_2D_local_hI]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:topleft)
	#scatter!([t_local_bif_sl_test_hI[ind_sl_ends_hI[1]]],[It_local_bif_sl_test_hI[ind_sl_ends_hI[1]]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends")

	xaxis!("Time (ms)",(t_local_bif_decr_hI[1]-2,t_local_bif_incr_hI[end]+2))
	yaxis!("Current",(-2.5,4))
end

# ╔═╡ 54f0a476-c944-4d29-901b-925c8386160e
begin 
	plot(plot_volt_time_full_hI,plot_It_time_full_hI,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_full_hI.pdf")
end

# ╔═╡ 28ad4f26-6477-4c32-a5ea-77a3afbe449c
begin 
	plot(plot_volt_time_roi_hI,plot_It_time_roi_hI,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_roi_hI.pdf")
end

# ╔═╡ 681e7068-1efb-4b4f-b36c-b23f8b1f865a
begin
	I_local_bif_hI = I_local_bif
	list_stable_hI = list_stable
	list_saddle_hI = list_saddle
	list_unstable_hI = list_unstable
	list_cycleV_hI = list_cycleV


	ind_bif_2D_I_hI = findfirst(x -> x>= Ibif_2D_local_hI,I_local_bif_hI)-1
	v_bif_2D_local_bif_hI = (list_saddle_hI[ind_bif_2D_I]+list_stable_hI[ind_bif_2D_I])/2

	plot(I_local_bif_hI, list_stable_hI,label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)

	ind_real_unstable_hI_ = findall(x -> x>=0 || x<0,list_unstable_hI)
	if length(ind_real_unstable_hI_)>0
		plot!(I_local_bif_hI, list_unstable_hI,label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	end
	plot!(I_local_bif_hI, list_saddle_hI,label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)

	#plot!(I_local_bif, list_cycleV[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	ind_real_hI_ = findall(x -> x>=0 || x<0,list_cycleV_hI[:,1])
	plot!(I_local_bif_hI[ind_real_hI_], (-35).*ones(size(list_cycleV_hI[ind_real_hI_,1])),label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5)

	plot!(I_local_bif_hI, list_cycleV_hI[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,legend=:outertopright)

	scatter!([Ibif_2D_local],[v_bif_2D_local_bif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	yaxis!("V")
	xaxis!("I")

	##Trajectory
	plot!(I_stable_local_bif_decr_hI,V_stable_local_bif_decr_hI,linecolor=RGB(0.0,0.6,0.2),linewidth=5,linealpha=0.5,label="Decreasing Vus, It <= Ibif")
	plot!(I_sl_cmin_local_bif_hI,V_sl_cmin_local_bif_hI,linecolor=RGB(0.1,0.6,0.6),linewidth=5,linealpha=0.5,label="Low Vus, It > Ibif")
	plot!(I_cycle_local_bif_incr_hI,(V_cyclemax_local_bif_incr_hI./maximum(V_cyclemax_local_bif_incr_hI)).*(-35),linecolor=RGB(0.8,0.4,0.8),linewidth=5,linealpha=0.5,label="Increasing Vus, It > 0")

	ind_cycle_for_min_I_hI = findfirst(x -> x==minimum(I_cycle_local_bif_incr_hI),I_cycle_local_bif_incr_hI)
	##arrow cycle to rest
	plot!([I_cycle_local_bif_incr_hI[ind_cycle_for_min_I_hI],I_stable_local_bif_decr_hI[1]],[(-35),V_stable_local_bif_decr_hI[1]],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,linestyle=:dashdotdot,label="High Vus, It <0")
	plot!([I_stable_local_bif_decr_hI[1]-0.15,I_stable_local_bif_decr_hI[1],I_stable_local_bif_decr_hI[1]+0.15],[V_stable_local_bif_decr_hI[1]+0.15,V_stable_local_bif_decr_hI[1],V_stable_local_bif_decr_hI[1]+0.15],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,label="High Vus, It <0")

	#savefig("bifdiag2D_traj_corresp_to_vus_3D.pdf")
end

# ╔═╡ 6c976673-88f2-48ab-9af8-fd84a40144c6
begin
	ind_v_local_bif_mod_hI = findall(x->x>0,v_local_bif_hI)
	v_local_bif_mod_hI = v_local_bif_hI
	for i=1:length(ind_v_local_bif_mod_hI)
		v_local_bif_mod_hI[ind_v_local_bif_mod_hI[i]]=0
	end
	md"""Modification of the pulse height for prettier display"""
end

# ╔═╡ 5cdae710-395a-4edd-a8a7-1f6bbf214e39
begin
	

	plot(I_local_bif_hI, list_stable_hI,label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)

	if length(ind_real_unstable_hI_)>0
		plot!(I_local_bif_hI, list_unstable_hI,label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	end
	plot!(I_local_bif_hI, list_saddle_hI,label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)

	#plot!(I_local_bif, list_cycleV[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	plot!(I_local_bif_hI[ind_real_hI_], (0).*ones(size(list_cycleV_hI[ind_real_hI_,1])),label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5)

	plot!(I_local_bif_hI, list_cycleV_hI[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,legend=:outertopright)

	scatter!([Ibif_2D_local],[v_bif_2D_local_bif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	yaxis!("V")
	xaxis!("I")

	##Trajectory
	plot!(I_stable_local_bif_decr_hI,V_stable_local_bif_decr_hI,linecolor=RGB(0.0,0.6,0.2),linewidth=5,linealpha=0.5,label="Decreasing Vus, It <= Ibif")
	plot!(I_sl_cmin_local_bif_hI,V_sl_cmin_local_bif_hI,linecolor=RGB(0.1,0.6,0.6),linewidth=5,linealpha=0.5,label="Low Vus, It > Ibif")
	plot!(I_cycle_local_bif_incr_hI,(V_cyclemax_local_bif_incr_hI./maximum(V_cyclemax_local_bif_incr_hI)).*(0),linecolor=RGB(0.8,0.4,0.8),linewidth=5,linealpha=0.5,label="Increasing Vus, It > 0")

	##arrow cycle to rest
	plot!([I_cycle_local_bif_incr_hI[ind_cycle_for_min_I_hI],I_stable_local_bif_decr_hI[1]],[(0),V_stable_local_bif_decr_hI[1]],linecolor=RGB(0.6,0.7,1),linewidth=2,linealpha=0.5,linestyle=:dashdotdot,label="High Vus, It <0",)
	
	##Trajectory V(It)
	#cool color maps
	#isoluminant_cgo_70_c39_n256
	#linear_ternary_blue_0_44_c57_n256
	#linear_bmy_10_95_c78_n256
	#neon
	indbob=length(It_local_bif_hI)
	indbob0=1
	t_local_bif_normed_hI=t_local_bif_hI.-t_local_bif_hI[1]
	#plot!(It_local_bif_hI[indbob0:indbob],v_local_bif_mod_hI[indbob0:indbob],line_z=t_local_bif_normed_hI,c=:matter,linewidth=2,linealpha=1,label="Trajectory V(It)",size=(800,500),colorbar=:bottom,legend=:outertopleft)
	plot!(It_local_bif_full_hI[10000:length(It_local_bif_full_hI)],sol_local_bif_hI[1,10000:length(It_local_bif_full_hI)],line_z=sol_local_bif_hI.t[10000:length(It_local_bif_full_hI)],c=:acton,linewidth=2,linealpha=1,label="Trajectory V(It)",size=(800,500),colorbar=:bottom,legend=:outertopleft)

	#savefig("bif_diag_2D_traj_IV_hI.pdf")
end

# ╔═╡ d64952ef-b0e2-41b0-bb97-171106c565d4
md"###### 2D phase planes +real-time 3D traj display for gif"

# ╔═╡ 61504dca-fe25-47e9-87ce-6c1f408c784e
begin
	##sampling
	n_sample_local_bif_decr_hI = 35
	n_sample_local_bif_incr_hI = 15
	
	ref_ind_hI = ind_Ibif_2D_local_hI+ind_sl_ends_hI[1]-1
	
	step_sample_local_bif_decr_hI = Int(round(length(t_local_bif_hI[1:ref_ind_hI])/n_sample_local_bif_decr_hI))
	step_sample_local_bif_incr_hI = Int(round(length(t_local_bif_hI[ref_ind_hI:end])/n_sample_local_bif_incr_hI))
	
	#add offset to be sure to take the last sample of It into account 
	offset_sample_local_bif_incr_hI = length(t_local_bif_hI[ref_ind_hI:end])-step_sample_local_bif_incr_hI*n_sample_local_bif_incr_hI
	
	sample_It_local_bif_hI = []
	sample_t_local_bif_hI = []
	sample_v_local_bif_hI = []
	sample_vs_local_bif_hI = []
	sample_index_local_bif_hI = []
	
	for i=0:n_sample_local_bif_decr_hI-1
		append!(sample_It_local_bif_hI,It_local_bif_hI[1+i*step_sample_local_bif_decr_hI])
		append!(sample_t_local_bif_hI,t_local_bif_hI[1+i*step_sample_local_bif_decr_hI])
		append!(sample_v_local_bif_hI,v_local_bif_hI[1+i*step_sample_local_bif_decr_hI])
		append!(sample_vs_local_bif_hI,vs_local_bif_hI[1+i*step_sample_local_bif_decr_hI])
		append!(sample_index_local_bif_hI,1+i*step_sample_local_bif_decr_hI)
	end
	
	for i=1:n_sample_local_bif_incr_hI
		append!(sample_It_local_bif_hI,It_local_bif_hI[ref_ind_hI-1+offset_sample_local_bif_incr_hI+i*step_sample_local_bif_incr_hI])
		append!(sample_t_local_bif_hI,t_local_bif_hI[ref_ind_hI-1+offset_sample_local_bif_incr_hI+i*step_sample_local_bif_incr_hI])
		append!(sample_v_local_bif_hI,v_local_bif_hI[ref_ind_hI-1+offset_sample_local_bif_incr_hI+i*step_sample_local_bif_incr_hI])
		append!(sample_vs_local_bif_hI,vs_local_bif_hI[ref_ind_hI-1+offset_sample_local_bif_incr_hI+i*step_sample_local_bif_incr_hI])
		append!(sample_index_local_bif_hI,ref_ind_hI-1+offset_sample_local_bif_incr_hI+i*step_sample_local_bif_incr_hI)
	end
	
	##reconstruction
	recon_It_local_bif_hI = reconstruc_sample_I(sample_It_local_bif_hI,sample_index_local_bif_hI,t_local_bif_hI,n_sample_local_bif_decr_hI,n_sample_local_bif_incr_hI)

	
	
	
	md"""Samples of v,vs,vus,t and It to show the evolution of the system in a 2D phase plane"""
end

# ╔═╡ 64c765ed-8e77-4592-a734-3cb17dc1cf9e
begin
	plot(t_local_bif_hI,It_local_bif_hI,label="It")
	plot!(t_local_bif_hI,recon_It_local_bif_hI,legend=:outertopright,label="Sampled It for gif")
	#plot!(t_local_bif_decr,It_local_bif_decr,line_z=t_local_bif,linewidth=2)
	yaxis!("Current",(-1.5,4))
	xaxis!("Time (ms)")
end

# ╔═╡ 60a82124-d2a0-4be8-9d91-9c31df91ebab
begin
	v_null_1_list_hI,v_null_2_list_hI,traj_v_list_hI,traj_vs_list_hI,fp1_list_hI,fp2_list_hI,stab_list_hI = give_me_phase_planes_l(v_potential,p_local_bif_hI,recon_It_local_bif_hI,v_local_bif_hI,vs_local_bif_hI,sample_index_local_bif_hI)
end

# ╔═╡ 213141c8-8d55-44b0-aab2-8617a42f280e
begin
	plot_list_hI = [make_my_phase_plane(v_potential,v_null_1_list_hI,v_null_2_list_hI,vs_null_2D,traj_v_list_hI,traj_vs_list_hI,fp1_list_hI,fp2_list_hI,stab_list_hI,recon_It_local_bif_hI,t_local_bif_hI,sample_index_local_bif_hI,i) for i=1:length(traj_v_list_hI)]
	

	#savefig("simu2D_for_spikelat_3D_phase.pdf")
end

# ╔═╡ 5e8fa693-bea2-4143-9a59-9daa9fef8d66
length(plot_list_hI)

# ╔═╡ a60a26cd-ebb3-4a53-b7a1-bc4982b1cfb7
begin
	in_hI = 50
	plt_hI = plot_list_hI[in_hI]
	#savefig("gif_part_hI_$(in).png")
end

# ╔═╡ 225c09b6-f2bd-42b4-a11b-6384d2fd147e
begin
	plot(plt_hI)
	savefig("gif_part_hI_$(in_hI).png")
	md"""save gif plot part"""
end

# ╔═╡ 1cbaf13b-c39e-48e3-8b0e-2e12d584b0a0
md"###### More restorative ultraslow feedback"

# ╔═╡ f60a5708-727b-4cf4-8554-2ec7b657d376
begin
	##Simulation of 3D model to get Vus(t) for these parameters
	p_local_bif_resto=(4,-40.0,-38.4,-50.0,1.0,1.0,0.5,0.015,10.0,1000.0)

	tspan_local_bif_resto =(0.0,10000.0)
	prob_local_bif_resto = ODEProblem(MQIF_3D!,u0_local_bif,tspan_local_bif_resto,p_local_bif_resto,callback=cb)
	sol_local_bif_resto = solve(prob_local_bif_resto,dense=false,dtmax=1,reltol=1e-6,abstol=1e-6)
	##Must retain only one cycle of Vus, chosen to be between the 2 last spikes observed in the simulation
	ind_tmin_resto = findfirst(x -> x>=20*maximum(tspan_local_bif_resto)/100, sol_local_bif_resto.t)
	ind_spike_resto_ = findall(x -> x>=maximum(sol_local_bif_resto[3,ind_tmin_resto:end])-0.01,sol_local_bif_resto[3,ind_tmin_resto:end])
	ind_spike_resto=[]
	for i=1:length(ind_spike_resto_)
		if i>1
			if ind_spike_resto_[i]-ind_spike_resto_[i-1]>10
				append!(ind_spike_resto,ind_spike_resto_[i-1])
			end
		end
	end
	append!(ind_spike_resto,ind_spike_resto_[end]) #add last spike

	vus_local_bif_resto = sol_local_bif_resto[3,(ind_spike_resto[end-1]-1)+ind_tmin_resto:(ind_spike_resto[end]-1)+ind_tmin_resto]
	vs_local_bif_resto = sol_local_bif_resto[2,(ind_spike_resto[end-1]-1)+ind_tmin_resto:(ind_spike_resto[end]-1)+ind_tmin_resto]
	v_local_bif_resto = sol_local_bif_resto[1,(ind_spike_resto[end-1]-1)+ind_tmin_resto:(ind_spike_resto[end]-1)+ind_tmin_resto]
	t_local_bif_resto = sol_local_bif_resto.t[(ind_spike_resto[end-1]-1)+ind_tmin_resto:(ind_spike_resto[end]-1)+ind_tmin_resto]
	
	It_local_bif_resto= zeros(length(t_local_bif_resto))
	for i=1:length(t_local_bif_resto)
		It_local_bif_resto[i] = compute_I_from_vus(vus_local_bif_resto[i],p_local_bif_resto)
	end

	##Separation of the Vus cycle into 2 parts : Vus is decreasing and Vus is increasing
	ind_min_resto = findfirst(x -> x == minimum(vus_local_bif_resto),vus_local_bif_resto)

	vus_local_bif_decr_resto = vus_local_bif_resto[1:ind_min_resto]
	vs_local_bif_decr_resto = vs_local_bif_resto[1:ind_min_resto]
	v_local_bif_decr_resto = v_local_bif_resto[1:ind_min_resto]
	t_local_bif_decr_resto = t_local_bif_resto[1:ind_min_resto]

	vus_local_bif_incr_resto = vus_local_bif_resto[ind_min_resto:end]
	vs_local_bif_incr_resto = vs_local_bif_resto[ind_min_resto:end]
	v_local_bif_incr_resto = v_local_bif_resto[ind_min_resto:end]
	t_local_bif_incr_resto = t_local_bif_resto[ind_min_resto:end]

	##Calculations to determine the 2D system convergence for It = I +Ius, depending on the value of Vus
	conv_local_bif_decr_resto = zeros(length(vus_local_bif_decr_resto),2)
	conv_local_bif_incr_resto = zeros(length(vus_local_bif_incr_resto),2)
	stab_fp_local_bif_decr_resto = zeros(length(vus_local_bif_decr_resto),2)
	stab_fp_local_bif_incr_resto = zeros(length(vus_local_bif_incr_resto),2)
	conv_type_local_bif_decr_resto = zeros(length(vus_local_bif_decr_resto))
	conv_type_local_bif_incr_resto = zeros(length(vus_local_bif_incr_resto))
	It_local_bif_decr_resto = zeros(length(vus_local_bif_decr_resto))
	It_local_bif_incr_resto = zeros(length(vus_local_bif_incr_resto))

	for i=1:length(vus_local_bif_decr_resto)
		#conv_local_bif_decr[i,:],conv_type_local_bif_decr[i],stab_fp_local_bif_decr[i,:],It_local_bif_decr[i]=local_I_bif_from2D_model(p_local_bif,vus_local_bif_decr[i])
		conv_local_bif_decr_resto[i,:],conv_type_local_bif_decr_resto[i],stab_fp_local_bif_decr_resto[i,:],It_local_bif_decr_resto[i]=local_I_bif_from2D_model_decr(p_local_bif_resto,vus_local_bif_decr_resto[i],vs_local_bif_decr_resto[i],v_local_bif_decr_resto[i])
	end
	for i=1:length(vus_local_bif_incr_resto)
		#conv_local_bif_incr[i,:],conv_type_local_bif_incr[i],stab_fp_local_bif_incr[i,:],It_local_bif_incr[i]=local_I_bif_from2D_model(p_local_bif,vus_local_bif_incr[i])

conv_local_bif_incr_resto[i,:],conv_type_local_bif_incr_resto[i],stab_fp_local_bif_incr_resto[i,:],It_local_bif_incr_resto[i]=local_I_bif_from2D_model_incr(p_local_bif_resto,vus_local_bif_incr_resto[i],vs_local_bif_incr_resto[i],v_local_bif_incr_resto[i])
	end
end

# ╔═╡ dab7a530-8240-4876-82d5-f350a73f2b44
begin
	V_stable_local_bif_decr_resto,I_stable_local_bif_decr_resto,t_stable_local_bif_decr_resto,V_unstable_local_bif_decr_resto,I_unstable_local_bif_decr_resto,t_unstable_local_bif_decr_resto,V_saddle_local_bif_decr_resto,I_saddle_local_bif_decr_resto,t_saddle_local_bif_decr_resto,V_cyclemax_local_bif_decr_resto,V_cyclemin_local_bif_decr_resto,I_cycle_local_bif_decr_resto,t_cycle_local_bif_decr_resto,V_spikelat_local_bif_decr_resto,V_sl_cmin_local_bif_decr_resto,I_spikelat_local_bif_decr_resto,t_spikelat_local_bif_decr_resto = separate_local_bif_conv_full(conv_local_bif_decr_resto,conv_type_local_bif_decr_resto,It_local_bif_decr_resto,t_local_bif_decr_resto,stab_fp_local_bif_decr_resto)

V_stable_local_bif_incr_resto,I_stable_local_bif_incr_resto,t_stable_local_bif_incr_resto,V_unstable_local_bif_incr_resto,I_unstable_local_bif_incr_resto,t_unstable_local_bif_incr_resto,V_saddle_local_bif_incr_resto,I_saddle_local_bif_incr_resto,t_saddle_local_bif_incr_resto,V_cyclemax_local_bif_incr_resto,V_cyclemin_local_bif_incr_resto,I_cycle_local_bif_incr_resto,t_cycle_local_bif_incr_resto,V_spikelat_local_bif_incr_resto,V_sl_cmin_local_bif_incr_resto,I_spikelat_local_bif_incr_resto,t_spikelat_local_bif_incr_resto = separate_local_bif_conv_full(conv_local_bif_incr_resto,conv_type_local_bif_incr_resto,It_local_bif_incr_resto,t_local_bif_incr_resto,stab_fp_local_bif_incr_resto)

	V_sl_cmin_local_bif_resto = zeros(length(V_sl_cmin_local_bif_decr_resto)+length(V_sl_cmin_local_bif_incr_resto))
	I_sl_cmin_local_bif_resto = zeros(size(V_sl_cmin_local_bif_resto))
	for i=1:length(V_sl_cmin_local_bif_resto)
		if i <= length(V_sl_cmin_local_bif_decr_resto)
			V_sl_cmin_local_bif_resto[i]=V_sl_cmin_local_bif_decr_resto[i]
			I_sl_cmin_local_bif_resto[i]=I_spikelat_local_bif_decr_resto[i]
		else
			V_sl_cmin_local_bif_resto[i]=V_sl_cmin_local_bif_incr_resto[i-length(V_sl_cmin_local_bif_decr_resto)]
			I_sl_cmin_local_bif_resto[i]=I_spikelat_local_bif_incr_resto[i-length(V_sl_cmin_local_bif_decr_resto)]
		end
	end
	md"""Separation of convergence value based on their type and stability"""
end

# ╔═╡ e5340f8f-ff7b-4c51-b780-9ac8832f8ad9
begin
	Ibif_2D_local_resto = bifI2D(p_3D_to_2D(p_local_bif_resto))
	ind_Ibif_2D_local_resto = findfirst(x -> x>= Ibif_2D_local,It_local_bif_decr_resto)
	ind_nbif_2D_local_resto = length(It_local_bif_incr_resto)-(findfirst(x -> x>= 0,reverse(It_local_bif_incr_resto))-2)
	md"""Location of the SN bifurcation in the 2D model and the bifurcation to rest in the 2D model for a negative current"""
end

# ╔═╡ 44a42273-955d-4918-b074-f15ccbee0922
begin
	It_local_bif_sl_test_resto = zeros(length(It_local_bif_decr_resto[ind_Ibif_2D_local_resto:end])+length(It_local_bif_incr_resto[1:ind_nbif_2D_local_resto-1]))
	t_local_bif_sl_test_resto = zeros(length(It_local_bif_sl_test_resto))
	Vs_max_low_V_null_2D_resto = zeros(length(It_local_bif_sl_test_resto))

	for i=1:length(It_local_bif_sl_test_resto)
		if i<=length(It_local_bif_decr_resto)-(ind_Ibif_2D_local_resto-1)
			It_local_bif_sl_test_resto[i] = It_local_bif_decr_resto[i+(ind_Ibif_2D_local_resto-1)]
			t_local_bif_sl_test_resto[i] = t_local_bif_decr_resto[i+(ind_Ibif_2D_local_resto-1)]
		else
			It_local_bif_sl_test_resto[i] = It_local_bif_incr_resto[i-(length(It_local_bif_decr_resto)-(ind_Ibif_2D_local_resto-1))]
			t_local_bif_sl_test_resto[i] = t_local_bif_incr_resto[i-(length(It_local_bif_decr_resto)-(ind_Ibif_2D_local_resto-1))]
		end
		Vs_max_low_V_null_2D_resto[i]=max_2D_low_V_null(It_local_bif_sl_test_resto[i])
	end
	md"""Find the maximum of the low part of the V nullcline in the 2D phase plane for the equivalent current modifiied by the third timescale It"""
end

# ╔═╡ 2023f73c-5b1b-4da9-9175-87f2111694af
begin
	ind_sl_ends_resto = [-1]
	for i=0:length(Vs_max_low_V_null_2D_resto)-1
		if Vs_max_low_V_null_2D_resto[i+1]<=sol_local_bif_resto[2,(ind_spike_resto[end-1]-1)+ind_tmin_resto+ind_Ibif_2D_local_resto+i] && ind_sl_ends_resto[1] <0
			ind_sl_ends_resto[1] = i+1
		end
	end
	md"""Location of the end of the spike latency"""
end

# ╔═╡ b8eae3d0-bef7-4ba1-ad84-47d5334422ed
begin
	plot_volt_time_full_resto = plot(sol_local_bif_resto.t,sol_local_bif_resto[1,:],label="V",linewidth=1,linecolor=:blue,background_color=:white)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[3,:],label="Vus",linewidth=1,linecolor=RGB(1.0,0.0,0.5))
	plot!(t_local_bif_resto,(maximum(vus_local_bif_resto)+1).*ones(size(t_local_bif_resto)),fill = ((minimum(v_local_bif_resto)-1).*ones(size(t_local_bif_resto)), 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest")

		#{
	#scatter!([t_local_bif_decr_resto[ind_Ibif_2D_local_resto]],[vus_local_bif_decr_resto[ind_Ibif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	#scatter!([t_local_bif_incr_resto[ind_nbif_2D_local_resto]],[vus_local_bif_incr_resto[ind_nbif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	#scatter!([t_local_bif_sl_test_resto[ind_sl_ends_resto]],[Vs_max_low_V_null_2D_resto[ind_sl_ends_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)
	

	yaxis!("Voltage",(-45,-20))
	xaxis!("Time (ms)")
end

# ╔═╡ 7cc5009e-4942-44dd-955f-dc0193744aeb
begin
	It_local_bif_full_resto = zeros(length(sol_local_bif_resto[3,:]))
	for i=1:length(sol_local_bif_resto[3,:])
		It_local_bif_full_resto[i] = compute_I_from_vus(sol_local_bif_resto[3,i],p_local_bif_resto)
	end
	
	plot_It_time_full_resto = plot(sol_local_bif_resto.t,It_local_bif_full_resto,label="It",linewidth=1.5,linecolor=RGB(0.65,0.29,0.65))	
	#plot!(sol_local_bif_resto.t,Ibif_2D_local_resto.*ones(size(It_local_bif_full_resto)),label="2D SN bifurcation",linewidth=1.5,linecolor=RGB(0.9,0.9,0),legend=:topleft)
	#plot!(sol_local_bif_resto.t,zeros(size(It_local_bif_full_resto)),label="Bifurcation to rest",linewidth=1.5,linecolor=RGB(0,0.5,1))
	
	plot!(sol_local_bif_resto.t,(4.5).*ones(size(sol_local_bif_resto.t)),fill = ([Ibif_2D_local_resto.*ones(size(sol_local_bif_resto.t))], 0.1, RGB(0,0.5,1)),linecolor=:purple,linealpha=0.1,label="Limit cycle")
	plot!(sol_local_bif_resto.t,Ibif_2D_local_resto.*ones(size(sol_local_bif_resto.t)),fill = ([zeros(size(It_local_bif_full_resto))], 0.1, RGB(1,1,0)),linecolor=:purple,linealpha=0.1,label="Bistable")
	plot!(sol_local_bif_resto.t,zeros(size(It_local_bif_full_resto)),fill = ([(-1.5).*ones(size(It_local_bif_full_resto))], 0.1, RGB(0,1,0)),linecolor=:purple,linealpha=0.1,label="Stable")	
	
	plot!(t_local_bif_resto,(maximum(It_local_bif_full_resto)+2).*ones(size(t_local_bif_resto)),fill = ((-1.5).*ones(size(t_local_bif_resto)), 0.6, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest",legend=:topleft)
	
	
	yaxis!("Current",(-1.5,4))
	xaxis!("Time (ms)")
	#savefig("simu3D_local_bif_fullIt.pdf")
end

# ╔═╡ d62d183e-1159-42c6-b9e2-a1b0e1d17b1d
begin
	plot_volt_time_roi_resto = plot(sol_local_bif_resto.t,sol_local_bif_resto[1,:],label="V",linewidth=1,linecolor=:blue,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[3,:],label="Vus",linewidth=1,linecolor=RGB(1,0,0.5))

	scatter!([t_local_bif_decr_resto[ind_Ibif_2D_local_resto]],[vus_local_bif_decr_resto[ind_Ibif_2D_local_resto]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	scatter!([t_local_bif_incr_resto[ind_nbif_2D_local_resto]],[vus_local_bif_incr_resto[ind_nbif_2D_local_resto]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	scatter!([t_local_bif_sl_test_resto[ind_sl_ends_resto]],[Vs_max_low_V_null_2D_resto[ind_sl_ends_resto]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)

	xaxis!("Time (ms)",(t_local_bif_decr_resto[1]-2,t_local_bif_incr_resto[end]+2))
	yaxis!("Voltage",((minimum(v_local_bif_resto)-1.5),maximum(vus_local_bif_resto)+1))

	#savefig("simu3D_local_bif_It.pdf")
end

# ╔═╡ 94e04adf-937e-411f-9cba-d1fe5e76f260
begin
	plot_It_time_roi_resto = plot(sol_local_bif_resto.t,It_local_bif_full_resto,linecolor=RGB(0.65,0.29,0.65),linewidth=2,background_color=RGB(1,0.87,0.9),background_color_outside=:white,background_color_legend=:white,label="It")

	scatter!([t_local_bif_decr_resto[ind_Ibif_2D_local_resto]],[Ibif_2D_local_resto],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation")
	scatter!([t_local_bif_incr_resto[ind_nbif_2D_local_resto]],[It_local_bif_incr_resto[ind_nbif_2D_local_resto]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:topleft)
	scatter!([t_local_bif_sl_test_resto[ind_sl_ends_resto[1]]],[It_local_bif_sl_test_resto[ind_sl_ends_resto[1]]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends")

	xaxis!("Time (ms)",(t_local_bif_decr_resto[1]-2,t_local_bif_incr_resto[end]+2))
	yaxis!("Current",(-1.5,4))
end

# ╔═╡ 094b1530-1a36-45ce-860b-21afcd6eceb2
begin
	plot(plot_volt_time_full_resto,plot_It_time_full_resto,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_full_uus_.pdf")
end

# ╔═╡ 4432f66e-ecb7-4552-a434-9cb9930ba907
begin
	plot(plot_volt_time_roi_resto,plot_It_time_roi_resto,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_local_bif_roi_uus.pdf")
end

# ╔═╡ 685a9c8b-bf8c-4937-8603-984f95d6f7df
begin
	plot_volt_time_full_resto_white = plot(sol_local_bif_resto.t,sol_local_bif_resto[1,:],label="V",linewidth=1,linecolor=:blue,background_color=:white)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[3,:],label="Vus",linewidth=1,linecolor=RGB(1.0,0.0,0.5))
	plot!(t_local_bif_resto,(maximum(vus_local_bif_resto)+1).*ones(size(t_local_bif_resto)),fill = ((minimum(v_local_bif_resto)-1).*ones(size(t_local_bif_resto)), 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest")

	scatter!([t_local_bif_decr_resto[ind_Ibif_2D_local_resto]],[vus_local_bif_decr_resto[ind_Ibif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	scatter!([t_local_bif_incr_resto[ind_nbif_2D_local_resto]],[vus_local_bif_incr_resto[ind_nbif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	scatter!([t_local_bif_sl_test_resto[ind_sl_ends_resto]],[Vs_max_low_V_null_2D_resto[ind_sl_ends_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)

	yaxis!("Voltage",(-45,-20))
	xaxis!("Time (ms)")
	#savefig("bob.pdf")
end

# ╔═╡ bbc09cd1-14af-4706-91cb-3b1dd84922a9
begin
	plot_volt_time_roi_resto__ = plot(sol_local_bif_resto.t,sol_local_bif_resto[1,:],label="V",linewidth=1,linecolor=:blue,background_color=:white)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[2,:],label="Vs",linewidth=1,linecolor = :orange)
	plot!(sol_local_bif_resto.t,sol_local_bif_resto[3,:],label="Vus",linewidth=1,linecolor=RGB(1.0,0.0,0.5))
	plot!(t_local_bif_resto,(maximum(vus_local_bif_resto)+1).*ones(size(t_local_bif_resto)),fill = ((minimum(v_local_bif_resto)-1).*ones(size(t_local_bif_resto)), 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Region of interest")

	scatter!([t_local_bif_decr_resto[ind_Ibif_2D_local_resto]],[vus_local_bif_decr_resto[ind_Ibif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="2D SN bifurcation",legend=:topright)
	scatter!([t_local_bif_incr_resto[ind_nbif_2D_local_resto]],[vus_local_bif_incr_resto[ind_nbif_2D_local_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.5,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Bifurcation to rest",legend=:outertopright)
	scatter!([t_local_bif_sl_test_resto[ind_sl_ends_resto]],[Vs_max_low_V_null_2D_resto[ind_sl_ends_resto]],markershape = :circle,markersize = 3,markeralpha = 0.6,markercolor = RGB(0,0.9,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Spike latency ends",legend=:topleft)
	xaxis!("Time (ms)",(t_local_bif_decr_resto[1]-2,t_local_bif_incr_resto[end]+2))
	yaxis!("Voltage",((minimum(v_local_bif_resto)-1.5),maximum(vus_local_bif_resto)+1))
end

# ╔═╡ ccc8d44b-1795-43a9-964b-808d7faa805c
begin
	bob2 = plot(plot_volt_time_full_resto_white,plot_volt_time_roi_resto__,layout=(2,1),legend=:outertopright)
	#savefig("simu3D_highervus0.pdf")
end

# ╔═╡ bd2c688f-e87c-4084-8143-a6e63b3e286b
gr()

# ╔═╡ 1e38d252-01b7-45d8-a5bf-9463a2e372cb
begin
	I_local_bif_resto = I_local_bif
	list_stable_resto = list_stable
	list_saddle_resto = list_saddle
	list_unstable_resto = list_unstable
	list_cycleV_resto = list_cycleV


	ind_bif_2D_I_resto = findfirst(x -> x>= Ibif_2D_local_hI,I_local_bif_hI)-1
	v_bif_2D_local_bif_resto = (list_saddle_resto[ind_bif_2D_I]+list_stable_resto[ind_bif_2D_I])/2

	plot(I_local_bif_resto, list_stable_resto,label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)

	ind_real_unstable_resto_ = findall(x -> x>=0 || x<0,list_unstable_resto)
	if length(ind_real_unstable_resto_)>0
		plot!(I_local_bif_resto, list_unstable_resto,label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	end
	plot!(I_local_bif_resto, list_saddle_resto,label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)


	ind_real_resto_ = findall(x -> x>=0 || x<0,list_cycleV_resto[:,1])
	plot!(I_local_bif_resto[ind_real_resto_], (0).*ones(size(list_cycleV_resto[ind_real_resto_,1])),label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5)

	plot!(I_local_bif_resto, list_cycleV_resto[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,legend=:outertopright)

	scatter!([Ibif_2D_local],[v_bif_2D_local_bif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	yaxis!("V")
	xaxis!("I")

	##Trajectory
	#plot!(I_stable_local_bif_decr_resto,V_stable_local_bif_decr_resto,linecolor=RGB(0.0,0.6,0.2),linewidth=5,linealpha=0.5,label="Decreasing Vus, It <= Ibif")
	#plot!(I_sl_cmin_local_bif_resto,V_sl_cmin_local_bif_resto,linecolor=RGB(0.1,0.6,0.6),linewidth=5,linealpha=0.5,label="Low Vus, It > Ibif")
	#plot!(I_cycle_local_bif_incr_resto,(V_cyclemax_local_bif_incr_resto./maximum(V_cyclemax_local_bif_incr_resto)).*(0),linecolor=RGB(0.8,0.4,0.8),linewidth=5,linealpha=0.5,label="Increasing Vus, It > 0")

	##arrow SN bif
	#plot!([I_stable_local_bif_decr_resto[end],I_sl_cmin_local_bif_resto[1]],[V_stable_local_bif_decr_resto[end],V_sl_cmin_local_bif_resto[1]],linecolor=RGB(1,0.8,0.0),linewidth=3,linealpha=0.5,linestyle=:dashdotdot,label="Low Vus, It -> Ibif")
	
	
	ind_v_local_bif_mod_resto = findall(x->x>0,v_local_bif_resto)
	v_local_bif_mod_resto = v_local_bif_resto
	for i=1:length(ind_v_local_bif_mod_resto)
		v_local_bif_mod_resto[ind_v_local_bif_mod_resto[i]]=0
	end

	
	t_local_bif_normed_resto=t_local_bif_resto.-t_local_bif_resto[1]
	plot!(It_local_bif_resto,v_local_bif_mod_resto,line_z=t_local_bif_normed_resto,c=:matter,linewidth=2,linealpha=1,label="Trajectory V(It)",size=(800,500),colorbar=:bottom,legend=:outertopleft)
	

	#savefig("bif_diag_2D_traj_IV_uus_.pdf")
end

# ╔═╡ bd9344f4-6cb3-4ead-9b7d-f826cff49eaa
md"#### Patterns simulations"

# ╔═╡ bc64c297-bff1-49ac-b80a-0e8b0eb37157
gr()

# ╔═╡ a55fb1de-3244-49a2-9eba-d878b976de15
function simulate_step(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc)
	tspan0=(t0,t0_step)
	
	if i_exc == false
		u00 =[-40,-40,-40]
	else
		u00 =[-10,-20,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_3D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
	
	tspan1=(t0_step,tf_step)
	u01 =[sol0[1,end],sol0[2,end],sol0[3,end]]
	p1=change_I_(Istep,p)
	
	prob1 = ODEProblem(MQIF_3D!,u01,tspan1,p1,callback=cb)
	sol1 = solve(prob1,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
	
	tspan2=(tf_step,tf)
	u02 =[sol1[1,end],sol1[2,end],sol1[3,end]]
	p2=change_I_(I0,p)
	
	prob2 = ODEProblem(MQIF_3D!,u02,tspan2,p2,callback=cb)
	sol2 = solve(prob2,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)

	I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))	
	t_I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))
	sim =zeros(3,length(sol0.t)+length(sol1.t)+length(sol2.t))
	
	for i=1:length(I)
		if i<=length(sol0.t)
			I[i] = I0
			t_I[i] = sol0.t[i]
			sim[1,i] = sol0[1,i]
			sim[2,i] = sol0[2,i]
			sim[3,i] = sol0[3,i]
		else
			if i<=length(sol0.t)+length(sol1.t)
				I[i] = Istep
				t_I[i] = sol1.t[i-length(sol0.t)]
				sim[1,i] = sol1[1,i-length(sol0.t)]
				sim[2,i] = sol1[2,i-length(sol0.t)]
				sim[3,i] = sol1[3,i-length(sol0.t)]
			else
				I[i] = I0
				t_I[i] = sol2.t[i-length(sol0.t)-length(sol1.t)]
				sim[1,i] = sol2[1,i-length(sol0.t)-length(sol1.t)]
				sim[2,i] = sol2[2,i-length(sol0.t)-length(sol1.t)]	
				sim[3,i] = sol2[3,i-length(sol0.t)-length(sol1.t)]	
			end
		end
	end
	return t_I,I,sim
end

# ╔═╡ 4a80d584-073b-48b0-a889-c22b17b9de0d
function simu_pulse(t0,tf,t_spike,n_spike,p,I0,Ispike,spike_duration,i_exc)
	tspan0=(t0,t_spike)
	
	if i_exc == false
		u00 =[-40,-40,-40]
	else
		u00 =[-10,-20,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_3D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
	u0s=[sol0[1,end],sol0[2,end],sol0[3,end]]
	
	v = []
	vs = []
	vus = []
	t = []
	I = []
	append!(v,sol0[1,:])
	append!(vs,sol0[2,:])
	append!(vus,sol0[3,:])
	append!(t,sol0.t)
	append!(I,I0.*ones(length(sol0.t)))
	
	for i=1:n_spike
		tspan_spike=(i*t_spike,i*t_spike+spike_duration)
		p_spike = change_I_(Ispike,p)
		
		probs = ODEProblem(MQIF_3D!,u0s,tspan_spike,p_spike,callback=cb)
		sols = solve(probs,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
		u0r = [sols[1,end],sols[2,end],sols[3,end]]
		
		append!(v,sols[1,:])
		append!(vs,sols[2,:])
		append!(vus,sols[3,:])
		append!(t,sols.t)
		append!(I,Ispike.*ones(length(sols.t)))
		
		if (i+1)*t_spike<=tf
			tspan_rest=(i*t_spike+spike_duration,(i+1)*t_spike)
		else
			tspan_rest=(i*t_spike+spike_duration,tf)
		end
		
		probr = ODEProblem(MQIF_3D!,u0r,tspan_rest,p0,callback=cb)
		solr = solve(probr,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
		
		u0s =[solr[1,end],solr[2,end],solr[3,end]]
		
		append!(v,solr[1,:])
		append!(vs,solr[2,:])
		append!(vus,solr[3,:])
		append!(t,solr.t)
		append!(I,I0.*ones(length(solr.t)))
	end
	
	sim=zeros(3,length(t))
	sim[1,:] = v
	sim[2,:] = vs
	sim[3,:] = vus
	
	return t,I,sim
end

# ╔═╡ ec7bc3f2-4ece-4a5f-9ada-bf926a607b49
function simu_ud_(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,i_exc)
	tspan0=(t0,t_u[1])

	if i_exc == false
		u00 =[-40,-40,-40]
	else
		u00 =[-10,-20,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_3D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
	u0u=[sol0[1,end],sol0[2,end],sol0[3,end]]
	
	v = []
	vs = []
	vus = []
	t = []
	I = []
	append!(v,sol0[1,:])
	append!(vs,sol0[2,:])
	append!(vus,sol0[3,:])
	append!(t,sol0.t)
	append!(I,I0.*ones(length(sol0.t)))
	
	for i=1:length(t_u)
		## up
		tspan_u=(t_u[i],t_u[i]+spike_duration)
		p_u = change_I_(Ispike,p)
		
		probu = ODEProblem(MQIF_3D!,u0u,tspan_u,p_u,callback=cb)
		solu = solve(probu,dtmax=0.1,DP5(),reltol=1e-5,abstol=1e-5)
		u0r = [solu[1,end],solu[2,end],solu[3,end]]
		
		append!(v,solu[1,:])
		append!(vs,solu[2,:])
		append!(vus,solu[3,:])
		append!(t,solu.t)
		append!(I,Ispike.*ones(length(solu.t)))
		
		## rest
		tspan_rest=(t_u[i]+spike_duration,t_d[i])
		
		probr = ODEProblem(MQIF_3D!,u0r,tspan_rest,p0,callback=cb)
		solr = solve(probr,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
		
		u0d =[solr[1,end],solr[2,end],solr[3,end]]
		
		append!(v,solr[1,:])
		append!(vs,solr[2,:])
		append!(vus,solr[3,:])
		append!(t,solr.t)
		append!(I,I0.*ones(length(solr.t)))
		
		## down
		tspan_d =(t_d[i],t_d[i]+spike_duration)
		p_d = change_I_(-Ispike,p)
		
		probd = ODEProblem(MQIF_3D!,u0d,tspan_d,p_d,callback=cb)
		sold = solve(probd,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)
		u0r = [sold[1,end],sold[2,end],sold[3,end]]
		
		append!(v,sold[1,:])
		append!(vs,sold[2,:])
		append!(vus,sold[3,:])
		append!(t,sold.t)
		append!(I,(-Ispike).*ones(length(sold.t)))
		
		
		if t_d[i]+spike_duration < tf
			if i<length(t_u)
				tspan_rest=(t_d[i]+spike_duration,t_u[i+1])
			else
				tspan_rest=(t_d[i]+spike_duration,tf)
			end

			probr = ODEProblem(MQIF_3D!,u0r,tspan_rest,p0,callback=cb)
			solr = solve(probr,dtmax=0.1,DP5(),reltol=1e-6,abstol=1e-6)

			u0u =[solr[1,end],solr[2,end],solr[3,end]]

			append!(v,solr[1,:])
			append!(vs,solr[2,:])
			append!(vus,solr[3,:])
			append!(t,solr.t)
			append!(I,I0.*ones(length(solr.t)))
		end
	end
	
	sim=zeros(3,length(t))
	sim[1,:] = v
	sim[2,:] = vs
	sim[3,:] = vus
	
	return t,I,sim
end

# ╔═╡ 7df52d7b-e145-462e-888b-a7e6a72e22d1
function subplot_simu_ud(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,both,i_exc)

	t,I,sim = simu_ud_(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,i_exc)
	
	plot_ud = plot(t,I,linecolor=RGB(0.7,0,0.1),label="I(t)     ")
	yaxis!("I")	
	xaxis!("Time (ms)")	
	
	if both==false
		plot_v = plot(t,sim[1,:],label="V(t)")
		lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:])])
		lim2 = maximum(sim[1,:])
		yaxis!("V",(lim1,lim2))	
		xaxis!("Time (ms)")	

		plot_vs = plot(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vs",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_vus = plot(t,sim[3,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_ud,plot_v,plot_vs,plot_vus,layout=@layout([a{0.2h};b ;c;d]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t,sim[1,:],label="V(t)")
				plot!(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	

		sub = plot(plot_ud,plot_v,layout=@layout([a{0.2h};b ]),linewidth = 1.5,legend=:outertopright,size=(800,500))
	end
	
	return sub
end

# ╔═╡ ffbc4063-57cd-46ba-81f3-5bb75edcba55
function subplot_simu_step(t0,tf,step_duration,p,I0,Istep,both,i_exc)
	t0_step = round((tf-t0)/2)-round(step_duration/2)
	tf_step = round((tf-t0)/2)+round(step_duration/2)
	
	t_I,I,sim = simulate_step(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc)
	
	plot_step = plot(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)     ")
	yaxis!("I")	
	xaxis!("Time (ms)")	
	
	plot_step_z = plot(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)     ")
	yaxis!("I")	
	xaxis!("Time (ms)",(t0_step -tf*0.01,tf_step +tf*0.01))	
	title!("------------ Zoom on pulse ------------ ",title_location=:centre,titlefontsize=16,titlefontcolor=RGB(0,0.4,0.95))
	
	if both==false
		plot_v = plot(t_I,sim[1,:],label="V(t)")
		lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:])])
		lim2 = maximum(sim[1,:])
		yaxis!("V",(lim1,lim2))	
		xaxis!("Time (ms)")	

		plot_vs = plot(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vs",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_vus = plot(t_I,sim[3,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_step,plot_v,plot_vs,plot_vus,layout=@layout([a{0.2h};b ;c;d]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	
		
		plot_zoom = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				#plot!(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)",(t0_step -tf*0.01,tf_step +tf*0.01))	
		
		title_zoom = plot(title = "Zoom", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:right,titlefontsize=14,titlefontcolor=RGB(0,0.4,0.95))

		sub = plot(plot_step,plot_v,plot_step_z,plot_zoom,layout=@layout([a{0.1h};b;c{0.1h};d ]),linewidth = 1.5,legend=:outertopright,size=(800,800))
		#sub = plot(plot_step,plot_v,title_zoom,plot_step_z,plot_zoom,layout=@layout([A{0.05h};B;C{0.01h};D{0.05h};E ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ c4275b0c-9ecc-45c6-9f5b-5f01e93ccb2f
function subplot_simu_spike(t0,tf,t_spike,p,I0,Ispike,spike_duration,both,i_exc)
	n_spike = floor((tf-t0)/t_spike)-1
	
	t,I,sim = simu_pulse(t0,tf,t_spike,n_spike,p,I0,Ispike,spike_duration,i_exc)
	
	plot_spike = plot(t,I,linecolor=RGB(0.7,0,0.1),label="I(t)     ")
	yaxis!("I")	
	xaxis!("Time (ms)")	
	
	if both==false
		plot_v = plot(t,sim[1,:],label="V(t)")
		lim1 = -60
		lim2 = 9.9
		yaxis!("V",(lim1,lim2))	
		xaxis!("Time (ms)")	

		plot_vs = plot(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vs",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_vus = plot(t,sim[3,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_spike,plot_v,plot_vs,plot_vus,layout=@layout([a{0.1h};b ;c;d]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t,sim[1,:],label="V(t)")
				plot!(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	

		sub = plot(plot_spike,plot_v,layout=@layout([a{0.1h};b ]),linewidth = 1.5,legend=:outertopright,size=(800,500))
	end
	
	return sub
end

# ╔═╡ d5b2ba0e-af1d-4c1a-baf0-918e02c0015b
begin
	t_u=[250,2000,7000]
	delta_ud = [600,1000,2500]
	t_d=t_u+delta_ud
	sub_ud =subplot_simu_ud(0.0,10000.0,t_u,t_d,p,0.4,20.0,3.0,true,false)
	#savefig("MQIF_3D-simuud.pdf")
end

# ╔═╡ 4a8fb301-09b5-4d64-b4fa-77da14e58304
begin
	sub =subplot_simu_step(0.0,10000.0,1200.0,p,0.0,20.0,true,false)	
	savefig("MQIF_3D-simustep_.pdf")
end

# ╔═╡ 324f30b4-742b-4407-b11c-fb4cf4a143e6
begin
	#t,I,sim = simu_pulse(0.0,100.0,20.0,5,p,0.0,30.0,0.5,true)
	subs = subplot_simu_spike(0.0,10000.0,1500.0,p,1.0,20.0,1,true,false)
	#savefig("MQIF_3D-simupulse.pdf")
end

# ╔═╡ 032578a3-4845-4df4-84f7-0e949c923b16


# ╔═╡ 67f78f73-b368-4d3c-aa53-74a03fa7b9aa
md"#### Additional gif // 4D"

# ╔═╡ b818ce21-6463-46a7-8357-71d8f63b92e9
function compute_Ius(vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	
	if length(vus)>1
		Ius = zeros(length(vus))
		for i=1:length(vus)
			Ius[i]=-gus*(vus[i]-vus0)^2
		end
	else
		Ius=-gus*(vus-vus0)^2
	end
	return Ius
end

# ╔═╡ 617ab2b0-ce6a-4432-ac4f-0c293bf3ee78
function compute_vus_null_It(p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	
	if I >=0
		vus = vus0 + sqrt((I)/gus)
	else
		vus = NaN
	end
	
	return vus
end

# ╔═╡ 00865f63-c16e-419a-85f6-0fb666016380
function compute_vs_extremas(vus,p)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	
	vus_bif_n = compute_vus_null_It(p)
	vsmax = vs0
	vsmin = vs0
	if vus<=vus_bif_n
		Ius = compute_Ius(vus,p)
		if I+Ius >=0
			vsmax = vs0 + sqrt((I+Ius)/gs)
			vsmin = vs0 - sqrt((I+Ius)/gs)
		end			
	end
	return vsmax,vsmin
end

# ╔═╡ 86841f5b-6281-486e-ba5c-87cc4fcdfc2e
function gif_well(v,vs,vus,t,p,i)
	I,v0,vs0,vus0,C,gf,gs,gus,ts,tus = p
	
	vus_ = collect(range(minimum(vus),stop=maximum(vus),length=100))
	vs_max = zeros(length(vus_))
	vs_min = zeros(length(vus_))
	for j=1:length(vus_)
		vs_max[j],vs_min[j] = compute_vs_extremas(I,vus_[j],p)
	end
	plt = plot(vus_,vs_max,linecolor=:red,label = "Extrema ordinate, vus = $(Int(round(vus_level)))")
	plot!(vus_,vs_min,linecolor=:red,label="Extrema ordinate, vus = $(Int(round(vus_level)))")
	plot!(collect(range(minimum(vus),stop=maximum(vus),length=length(vus_))),vs_max,fill=(vs_min, 0.5, :orange),linealpha=0,label="Excitability well, vus = $(Int(round(vus_level)))",size=(500,400))
	plot!(vus[1:i],vs[1:i],line_z=t,c=:cyclic_mygbm_30_95_c78_n256)	
	yaxis!("Vs",(-80,0))
	xaxis!("Vus")
	title!("  C. ",titlelocation=:left,titlefontsize=18)
	return plt
end

# ╔═╡ 34f78910-4e5c-4589-8036-72d93f4e7e95
function gif_It_time(It,t,p,i)
	Ilim1 = maximum(It)
	Ilim2 = minimum(It)
	tlim1 = maximum(t)
	tlim2 = minimum(t)
	plt = plot(t[1:i],It[1:i],linecolor="red",label="It(t)")
	yaxis!("Current",(Ilim2,Ilim1))
	xaxis!("Time")
	title!("  D. ",titlelocation=:left,titlefontsize=18)
	return plt
end

# ╔═╡ eee076ba-8752-4445-95af-35598194825c
function make_my_time_plot(v_potential,v_null_1_list,v_null_2_list,vs_null_2D, v,vs,vus,vuus,fp1_list,fp2_list,stab_list,recon_It,t,sample_index,p3D,i)
	
	p2D = p_3D_to_2D(p3D)
	
	if i<length(sample_index)
		plot_test = plot(t[1:sample_index[i+1]],v[1:sample_index[i+1]],linecolor="blue",label="V(t)")
	plot!(t[1:sample_index[i+1]],vs[1:sample_index[i+1]],linecolor="orange",label="Vs(t)")
	plot!(t[1:sample_index[i+1]],vus[1:sample_index[i+1]],linecolor="black",label="Vus(t)")

	else
		t_ = zeros(sample_index[i]+1)
		v_ = zeros(sample_index[i]+1)
		vs_ = zeros(sample_index[i]+1)
		vus_ = zeros(sample_index[i]+1)
		t_[1:end-1]= t[1:sample_index[i]]
		v_[1:end-1]= v[1:sample_index[i]]
		vs_[1:end-1]= vs[1:sample_index[i]]
		vus_[1:end-1]= vus[1:sample_index[i]]
		t_[end] = t[end]  
		v_[end] = v[end]  
		vs_[end] = vs[end]  
		vus_[end] = vus[end]   
		plot_test = plot(t_,v_,linecolor="blue",label="V(t)")
		plot!(t_,vs_,linecolor="orange",label="Vs(t)")
		plot!(t_,vus_,linecolor="black",label="Vus(t)")
	end
	
	xaxis!("Time")
	yaxis!("Voltage",(-80,0))
	title!("  B. ",titlelocation=:left,titlefontsize=18)
	return plot_test
end

# ╔═╡ 2c87d0d2-c4bf-40c7-94c8-d96f79a31ab4
function find_period(vus,t)
	ind_t_min = findfirst(x -> x>=maximum(t)/10,t)
	ind_vus_spike_ = findall(x -> x>=(maximum(vus[ind_t_min:end])-0.03),vus[ind_t_min:end])
	
	delta=zeros(length(ind_vus_spike_)-1)
	for i=1:length(ind_vus_spike_)-1
		delta[i]=ind_vus_spike_[i+1]-ind_vus_spike_[i]
	end
	
	ind_delta=findfirst(x -> x==maximum(delta),delta)
	delta_=reverse(delta[1:ind_delta])
	to_add_=[]
	for i=2;ind_delta
		if delta_[i] >10
			break
		end
		append!(to_add_,delta_[i])
	end
	delta_period = Int(maximum(delta)+sum(to_add_))
	return delta_period,ind_t_min
end

# ╔═╡ 9d68f4bc-5cf2-4689-8930-3322b51885a0
function find_ind_period(v,vus,t)
	delta_period,ind_t_min = find_period(vus,t)
	ind_spike = findfirst(x -> x>=maximum(v[ind_t_min:end])*0.99,v[ind_t_min:end])
	ind_period1 =ind_spike+ind_t_min-1
	ind_period2 = findfirst(x -> x>=maximum(v[ind_period1+delta_period-10:end])*0.99,v[ind_period1+delta_period-10:end])+ind_period1+delta_period-10-1
	#ind_period2 =ind_period1+delta_period
	return ind_period1,ind_period2
end

# ╔═╡ cf21811a-11ac-49dd-bbf0-2eae13542cff
begin
	ind1,ind2 = find_ind_period(sol_local_bif_resto[1,:],sol_local_bif_resto[3,:],sol_local_bif_resto.t)
	t_local_resto = sol_local_bif_resto.t[ind1:ind2]
	v_local_resto = sol_local_bif_resto[1,ind1:ind2]
	vs_local_resto = sol_local_bif_resto[2,ind1:ind2]
	vus_local_resto = sol_local_bif_resto[3,ind1:ind2]
	
	Ius_local_resto = compute_Ius(vus_local_resto,p_local_bif_resto)
	It_local_resto = p_local_bif_resto[1].+Ius_local_resto
	
	plot(t_local_resto,It_local_resto)
end

# ╔═╡ 2284d700-0ff4-4056-86ed-f8842e027f3b


# ╔═╡ af2290b0-ee41-4709-b1c2-b5455c9ae93c
function sample_It(n_sample_decr,n_sample_incr,t,v,vs,vus,It)
	
	ref_ind = findfirst(x -> x == maximum(It),It)
	
	step_sample_decr = Int(round(length(t[1:ref_ind])/n_sample_decr))
	step_sample_incr = Int(round(length(t[ref_ind:end])/n_sample_incr))
	
	#add offset to be sure to take the last sample of It into account 
	offset_sample_incr = length(t[ref_ind:end])-step_sample_incr*n_sample_incr
	
	sample_It = []
	sample_t = []
	sample_v = []
	sample_vs = []
	sample_index = []
	
	for i=0:n_sample_decr-1
		append!(sample_It,It[1+i*step_sample_decr])
		append!(sample_t,t[1+i*step_sample_decr])
		append!(sample_v,v[1+i*step_sample_decr])
		append!(sample_vs,vs[1+i*step_sample_decr])
		append!(sample_index,1+i*step_sample_decr)
	end
	
	for i=1:n_sample_incr
		append!(sample_It,It[ref_ind-1+offset_sample_incr+i*step_sample_incr])
		append!(sample_t,t[ref_ind-1+offset_sample_incr_hI+i*step_sample_incr])
		append!(sample_v,v[ref_ind-1+offset_sample_incr+i*step_sample_incr])
		append!(sample_vs,vs[ref_ind-1+offset_sample_incr+i*step_sample_incr])
		append!(sample_index,ref_ind-1+offset_sample_incr_hI+i*step_sample_incr)
	end
	
	##reconstruction
	recon_It = reconstruc_sample_I(sample_It,sample_index,t,n_sample_decr,n_sample_incr)

	
	
	
	md"""Samples of v,vs,vus,t and It to show the evolution of the system in a 2D phase plane"""
end

# ╔═╡ Cell order:
# ╟─792f7770-95e7-11eb-23a5-45e952a738da
# ╠═bd35d5de-7544-11eb-2023-9345c9da0b0f
# ╠═05913e12-7545-11eb-2eec-875820f00015
# ╠═8da532ce-95e7-11eb-1d64-a37503684718
# ╟─97c1b30e-95e7-11eb-2047-27179c477642
# ╠═009dbc70-7546-11eb-1ae1-8dbac638d9e9
# ╠═79a51975-67f2-4992-8b7c-d12176075249
# ╠═9967a94d-ef13-45a0-b936-48bca8a64526
# ╠═262a2b7f-7f42-4002-b91b-4bea4dea074c
# ╟─da7703c2-95ee-11eb-2ce4-399430fe71e5
# ╟─2b0581dc-f37c-4a44-bbe1-a16944a5e1ab
# ╠═756351fa-5648-4cfb-a804-3fe5bda48cb9
# ╟─616560f0-95e7-11eb-0ae3-ddcbdbe0be8a
# ╠═63a575ae-7555-11eb-0cb8-43a1d04dabde
# ╟─4fd67ae0-79e5-11eb-0278-838411a5c092
# ╟─20f2eb20-79e3-11eb-268e-bbd0237d7af2
# ╟─a9a5bb00-79e3-11eb-0ecf-2397fdc5fa21
# ╟─939c4bf2-7ac2-11eb-1bc3-1bd1119508d6
# ╟─d60571d0-7b60-11eb-0e78-93ffbb86bc10
# ╟─e6f8e440-7b60-11eb-07b0-cba945d97277
# ╟─310f1bd2-7b61-11eb-3e3e-65e1334b589f
# ╟─6d492f43-e36b-4725-82da-905d02e83d73
# ╟─e4594047-eb27-44ec-aea2-c2df1defba77
# ╟─3444011c-122f-44ce-b53a-ddf04e46f1dc
# ╟─841b5bf9-8b23-469d-bf4d-31764ef5d900
# ╟─e2e956e1-e100-4315-83f1-557a2634c3f9
# ╟─b2850211-edf0-4f48-aeff-942b23948ad8
# ╟─1abb4876-62be-4df3-bd98-555ef0c3fd2a
# ╠═5d5293ed-0060-4615-b068-3acf16b2c1b5
# ╠═6ffd59cc-13e2-41c7-8ac6-4031438f8780
# ╟─d5115900-95e7-11eb-34e4-ff02328e7981
# ╟─5702cbc0-948e-11eb-2b59-1d53d6a15304
# ╟─5f26b7f2-9b35-4390-aa0e-29926da44b82
# ╠═e5b70336-c1c2-4ec4-858f-33a8d481599b
# ╟─6bdb76b4-7b9e-4657-8b00-213e8453c8e4
# ╟─14e62c79-d9b3-4f2d-9a8c-f16e2f0fd6b1
# ╟─e1037f19-64a6-42f8-bef5-5a470b8e438d
# ╟─3a8758c2-84e3-47dc-bdac-52fe7bec5800
# ╟─af7f680e-5209-46a3-a6b3-5323d1b66729
# ╟─df44c90b-868a-4da5-833d-2a508a3b8e91
# ╟─35b7b20d-7a32-4bdf-a682-9e2f72eba110
# ╟─d60d0341-703b-4908-9017-996b4d222c37
# ╟─7736932f-d9be-49e3-be14-888b780f3cca
# ╟─a3edd5c2-00bc-4925-a18d-e4b51b54f96e
# ╟─19188423-c2e3-403e-9e16-7ccfeaca4f45
# ╟─2131d4e4-5c7a-4b6e-b6be-05d5e78ce8f9
# ╟─4fcc4093-0b2b-4645-82ba-b6bee0213695
# ╟─8aa1c774-8026-4cac-bd58-f9ac89c5390c
# ╟─fb9b0595-c11c-4be4-89df-6c07e0f9a02d
# ╟─56d892b0-95e7-11eb-1b32-6db8dc815056
# ╟─78f557b0-7545-11eb-2810-9bd1d2f85b33
# ╟─8a1d2b20-7546-11eb-2d7a-494ae827de36
# ╠═9173ad8e-7546-11eb-1c3e-51967fc5a4ce
# ╟─bb367180-7546-11eb-361d-53889265d7be
# ╟─c0b88210-7546-11eb-1343-31a62c6e60b6
# ╠═ddc72af0-7546-11eb-3d4e-f9ff68884993
# ╠═e83ae550-95e7-11eb-11d3-0bb50817f70e
# ╠═0dc4c81e-7547-11eb-3439-3d4a5750c2cb
# ╠═93032350-9618-11eb-0a3f-95b734f703bb
# ╠═56bf7d90-7547-11eb-054a-6ff0654fd396
# ╟─dc5c0f90-79e3-11eb-34f6-bfda108d916f
# ╠═3cce0240-960c-11eb-205e-bb2765422b75
# ╠═e07ed7c0-7b63-11eb-3b80-8d22fcd902e3
# ╟─a7ae4c37-3547-44ee-a9cc-8d3d3f5df221
# ╟─c28bddf0-79e5-11eb-1c82-014cf00fadd9
# ╠═4cb043a9-20aa-488a-8020-d9383d83b774
# ╠═55b08c31-a827-4ea1-bea6-e5e9ecf7d639
# ╠═b598cd61-be3c-4e09-85b2-0518a1b85054
# ╠═25c6a491-28a0-4e19-829f-7b790fc97d5f
# ╠═ae0b4f6d-7269-44e5-9acc-4c618875b084
# ╠═f1111e66-63d0-4a3e-bb59-49b729d0281c
# ╠═31f42cfb-1e08-4120-b48d-d195788eadcb
# ╠═5f8af83d-577f-4a6b-acb4-0b67dbbdf113
# ╠═362b5e4f-f9ba-4e41-9164-e22dcb1a6731
# ╟─81513370-79eb-11eb-2257-f39b0b5789af
# ╟─eb8f7590-9489-11eb-353d-dff43f1b970c
# ╟─d8db632e-9490-11eb-27d0-73a81a9eeab6
# ╟─8e4240b0-95e8-11eb-1d42-cdae28c3935c
# ╟─0a06f401-42db-430d-acce-6cb609d654d9
# ╠═12d22342-f5ec-4a88-ab9e-5c4d76f0f92a
# ╠═7ba8807a-3f9f-4919-a267-70539fa955a4
# ╠═37e4ba5e-33c3-4083-a483-f20b111c217d
# ╠═3d22c39f-d723-4c99-87b6-6b60f861d816
# ╠═ded847d4-16e7-4e77-ab39-8f2ab65d66f1
# ╟─381d77d7-0723-47e9-b724-fed87d2f93a9
# ╟─8aae255d-a738-4236-a4a9-15ad137eff0a
# ╟─f30c88e3-1bc0-4eaa-96a0-37aa7ce9a25c
# ╟─b19620bd-f7eb-4e53-b4e7-4c7b42675d31
# ╟─20fa3ea1-7a0e-43ea-947b-ce78a8929e83
# ╟─13feea03-5093-4c99-bae0-d88da15a2925
# ╟─ec0e5f62-2114-43b3-931f-324cb45046a5
# ╟─781d46f6-0a91-408c-a2bb-4c78cbb91150
# ╟─2f225060-6402-4fae-bced-192328ecf4f7
# ╟─e2726064-c076-4b44-98d7-f80ce2d572aa
# ╠═7ebed900-95e8-11eb-0e9d-d7bd5f6748b1
# ╟─82582c60-95e8-11eb-1a32-3998de2baa18
# ╟─2497dbf0-95ea-11eb-2f36-dfcad08b86da
# ╟─f512e270-95ef-11eb-0ba6-91e5e91857c2
# ╟─357243f0-95f1-11eb-37be-b58f6942268f
# ╟─243d6310-95ee-11eb-3a7e-c19eb0280760
# ╟─7938db90-960e-11eb-1815-b7ddd5a1853c
# ╟─342d6f3e-95ee-11eb-1e02-f15eeda830f5
# ╟─63230109-92cf-44bc-8b95-4489053cb92a
# ╟─ff4c958c-1a63-4a3f-a3b8-5e8aad4e5992
# ╟─ce60c8d0-960e-11eb-17d1-a1855ad6a7c4
# ╠═84475afd-4453-4c85-83b3-951198d85fba
# ╟─ed9077c2-b76b-40e7-8d11-7c1a0b144df6
# ╟─e0ff370f-d23e-44bb-b14a-9b6f1972421d
# ╟─e914d015-7dd2-4df0-af05-d298fbd262f5
# ╠═0f28b311-cc76-40a8-8cea-cb357e450834
# ╟─a6cfb516-8f6c-43a2-a25d-672994a754eb
# ╟─37ebe428-4359-4203-a2e9-93cdcc821c58
# ╟─03e6717a-4959-4c54-bd2b-50c78e7282c7
# ╟─b1fb80f0-f9b2-470e-a628-c14a0b6a4dbe
# ╟─c550a58f-615e-442e-862a-baca5f741227
# ╟─cb6382ce-9f1d-40b6-99b9-e9f3f5133eed
# ╟─7c3fb682-0460-4014-934a-be1d2c0f7acc
# ╟─3b8c3790-aec6-45ca-90b9-5e48758c70e1
# ╟─40097d95-9548-4291-8101-a74c6f405440
# ╟─cbff4346-cef2-429e-86b9-d79e92f1389f
# ╟─f6388d7c-b102-4cc8-a7a9-fc9aa24e13bf
# ╟─b0bf9322-dba5-4084-99df-f6b77e39273f
# ╟─3a3c7729-d8c3-4b4a-92b9-020da890792e
# ╟─05107ac0-69f7-493a-a78a-64f84167f547
# ╟─5dfa49e0-b3ec-4563-ba63-dc945cdff46d
# ╟─8726355e-4969-4eb4-b435-511b507bd6bd
# ╟─2c643ebe-95f7-11eb-3174-59c327320e07
# ╟─ee87bfb5-68b0-40e7-8b8a-db2c2a4130a0
# ╠═b3318e49-5cee-4b9e-82b2-d9daacf0f480
# ╠═20175480-95f8-11eb-025f-ed13fb0faec1
# ╟─afd86ff3-3457-4e1d-964b-b5d20ab65da9
# ╠═83f85190-9619-11eb-2b29-cb5032c8b81d
# ╟─b3a6d3be-4a94-40fb-aedf-a0fa3b5a8bf3
# ╟─37486ae0-95ee-11eb-1912-b1e2ea244526
# ╟─243e2330-95ea-11eb-10c8-7707ea0c7402
# ╟─87273901-f5b8-4db1-8810-28c7f5e8d60a
# ╟─d5a7342d-d27d-4fef-b2ff-3a768a4a5b7c
# ╟─64cac8b3-1d84-4152-bef7-aa4f7ee3052f
# ╟─71b19713-3e2c-463d-a44a-fcd7c46a12a1
# ╟─09266a70-9c17-4f5f-a1e9-8917bbbcc0de
# ╟─835e25e2-a335-40e5-af56-d75f0db67c8a
# ╟─7f469400-632a-40c6-96c7-2b674324ed60
# ╟─9f0566ed-09a1-4c4f-ac7c-dbc1b2452f68
# ╠═9c5553f3-751b-4b28-abd1-c06ca67a80aa
# ╟─17656593-6727-447e-97b0-d7f67faae4e2
# ╠═86854a97-660b-4b10-a769-341a17ac71ae
# ╟─a9cb8c87-dfd9-46bf-b4b9-ab9c004c2667
# ╠═378a31a7-4c78-4f56-b225-73fa47921fa4
# ╟─1788b635-ae11-4617-958e-975a4d106c80
# ╟─7d998b8e-6eae-43f6-8e67-684d6b6df2f8
# ╟─9cd94904-cb24-48b6-9576-b73c2d3aa391
# ╟─997f5fe7-d3e6-4968-89f3-ab1348ffa1d4
# ╟─c5a958b2-434b-4d66-b5a6-63fb293fa295
# ╟─6824efca-1041-4e55-a660-db38962f330b
# ╟─3ba6624a-2219-498b-a435-a9d49fc6db14
# ╟─465e64e3-ccc0-4951-b3c9-37710919e8a7
# ╠═0c2c6410-3d15-4ed9-9adb-568d986d0cad
# ╟─f46c788e-dd72-4dde-95bc-f5ff15ded9bc
# ╟─579383f0-ef43-42e7-b771-9a0b2dbf1f66
# ╟─933002ad-57f5-4573-8f18-fa6c77bef9d3
# ╟─db7220d3-d245-44e2-86a4-b5e6f8c9ef90
# ╟─298033d0-6e82-4be3-9d87-70a10f734d2b
# ╟─5ea20fe5-3c23-4cb0-95bb-04dd8672a08c
# ╟─3daae77a-ef7a-4766-b354-cb06598e2047
# ╟─a516dbb5-57b5-42ef-a10b-1358694b9daa
# ╟─abdf6c7f-e761-46e7-8f5b-d0d2d7a2d285
# ╠═771927b2-814a-4422-8d78-c4020cc89cda
# ╠═d3db0388-7f11-41ef-887e-46187b3e4b53
# ╟─733c2729-3558-4660-a059-93d9d8220ad4
# ╟─053be29e-3c94-4ad2-b89f-7dafd86b06f9
# ╟─b88a07e8-97ce-476b-8e0d-f5de53fda620
# ╟─2517434c-1aa9-4cb5-8ed5-9b5943632538
# ╟─9064707b-c621-44be-94f8-e4e15242e447
# ╠═ef6bdd3f-00bc-43ff-b8b6-0d266abdfd57
# ╠═f61a1af4-5f02-4fe8-9c51-912d8ff2de11
# ╟─4802d83b-9051-46ef-94aa-b60b83a0fcde
# ╠═7358098d-88a6-42a7-b401-a3088c98b8fd
# ╟─bd4faf19-cbee-46fd-9fbf-84d30da07cbb
# ╟─9c063b76-e3ef-40b1-9691-aadcdd46b82b
# ╟─88c95811-b6cc-4512-9756-c796ddcb6ea2
# ╟─676cd0a8-4f10-407d-88e5-d777acce0758
# ╟─c03e6b55-6d09-4c71-987e-99227614b2d1
# ╟─dd3e2e91-9dbc-45ad-88db-2fee49fa0ece
# ╟─7b1ec15f-f926-49f5-bbd8-fedaec8432c9
# ╟─262f3030-053a-457f-bebc-32ddaa19b154
# ╟─e2fea39c-6a9f-4e68-9be2-ffae78c9a807
# ╟─840c7361-a9b1-4fdb-ab67-8116428da226
# ╟─b8fd420a-da3d-485a-a98d-c2271b24b6fa
# ╠═e30ec005-016a-45a6-abf8-3d6eda1de911
# ╠═d13d81a7-c7e8-4096-833a-cda6f5f99c32
# ╠═f2ce1d7b-8392-4f20-b567-3593889b73c8
# ╠═3d6b6513-397b-4369-895b-7eee2678ec37
# ╟─8b41be33-fad7-4f3f-94b2-dc4be766102c
# ╟─25e64d32-8603-44d4-8a9a-e3396c9ed17a
# ╟─82ee2ffb-5375-4271-ad7b-f15c5f520d4a
# ╟─73cd72e6-c389-4c40-803e-c4513995d438
# ╟─e9b8a4ca-05cb-4aad-b2fd-3b22c3fc23c2
# ╟─e7637146-0496-495b-892a-08668620e00b
# ╟─07b5937f-f95a-49e4-8c66-8c7ae791b7d9
# ╟─5dcaee47-e015-4592-8144-979d2616804e
# ╟─0f0ae675-3b24-4b60-aabc-6d8dff897a2e
# ╟─54f0a476-c944-4d29-901b-925c8386160e
# ╟─28ad4f26-6477-4c32-a5ea-77a3afbe449c
# ╟─681e7068-1efb-4b4f-b36c-b23f8b1f865a
# ╟─6c976673-88f2-48ab-9af8-fd84a40144c6
# ╟─5cdae710-395a-4edd-a8a7-1f6bbf214e39
# ╟─d64952ef-b0e2-41b0-bb97-171106c565d4
# ╠═61504dca-fe25-47e9-87ce-6c1f408c784e
# ╟─64c765ed-8e77-4592-a734-3cb17dc1cf9e
# ╠═60a82124-d2a0-4be8-9d91-9c31df91ebab
# ╠═213141c8-8d55-44b0-aab2-8617a42f280e
# ╠═5e8fa693-bea2-4143-9a59-9daa9fef8d66
# ╠═a60a26cd-ebb3-4a53-b7a1-bc4982b1cfb7
# ╠═225c09b6-f2bd-42b4-a11b-6384d2fd147e
# ╟─1cbaf13b-c39e-48e3-8b0e-2e12d584b0a0
# ╠═f60a5708-727b-4cf4-8554-2ec7b657d376
# ╟─dab7a530-8240-4876-82d5-f350a73f2b44
# ╟─e5340f8f-ff7b-4c51-b780-9ac8832f8ad9
# ╟─44a42273-955d-4918-b074-f15ccbee0922
# ╟─2023f73c-5b1b-4da9-9175-87f2111694af
# ╠═b8eae3d0-bef7-4ba1-ad84-47d5334422ed
# ╟─7cc5009e-4942-44dd-955f-dc0193744aeb
# ╟─d62d183e-1159-42c6-b9e2-a1b0e1d17b1d
# ╟─94e04adf-937e-411f-9cba-d1fe5e76f260
# ╠═094b1530-1a36-45ce-860b-21afcd6eceb2
# ╟─4432f66e-ecb7-4552-a434-9cb9930ba907
# ╟─685a9c8b-bf8c-4937-8603-984f95d6f7df
# ╟─bbc09cd1-14af-4706-91cb-3b1dd84922a9
# ╠═ccc8d44b-1795-43a9-964b-808d7faa805c
# ╠═bd2c688f-e87c-4084-8143-a6e63b3e286b
# ╠═1e38d252-01b7-45d8-a5bf-9463a2e372cb
# ╟─bd9344f4-6cb3-4ead-9b7d-f826cff49eaa
# ╠═bc64c297-bff1-49ac-b80a-0e8b0eb37157
# ╟─a55fb1de-3244-49a2-9eba-d878b976de15
# ╟─4a80d584-073b-48b0-a889-c22b17b9de0d
# ╟─ec7bc3f2-4ece-4a5f-9ada-bf926a607b49
# ╟─7df52d7b-e145-462e-888b-a7e6a72e22d1
# ╟─ffbc4063-57cd-46ba-81f3-5bb75edcba55
# ╟─c4275b0c-9ecc-45c6-9f5b-5f01e93ccb2f
# ╠═d5b2ba0e-af1d-4c1a-baf0-918e02c0015b
# ╠═4a8fb301-09b5-4d64-b4fa-77da14e58304
# ╠═324f30b4-742b-4407-b11c-fb4cf4a143e6
# ╟─032578a3-4845-4df4-84f7-0e949c923b16
# ╟─67f78f73-b368-4d3c-aa53-74a03fa7b9aa
# ╠═b818ce21-6463-46a7-8357-71d8f63b92e9
# ╟─617ab2b0-ce6a-4432-ac4f-0c293bf3ee78
# ╟─00865f63-c16e-419a-85f6-0fb666016380
# ╟─86841f5b-6281-486e-ba5c-87cc4fcdfc2e
# ╟─34f78910-4e5c-4589-8036-72d93f4e7e95
# ╟─eee076ba-8752-4445-95af-35598194825c
# ╟─2c87d0d2-c4bf-40c7-94c8-d96f79a31ab4
# ╠═9d68f4bc-5cf2-4689-8930-3322b51885a0
# ╠═cf21811a-11ac-49dd-bbf0-2eae13542cff
# ╠═2284d700-0ff4-4056-86ed-f8842e027f3b
# ╟─af2290b0-ee41-4709-b1c2-b5455c9ae93c
