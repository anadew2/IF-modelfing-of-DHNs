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

# ╔═╡ eb47cf0d-f73a-4719-afda-3188b2998e52
begin
	using DifferentialEquations
	using LinearAlgebra
	using Plots
	using Statistics
	using ProgressBars
end

# ╔═╡ eb93ed6e-a1f2-11eb-00a6-d5bbd921a4be
md" #### Packages"

# ╔═╡ b9b842e4-7861-4bc5-9956-2db1d09c1c7d
md" #### Problem parameters"

# ╔═╡ d49d2169-a568-4319-9287-d13044c30bd4
p=(70,-40.0,-38.4,-19.0,-50,1.0,1.0,0.5,0.1,0.01,10.0,100.0,1000.0)

# ╔═╡ 97629afb-5289-422a-bb17-5c1a0afebc00
function compute_Iuus(vuus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	if length(vuus)>1
		Iuus = zeros(length(vuus))
		for i=1:length(vuus)
			Iuus[i]=-guus*(vuus[i]-vuus0)^2
		end
	else
		Iuus=-guus*(vuus-vuus0)^2
	end
	return Iuus
end

# ╔═╡ 94979c24-7cc6-4c67-a0fa-13e0ea2c51f7
function compute_Ius(vus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
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

# ╔═╡ 3ae16e48-9afb-4d31-9123-cf6edb092b51
function compute_It(vus,vuus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	Iuus = compute_Iuus(vuus,p)
	Ius = compute_Ius(vus,p)
	
	It = I.+Ius.+Iuus
	
	return It	
end

# ╔═╡ 9a5b6a2b-4080-4128-aabb-c1846eab7a40
function compute_vuus_null_It(vus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	Ius = compute_Ius(vus,p)
	if I+Ius >=0
		vuus = vuus0 + sqrt((I+Ius)/guus)
	else
		vuus = NaN
	end
	
	return vuus
end

# ╔═╡ c7930bda-a02c-4ec1-a975-83a61de2592a
function compute_vs_extremas(vus,vuus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	vuus_bif_n = compute_vuus_null_It(vus,p)
	vsmax = vs0
	vsmin = vs0
	if vuus<=vuus_bif_n
		Iuus = compute_Iuus(vuus,p)
		Ius = compute_Ius(vus,p)
		if I+Ius+Iuus >=0
			vsmax = vs0 + sqrt((I+Ius+Iuus)/gs)
			vsmin = vs0 - sqrt((I+Ius+Iuus)/gs)
		end			
	end
	return vsmax,vsmin
end

# ╔═╡ 042c10bf-0755-48f0-be92-f7451765a816
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

# ╔═╡ d21074d8-6462-4ec1-b2a2-dd46028f675c
function change_gs_(gs,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==8
			p2[i]=gs
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 139362f6-4a22-48a0-8d35-219deb803372
function change_gus_(gus,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==9
			p2[i]=gus
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 58cbec76-b9bd-4e3c-bd7c-295bfd036f5b
function change_guus_(guus,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==10
			p2[i]=guus
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 25553d60-2155-45ad-a6cd-d7bfb247e7c8
function change_vuus0_(vuus0,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==5
			p2[i]=vuus0
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ be4c0e03-30c1-4b24-b491-4c8217ee27c8
function change_vus0_(vus0,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==4
			p2[i]=vus0
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ c5bc2248-002b-4da9-9640-42db8bca0ee6
function change_vs0_(vs0,p)
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

# ╔═╡ 774f5bdf-863d-48f7-8115-625a06c171a3
function find_global_bif(p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	numerator = gf*gs*(v0-vs0)^2 +gf*gus*(v0-vus0)^2 +gf*guus*(v0-vuus0)^2 -gs*gus*(vs0-vus0)^2 -gs*guus*(vs0-vuus0)^2 -gus*guus*(vus0-vuus0)^2
	Ibif = numerator /(gf-gs-gus-guus)
end

# ╔═╡ 78d2b9c6-9953-49dd-81a2-e1f22c887d30
Ibif = find_global_bif(p)

# ╔═╡ 1b9fd190-0614-4730-b9f4-b48619abd92d
function find_local_bif(p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	numerator = gf*gs*(v0-vs0)^2 
	Ibif = numerator /(gf-gs)
end

# ╔═╡ 4d024409-a4f4-4edd-9039-988e3c37c09d
md" #### Nullclines functions"

# ╔═╡ 7a533a24-1a7d-4502-bfc9-075ed377866d
function testplane(x,y)
	z=x
	return z
end

# ╔═╡ 455c242e-c3fc-4b2f-ae99-a15897ff4613
function testplane(x,y,w)
	z=x.+w
	return z
end

# ╔═╡ 9d5d2669-dba9-4c00-9b5c-5e27428d05b1
function Vnullcline1(v,vus,vuus,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	if (gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(vuus-vuus0)^2 + I) >= 0
		vs1 = vs0 + sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(vuus-vuus0)^2 + I)/gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ acec6140-f6b1-446f-8401-777d78c648e5
function Vnullcline2(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	if (gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(vuus-vuus0)^2 + I) >= 0
		vs2 = vs0 - sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(vuus-vuus0)^2 + I)/gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 507962d3-b77c-478a-b39b-ebb75c0e98c5
function V_Vs_nullcline1(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	#compute US component
	US = gus*(vus-vus0)^2
	UUS = guus*(vuus-vuus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US-UUS)

	if delta >= 0 #there is 1 (2) fixed point
		vs1 = ((gf*v0)-(gs*vs0)+sqrt(delta))/(gf-gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 4e0b1fd0-107d-42c1-9e9c-c0fb0e7de3df
function V_Vs_nullcline2(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	#compute US component
	US = gus*(vus-vus0)^2
	UUS = guus*(vuus-vuus0)^2
	delta = gf*gs*(v0-vs0)^2 - (gf-gs)*(I-US-UUS)

	if delta >= 0 #there is 1 (2) fixed point
		vs2 = ((gf*v0)-(gs*vs0)-sqrt(delta))/(gf-gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ cbae9324-dd13-4769-a381-8089c1114a98
function V_Vuus_nullcline1(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	if v == vuus
		if gf*(v-v0)^2 -gus*(vus-vus0)^2 - guus*(v-vuus0)^2 +I >=0
			vs1 = vs0 + sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(v-vuus0)^2 +I)/gs)
		else
			vs1 = NaN
		end
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 50a836a6-c023-49ab-8aad-79983dd7c2b8
function V_Vuus_nullcline2(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	if v == vuus
		if gf*(v-v0)^2 -gus*(vus-vus0)^2 - guus*(v-vuus0)^2 +I >=0
			vs2 = vs0 - sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(v-vuus0)^2 +I)/gs)
		else
			vs2 = NaN
		end
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 1cfc2491-8e71-42ac-a8c3-5001adbe3597
function V_Vus_nullcline1(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	if v==vus
		if gf*(v-v0)^2 -gus*(vus-vus0)^2 - guus*(v-vuus0)^2 +I >=0
			vs1 = vs0 + sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(v-vuus0)^2 +I)/gs)
		else
			vs1 = NaN
		end
	else
		vs1 = NaN
	end
	
	return vs1
end

# ╔═╡ 5f609d44-ea73-4c26-95e8-1158e0d55e44
function V_Vus_nullcline2(v,vus,vuus,p)
	#fixed vus for 3D display 
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	if v==vus
		if gf*(v-v0)^2 -gus*(vus-vus0)^2 - guus*(v-vuus0)^2 +I >=0
			vs2 = vs0 - sqrt((gf*(v-v0)^2 -gus*(vus-vus0)^2 -guus*(v-vuus0)^2 +I)/gs)
		else
			vs2 = NaN
		end
	else
		vs2 = NaN
	end
	
	return vs2
end

# ╔═╡ d1c422a3-82ab-42b9-8a17-be40b3f8664b
md"##### 2D model"

# ╔═╡ 1c92ffd6-3284-444e-bbe6-22103df0dbff
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

# ╔═╡ 9fa7fb2b-2d9c-4775-9a02-f98bb27133d9
function p_4D_to_2D(p4D)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p4D
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

# ╔═╡ 8d571a7a-78ac-4841-8ae1-6c4a6adcc650
function jacobian_2D(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	J = zeros(2,2)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[2,1] = 1
	J[2,2] = -1

	return J
end

# ╔═╡ fe3f6a22-5450-4983-bfa2-d9440f72e4f3
function Vnullcline1_2D(v,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gf*(v-v0)^2 + I) >= 0
		vs1 = vs0 + sqrt((gf*(v-v0)^2 + I)/gs)
	else
		vs1 = NaN
	end

	return vs1
end

# ╔═╡ 8172c9ee-117c-468a-9c37-71f1de0708dd
function Vnullcline2_2D(v,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gf*(v-v0)^2 + I) >= 0
		vs2 = vs0 - sqrt((gf*(v-v0)^2 + I)/gs)
	else
		vs2 = NaN
	end

	return vs2
end

# ╔═╡ 07633dbf-9633-45ef-84a9-3f839122386f
function Vnullcline1_2D_(vs,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gs*(vs-vs0)^2 - I) >= 0
		v1 = v0 + sqrt((gs*(vs-vs0)^2 - I)/gf)
	else
		v1 = NaN
	end

	return v1
end

# ╔═╡ 120b0cc1-d150-4720-97dd-9e473117bd38
function Vnullcline2_2D_(vs,p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D

	if (gs*(vs-vs0)^2 - I) >= 0
		v2 = v0 - sqrt((gs*(vs-vs0)^2 - I)/gf)
	else
		v2 = NaN
	end

	return v2
end

# ╔═╡ 02a4331e-ad81-49fa-bcc5-b2131af42691
function Vsnullcline_2D(v,p_2D)
	return v
end

# ╔═╡ f1226b78-526e-46fb-96d4-a7da54c6126e
function bifI2D(p_2D)
	I,v0,vs0,C,gf,gs,ts = p_2D
	Ibif = (gf*gs*(v0-vs0)^2)/(gf-gs)
	return Ibif
end

# ╔═╡ 9d241efb-d465-4ba1-987a-f509ba0b81aa
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

# ╔═╡ 41fcfdfe-0f90-4362-9c95-c8ff6201f643
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

# ╔═╡ cb1e1b89-f2c3-4f3a-ba56-dd1c4f6d2c98
md" #### ODE problem functions"

# ╔═╡ 7530943d-039a-4467-9c7d-8f565c5e5298
function MQIF_4D!(du,u,p,t)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 - gus*(u[3]-vus0)^2 - guus*(u[4]-vuus0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
	du[3] = (u[1]-u[3])/tus
	du[4] = (u[1]-u[4])/tuus
end

# ╔═╡ ae088c53-ab45-4dd8-8e01-384940acb491
begin
	Vmax = 30 
	Vr = -40
	Vsr = -25
	DVusr = 3
	DVUusr = 3

	md"""Voltages used for reset"""
end

# ╔═╡ 0d58c293-3a2b-4518-ac0f-f26620ca4e44
function spike(x)  # spikes when spike(x) goes from negative to positive
    (x[1] - Vmax)
end

# ╔═╡ 411e46b8-8ebc-4ef7-bb0c-688020416e5b
# event when event_f(u,t) == 0
function condition(x,t,integrator) #
    spike(x)
end

# ╔═╡ 36f2a546-2716-494c-ad21-09c21f78cb8a
function reset!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
	x[3] = x[3] + DVusr
	x[4] = x[4] + DVUusr
end

# ╔═╡ 44d2c44b-44a3-490d-a933-12f74d0ef5b1
# when condition == 0 and upcrossing (from negative to positive)
function affect!(integrator)
    reset!(integrator.u)
end

# ╔═╡ 47bccea1-4cb9-40ac-a5dc-62322512ce1a
cb   = ContinuousCallback(condition,affect!,nothing)

# ╔═╡ be306caa-d310-4a27-bfcf-16b5cff5590b
md" #### Simulations"

# ╔═╡ 006de583-9601-432c-95a7-a3f786b0c937
begin
	u0=[-40,-40,-40,-40]
	tspan = (0.0,12000.0)
	prob = ODEProblem(MQIF_4D!,u0,tspan,p,callback=cb)
	sol = solve(prob,dense=false)
end

# ╔═╡ 61645246-ac2d-4d11-9252-ab4929841135
begin
	v_=sol[1,:]
	Iuus = compute_Iuus(sol[4,:],p)
	Ius = compute_Ius(sol[3,:],p)
	md"""Compute Iuus"""
	It=Iuus.+p[1].+Ius
	
	ind_I_0_ = findall(x -> x<-1 && x>-5, It)
	ind_I_0 = []
	for i=1:length(ind_I_0_)-1
		if ind_I_0_[i+1]-ind_I_0_[i]>30 
			append!(ind_I_0,ind_I_0_[i])
		end
	end
	append!(ind_I_0,ind_I_0_[end])
	md"""Find when the current cross the axis I=0"""
end

# ╔═╡ 20c01555-448a-48b8-ad43-345ad4c82c34
function find_period(vuus,t)
	ind_t_min = findfirst(x -> x>=maximum(t)/10,t)
	ind_vuus_spike_ = findall(x -> x>=(maximum(vuus[ind_t_min:end])-0.03),vuus[ind_t_min:end])
	
	delta=zeros(length(ind_vuus_spike_)-1)
	for i=1:length(ind_vuus_spike_)-1
		delta[i]=ind_vuus_spike_[i+1]-ind_vuus_spike_[i]
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

# ╔═╡ a7234266-83ee-4f54-b3a8-3d2ada1c2a00
function find_ind_period(v,vuus,t)
	delta_period,ind_t_min = find_period(vuus,t)
	ind_spike = findfirst(x -> x>=maximum(v[ind_t_min:end])*0.99,v[ind_t_min:end])
	ind_period1 =ind_spike+ind_t_min-1
	ind_period2 = findfirst(x -> x>=maximum(v[ind_period1+delta_period-10:end])*0.99,v[ind_period1+delta_period-10:end])+ind_period1+delta_period-10-1
	#ind_period2 =ind_period1+delta_period
	return ind_period1,ind_period2
end

# ╔═╡ 7ab02730-0c95-445f-a1f3-f3d645776f9b
md" #### Phase Plane"

# ╔═╡ 3e13d7ce-1056-4562-86ad-412598a9eadf
function compute_nullclines(v,vus_level,vuus,p)
	V1_null_(a,b) = Vnullcline1(a,vus_level,b,p)
	V2_null_(a,b) = Vnullcline2(a,vus_level,b,p)
	
	inter_V_Vs1 = zeros(length(vuus))
	inter_V_Vs2 = zeros(length(vuus))
	
	inter_V_Vus1 = zeros(size(vuus))
	inter_V_Vus2 = zeros(size(vuus))
	for i=1:length(vuus)
		inter_V_Vs1[i] = V_Vs_nullcline1(NaN,vus_level,vuus[i],p)
		inter_V_Vs2[i] = V_Vs_nullcline2(NaN,vus_level,vuus[i],p)
		inter_V_Vus1[i] = V_Vus_nullcline1(vus_level,vus_level,vuus[i],p)
		inter_V_Vus2[i] = V_Vus_nullcline2(vus_level,vus_level,vuus[i],p)
	end

	inter_V_Vuus1 = zeros(size(v))
	inter_V_Vuus2 = zeros(size(v))
	for i=1:length(v)
		inter_V_Vuus1[i] = V_Vuus_nullcline1(v[i],vus_level,v[i],p)
		inter_V_Vuus2[i] = V_Vuus_nullcline2(v[i],vus_level,v[i],p)
	end
	
	return V1_null_,V2_null_,inter_V_Vs1,inter_V_Vs2,inter_V_Vuus1,inter_V_Vuus2, inter_V_Vus1,inter_V_Vus2
end

# ╔═╡ 0ef6fa1a-efe7-441c-8360-9afba715e066
function plot_3D_pp(v,vus_level,vuus,p)
		V1_null_,V2_null_,inter_V_Vs1,inter_V_Vs2,inter_V_Vuus1,inter_V_Vuus2, inter_V_Vus1,inter_V_Vus2 = compute_nullclines(v,vus_level,vuus,p)
		
	pp1 =plot(v,vuus,V1_null_,st=:surface,c=:sun,colorbar_entry=false,xlabel="V", ylabel="Vuus",zlabel="Vs",camera=(60,20)) #yellow-orange
	plot!(v,vuus,V2_null_,st=:surface,c=:sun,colorbar_entry=false,label="dV/dt = 0") #yellow-orange
	plot!([inter_V_Vs1,inter_V_Vs2],[vuus,vuus],[inter_V_Vs1,inter_V_Vs2],linecolor=RGB(1,0.5,0.7),linewidth=3,label="dVs/dt = 0") #red 
	plot!([v,v],[v,v],[inter_V_Vuus1,inter_V_Vuus2],linecolor=RGB(0.31,0.66,1),linewidth=3,label="dVuus/dt = 0") #blue
	plot!([vus_level.*ones(size(vuus)),vus_level.*ones(size(vuus))],[vuus,vuus],[inter_V_Vus1,inter_V_Vus2],linecolor=:green,linewidth=3,label="dVus/dt = 0") #blue
	title!("V nullcline and its intersections")
	
	pp1_side = plot(pp1,camera=(0,0),legend=:outertopright)
	pp1_oside = plot(pp1,camera=(30,20),legend=:outertopright,colorbar=true)
	#pp1_supp = plot()
	
	#plot(pp1,pp1_side,pp1_oside,pp1_supp,layout=(2,2))
	
	vs_lim1 = maximum(inter_V_Vuus1)
	vs_lim2 = minimum(inter_V_Vuus2)
	return pp1,vs_lim1,vs_lim2
end

# ╔═╡ db5df20a-fea3-4f51-827f-5ff657bd21a4
function plot_3D_pp_gif(v,vus_level,vuus,p)
		V1_null_,V2_null_,inter_V_Vs1,inter_V_Vs2,inter_V_Vuus1,inter_V_Vuus2, inter_V_Vus1,inter_V_Vus2 = compute_nullclines(v,vus_level,vuus,p)
		
	pp1 =plot(v,vuus,V1_null_,st=:surface,c=:sun,colorbar_entry=false,xlabel="V", ylabel="Vuus",zlabel="Vs",camera=(60,20),size=(500,400)) #yellow-orange
	plot!(v,vuus,V2_null_,st=:surface,c=:sun,colorbar_entry=false,label="dV/dt = 0") #yellow-orange
	plot!([inter_V_Vs1,inter_V_Vs2],[vuus,vuus],[inter_V_Vs1,inter_V_Vs2],linecolor=RGB(1,0.5,0.7),linewidth=3,label="dVs/dt = 0") #red 
	plot!([v,v],[v,v],[inter_V_Vuus1,inter_V_Vuus2],linecolor=RGB(0.31,0.66,1),linewidth=3,label="dVuus/dt = 0") #blue
	plot!([vus_level.*ones(size(vuus)),vus_level.*ones(size(vuus))],[vuus,vuus],[inter_V_Vus1,inter_V_Vus2],linecolor=:green,linewidth=3,label="dVus/dt = 0") #blue
	title!("V nullcline and its intersections")
	
	pp1_side = plot(pp1,camera=(0,0),legend=:outertopright)
	pp1_oside = plot(pp1,camera=(30,20),legend=:outertopright,colorbar=true)
	#pp1_supp = plot()
	
	#plot(pp1,pp1_side,pp1_oside,pp1_supp,layout=(2,2))
	
	vs_lim1 = maximum(inter_V_Vuus1)
	vs_lim2 = minimum(inter_V_Vuus2)
	return pp1
end

# ╔═╡ 17856e2b-0541-42fb-b393-9858be2411be
md" #### Excitability well Plane"

# ╔═╡ 35a023dc-16c4-4054-903a-1659474ac2a5
md"""Study of one period of the reference simulation"""

# ╔═╡ e784cdff-fab4-4344-883f-6fcc17af8ad0
function gif_well(v,vs,vus,vuus,t,p,i)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	vus_level = vus[i]
	vuus_ = collect(range(minimum(vuus),stop=maximum(vuus),length=100))
	vs_max = zeros(length(vuus_))
	vs_min = zeros(length(vuus_))
	for j=1:length(vuus_)
		vs_max[j],vs_min[j] = compute_vs_extremas(vus_level,vuus_[j],p)
	end
	plt = plot(vuus_,vs_max,linecolor=:red,label = "Extrema ordinate, vus = $(Int(round(vus_level)))")
	plot!(vuus_,vs_min,linecolor=:red,label="Extrema ordinate, vus = $(Int(round(vus_level)))")
	plot!(collect(range(minimum(vuus),stop=maximum(vuus),length=length(vuus_))),vs_max,fill=(vs_min, 0.5, :orange),linealpha=0,label="Excitability well, vus = $(Int(round(vus_level)))",size=(500,400))
	plot!(vuus[1:i],vs[1:i],line_z=t,c=:cyclic_mygbm_30_95_c78_n256)	
	yaxis!("Vs",(vslim21,vslim11))
	xaxis!("Vuus")
	title!("  C. ",titlelocation=:left,titlefontsize=18)
	return plt
end

# ╔═╡ 80acbb79-335f-4f40-90a5-fe3a03c937aa
function gif_time(v,vs,vus,vuus,t,p,i)
	vlim1 = maximum([maximum(v),maximum(vs),maximum(vus),maximum(vuus)])
	vlim2 = minimum([minimum(v),minimum(vs),minimum(vus),minimum(vuus)])
	tlim1 = maximum(t)
	tlim2 = minimum(t)
	plt = plot(t[1:i],v[1:i],linecolor="blue",label="V(t)")
	plot!(t[1:i],vs[1:i],linecolor="orange",label="Vs(t)")
	plot!(t[1:i],vus[1:i],linecolor="pink",label="Vus(t)")
	plot!(t[1:i],vuus[1:i],linecolor="black",label="Vuus(t)")
	plot!(t[1:i],p[3]ones(size(t[1:i])),linecolor="red",label="Vs0")
	plot!(t[1:i],p[4]ones(size(t[1:i])),linecolor="magenta",label="Vus0")
	plot!(t[1:i],p[5]ones(size(t[1:i])),linecolor="green",label="Vuus0",legend=:outertopright)
	yaxis!("Voltage",(vlim2,vlim1))
	xaxis!("Time")
	return plt
end

# ╔═╡ 14e0a258-d957-45cc-8543-9203fea533be
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

# ╔═╡ 67122038-ec98-4c9f-9d18-51f02e28a183
function last_max_It(v,vs,vus,vuus,p)
	Iuus = compute_Iuus(vuus,p)
	Ius = compute_Ius(vus,p)
	
	It=Iuus.+p[1].+Ius
	ind_it_min = findfirst(x->x==minimum(It),It)
	ind_I_0_ = findall(x -> x<-1 && x>-5, It[1:ind_it_min])
	
	ind_I_0=[]
	append!(ind_I_0,ind_I_0_[end])
	return ind_I_0[1]
end

# ╔═╡ 4f69ac27-61e1-4ec2-92ce-4bd3b9d49bf8
function reconstruc_sample_I(sample_It,sample_index,t,n_sample_acc,n_sample_slo,n_sample_neg)
	
	recon_It = ones(size(t))
	for i=2:n_sample_acc
		ind1 = sample_index[i-1]
		ind2 = sample_index[i]
		recon_It[ind1:ind2].=sample_It[i-1]
	end
	
	indtrans1 = sample_index[n_sample_acc]
	indtrand2 = sample_index[1+n_sample_acc]
	recon_It[indtrans1:indtrand2].=sample_It[n_sample_acc]
	
	for i=2:n_sample_slo
		ind1 = sample_index[i-1+n_sample_acc]
		ind2 = sample_index[i+n_sample_acc]
		recon_It[ind1:ind2].=sample_It[i-1+n_sample_acc]
	end
	
	indtrans1 = sample_index[n_sample_acc+n_sample_slo]
	indtrand2 = sample_index[1+n_sample_acc+n_sample_slo]
	recon_It[indtrans1:indtrand2].=sample_It[n_sample_acc+n_sample_slo]
	
	for i=2:n_sample_neg
		ind1 = sample_index[i-1+n_sample_acc+n_sample_slo]
		ind2 = sample_index[i+n_sample_acc+n_sample_slo]
		recon_It[ind1:ind2].=sample_It[i-1+n_sample_acc+n_sample_slo]
	end
	return recon_It
end

# ╔═╡ 33117e03-5d89-47fa-90a5-e89dc2071512
function sample_It(n_sample_acc,n_sample_slo,n_sample_neg,t,v,vs,vus,vuus,It)
	##sampling
	#n_sample_decr = 30
	#n_sample_incr = 7
	
	#ref_ind = ind_Ibif_2D_local+ind_sl_ends[1]-1
	
	ref_ind1 = findfirst(x -> x == maximum(It),It)
	ref_ind2 = last_max_It(v,vs,vus,vuus,p)[1]
	
	step_sample_acc = Int(round(length(t[1:ref_ind1])/n_sample_acc))
	step_sample_slo = Int(round(length(t[ref_ind1:ref_ind2])/n_sample_slo))
	step_sample_neg = Int(round(length(t[ref_ind2:end])/n_sample_neg))
	
	#add offset to be sure to take the last sample of It into account 
	offset_sample_neg = length(t[ref_ind2:end])-step_sample_neg*n_sample_neg
	
	sample_It = []
	sample_t = []
	sample_v = []
	sample_vs = []
	sample_vus = []
	sample_vuus = []
	sample_index = []
	
	for i=0:n_sample_acc-1
		append!(sample_It,It[1+i*step_sample_acc])
		append!(sample_t,t[1+i*step_sample_acc])
		append!(sample_v,v[1+i*step_sample_acc])
		append!(sample_vs,vs[1+i*step_sample_acc])
		append!(sample_vus,vus[1+i*step_sample_acc])
		append!(sample_vuus,vuus[1+i*step_sample_acc])
		append!(sample_index,1+i*step_sample_acc)
	end
	
	for i=0:n_sample_slo-1
		append!(sample_It,It[ref_ind1+i*step_sample_slo])
		append!(sample_t,t[ref_ind1+i*step_sample_slo])
		append!(sample_v,v[ref_ind1+i*step_sample_slo])
		append!(sample_vs,vs[ref_ind1+i*step_sample_slo])
		append!(sample_vus,vus[ref_ind1+i*step_sample_slo])
		append!(sample_vuus,vuus[ref_ind1+i*step_sample_slo])
		append!(sample_index,ref_ind1+i*step_sample_slo)
	end
	
	for i=1:n_sample_neg
		append!(sample_It,It[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_t,t[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_v,v[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_vs,vs[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_vus,vus[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_vuus,vuus[ref_ind2-1+offset_sample_neg+i*step_sample_neg])
		append!(sample_index,ref_ind2-1+offset_sample_neg+i*step_sample_neg)
	end
	
	##reconstruction
	recon_It = reconstruc_sample_I(sample_It,sample_index,t,n_sample_acc,n_sample_slo,n_sample_neg)
	
	return recon_It,sample_index
end

# ╔═╡ 254f1c1d-5f79-4078-b686-51f2d1f6f1a3
function give_me_v_null_1(v_potential,p4D,It,sample_index,i)
	v_null_2D_1 = zeros(size(v_potential))
	p_2D = change_I_(It[sample_index[i]],p_4D_to_2D(p4D))
	for j=1:length(v_potential)
		if It[sample_index[i]] >=0
			v_null_2D_1[j] = Vnullcline1_2D(v_potential[j],p_2D)
		else
			v_null_2D_1[j] = Vnullcline1_2D_(v_potential[j],p_2D)
		end
	end
	return v_null_2D_1
end

# ╔═╡ 316a370c-e8bf-4da6-8782-381fba2e0d61
function give_me_v_null_2(v_potential,p4D,It,sample_index,i)
	v_null_2D_2 = zeros(size(v_potential))
	p_2D = change_I_(It[sample_index[i]],p_4D_to_2D(p4D))
	for j=1:length(v_potential)
		if It[sample_index[i]] >=0
			v_null_2D_2[j] = Vnullcline2_2D(v_potential[j],p_2D)
		else
			v_null_2D_2[j] = Vnullcline2_2D_(v_potential[j],p_2D)
		end
	end
	return v_null_2D_2
end

# ╔═╡ 4f6c0eb5-819d-4416-9449-e74f95c373f0
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

# ╔═╡ 2a6eb3be-ce0f-4e61-ae4f-d3e4abfa670e
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

# ╔═╡ 7f690094-96ac-439a-bf4a-8023cfdff5a0
function give_me_fp1_2D(p4D,It,sample_index,i)
	p_2D = change_I_(It[sample_index[i]],p_4D_to_2D(p4D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return fp1	
end

# ╔═╡ de46282f-4ce1-44f8-b645-c23ee550be02
function give_me_fp2_2D(p4D,It,sample_index,i)
	p_2D = change_I_(It[sample_index[i]],p_4D_to_2D(p4D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return fp2	
end

# ╔═╡ 20bfa55e-0ce0-46ae-8e7b-dc28004d5a38
function give_me_stab_2D(p4D,It,sample_index,i)
	p_2D = change_I_(It[sample_index[i]],p_4D_to_2D(p4D))
	fp1,fp2,stability = fixedpoints_2D(p_2D)
	return stability	
end

# ╔═╡ 0d082a4b-9fbf-4045-b710-1c2930c1c530
function give_me_phase_planes_l(v_potential,p4D,It,v,vs,sample_index)
	
	v_null_1_list = [give_me_v_null_1(v_potential,p4D,It,sample_index,i) for i=1:length(sample_index)]
	v_null_2_list = [give_me_v_null_2(v_potential,p4D,It,sample_index,i) for i=1:length(sample_index)]
	traj_v_list = [give_me_traj_v(v,sample_index,i) for i=1:length(sample_index)]
	traj_vs_list = [give_me_traj_vs(vs,sample_index,i) for i=1:length(sample_index)]
	fp1_list = [give_me_fp1_2D(p4D,It,sample_index,i) for i=1:length(sample_index)]
	fp2_list = [give_me_fp2_2D(p4D,It,sample_index,i) for i=1:length(sample_index)]
	stab_list = [give_me_stab_2D(p4D,It,sample_index,i) for i=1:length(sample_index)]
	
	return v_null_1_list,v_null_2_list,traj_v_list,traj_vs_list,fp1_list,fp2_list, stab_list
end

# ╔═╡ a2b418be-4fb6-4abf-95aa-f914d29ee325
function make_my_phase_plane(v_potential,v_null_1_list,v_null_2_list,vs_null_2D, traj_v_list,traj_vs_list,fp1_list,fp2_list,stab_list,recon_It,t,sample_index,p4D,i)
	
	p2D = p_4D_to_2D(p)
	
	plot_test = plot(v_potential,vs_null_2D,linewidth=2,linecolor=RGB(0,0.6,0.7),label="dVs/dt=0")
	
	if recon_It[sample_index[i]]>=0
		plot!(v_potential,v_null_1_list[i],linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
		plot!(v_potential,v_null_2_list[i],linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
	else
		plot!(v_null_1_list[i],v_potential,linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
		plot!(v_null_2_list[i],v_potential,linewidth=2,linecolor=RGB(1,0.4,0.1), label="dV/dt=0")
	end

	plot!(traj_v_list[i],traj_vs_list[i],linewidth=2,line_z=t,c=:cyclic_mygbm_30_95_c78_n256,label="Trajectory",linealpha=	0.9,legend=:topleft)
	
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
	
	
	xaxis!("V",(-80,0))
	yaxis!("Vs",(-80,0))
	#title!("t =$(Int((round(t[sample_index[i]])))),V0 = $(Int(p2D[2])), Vs0 = $(Int(p2D[3]))")
	title!("  A. ",titlelocation=:left,titlefontsize=18)
	return plot_test
end

# ╔═╡ 04ed6d51-354a-4432-904d-76a72d7fa322
function make_my_time_plot(v_potential,v_null_1_list,v_null_2_list,vs_null_2D, v,vs,vus,vuus,fp1_list,fp2_list,stab_list,recon_It,t,sample_index,p4D,i)
	
	p2D = p_4D_to_2D(p)
	
	if i<length(sample_index)
		plot_test = plot(t[1:sample_index[i+1]],v[1:sample_index[i+1]],linecolor="blue",label="V(t)")
	plot!(t[1:sample_index[i+1]],vs[1:sample_index[i+1]],linecolor="orange",label="Vs(t)")
	plot!(t[1:sample_index[i+1]],vus[1:sample_index[i+1]],linecolor="pink",label="Vus(t)")
	plot!(t[1:sample_index[i+1]],vuus[1:sample_index[i+1]],linecolor="black",label="Vuus(t)")
	else
		t_ = zeros(sample_index[i]+1)
		v_ = zeros(sample_index[i]+1)
		vs_ = zeros(sample_index[i]+1)
		vus_ = zeros(sample_index[i]+1)
		vuus_ = zeros(sample_index[i]+1)
		t_[1:end-1]= t[1:sample_index[i]]
		v_[1:end-1]= v[1:sample_index[i]]
		vs_[1:end-1]= vs[1:sample_index[i]]
		vus_[1:end-1]= vus[1:sample_index[i]]
		vuus_[1:end-1]= vuus[1:sample_index[i]]
		t_[end] = t[end]  
		v_[end] = v[end]  
		vs_[end] = vs[end]  
		vus_[end] = vus[end]  
		vuus_[end] = vuus[end]  
		plot_test = plot(t_,v_,linecolor="blue",label="V(t)")
		plot!(t_,vs_,linecolor="orange",label="Vs(t)")
		plot!(t_,vus_,linecolor="pink",label="Vus(t)")
		plot!(t_,vuus_,linecolor="black",label="Vuus(t)")
	end
	
	xaxis!("Time")
	yaxis!("Voltage",(-80,0))
	title!("  B. ",titlelocation=:left,titlefontsize=18)
	return plot_test
end

# ╔═╡ 57202474-8e51-43b3-8a71-17cc786aead6
md"#### Frenquency study" 

# ╔═╡ 61bff036-322e-4e53-9114-c99c6457e115
function frequency_from_period(T)
	F=zeros(length(T))
	for i=1:length(T)
		F[i]=1/T[i]
	end
	return F
end

# ╔═╡ 22a4c6d3-ab80-4f99-9b04-8f467edaa965
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

# ╔═╡ 7d63576b-39f7-41e0-a82b-ad570d31f8e5
function find_Vspikes_(solt,solv,tspan_max,mintspan)
	if mintspan <1
		i_tmin = findfirst(x -> x>=mintspan*tspan_max, solt)
	else
		i_tmin = findfirst(x -> x>=mintspan, solt)
	end
	solv_=solv[i_tmin:end]
	i_spike_ = findall(x -> x >= (Vmax-1), solv_)
	i_spike=[]
	for i=1:length(i_spike_)
		if length(solv_)>= i_spike_[i]+1
			if solv_[i_spike_[i]+1]<solv_[i_spike_[i]]
				append!(i_spike,i_spike_[i])
			end
		end
	end
	return i_spike.+ (i_tmin-1)
end

# ╔═╡ 6146153e-49d6-4b51-9b03-00c56c5c08c8
function find_Vuusspikes(solt,solvuus,tspan_max)
	i_tmin = findfirst(x -> x>=50*tspan_max/100, solt)
	solvuus_=solvuus[i_tmin:end]
	i_spike_ = findall(x -> x <= minimum(solvuus_)+1, solvuus_)
	i_spike=[]
	for i=1:length(i_spike_)
		if i_spike_[i]<length(solvuus_)
			if i_spike_[i]>1
				if solvuus_[i_spike_[i]+1]>solvuus_[i_spike_[i]] && solvuus_[i_spike_[i]-1]>solvuus_[i_spike_[i]]
					append!(i_spike,i_spike_[i])
				end
			end
		end
	end
	return i_spike.+ (i_tmin-1)
end

# ╔═╡ bcd46d0c-4738-4885-b6d4-ffe13143fce1
function burst_frequency_(cycle,i,sol,tspan)
	if cycle==true 
			ind_spikeV = find_Vspikes(sol.t,sol[1,:],maximum(tspan)) 
			ind_spikeVuus = find_Vuusspikes(sol.t,sol[4,:],maximum(tspan))
			
		delta_t_Vuus = sol.t[ind_spikeVuus[2:end]] - sol.t[ind_spikeVuus[1:end-1]]
		delta_i_Vuus = ind_spikeVuus[end]-ind_spikeVuus[end-1]
			t_period_inter= mean(delta_t_Vuus)
			
			
			delta_t_V = sol.t[ind_spikeV[2:end]] - sol.t[ind_spikeV[1:end-1]]
			i_p_n_1 =findfirst(x-> x>=ind_spikeVuus[end-1],ind_spikeV)
			i_p_n_2 =length(ind_spikeV) - (findfirst(x-> x<=ind_spikeVuus[end],reverse(ind_spikeV))-1)
			
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

# ╔═╡ 26dd3baf-c6ed-47c1-8337-69a0589a237a
function burst_frequency(cycle,i,sol,tspan)
	if cycle==true 
			ind_spikeV = find_Vspikes(sol.t,sol[1,:],maximum(tspan)) 
			ind_spikeVuus = find_Vuusspikes(sol.t,sol[4,:],maximum(tspan))
			
		delta_t_Vuus = sol.t[ind_spikeVuus[2:end]] - sol.t[ind_spikeVuus[1:end-1]]
		delta_i_Vuus = ind_spikeVuus[end]-ind_spikeVuus[end-1]
			t_period_inter= mean(delta_t_Vuus)
			
			
			delta_t_V = sol.t[ind_spikeV[2:end]] - sol.t[ind_spikeV[1:end-1]]
			i_p_n_1 =findfirst(x-> x>=ind_spikeVuus[end-1],ind_spikeV)
			i_p_n_2 =length(ind_spikeV) - (findfirst(x-> x<=ind_spikeVuus[end],reverse(ind_spikeV))-1)
			
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

# ╔═╡ b77a655e-4d54-4e3c-8bdb-a557141aa774
begin
	t_period_inter_ref,t_period_intra1_ref,t_period_intra2_ref,t_period_intra_av_ref, n_spike_ref,t_period_intra_full = burst_frequency_(true,NaN,sol,tspan)
	
	f_intra_full=frequency_from_period(t_period_intra_full[:]).*1000
	
	t_burst = ceil(sum(t_period_intra_full))
	t_period = collect(range(0.0,stop=t_burst,step=0.1))
	f_intra = zeros(size(t_period))
	ind_t_period_s = [1]
	for j=1:length(t_period_intra_full)
		ind_t_period_s_ = findfirst(x -> x>=sum(t_period_intra_full[1:j]),t_period)
		f_intra[ind_t_period_s[end]:ind_t_period_s_].= f_intra_full[j]
		append!(ind_t_period_s,ind_t_period_s_)
	end
	
	f_inter_ref = frequency_from_period(t_period_inter_ref).*1000
	f_intra_min_ref = frequency_from_period(t_period_intra1_ref).*1000
	f_intra_max_ref = frequency_from_period(t_period_intra2_ref).*1000
	f_intra_av_ref = frequency_from_period(t_period_intra_av_ref).*1000
end

# ╔═╡ f1789d85-45bc-4544-a175-197225e8ed6f
begin
	I_global_bif= collect(range(0,stop=100,step=0.5))
	V_cycle_global_bif = zeros(length(I_global_bif),2) #[max,min]
	t_period_inter = zeros(length(I_global_bif)) #[max,min]
	
	t_period_intra = zeros(length(I_global_bif),2) #[max,min]
	t_period_intra_av = zeros(length(I_global_bif))
	n_spike = zeros(length(I_global_bif))
	for i=1:length(I_global_bif)
		p_global_bif=change_I_(I_global_bif[i],p)
		tspan = (0.0,18000.0)
		u0=[-40.0,-40.0,-40.0,-40.0]
		prob = ODEProblem(MQIF_4D!,u0,tspan,p_global_bif,callback=cb)
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
	
	f_inter = frequency_from_period(t_period_inter).*1000
	f_intra_min = frequency_from_period(t_period_intra[:,1]).*1000
	f_intra_max = frequency_from_period(t_period_intra[:,2]).*1000
	f_intra_av = frequency_from_period(t_period_intra_av).*1000

	md"""Try to calculate max and min of the cycle in the 3D model after the SN bif but not smooth"""
end

# ╔═╡ 481602fa-f9e9-450f-a831-f0caf68eca54
begin
	prob_test = ODEProblem(MQIF_4D!,u0,tspan,change_I_(82.0,p),callback=cb)
	sol_test = solve(prob_test,dense=false)
	ind1_test,ind2_test = find_ind_period(sol_test[1,:],sol_test[4,:],sol_test.t)
	
	plot(sol_test.t[ind1_test:ind2_test],sol_test[1,ind1_test:ind2_test],linecolor="blue",label="V(t)")
	plot!(sol_test.t[ind1_test:ind2_test],sol_test[2,ind1_test:ind2_test],linecolor="orange",label="Vs(t)")
	plot!(sol_test.t[ind1_test:ind2_test],sol_test[3,ind1_test:ind2_test],linecolor="pink",label="Vus(t)")
	plot!(sol_test.t[ind1_test:ind2_test],sol_test[4,ind1_test:ind2_test],linecolor="black",label="Vuus(t)")
	plot!(sol_test.t[ind1_test:ind2_test],p[3]ones(size(sol_test.t[ind1_test:ind2_test])),linecolor="red",label="Vs0")
	plot!(sol_test.t[ind1_test:ind2_test],p[4]ones(size(sol_test.t[ind1_test:ind2_test])),linecolor="magenta",label="Vus0")
	plot!(sol_test.t[ind1_test:ind2_test],p[5]ones(size(sol_test.t[ind1_test:ind2_test])),linecolor="green",label="Vuus0",legend=:outertopright)
	md"""Test computation and plot"""
end

# ╔═╡ 4c197fd0-9ae6-447f-90d4-4e6f395613fd
function jacobian(v,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	J = zeros(4,4)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[1,3] = -2*gus*(v-vus0)
	J[1,4] = -2*guus*(v-vuus0)
	J[2,1] = 1
	J[2,2] = -1
	J[2,3] = 0
	J[2,4] = 0
	J[3,1] = 1
	J[3,2] = 0
	J[3,3] = -1
	J[3,4] = 0
	J[4,1] = 1
	J[4,2] = 0
	J[4,3] = 0
	J[4,4] = -1

	return J
end

# ╔═╡ 7307d524-5bd4-4a24-9fd0-99313af51db9
function fp_stability_4D(fp,p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	stability = ones(1,2)*66 #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable
	for i=1:length(fp)
		J = jacobian(fp[i],p)
		lambda1,lambda2,lambda3,lambda4 = eigvals(J)
		if real(lambda1)>0 && real(lambda2)>0 && real(lambda3)>0 && real(lambda4)>0
			stability[i]=1.0
		end
		if (real(lambda1)>0 && (real(lambda2)<0 ||( real(lambda3)<0 || real(lambda4)<0))) || (real(lambda1)<0 && (real(lambda2)>0 || (real(lambda3)>0 || real(lambda4)>0)))
			stability[i]=0.0
		end
		if real(lambda1)<0 && real(lambda2)<0 && real(lambda3)<0 && real(lambda4)<0
			stability[i]=-1.0
		end
	end
	return stability
end

# ╔═╡ 621533e4-1c67-4436-baf8-8d55598c1184
function global_I_bifurcation(p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p

	common_term = gf*v0 - gs*vs0 - gus*vus0 - guus*vuus0
	delta = gf*gs*(v0-vs0)^2 + gf*gus*(v0-vus0)^2 + gf*guus*(v0-vuus0)^2 - gs*gus*(vs0-vus0)^2 - gs*guus*(vs0-vuus0)^2 - gus*guus*(vus0-vuus0)^2 -(gf-gs-gus-guus)*I

	if delta >= 0
		numerator1 = common_term + sqrt(delta)
		numerator2 = common_term - sqrt(delta)
		denominator = gf-gs-gus-guus
		V_fp_1 = numerator1/denominator
		V_fp_2 = numerator2/denominator
		V_fp = [V_fp_1,V_fp_2]

		stability = fp_stability_4D(V_fp,p)
	else
		V_fp_1 = NaN
		V_fp_2 = NaN
		V_fp = [V_fp_1,V_fp_2]

		stability = [NaN,NaN]
	end
	return V_fp,stability
end

# ╔═╡ 90318732-189c-49e0-b24a-cc499b19f152
begin
	V_fp_global_bif = zeros(length(I_global_bif),2)
	stab_fp_global_bif = zeros(length(I_global_bif),2)
	
	V_fp_global_bif_r = zeros(length(I_global_bif),2)
	stab_fp_global_bif_r = zeros(length(I_global_bif),2)

	for i=1:length(I_global_bif)
		p_global_bif = change_I_(I_global_bif[i],p)
		V_fp_global_bif[i,:],stab_fp_global_bif[i,:]=global_I_bifurcation(p_global_bif)
		#V_fp_global_bif_r[i,:],stab_fp_global_bif_r[i,:]=global_I_bifurcation(p_global_bif_r)
	end

	ind_bif = findfirst(x-> x>Ibif,I_global_bif)-1
	Vbif = mean(V_fp_global_bif[ind_bif,:])
	
	#Ibif_r = find_global_bif(p_regen)
	#ind_bif_r = findfirst(x-> x>Ibif_r,I_global_bif)-1
	#Vbif_r = mean(V_fp_global_bif_r[ind_bif_r,:])
	md"""Fixed point computation for a range of current values"""
end

# ╔═╡ e02a4195-9b39-4388-84e6-fad12586ac5e
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

# ╔═╡ f2a9315d-b76d-4017-9bc2-5280f1dc8a97
md"#### Deep DHN behaviors simulation"

# ╔═╡ eae50c5b-81bd-41f4-995c-4df44d541c05
begin
	md"""code for levels"""
	#@bind gus_level html"<input type=range min=0 max=0.5 step=0.05>"
	#@bind guus_level html"<input type=range min=0 max=0.01 step=0.0002>"
	#gus_level,guus_level
end

# ╔═╡ 8dfb9cd4-d16b-49a7-81e1-5ba4085c614a
md"#### Dynamic input conductances"

# ╔═╡ 0d3b4e72-8498-4c99-b380-cc269dd74b1c
function dic(p,norm)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	

	v=collect(range(-70,Vmax,step=0.1))
	dic_gf = []
	dic_gs = []
	dic_gus = []
	dic_guus = []
	for i=1:length(v)
		J = jacobian(v[i],p)
		append!(dic_gf,J[1,1])
		append!(dic_gs,J[1,2])
		append!(dic_gus,J[1,3])
		append!(dic_guus,J[1,4])
	end
	
	dic_plot = plot(v,dic_gf,linecolor="blue",label="gf")
	plot!(v,dic_gs,linecolor="orange",label="gs")
	plot!(v,dic_gus,linecolor=RGB(1,0.5,0.8),label="gus")
	plot!(v,dic_guus,linecolor="black",label="guus")
	scatter!([v0],[jacobian(v0,p)[1,1]],markershape = :circle,markersize = 4,markeralpha = 0.6,markercolor = "blue",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="V0")
	scatter!([vs0],[jacobian(vs0,p)[1,2]],markershape = :circle,markersize = 4,markeralpha = 0.6,markercolor = "orange",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Vs0")
	scatter!([vus0],[jacobian(vus0,p)[1,3]],markershape = :circle,markersize = 4,markeralpha = 0.6,markercolor = RGB(1,0.5,0.8),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Vus0")
	scatter!([vuus0],[jacobian(vuus0,p)[1,4]],markershape = :circle,markersize = 4,markeralpha = 0.6,markercolor = "black",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Vuus0")
	xaxis!("V")
	yaxis!("DIC",(-50,50))
	
	fp_l,stab_l=global_I_bifurcation(change_I_(Ibif-30.,p))
	rest = fp_l[2]
	rest=-50
	ind_rest = findfirst(x -> x>=rest,v)
	dic_ = collect(range(-50,stop=50,step=5))
	if length(fp_l)>1
		plot!(dic_plot,ones(size(dic_)).*fp_l[2],dic_,linestyle=:dot,linecolor="green",label="Stable")
		plot!(dic_plot,ones(size(dic_)).*fp_l[1],dic_,linestyle=:dash,linecolor="purple",label="Saddle")
	end
	
	
	s = ["gf", "gs", "gus", "guus"]
	max_cond_plt = plot()
	bar!([dic_gf[end],0,0,0], xticks=(1:4, s),linealpha=0,c = [:blue])
	bar!([0,dic_gs[end],0,0], xticks=(1:4, s),linealpha=0,c = [:orange])
	bar!([0,0,dic_gus[end],0], xticks=(1:4, s),linealpha=0,c = RGB(1,0.5,0.8))
	bar!([0,0,0,dic_guus[end]], xticks=(1:4, s),linealpha=0,c = [:black])
	title!("Near max ($(round(Vmax)))")
	yaxis!((-150,150))
	
	rest_cond_plt = plot()
	bar!([dic_gf[ind_rest],0,0,0], xticks=(1:4, s),linealpha=0,c = [:blue])
	bar!([0,dic_gs[ind_rest],0,0], xticks=(1:4, s),linealpha=0,c = [:orange])
	bar!([0,0,dic_gus[ind_rest],0], xticks=(1:4, s),linealpha=0,c = RGB(1,0.5,0.8))
	bar!([0,0,0,dic_guus[ind_rest]], xticks=(1:4, s),linealpha=0,c = [:black])
	title!("Near rest ($(round(rest)))")
	
	cond_plt = plot(max_cond_plt,rest_cond_plt,layout=(1,2),legend=:none)
	
	
	return dic_plot,cond_plt
end

# ╔═╡ 33915b18-1f42-4b0e-8f15-2d35dea8451f
function tune_bif(p)
	vs0 = collect(range(-40,stop=-30,step=0.1))
	vus0 = collect(range(-40,stop=30,step=0.2))
	vuus0 = collect(range(-60,stop=-40,step=0.1))
	
	gs = collect(range(0,stop=2,step=0.01))
	gus = collect(range(0,stop=0.5,step=0.001))
	guus = collect(range(0,stop=0.5,step=0.0001))
	
	B_vs0 =[]
	B_vus0 =[]
	B_vuus0 =[]
	for i=1:length(vs0)
		append!(B_vs0,find_global_bif(change_vs0_(vs0[i],p)))
	end
	for i=1:length(vus0)
		append!(B_vus0,find_global_bif(change_vus0_(vus0[i],p)))
	end
	for i=1:length(vuus0)
		append!(B_vuus0,find_global_bif(change_vuus0_(vuus0[i],p)))
	end
	
	B_gs =[]
	B_gus =[]
	B_guus =[]
	for i=1:length(gs)
		append!(B_gs,find_global_bif(change_gs_(gs[i],p)))
	end
	for i=1:length(gus)
		append!(B_gus,find_global_bif(change_gus_(gus[i],p)))
	end
	for i=1:length(guus)
		append!(B_guus,find_global_bif(change_guus_(guus[i],p)))
	end
	
	plt_vs0 = plot(vs0,B_vs0)
	xaxis!("Vs0")
	plt_vus0 = plot(vus0,B_vus0)
	xaxis!("Vus0")
	plt_vuus0 = plot(vuus0,B_vuus0)
	xaxis!("Vuus0")
	
	plt_gs = plot(gs,B_gs)
	xaxis!("gs")
	yaxis!("I",(-500,500))
	plt_gus = plot(gus,B_gus)
	xaxis!("gss",(0.0,0.1))
	yaxis!("I",(-500,500))
	plt_guus = plot(guus,B_guus)
	xaxis!("gus",(0.0,0.1))
	yaxis!("I",(-500,500))
	
	plt_v_ = plot(plt_vs0,plt_vus0,plt_vuus0,layout=(1,3))
	plt_g_ = plot(plt_gs,plt_gus,plt_guus,layout=(1,3))
	#plt = plot(plt_v_,plt_g_,layout=(2,1))
	plt = plot(plt_g_,legend=:none)
	
	return plt,B_gs,gs
end

# ╔═╡ 1fa67863-e6c0-48a3-a937-a8c42cd49c29
plt_bif,b_gs,gs = tune_bif(change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p)))))))

# ╔═╡ 672f4668-856f-4da7-977c-a1f2f28bd68e
begin
	plot(plt_bif,size=(800,500))
	#savefig("4Dbifcurrent.pdf")
end

# ╔═╡ f75aee30-303c-48f4-a554-9d144ff3611d
begin
	p_plateau_pot_no_burst_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00001,change_gus_(0.002,p))))))
	#change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.000001,change_gus_(0.0015,p)))))) #with more spikes after pulse
	
	p_tonic_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.0004,change_gus_(0.0002,p)))))) #Validated 23/05,decreasing but delta around 3Hz, mean at 167Hz
	
	p_decel_reg_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.0009,change_gus_(0.000001,p))))))#change of 10Hz, considered as tonic 
	#change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.25,change_guus_(0.001,change_gus_(0.000001,p)))))) #change pf 30Hz (decreasing)
	
	p_decel_plateau_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.75,change_guus_(0.0001,change_gus_(0.000001,p))))))#change <1Hz during pulse, tonic during pulse 
	
	#p_stable_sh = change_I_(Ibif_l2+5.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.99,change_guus_(0.0005,change_gus_(0.0000001,p))))))
	
	p_unstable_full_sh = 0 #not found, low current may be too high 
	
	p_unstable_after_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.000002,change_gus_(0.00001,p)))))) #delta f during spike of 0.3Hz
	
	p_plateau_pot_nb_too_long_pulse_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.2,change_guus_(0.00005,change_gus_(0.002,p))))))#f increases and decreases during the pulse ; Final f at end of pulse higher that initial f start of pulse 
	
	p_accelerating_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00001,change_gus_(0.008,p)))))) #f increases, plateau pot early stopped
	
	p_burst_sh = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.001,change_gus_(0.008,p))))))
	
	dur_pulse_sh = 100.0
	Il_sh = 0.5
	Ih_sh = 20.0
	
	md"""Parameters used to simulate various responses to a pulse of 100ms"""
end

# ╔═╡ e7e0b571-c282-4db7-8a0c-842b3fe04f29
p_level4 = p_tonic_sh#change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.1,change_gus_(0.000001,p))))))

# ╔═╡ 8d5b0ab0-0825-446e-a6a3-132398e2ec04
begin
	plt_dic,plt_c_max = dic(p_level4,true)
	plt_c_max
end

# ╔═╡ 12043a36-7ae2-4f1d-a5a6-21827356335e
begin
	u0_l4=[-50,-50,-50,-50]
	tspan_l4 = (0.0,7500.0)
	prob_l4 = ODEProblem(MQIF_4D!,u0_l4,tspan_l4,p_level4,callback=cb)
	sol_l4 = solve(prob_l4,DP5(),reltol=1e-12,abstol=1e-12)
	
	plot(sol_l4.t,sol_l4[1,:],linecolor="blue",label="V(t)")
	plot!(sol_l4.t,sol_l4[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_l4.t,sol_l4[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_l4.t,sol_l4[4,:],linecolor="black",label="Vuus(t)",size=(300,200))
	
	#scatter!(sol.t[ind_I_0],v[ind_I_0])
	#scatter!([sol.t[ind1],sol.t[ind2]],[v[ind1],v[ind2]])
	yaxis!("Voltage")
	xaxis!("Time [ms]")
end

# ╔═╡ 62899eda-55e0-4ebf-bfdf-86835869e528
begin
	p_list_sh = [p_plateau_pot_no_burst_sh,p_tonic_sh,p_decel_reg_sh,p_unstable_after_sh,p_accelerating_sh,p_burst_sh]
	p_list = [p_plateau_pot_no_burst,p_tonic,p_decel_reg,p_unstable_after,p_accelerating]
	plt_p_ = plot()
	#plt_p = plot_cond_space_uus_us(p_list,plt_p_,0.0,10000.0,250.0,5.0,30.0)
	md"""Plot set of chosen parameters"""
end

# ╔═╡ 0f7db3d9-4037-4aea-9314-2d83c7e32f14
 #a,b = p_list_mesh_I(p_plateau_pot_no_burst_sh,true,false,true,false,true,0.5,NaN) #big simu

# ╔═╡ bdb2bf03-a44e-4c1a-b303-8533ef3907af
a[81][9]

# ╔═╡ fdc42466-d089-4f44-8fd1-0d71de8e9514
#plt_full_mesh = plot_cond_space_uus_us_(p_list_sh_mesh[1:end],plot(),0.0,10000.0,250.0,0.5,10.0,false,true,true)
#plt_full_mesh = plot_cond_space_uus_us_(p_list_sh_mesh[1:end],plot(),0.0,10000.0,100.0,0.5,20.0,true,false,true)

# ╔═╡ 6b0d2310-8252-4790-8fea-4a542ab0cf78
begin
	plot(plt_full_mesh,legend=:outertopright,yaxis=:none)
	#savefig("MQIF4D-condvar-result_gs_guus_sh.pdf")
end

# ╔═╡ 3f7ff72d-971a-492c-8d1c-80359c573c15
begin 
	gus = []
	gus_1_ = collect(10 .^(range(-6,stop=-2,length=20))) #gus
	gus_2_ = collect((range(10 .^(-2),stop=10 .^(-1),length=15))) #guus
	append!(gus,gus_1_)
	append!(gus,gus_2_)
end

# ╔═╡ 1a665a69-8122-4547-ba84-223bc4aba2ab
begin 
	guus = []
	guus_1_ = collect(10 .^(range(-6,stop=-2,length=20))) #guus
	guus_2_ = collect((range(10 .^(-2),stop=10 .^(-1),length=10))) #guus
	append!(guus,guus_1_)
	append!(guus,guus_2_)
end

# ╔═╡ ee8900b4-b0a1-4f2d-9f71-53e9ce6edf72
RGB(0.9,0.2,0.0), RGB(0.0,0.8,0.3),RGB(0.4,0.4,0.4),RGB(0.9,0.5,0.2) ,RGB(0.2,0.6,1), RGB(0.2,0.2,0.9),RGB(0.9,0.7,0.75),RGB(0.7,0.5,0.6),RGB(0.75,0.7,0.9), RGB(0.6,0.5,0.7) ,RGB(0.9,0.7,0.2),RGB(0.5,0.99,0.5),RGB(0.5,0.5,0.0) 

# ╔═╡ 932f0f93-1cb5-44ef-97be-cd36f352805d
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# ╔═╡ 047609ef-d362-431e-a16c-50570f2eaf92
function p_list_mesh(p,check_gs,check_gus,check_guus)
	extendedgus = false 
	extendedguus = false 
	
	if check_gus==true
		if extendedgus
			gus_range = collect(10 .^(range(-5,stop=-1,step=0.1)))
		else
			gus_range = []
			gus_1_ = collect(10 .^(range(-6,stop=-2,length=20))) #gus
			gus_2_ = collect((range(10 .^(-2),stop=10 .^(-1),length=15))) #guus
			append!(gus_range,gus_1_)
			append!(gus_range,gus_2_)
		end
	end
	
	if check_guus==true
		if extendedguus
			guus_range = collect(10 .^(range(-5,stop=-1,step=0.1)))
		else
			guus_range = []
			guus_1_ = collect(10 .^(range(-6,stop=-2,length=20))) #guus
			guus_2_ = collect((range(10 .^(-2),stop=10 .^(-1),length=10))) #guus
			append!(guus_range,guus_1_)
			append!(guus_range,guus_2_)
		end
	end
	
	if check_gs==true
		gs_range = collect(range(0.05,stop=0.95,step=0.05))
	end
	
	if check_guus==true && check_gus==true
		guus_grid,gus_grid = meshgrid(guus_range,gus_range)
		p_list = []
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i],change_gus_(gus_grid[i],p))
			append!(p_list,[p_grid])
		end
	end
	
	if check_guus==true && check_gs==true
		guus_grid,gs_grid = meshgrid(guus_range,gs_range)
		p_list = []
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i],change_gs_(gs_grid[i],p))
			append!(p_list,[p_grid])
		end
	end
	
	if check_gus==true && check_gs==true
		gus_grid,gs_grid = meshgrid(gus_range,gs_range)
		p_list = []
		for i=1:length(gus_grid)
			p_grid = change_gus_(gus_grid[i],change_gs_(gs_grid[i],p))
			append!(p_list,[p_grid])
		end
	end
	
	
	
	return p_list 
end

# ╔═╡ 78ded07f-ceac-4ac4-82bd-991da87d8d78
begin
	#p_list_sh_mesh = p_list_mesh(p_plateau_pot_no_burst_sh,false,true,true)
	p_list_sh_mesh = p_list_mesh(p_plateau_pot_no_burst_sh,true,false,true)
end

# ╔═╡ fc7cb2be-9810-4dbc-bce3-52f1aa45d4ce
xx,yy = meshgrid(guus,collect(range(0.05,stop=0.95,step=0.05)))

# ╔═╡ d5f0e41d-9b82-4c5a-af89-3560464895d0
xx2,yy2 = meshgrid(guus,gus)

# ╔═╡ 46a50de2-2aa3-462e-9fc4-9dd96b3b52bc
scatter(xx2,yy2)

# ╔═╡ a6632203-f59c-4793-ad46-4b508214d08d
md"###### Effect of applied current"

# ╔═╡ d99f7e3c-3de6-4a6d-86d6-2d63755e8e28
length(collect(range(0.0,stop=20.0,step=1))),length(collect(range(5.0,stop=25.0,step=1)))

# ╔═╡ 34ed2f38-0376-454f-877f-2e6ac7fb34be
begin
	u0_t=[-70,-70,-70,-70]
	tspan_t = (0.0,7500.0)
	prob_t = ODEProblem(MQIF_4D!,u0_t,tspan_t,change_I_(5.0,change_guus_(0.00005,p_tonic_sh)),callback=cb)
	sol_t = solve(prob_t,DP5(),reltol=1e-12,abstol=1e-12)
	
	plot(sol_t.t,sol_t[1,:],linecolor="blue",label="V(t)")
	plot!(sol_t.t,sol_t[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol_t.t,sol_t[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol_t.t,sol_t[4,:],linecolor="black",label="Vuus(t)",size=(300,200))
	
	#scatter!(sol.t[ind_I_0],v[ind_I_0])
	#scatter!([sol.t[ind1],sol.t[ind2]],[v[ind1],v[ind2]])
	yaxis!("Voltage")
	xaxis!("Time [ms]")
end

# ╔═╡ 7dd65341-de1f-49a0-bd39-d374e287f23a
gr()

# ╔═╡ 6e99ef78-e4cf-4c91-9e26-09a6cf31b165
#plt_I,class_I = Il_Ih_effect(change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p))))))) 

# ╔═╡ 915b6ea6-8ec6-4a34-93d0-5c8d2596980e
plt_I

# ╔═╡ 4bf33fbd-e7db-4dd1-a516-a35657f4055a
begin
	plot(plt_I,legend=:outertopright)
	#savefig("choiceofIlandIh.pdf")
end

# ╔═╡ 833e6454-c96c-4d40-838a-8b631f7d171b
begin
	#==
	title_test = plot(title = translate_classification(conv_I), grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=20,titlefontcolor=color_classification(conv_I))
	full_test = subplot_simu_step_full(0.0,10000.0,pulse_dur,p_plateau_pot_no_burst_sh,0.5,8.0,false,false,true)
	plot(title_test,full_test,
	layout = @layout([A{0.03h} ;B ]),size=(800,800))
	xaxis!((5050,5150))
	==#
end

# ╔═╡ dce92cf6-b5f0-40e7-b0b9-79eed0a74ffe
begin
	#==
	I_l_list_ = collect(range(0.0,stop=5.0,step=0.25))
	p_list_mesh(p,check_gs,check_gus,check_guus,check_Il,check_Ih,Il,Ih)
	==#
end

# ╔═╡ 2cedbac4-9944-4299-a6ec-331f7cbf1e2c
b[81],a[81]

# ╔═╡ bf8883c8-1387-4b4c-aa0a-dfe184d48816
limits_I_

# ╔═╡ 3ee5ac17-abd8-4dd7-8ec5-fa94dd93f013


# ╔═╡ 89328195-8ed9-47ab-912c-16fbb12bf5eb
function translate_classification(conv)
	if conv == -1
		return "Unstable"
	end
	if conv == 0
		return "Stable"
	end
	if conv == 1
		return "Burst"
	end
	if conv == 2
		return "Tonic"
	end
	if conv == 3
		return "Accelerating"
	end
	if conv == 4
		return "Plateau potentials"
	end
	if conv == 5
		return "Decelerating"
	end
	if conv == 6
		return "Decelerating with plateau"
	end
	if conv == 7
		return "Curved"
	end
	if conv == 8
		return "Curved with plateau"
	end
	if conv == 9
		return "Tonic with plateau"
	end
	if conv == 10
		return "Latent"
	end
	if conv == NaN
		return "Undefined"
	end
end

# ╔═╡ b459c6fe-696a-4af2-80b4-a7c81f7136ea
function color_classification(conv)
	if conv == -1
		#unstable
		return RGB(0.9,0.2,0.0)
	end
	if conv == 0
		#stable
		return RGB(0.0,0.8,0.3)
	end
	if conv == 1
		#burst
		return RGB(0.4,0.4,0.4)
	end
	if conv == 2
		#tonic
		return RGB(0.9,0.5,0.2) 
	end
	if conv == 3
		#accelerating 
		return RGB(0.2,0.6,1)
	end
	if conv == 4
		#plateau potentials 
		return RGB(0.2,0.2,0.9)
	end
	if conv == 5
		#deceleraing 
		return RGB(0.9,0.7,0.75)
	end
	if conv == 6
		#deceleraing with plateau 
		return RGB(0.7,0.5,0.6)
	end
	if conv == 7
		#curved
		return RGB(0.75,0.7,0.9)
	end
	if conv == 8
		#curved with plateau 
		return RGB(0.6,0.5,0.7)
	end
	if conv == 9
		#tonic with plateau 
		return RGB(0.9,0.7,0.2)
	end
	if conv == 10
		#stable due to spike latency
		return RGB(0.5,0.99,0.5)
	end
	if conv == NaN
		#undefined
		return RGB(0.5,0.5,0.0)
	end
end

# ╔═╡ 4d7c6062-b61e-44bc-a90d-7fd8f62b178e
function freq_during_pulse(t,tspan_pulse,f_ins,t_ins)
	ind_tmin = findfirst(x -> x >=minimum(tspan_pulse),t)
	ind_tmax = findfirst(x -> x >=maximum(tspan_pulse),t)
	
	tmin = t[ind_tmin]
	tmax = t[ind_tmax]
	
	ind_spike_pulse_0 = findfirst(x -> x >= tmin, t_ins)
	if length(t_ins)>1
		if t_ins[end]>tmax
			ind_spike_pulse_f = findfirst(x -> x >= tmax, t_ins)-1
		else
			ind_spike_pulse_f = length(t_ins)
		end
			ind_f_max_pulse = findfirst(x -> x == maximum(f_ins[ind_spike_pulse_0:ind_spike_pulse_f]),f_ins[ind_spike_pulse_0:ind_spike_pulse_f])+ind_spike_pulse_0-1

		if abs(maximum(f_ins)-minimum([f_ins[ind_spike_pulse_0],f_ins[ind_spike_pulse_f]])) <=(abs(maximum(f_ins))*0.02)
			#constant 
			f_evol = 0
		else
			if ind_f_max_pulse <= ind_spike_pulse_0 + (ind_spike_pulse_f-ind_spike_pulse_0)/3
				#decelerating 
				f_evol = -1
			else
				if ind_f_max_pulse >= ind_spike_pulse_0 + 2*(ind_spike_pulse_f-ind_spike_pulse_0)/3
					#accelerating 
					f_evol =1
				else
					#curved 
					f_evol = 2
				end
			end
		end
	else
		f_evol=NaN
		return f_evol 
	end
	
	return f_evol 
end

# ╔═╡ a35ef1ed-b9b5-4702-aa91-555a28911578
function check_if_tonic(sol,tspan_pulse,f_ins,t_ins)
	ind_tmin = findfirst(x -> x >=minimum(tspan_pulse)/2,sol.t)
	ind_tmax = findfirst(x -> x >=maximum(tspan_pulse)/2,sol.t)
	
	tmin = sol.t[ind_tmin]
	tmax = sol.t[ind_tmax]
	
	ind_spike_pulse_0 = findfirst(x -> x >= tmin, t_ins)
	ind_spike_pulse_f = findfirst(x -> x >= tmax, t_ins)-1
	
	if abs(f_ins[ind_spike_pulse_0]-f_ins[ind_spike_pulse_f])<=5 && (abs(maximum(f_ins)-minimum(f_ins[ind_spike_pulse_0],f_ins[ind_spike_pulse_f])) <=5 )
		return true
	else
		return false
	end
end

# ╔═╡ 842a06d5-72ee-44c6-b38c-f549e439b95c
function check_if_plateau(t,sim,tspan,tspan_pulse)
	pulse_dur = maximum(tspan_pulse)-minimum(tspan_pulse)
	ind_tmin = findfirst(x -> x >maximum(tspan_pulse)+(0.01*pulse_dur),t)
	ind_spike = findall(x -> x>=Vmax-1,sim[1,ind_tmin:end]).+(ind_tmin-1)
	
	if length(ind_spike)>2
		return true
	else
		return false
	end
end

# ╔═╡ 18f4d337-2549-4dc9-97f0-abb3130fef8f
function check_if_unstable(t,sim,tspan)
	ind_tmin = findfirst(x -> x >=t[end]-1000,t)
	ind_spike = findall(x -> x>=Vmax-1,sim[1,ind_tmin:end]).+(ind_tmin-1)
	
	if length(ind_spike)>0
		return t[ind_spike[end]]>=t[end]-100
	else
		return false
	end
end

# ╔═╡ 9688f176-b51e-4dcd-bf21-fbe5af0c366c
function check_if_stable(sol,tspan)
	ind_tmin = findfirst(x -> x >=sol.t[end]/2,sol.t)
	ind_spike = findall(x -> x>=Vmax-1,sol[1,ind_tmin:end])
	
	if length(ind_spike)<=1
		return true
	else
		return false
	end
end

# ╔═╡ a7e59710-2ee0-47fa-a7ba-63054e8cc708
function full_frequency(t,sim,tspan,tspan_pulse)
	#to use if we are sure we do not have burst 
	
	ind_spikeV = find_Vspikes_(t,sim[1,:],maximum(tspan),0) #on the full period
	delta_t_V = t[ind_spikeV[2:end]] - t[ind_spikeV[1:end-1]]
	
	t_period_full =[]
	t_sampled_f_ins = []

	for j=1:length(delta_t_V)
		append!(t_period_full,delta_t_V[j])
		append!(t_sampled_f_ins,sum(t_period_full)+t[ind_spikeV[1]])
	end
		
	f_ins_full=frequency_from_period(t_period_full[:]).*1000
	
	return f_ins_full,t_sampled_f_ins
end

# ╔═╡ adb918d6-6614-4749-8f50-222138ac2518
begin
	sim=zeros(4,length(sol_t.t))
	sim[1,:]=sol_t[1,:]
	sim[2,:]=sol_t[2,:]
	sim[3,:]=sol_t[3,:]
	sim[4,:]=sol_t[4,:]
	full_frequency(sol_t.t,sim,tspan,(NaN,NaN))
end

# ╔═╡ 953b5d8c-5fbd-4ed0-8b46-21e4bee0e565
function check_if_cycle(sol,tspan)
	sim = zeros(4,length(sol.t))
	sim[1,:].= sol[1,:]
	sim[2,:].= sol[2,:]
	sim[3,:].= sol[3,:]
	sim[4,:].= sol[4,:]
	f,t = full_frequency(sol.t,sim,tspan,(0.0,0.0))
	
	i_t_min = findall(x -> x >=maximum(tspan)/2,t)
	if length(i_t_min)>=1
		ind_min_f = findall(x -> x <= minimum(f[i_t_min[1]:end])+minimum(f[i_t_min[1]:end])*0.1,f[i_t_min[1]:end])
	else
		ind_min_f = findall(x -> x <= minimum(f)+minimum(f)*0.1,f)
	end
	
	if length(ind_min_f)>1
		d_ind_min_f = ind_min_f[2:end]-ind_min_f[1:end-1]
		check_cycle =  minimum(d_ind_min_f) ==1
	else
		check_cycle = false
	end
	return check_cycle
end

# ╔═╡ 3d9b0adc-389b-4a58-869c-b1a5c1af22a6
function check_if_burst(sol,tspan)
	sim = zeros(4,length(sol.t))
	sim[1,:].= sol[1,:]
	sim[2,:].= sol[2,:]
	sim[3,:].= sol[3,:]
	sim[4,:].= sol[4,:]
	f,t = full_frequency(sol.t,sim,tspan,(0.0,0.0))
	
	i_t_min = findall(x -> x >=maximum(tspan)/2,t)
	if length(i_t_min)>=1
		ind_min_f = findall(x -> x <= minimum(f[i_t_min[1]:end])+minimum(f[i_t_min[1]:end])*0.1,f[i_t_min[1]:end])
	else
		ind_min_f = findall(x -> x <= minimum(f)+minimum(f)*0.1,f)
	end
	
	if length(ind_min_f)>1
		d_ind_min_f = ind_min_f[2:end]-ind_min_f[1:end-1]
		check_burst = maximum(d_ind_min_f)-minimum(d_ind_min_f)<5 && minimum(d_ind_min_f) >1
	else
		check_burst = false
	end
	return check_burst 
end

# ╔═╡ d2e58edd-5759-42f3-920c-0ab640c063fd
function plot_I_effect(p_root,Iapp)
	Ibif_2D_model_to_cycle = find_local_bif(p_root)
	Ibif_2D_model_to_stable = 0 
	conv_4D=[]
	It = zeros(2,length(Iapp))
	f_i = zeros(2,length(Iapp))
	
	for i=1:length(Iapp)
		p = change_I_(Iapp[i], p_root)
		
		##4D model simulation 
		u0=[-60,-60,-60,-60]
		tspan = (0.0,7500.0)
		prob = ODEProblem(MQIF_4D!,u0,tspan,p,callback=cb)
		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		
		##minimum time ensuring convergence 
		t_min = findfirst(x -> x >=maximum(tspan)*0.5,sol.t)
		
		##compute It in 2D model 
		Iuus = compute_Iuus(sol[4,:],p)
		Ius = compute_Ius(sol[3,:],p)
		It_ = p[1].+Iuus.+Ius
		It[1,i] = maximum(It_[t_min:end])
		It[2,i] = minimum(It_[t_min:end])
		
		##convergence of the 4D model in response to a constant current
		conv_4D_ = [NaN,NaN,NaN]
		conv_4D_[1]=check_if_stable(sol,tspan)
		conv_4D_[2]=check_if_burst(sol,tspan)
		conv_4D_[3]=check_if_cycle(sol,tspan)
		
		ind_conv_4D_ = findall(x -> x == true, conv_4D_)
		if length(ind_conv_4D_)==1
			if ind_conv_4D_[1]==1
				push!(conv_4D,["stable"])
			else
				if ind_conv_4D_[1]==2
					push!(conv_4D,["burst"])
				else
					if ind_conv_4D_[1]==3
						push!(conv_4D,["cycle"])
					end
				end
			end
		else
			push!(conv_4D,["undefined"])
		end
		
		##frequency compuattion 
		sim=zeros(4,length(sol.t))
		sim[1,:]=sol[1,:]
		sim[2,:]=sol[2,:]
		sim[3,:]=sol[3,:]
		sim[4,:]=sol[4,:]
		f_ins,t_sampled_f = full_frequency(sol.t,sim,tspan,(NaN,NaN))
		if conv_4D_[1]==false
			t_sampled_min = findfirst(x -> x>= sol.t[t_min],t_sampled_f)
			f_i[1,i] = maximum(f_ins[t_sampled_min:end])
			f_i[2,i] = minimum(f_ins[t_sampled_min:end])
		else
			f_i[1,i] = NaN
			f_i[2,i] = NaN
		end
		
		
		#plot I,It based on culculation of It and convergence in 2D model 
		#plot I,f  and add colors to remind behavior of 4D model (burst to tonic?)
	end
	
	if Iapp[1]>0.0
		step_I = Iapp[1]
	else
		step_I = Iapp[2]
	end
	
	ind_l = findall(x -> abs(x-Ibif_2D_model_to_stable)<=step_I/2,It[1,:])
	target_diff = abs.(It[1,:].-Ibif_2D_model_to_cycle)
	ind_h = findall(x -> x <=step_I/2 ,target_diff)
	
	if length(ind_h)>0
		Ih=Iapp[ind_h[1]]
		Ith = It[1,ind_h[1]]
	else
		ind_h_ = findall(x -> x >=Ibif_2D_model_to_cycle ,It[1,:])
		if length(ind_h_)>0
			Ih=Iapp[ind_h_[1]]
			Ith = It[1,ind_h_[1]]
		else
			Ih=NaN
			Ith=NaN
		end
	end
	
	if length(ind_l)>0
		Il=Iapp[ind_l[1]]
		Itl = It[1,ind_l[1]]
	else
		Il=NaN
		Itl=NaN
	end
	
	##Plot I;It
	plot_I_It = plot(Iapp,It[1,:],linecolor=RGB(0.9,0.2,0.4),linewidth=2,label="Maximum")
				plot!(Iapp,It[2,:],linecolor=RGB(1,0.4,0.4),linewidth=2,linestyle=:dot,label="Minimum")
				#plot!(Iapp,Ibif_2D_model_to_cycle.*ones(size(Iapp)),label="SN bifurcation",linewidth=1.5,linecolor=RGB(0.9,0.9,0),legend=:topleft)
				#plot!(Iapp,zeros(size(Iapp)),label="Bifurcation to rest",linewidth=1.5,linecolor=RGB(0,0.5,0.5))
	
				plot!(Iapp,maximum([maximum(It[1,:]),Ibif_2D_model_to_cycle]).*ones(size(Iapp)),fill = ([Ibif_2D_model_to_cycle.*ones(size(Iapp))], 0.1, RGB(0,0.5,1)),linecolor=:purple,linealpha=0.1,label="Limit cycle")
				plot!(Iapp,Ibif_2D_model_to_cycle.*ones(size(Iapp)),fill = ([Ibif_2D_model_to_stable.*ones(size(Iapp))], 0.1, RGB(1,1,0)),linecolor=:purple,linealpha=0.1,label="Bistable")
				plot!(Iapp,Ibif_2D_model_to_stable.*ones(size(Iapp)),fill = ([minimum([minimum(It[2,:]),Ibif_2D_model_to_stable]).*ones(size(Iapp))], 0.1, RGB(0,1,0)),linecolor=:purple,linealpha=0.1,label="Stable")	
				scatter!([Ih],[Ith],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.0,0.0,1),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Ih = $(Ih)")
				scatter!([Il],[Itl],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,0.5,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Il = $(Il)",legend=:outertopright)
				xaxis!("I")
				yaxis!("It")
	
	##Find points with the same convergence
	i_stable = findall(x-> x == ["stable"],conv_4D)
	i_burst = findall(x-> x == ["burst"],conv_4D)
	i_cycle = findall(x-> x == ["cycle"],conv_4D)
	
	f_i_exist = findall(x -> x>=0,f_i[1,:])
	
	plot_I_fi = plot()
				plot!(Iapp,f_i[1,:],linecolor=RGB(0.9,0.2,0.4),linewidth=2)
				plot!(Iapp,f_i[2,:],linecolor=RGB(1,0.4,0.4),linewidth=2,linestyle=:dot)
				plot!(Iapp[i_stable],maximum(f_i[1,f_i_exist]).*ones(size(Iapp[i_stable])),fill = ([minimum(f_i[2,f_i_exist]).*ones(size(Iapp[i_stable]))], 0.2, RGB(0.5,1,0.5)),linecolor=:purple,linealpha=0.1,label="Stable")
		if length(i_cycle)>0
				plot!(Iapp[i_cycle],maximum(f_i[1,f_i_exist]).*ones(size(Iapp[i_cycle])),fill = ([minimum(f_i[2,f_i_exist]).*ones(size(Iapp[i_cycle]))], 0.2, RGB(0.65,0.3,0.65)),linecolor=:purple,linealpha=0.1,label="Cycle")
		end
		if length(i_burst)>0
				plot!(Iapp[i_burst],maximum(f_i[1,f_i_exist]).*ones(size(Iapp[i_burst])),fill = ([minimum(f_i[2,f_i_exist]).*ones(size(Iapp[i_burst]))], 0.2, RGB(1,0.5,0)),linecolor=:purple,linealpha=0.1,label="Burst")
		end
				xaxis!("I")
				yaxis!("f")
				
	
	#plot_I = plot(plot_I_It,plot_I_fi,layout=(2,1),size=(500,800))
	plot_I=plot_I_It
	
	#return plot_I,[Iapp[ind_h],Iapp[ind_l]]
	return plot_I,[Iapp,It[1,:]]
end

# ╔═╡ cbe7686f-4d76-4bd0-8983-c0a5633c22d6
begin
	Iapp=collect(range(0.0,stop=30,step=0.25))
	plot_I,I_limits=plot_I_effect(change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p)))))),Iapp)
end

# ╔═╡ 497f8d3c-afcb-42ce-a50f-00b5bc6ecbb1
begin
	Il = floor(I_limits[2])
	Ih = ceil(I_limits[1])+2.5
	pulse_dur = 100.0
end

# ╔═╡ 5c8ec9ad-23c0-4cb7-a5cc-fa95d4a1368d
Il,Ih

# ╔═╡ cee6d373-c5f9-441f-ba98-b7c1cd2ff394
begin
	plot(plot_I,legend=:outertopright)
	#savefig("4DIl_It_def.pdf")
end

# ╔═╡ 990efd90-eb3c-4f1d-8f2b-7de99916cbf4
function limits_I_effect(p_root,Iapp)
	Ibif_2D_model_to_cycle = find_local_bif(p_root)
	Ibif_2D_model_to_stable = 0 
	conv_4D=[]
	It = zeros(2,length(Iapp))
	f_i = zeros(2,length(Iapp))
	
	for i=1:length(Iapp)
		p = change_I_(Iapp[i], p_root)
		
		##4D model simulation 
		u0=[-60,-60,-60,-60]
		tspan = (0.0,7500.0)
		prob = ODEProblem(MQIF_4D!,u0,tspan,p,callback=cb)
		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		
		##minimum time ensuring convergence 
		t_min = findfirst(x -> x >=maximum(tspan)*0.5,sol.t)
		
		##compute It in 2D model 
		Iuus = compute_Iuus(sol[4,:],p)
		Ius = compute_Ius(sol[3,:],p)
		It_ = p[1].+Iuus.+Ius
		It[1,i] = maximum(It_[t_min:end])
		It[2,i] = minimum(It_[t_min:end])
		
		##convergence of the 4D model in response to a constant current
		conv_4D_ = [NaN,NaN,NaN]
		conv_4D_[1]=check_if_stable(sol,tspan)
		conv_4D_[2]=check_if_burst(sol,tspan)
		conv_4D_[3]=check_if_cycle(sol,tspan)
		
		ind_conv_4D_ = findall(x -> x == true, conv_4D_)
		if length(ind_conv_4D_)==1
			if ind_conv_4D_[1]==1
				push!(conv_4D,["stable"])
			else
				if ind_conv_4D_[1]==2
					push!(conv_4D,["burst"])
				else
					if ind_conv_4D_[1]==3
						push!(conv_4D,["cycle"])
					end
				end
			end
		else
			push!(conv_4D,["undefined"])
		end
		
		##frequency compuattion 
		sim=zeros(4,length(sol.t))
		sim[1,:]=sol[1,:]
		sim[2,:]=sol[2,:]
		sim[3,:]=sol[3,:]
		sim[4,:]=sol[4,:]
		f_ins,t_sampled_f = full_frequency(sol.t,sim,tspan,(NaN,NaN))
		if conv_4D_[1]==false
			t_sampled_min = findfirst(x -> x>= sol.t[t_min],t_sampled_f)
			f_i[1,i] = maximum(f_ins[t_sampled_min:end])
			f_i[2,i] = minimum(f_ins[t_sampled_min:end])
		else
			f_i[1,i] = NaN
			f_i[2,i] = NaN
		end
		
		
		#plot I,It based on culculation of It and convergence in 2D model 
		#plot I,f  and add colors to remind behavior of 4D model (burst to tonic?)
	end
	
	if Iapp[1]>0.0
		step_I = Iapp[1]
	else
		step_I = Iapp[2]
	end
	
	ind_l = findall(x -> abs(x-Ibif_2D_model_to_stable)<=step_I/2,It[1,:])
	target_diff = abs.(It[1,:].-Ibif_2D_model_to_cycle)
	ind_h = findall(x -> x <=step_I/2 ,target_diff)
	
	if length(ind_h)>0
		Ih=Iapp[ind_h[1]]
	else
		ind_h_ = findall(x -> x >=Ibif_2D_model_to_cycle ,It[1,:])
		if length(ind_h_)>0
			Ih=Iapp[ind_h_[1]]
		else
			Ih=NaN
		end
	end
	
	if length(ind_l)>0
		Il=Iapp[ind_l[1]]
	else
		Il=NaN
	end
	
	return [Ih,Il]
end

# ╔═╡ 9d4c91d7-2e51-4cb2-9ed6-0f1f67c694c7
function p_list_mesh_I(p,check_gs,check_gus,check_guus,check_Il,check_Ih,Il,Ih)
	extendedgus = true 
	extendedguus = true 
	
	if check_gus==true
		if extendedgus
			#gus_range = collect(10 .^(range(-5,stop=-1.5,step=0.1)))
			gus_range = collect(10 .^(range(-6,stop=-2,length=13)))
		else
			gus_range = []
			gus_1_ = collect(10 .^(range(-6,stop=-2,length=10))) #gus
			gus_2_ = collect((range(10 .^(-2),stop=10 .^(-1.5),length=3))) #guus
			append!(gus_range,gus_1_)
			append!(gus_range,gus_2_)
		end
	end
	
	if check_guus==true
		if extendedguus
			#guus_range = collect(10 .^(range(-5,stop=-2,step=0.1)))
			guus_range = collect(10 .^(range(-5,stop=-1,length=13)))
		else
			guus_range = []
			guus_1_ = collect(10 .^(range(-6,stop=-2,length=10))) #guus
			guus_2_ = collect((range(10 .^(-2),stop=10 .^(-1),length=3))) #guus
			append!(guus_range,guus_1_)
			append!(guus_range,guus_2_)
		end
	end
	
	if check_gs==true
		gs_range = collect(range(0.05,stop=0.95,length=10))
	end
	
	if check_Il==true
		Il_range = collect(range(0.0,stop=Ih-1,length=21))
	end
	
	
	I_list = []
	
	
	if check_guus==true && check_gus==true
		guus_grid,gus_grid = meshgrid(guus_range,gus_range)
		p_list = []
		if check_Ih==true
			Iapp=collect(range(0.0,stop=10,step=0.05))
			Iapp_max = maximum(Iapp)
			Iapp_min = minimum(Iapp)
			Ih_range = []
			Il_mem = []
		end
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i],change_gus_(gus_grid[i],p))
			append!(p_list,[p_grid])
			if check_Ih==true
				if i>1
					if Ih_range[i-1]>Iapp_max*2/3
						if Il_mem[i-1]>Iapp_max*1/2
							Iapp_=collect(range(Iapp_max*1/2,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						else
							Iapp_=collect(range(Iapp_min,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						end
					else
						Ih_,Il_=limits_I_effect(p_grid,Iapp)
					end
				else
					Ih_,Il_=limits_I_effect(p_grid,Iapp)
				end
				append!(Ih_range,Ih_)
				append!(Il_mem,Il_)
			end
		end
		if check_Ih==true
			I_list = Ih_range
		end
	end
	
	if check_guus==true && check_gs==true
		guus_grid,gs_grid = meshgrid(guus_range,gs_range)
		p_list = []
		if check_Ih==true
			Iapp=collect(range(0.0,stop=10,step=0.05))
			Iapp_max = maximum(Iapp)
			Iapp_min = minimum(Iapp)
			Ih_range = []
			Il_mem = []
		end
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i],change_gs_(gs_grid[i],p))
			append!(p_list,[p_grid])
			if check_Ih==true
				if i>1
					if Ih_range[i-1]>Iapp_max*2/3
						if Il_mem[i-1]>Iapp_max*1/2
							Iapp_=collect(range(Iapp_max*1/2,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						else
							Iapp_=collect(range(Iapp_min,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						end
					else
						Ih_,Il_=limits_I_effect(p_grid,Iapp)
					end
				else
					Ih_,Il_=limits_I_effect(p_grid,Iapp)
				end
				append!(Ih_range,Ih_)
				append!(Il_mem,Il_)
			end
		end
		if check_Ih==true
			I_list = Ih_range
		end
	end
	
	if check_gus==true && check_gs==true
		gus_grid,gs_grid = meshgrid(gus_range,gs_range)
		p_list = []
		if check_Ih==true
			Iapp=collect(range(0.0,stop=10,step=0.05))
			Iapp_max = maximum(Iapp)
			Iapp_min = minimum(Iapp)
			Ih_range = []
			Il_mem = []
		end
		for i=1:length(gus_grid)
			p_grid = change_gus_(gus_grid[i],change_gs_(gs_grid[i],p))
			append!(p_list,[p_grid])
			if check_Ih==true
				if i>1
					if Ih_range[i-1]>Iapp_max*2/3
						if Il_mem[i-1]>Iapp_max*1/2
							Iapp_=collect(range(Iapp_max*1/2,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						else
							Iapp_=collect(range(Iapp_min,stop=Iapp_max*3, step=0.125))
							Ih_,Il_=limits_I_effect(p_grid,Iapp_)
						end
					else
						Ih_,Il_=limits_I_effect(p_grid,Iapp)
					end
				else
					Ih_,Il_=limits_I_effect(p_grid,Iapp)
				end
				append!(Ih_range,Ih_)
				append!(Il_mem,Il_)
			end
		end
	end
	
	if check_guus==true && check_Il==true
		guus_grid,Il_grid = meshgrid(guus_range,Il_range)
		p_list = []
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i])
			append!(p_list,[p_grid])
			append!(I_list,[Il_grid])
		end
	end
	
	#==
	if check_guus==true && check_Ih==true
		guus_grid,Ih_grid = meshgrid(guus_range,Ih_range)
		p_list = []
		for i=1:length(guus_grid)
			p_grid = change_guus_(guus_grid[i])
			append!(p_list,[p_grid])
			append!(I_list,[Ih_grid])
		end
	end
	==#
	
	if check_gus==true && check_Il==true
		gus_grid,Il_grid = meshgrid(gus_range,Il_range)
		p_list = []
		for i=1:length(gus_grid)
			p_grid = change_gus_(gus_grid[i])
			append!(p_list,[p_grid])
			append!(I_list,[Il_grid])
		end
	end
	
	#==
	if check_gus==true && check_Ih==true
		gus_grid,Ih_grid = meshgrid(gus_range,Ih_range)
		p_list = []
		for i=1:length(gus_grid)
			p_grid = change_gus_(gus_grid[i])
			append!(p_list,[p_grid])
			append!(I_list,[Ih_grid])
		end
	end
	==#
	
	
	return p_list,I_list
end

# ╔═╡ 2918e0cb-8d0a-403d-a04d-679df9660b5c
function Il_Ih_limits_variation(guus,gus,p_root) 
	Iapp=collect(range(0.0,stop=10,step=0.125))
	Ih_guus = []
	Il_guus = []
	Ih_gus = []
	Il_gus = []
	
	guus_root = p_root[10]
	gus_root = p_root[9]
	
	for i=1:length(guus)
		p = change_guus_(guus[i],p_root)
		Ih,Il = limits_I_effect(p,Iapp)
		append!(Ih_guus,Ih)
		append!(Il_guus,Il)
	end
	Iapp_=collect(range(0.0,stop=15,step=0.2))
	for i=1:length(gus)
		p = change_gus_(gus[i],p_root)
		Ih,Il = limits_I_effect(p,Iapp_)
		append!(Ih_gus,Ih)
		append!(Il_gus,Il)
	end
	
	d_guus_root = abs.(guus.- guus_root)
	i_approx_guus_root = findfirst(x -> x == minimum(d_guus_root),d_guus_root)
	
	#==
	plot_guus = plot(guus,Ih_guus,label="Ih")
	plot!(guus,Il_guus,label="Il")
	scatter!([guus(i_approx_guus_root),guus(i_approx_guus_root)],[Ih_guus(i_approx_guus_root),Il_guus(i_approx_guus_root)],label = "guus root")
	xaxis!("guus")
	yaxis!("I")
	title!("gus = $(gus_root)")
	
	d_gus_root = abs.(gus.- gus_root)
	i_approx_gus_root = findfirst(x -> x == minimum(d_gus_root),d_gus_root)
	
	plot_gus = plot(gus,Ih_gus,label="Ih")
	plot!(gus,Il_gus,label="Il")
	scatter!([gus(i_approx_gus_root),gus(i_approx_gus_root)],[Ih_gus(i_approx_gus_root),Il_gus(i_approx_gus_root)],label = "gus root")
	xaxis!("gus")
	yaxis!("I")
	title!("guus = $(guus_root)")
	
	plot_full = plot(plot_guus,plot_gus,layout=(1,2))
	return plot_full
	==#
	return gus,guus,Ih_gus,Il_gus,Ih_guus,Il_guus
end

# ╔═╡ 20f0f124-055c-44af-ae46-26507b04e3a7
gus_I,guus_I,Ih_gus,Il_gus,Ih_guus,Il_guus = Il_Ih_limits_variation(guus,gus,p_tonic_sh) 

# ╔═╡ 525bfa60-4254-4128-a090-ac63da977919
begin
	guus_root = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p))))))[10]
	gus_root = change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p))))))[9]
	
	d_guus_root = abs.(guus_I.- guus_root)
	i_approx_guus_root = findfirst(x -> x == minimum(d_guus_root),d_guus_root)
	
	
	plot_guus = plot(guus_I,Ih_guus,label="Ih")
	plot!(guus_I,Il_guus,label="Il")
	#scatter!([guus_I[i_approx_guus_root],guus_I[i_approx_guus_root]],[Ih_guus[i_approx_guus_root],Il_guus[i_approx_guus_root]],label = "gus root")
	xaxis!("gus")
	yaxis!("I")
	title!("gss = $(gus_root)")
	
	d_gus_root = abs.(gus_I.- gus_root)
	i_approx_gus_root = findfirst(x -> x == minimum(d_gus_root),d_gus_root)
	
	plot_gus = plot(gus,Ih_gus,label="Ih")
	plot!(gus_I,Il_gus,label="Il")
	#scatter!([gus_I[i_approx_gus_root],gus_I[i_approx_gus_root]],[Ih_gus[i_approx_gus_root],Il_gus[i_approx_gus_root]],label = "gss root")
	xaxis!("gss")
	yaxis!("I")
	title!("gus = $(guus_root)")
	
	plot_full = plot(plot_guus,plot_gus,layout=(1,2),xaxis=:log)
	#savefig("4DIhIlevolcond.pdf")
end

# ╔═╡ 3b5c89f8-f251-4016-943b-7cdb4b9cddb5
limits_I = limits_I_effect(change_gus_(gus_I[18],p_tonic_sh),Iapp)

# ╔═╡ bf971085-7df9-48e0-a126-f0fd527248c0
 limits_I_effect(change_I_(20.0,change_vs0_(-38.4,change_vus0_(-10,change_gs_(0.5,change_guus_(0.00015,change_gus_(0.015,p)))))),Iapp)

# ╔═╡ 1b4a0e49-b64c-4dd7-a03c-94e2611613f1
function full_frequency_(t,sim,tspan,tspan_pulse)
	#to use if we are sure we do not have burst 
	
	ind_spikeV = find_Vspikes_(t,sim[1,:],maximum(tspan),0) #on the full period
	delta_t_V = t[ind_spikeV[2:end]] - t[ind_spikeV[1:end-1]]
	
	ind_pulse_0 =findall(x-> x>=tspan_pulse[1],t)[1]
	ind_pulse_f =findall(x-> x>=tspan_pulse[2],t)[1]
	
	i_high_1_ = findall(x-> x>=ind_pulse_0,ind_spikeV)
	i_high_2_ = findall(x-> x<=ind_pulse_f,reverse(ind_spikeV))
	
	if length(i_high_1_)>=1
		i_high_1 = i_high_1_[1]
	else
		i_high_1 = -1
	end
	if length(i_high_2_)>=1
		i_high_2 =length(ind_spikeV) - (i_high_2_[1]-1)
	else
		i_high_2 = -1
	end	
	
	if i_high_1==-1 && i_high_2 ==-1
		n_spike = 0
	else
		n_spike= i_high_2-i_high_1+1 #to put during pulse
	end
	
	t_period_full =[]
	t_sampled_f_ins = []

	for j=1:length(delta_t_V)
		append!(t_period_full,delta_t_V[j])
		append!(t_sampled_f_ins,sum(t_period_full)+t[ind_spikeV[1]])
	end
		
	f_ins_full=frequency_from_period(t_period_full[:]).*1000
		
	#==f_ins = zeros(size(t))
	ind_t_period = [ind_spikeV[1]]
	for j=1:length(t_period_full)
		ind_t_period_ = findfirst(x -> x>(sum(t_period_full[1:j])+t[ind_spikeV[1]]),t)
	
		f_ins[ind_t_period[end]:ind_t_period_].= f_ins_full[j]
		append!(ind_t_period,ind_t_period_)
	end
	==#
	
	return f_ins_full,t_sampled_f_ins
end

# ╔═╡ 4471dd69-4a70-4f8c-aacf-589d26a53202
md"#### Pattern simulations" 

# ╔═╡ 97d078a4-2005-443b-ac35-3802c85dfa21
function simulate_cst(t0,tf,p,I0,i_exc)
	tspan0=(t0,t0_step)
	
	if i_exc == false
		u00 =[-60,-60,-60,-60]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-12,abstol=1e-12)

	I = ones(length(sol0.t)).*I0	
	t_I = sol0.t
	sim =zeros(4,length(sol0.t))
	sim[1,:] .= sol0[1,:]
	sim[2,:] .= sol0[2,:]
	sim[3,:] .= sol0[3,:]
	sim[4,:] .= sol0[4,:]

	return t_I,I,sim
end

# ╔═╡ 13b1f81a-cec8-4ba3-b2d7-4beaf3103c39
function simulate_step(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc)
	tspan0=(t0,t0_step)
	
	if i_exc == false
		u00 =[-40,-40,-40,-40]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan1=(t0_step,tf_step)
	u01 =[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	p1=change_I_(Istep,p)
	
	prob1 = ODEProblem(MQIF_4D!,u01,tspan1,p1,callback=cb)
	sol1 = solve(prob1,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan2=(tf_step,tf)
	u02 =[sol1[1,end],sol1[2,end],sol1[3,end],sol1[4,end]]
	p2=change_I_(I0,p)
	
	prob2 = ODEProblem(MQIF_4D!,u02,tspan2,p2,callback=cb)
	sol2 = solve(prob2,DP5(),reltol=1e-12,abstol=1e-12)

	I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))	
	t_I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))
	sim =zeros(4,length(sol0.t)+length(sol1.t)+length(sol2.t))
	
	for i=1:length(I)
		if i<=length(sol0.t)
			I[i] = I0
			t_I[i] = sol0.t[i]
			sim[1,i] = sol0[1,i]
			sim[2,i] = sol0[2,i]
			sim[3,i] = sol0[3,i]
			sim[4,i] = sol0[4,i]
		else
			if i<=length(sol0.t)+length(sol1.t)
				I[i] = Istep
				t_I[i] = sol1.t[i-length(sol0.t)]
				sim[1,i] = sol1[1,i-length(sol0.t)]
				sim[2,i] = sol1[2,i-length(sol0.t)]
				sim[3,i] = sol1[3,i-length(sol0.t)]
				sim[4,i] = sol1[4,i-length(sol0.t)]
			else
				I[i] = I0
				t_I[i] = sol2.t[i-length(sol0.t)-length(sol1.t)]
				sim[1,i] = sol2[1,i-length(sol0.t)-length(sol1.t)]
				sim[2,i] = sol2[2,i-length(sol0.t)-length(sol1.t)]	
				sim[3,i] = sol2[3,i-length(sol0.t)-length(sol1.t)]	
				sim[4,i] = sol2[4,i-length(sol0.t)-length(sol1.t)]	
			end
		end
	end
	return t_I,I,sim
end

# ╔═╡ b66d7b10-7426-4f2b-a4bf-b55772244507
function simulate_step_(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc,mod_ci)
	tspan0=(t0,t0_step)
	
	if i_exc == false
		u00 =[-40,-40,-40,-40]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	if mod_ci == true
		u00 =[-60,-60,-60,-60]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan1=(t0_step,tf_step)
	u01 =[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	p1=change_I_(Istep,p)
	
	prob1 = ODEProblem(MQIF_4D!,u01,tspan1,p1,callback=cb)
	sol1 = solve(prob1,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan2=(tf_step,tf)
	u02 =[sol1[1,end],sol1[2,end],sol1[3,end],sol1[4,end]]
	p2=change_I_(I0,p)
	
	prob2 = ODEProblem(MQIF_4D!,u02,tspan2,p2,callback=cb)
	sol2 = solve(prob2,DP5(),reltol=1e-12,abstol=1e-12)

	I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))	
	t_I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))
	sim =zeros(4,length(sol0.t)+length(sol1.t)+length(sol2.t))
	
	for i=1:length(I)
		if i<=length(sol0.t)
			I[i] = I0
			t_I[i] = sol0.t[i]
			sim[1,i] = sol0[1,i]
			sim[2,i] = sol0[2,i]
			sim[3,i] = sol0[3,i]
			sim[4,i] = sol0[4,i]
		else
			if i<=length(sol0.t)+length(sol1.t)
				I[i] = Istep
				t_I[i] = sol1.t[i-length(sol0.t)]
				sim[1,i] = sol1[1,i-length(sol0.t)]
				sim[2,i] = sol1[2,i-length(sol0.t)]
				sim[3,i] = sol1[3,i-length(sol0.t)]
				sim[4,i] = sol1[4,i-length(sol0.t)]
			else
				I[i] = I0
				t_I[i] = sol2.t[i-length(sol0.t)-length(sol1.t)]
				sim[1,i] = sol2[1,i-length(sol0.t)-length(sol1.t)]
				sim[2,i] = sol2[2,i-length(sol0.t)-length(sol1.t)]	
				sim[3,i] = sol2[3,i-length(sol0.t)-length(sol1.t)]	
				sim[4,i] = sol2[4,i-length(sol0.t)-length(sol1.t)]	
			end
		end
	end
	return t_I,I,sim
end

# ╔═╡ 5bc8263f-184c-4549-8f2b-681bdf697ffd
function simulate_step__(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc,mod_ci)
	tspan0=(t0,t0_step)
	
	if i_exc == false
		u00 =[-40,-40,-40,-40]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	if mod_ci == true
		u00 =[-60,-60,-60,-60]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan1=(t0_step,tf_step)
	u01 =[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	p1=change_I_(Istep,p)
	
	prob1 = ODEProblem(MQIF_4D!,u01,tspan1,p1,callback=cb)
	sol1 = solve(prob1,DP5(),reltol=1e-12,abstol=1e-12)
	
	tspan2=(tf_step,tf)
	u02 =[sol1[1,end],sol1[2,end],sol1[3,end],sol1[4,end]]
	p2=change_I_(I0,p)
	
	prob2 = ODEProblem(MQIF_4D!,u02,tspan2,p2,callback=cb)
	sol2 = solve(prob2,DP5(),reltol=1e-12,abstol=1e-12)

	I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))	
	t_I = zeros(length(sol0.t)+length(sol1.t)+length(sol2.t))
	sim =zeros(4,length(sol0.t)+length(sol1.t)+length(sol2.t))
	
	sim[1,1:length(sol0.t)].=sol0[1,:]
	sim[2,1:length(sol0.t)].=sol0[2,:]
	sim[3,1:length(sol0.t)].=sol0[3,:]
	sim[4,1:length(sol0.t)].=sol0[4,:]
	t_I[1:length(sol0.t)].=sol0.t[:]
	I[1:length(sol0.t)].=I0
	
	sim[1,1+length(sol0.t):length(sol0.t)+length(sol1.t)].=sol1[1,:]
	sim[2,1+length(sol0.t):length(sol0.t)+length(sol1.t)].=sol1[2,:]
	sim[3,1+length(sol0.t):length(sol0.t)+length(sol1.t)].=sol1[3,:]
	sim[4,1+length(sol0.t):length(sol0.t)+length(sol1.t)].=sol1[4,:]
	t_I[1+length(sol0.t):length(sol0.t)+length(sol1.t)].=sol1.t[:]
	I[1+length(sol0.t):length(sol0.t)+length(sol1.t)].=Istep
	
	sim[1,1+length(sol0.t)+length(sol1.t):end] .=sol2[1,:]
	sim[2,1+length(sol0.t)+length(sol1.t):end] .=sol2[2,:]
	sim[3,1+length(sol0.t)+length(sol1.t):end] .=sol2[3,:]
	sim[4,1+length(sol0.t)+length(sol1.t):end] .=sol2[4,:]
	t_I[1+length(sol0.t)+length(sol1.t):end].=sol2.t[:]
	I[1+length(sol0.t)+length(sol1.t):end].=I0
	
	return t_I,I,sim
end

# ╔═╡ 212a5853-8072-4614-bf6b-15c3029e26a9
function classify_pattern(t0,tf,step_duration,p,Il,Ih)
	t0_step = 5000
	tf_step = t0_step+step_duration
	
	tspan = (t0,tf)
	tspan_pulse = (t0_step,tf_step)
	
	u0 = [-60,-60,-60,-60]
	
	V_fp,stab_fp=global_I_bifurcation(change_I_(Il,p))
	#check for NaNs (no fixed points)
	ind_check = findall(x -> x >=0 || x<0, stab_fp)
	if length(ind_check)>0
		ind_stable = findfirst(x -> x <0, stab_fp)#cartesian index
		if length(ind_stable)>0
			u0=V_fp[ind_stable[2]].*[1,1,1,1]
		end
	end
	
	##stable ? 
	ph = change_I_(Ih,p)
	tspan=(0.0,7500.0)
	
	probh = ODEProblem(MQIF_4D!,u0,tspan,ph,callback=cb)
	solh = solve(probh,DP5(),reltol=1e-12,abstol=1e-12)
	
	conv = NaN
	
	if check_if_stable(solh,tspan)
		conv = 0 #stable
	else
		##burst ? 
		if check_if_burst(solh,tspan)
			conv = 1 #burst
		else
			pl = change_I_(Il,p)
			probl = ODEProblem(MQIF_4D!,u0,tspan,pl,callback=cb)
			soll = solve(probl,DP5())
			
			if check_if_burst(soll,tspan)
				conv = 1 #burst
			else
				t,I,sim = simulate_step__(t0,tf,t0_step,tf_step,p,Il,Ih,false,true)
				f_ins,t_s_f_ins = full_frequency(t,sim,tspan,tspan_pulse)
				
				if length(f_ins)<=1
					conv = 10 #stable due to spike latency
				else
					if check_if_unstable(t,sim,tspan)
						conv=-1 #unstable
					else
						plateau = check_if_plateau(t,sim,tspan,tspan_pulse)
						f_evol = freq_during_pulse(t,tspan_pulse,f_ins,t_s_f_ins)

						if plateau == false && f_evol==0
							conv=2 #tonic
						end
						if plateau == false && f_evol==1
							conv=3 #accelerating 
						end
						if plateau == true && f_evol==1
							conv=4 #plateau potentials 
						end
						if plateau == false && f_evol==-1
							conv=5 #decelerating 
						end
						if plateau == true && f_evol==-1
							conv=6 #decelerating with plateau 
						end
						if plateau == false && f_evol==2
							conv=7 #curved without plateau 
						end
						if plateau == true && f_evol==2
							conv=8 #curved with plateau 
						end
						if plateau == true && f_evol==0
							conv=9 #tonic with plateau 
						end
					end
				end
			end
		end
	end
	#later : return spike latency and delta frequency to have the magnitude of the acceleration/deceleration 
	return conv
end

# ╔═╡ cef9fb8f-8824-452a-b335-2e9f31471c92
begin
	conv = classify_pattern(0.0,10000.0,pulse_dur,p_level4,Il,Ih)
	translate_classification(conv)
end

# ╔═╡ 5606e351-b44c-46aa-b5d9-b7d5b43928da
function plot_cond_space_uus_us(p_list,plt,t0,tf,step_duration,Il,Ih)
	conv_mem = []
	
	for i=1:length(p_list)
		I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p_list[i]
		conv = classify_pattern(t0,tf,step_duration,p_list[i],Il,Ih)
		scatter!(plt,[guus],[gus],markercolor = color_classification(conv),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv))
		append!(conv_mem,conv)
	end
	xaxis!("guus")
	yaxis!("gus")
	
	
	
	plt_gs_guus = plot()
	for i=1:length(p_list)
		I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p_list[i]
		conv = conv_mem[i]
		scatter!(plt_gs_guus,[guus],[gs],markercolor = color_classification(conv),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv))
	end
	xaxis!("guus")
	yaxis!("gs")
	
	#sub = plot(plt,plt_gs_guus,layout=(1,2),legend=:outertopright,size=(1000,650))
	sub=plt
	
	return sub
end

# ╔═╡ 048a13f5-b5e2-46b6-b370-6c388dcf2d4d
function plot_cond_space_uus_us_(p_list,plt,t0,tf,step_duration,Il,Ih,check_gs,check_gus,check_guus)
	conv_mem = []
	class = [[],[],[],[],[],[],[],[],[],[],[],[]]
	
	for i=1:length(p_list)
		I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p_list[i]
		conv = classify_pattern(t0,tf,step_duration,p_list[i],Il,Ih)
		append!(conv_mem,conv)
		if conv != NaN
			append!(class[conv+2],i)
		end
	end
	
	
	for i=1:length(class)
		if length(class[i])>1
			X= []
			Y= []
			for j=1:length(class[i])
				if check_gus == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][9])
				end
				if check_gs == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][8])
				end
				if check_gs == true && check_gus ==true
					append!(X,p_list[class[i]][j][9])
					append!(Y,p_list[class[i]][j][8])
				end
			end
			scatter!(plt,X,Y,markercolor = color_classification(conv_mem[class[i][1]]),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv_mem[class[i][1]]))
			if check_gus == true && check_guus ==true
				xaxis!("guus")
				yaxis!("gus")
			end
			if check_gs == true && check_guus ==true
				xaxis!("guus")
				yaxis!("gs")
			end
			if check_gs == true && check_gus ==true
				xaxis!("gus")
				yaxis!("gs")
			end
		end
	end
	
	#sub = plot(plt,plt_gs_guus,layout=(1,2),legend=:outertopright,size=(1000,650))
	sub=plt
	
	return sub
end

# ╔═╡ 3a8f030f-4da3-4d3b-8fb2-000dd162823b
function Il_Ih_effect(p_root) 
	#Iapp=collect(range(0.0,stop=7,step=0.5))
	Iapp=collect(range(0.0,stop=20,step=0.25))
	I_limits = limits_I_effect(p_root,Iapp)
	Ih = ceil(I_limits[1])+1
	Il = floor(I_limits[2])
	
	I_test_h_ = collect(range(5.0,stop=25.0,step=1))
	I_test_l_ = collect(range(0.0,stop=20.0,step=1))
	
	I_test_h,I_test_l=meshgrid(I_test_h_,I_test_l_)
	conv_mem = []
	class = [[],[],[],[],[],[],[],[],[],[],[],[]]
	class_skip = [[],[],[],[],[],[],[],[],[],[],[],[]]
	test = []
	
	counter = zeros(size(I_test_h))
	plt =  plot()
	
	for i=1:length(I_test_h)
		if I_test_h[i]>I_test_l[i]
			conv = classify_pattern(0.0,10000.0,300.0,p_root,I_test_l[i],I_test_h[i])
			append!(conv_mem,conv)
			if conv != NaN
				if i>1
					append!(class_skip[conv+2],i-counter[i-1])
				else
					append!(class_skip[conv+2],i)
				end
				append!(class[conv+2],i)	
			end
		else
			if i>1
				counter[i] = counter[i-1]+1
			else
				counter[i]=i
			end
				
		end
		if I_test_h[i]==8.0 && I_test_l[i]==0.5
				append!(test,conv)
		end
	end
	
	plt = plot()
	for i=1:length(class)
		if length(class[i])>=1
			X= []
			Y= []
			for j=1:length(class[i])
				append!(X,I_test_l[class[i]][j])
				append!(Y,I_test_h[class[i]][j])
			end
			scatter!(X,Y,markercolor = color_classification(i-2),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,label=translate_classification(i-2))
			xaxis!("Il")
			yaxis!("Ih")
		end
	end
	
	conv = classify_pattern(0.0,10000.0,300.0,p_root,Il,Ih)
	#scatter!(plt,[Il],[Ih],markercolor = color_classification(conv),markershape = :star5,markersize = 10,markeralpha = 0.9,markerstrokewidth = 1,markerstrokealpha = 0.9, markerstrokecolor = :black,markerstrokestyle = :dot,label=translate_classification(conv)*", limit")
	==#
	xaxis!("Il")
	yaxis!("Ih")
	
	return plt,I_limits
end

# ╔═╡ 549bfbf4-db61-4887-9334-aa163a0153da
function plot_cond_space_(p_list,I_list,t0,tf,step_duration,Il,Ih,check_gs,check_gus,check_guus,check_Il,check_Ih)
	
	##to modify, calculate limits with each p 
	conv_mem = []
	class = [[],[],[],[],[],[],[],[],[],[],[],[]]
	
	for i=1:length(p_list)
		I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p_list[i]
		if check_Ih
			Ih = I_list[i]
		end
		if check_Il
			Il = I_list[i]
		end
		conv = classify_pattern(t0,tf,step_duration,p_list[i],Il,Ih)
		append!(conv_mem,conv)
		if conv != NaN
			append!(class[conv+2],i)
		end
	end
	
	
	for i=1:length(class)
		if length(class[i])>1
			X= []
			Y= []
			for j=1:length(class[i])
				if check_gus == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][9])
				end
				if check_gs == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][8])
				end
				if check_gs == true && check_gus ==true
					append!(X,p_list[class[i]][j][9])
					append!(Y,p_list[class[i]][j][8])
				end
				if check_Il == true && check_gus ==true
					append!(Y,p_list[class[i]][j][9])
					append!(X,I_list[class[i]][j])
				end
				if check_Ih == true && check_gus ==true
					append!(Y,p_list[class[i]][j][9])
					append!(X,I_list[class[i]][j])
				end
				if check_Il == true && check_guus ==true
					append!(Y,p_list[class[i]][j][10])
					append!(X,I_list[class[i]][j])
				end
				if check_Ih == true && check_guus ==true
					append!(Y,p_list[class[i]][j][10])
					append!(X,I_list[class[i]][j])
				end
			end
			scatter!(plt,X,Y,markercolor = color_classification(conv_mem[class[i][1]]),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv_mem[class[i][1]]))
			if check_gus == true && check_guus ==true
				xaxis!("guus")
				yaxis!("gus")
			end
			if check_gs == true && check_guus ==true
				xaxis!("guus")
				yaxis!("gs")
			end
			if check_gs == true && check_gus ==true
				xaxis!("gus")
				yaxis!("gs")
			end
			if check_Il == true && check_gus ==true
				xaxis!("Il")
				yaxis!("gus")
			end
			if check_Ih == true && check_gus ==true
				xaxis!("Ih")
				yaxis!("gus")
			end
			if check_Il == true && check_guus ==true
				xaxis!("Il")
				yaxis!("guus")
			end
			if check_Ih == true && check_guus ==true
				xaxis!("Ih")
				yaxis!("guus")
			end
		end
	end
	
	#sub = plot(plt,plt_gs_guus,layout=(1,2),legend=:outertopright,size=(1000,650))
	sub=plt
	
	return sub
end

# ╔═╡ 00505f0c-12a6-4144-b4a8-46c3d9ead091
function classify_pattern_(t0,tf,step_duration,p,Il,Ih)
	t0_step = 5000
	tf_step = t0_step+step_duration
	
	tspan = (t0,tf)
	tspan_pulse = (t0_step,tf_step)
	
	u0 = [-60,-60,-60,-60]
	
	V_fp,stab_fp=global_I_bifurcation(change_I_(Il,p))
	#check for NaNs (no fixed points)
	ind_check = findall(x -> x >=0 || x<0, stab_fp)
	if length(ind_check)>0
		ind_stable_ = findall(x -> x <0, stab_fp)#cartesian index
		if length(ind_stable_)>0
			ind_stable = findfirst(x -> x <0, stab_fp)#cartesian index
			u0=V_fp[ind_stable[2]].*[1,1,1,1]
		end
	end
	
	##stable ? 
	ph = change_I_(Ih,p)
	tspan=(0.0,7500.0)
	
	probh = ODEProblem(MQIF_4D!,u0,tspan,ph,callback=cb)
	solh = solve(probh,DP5(),reltol=1e-12,abstol=1e-12)
	
	conv = NaN
	
	if check_if_stable(solh,tspan)
		conv = 0 #stable
	else
		##burst ? 
		if check_if_burst(solh,tspan)
			conv = 1 #burst
		else
			pl = change_I_(Il,p)
			probl = ODEProblem(MQIF_4D!,u0,tspan,pl,callback=cb)
			soll = solve(probl,DP5())
			
			if check_if_burst(soll,tspan)
				conv = 1 #burst
			else
				t,I,sim = simulate_step__(t0,tf,t0_step,tf_step,p,Il,Ih,false,true)
				f_ins,t_s_f_ins = full_frequency(t,sim,tspan,tspan_pulse)
				
				if length(f_ins)<=1
					conv = 10 #stable due to spike latency
				else
					if check_if_unstable(t,sim,tspan)
						conv=-1 #unstable
					else
						plateau = check_if_plateau(t,sim,tspan,tspan_pulse)
						f_evol = freq_during_pulse(t,tspan_pulse,f_ins,t_s_f_ins)

						if plateau == false && f_evol==0
							conv=2 #tonic
						end
						if plateau == false && f_evol==1
							conv=3 #accelerating 
						end
						if plateau == true && f_evol==1
							conv=4 #plateau potentials 
						end
						if plateau == false && f_evol==-1
							conv=5 #decelerating 
						end
						if plateau == true && f_evol==-1
							conv=6 #decelerating with plateau 
						end
						if plateau == false && f_evol==2
							conv=7 #curved without plateau 
						end
						if plateau == true && f_evol==2
							conv=8 #curved with plateau 
						end
						if plateau == true && f_evol==0
							conv=9 #tonic with plateau 
						end
					end
				end
			end
		end
	end
	#later : return spike latency and delta frequency to have the magnitude of the acceleration/deceleration 
	return conv
end

# ╔═╡ b97e4725-8b6d-46d5-ad95-f2ae8180079c
function plot_cond_space_I_(p_list,I_list,plt,t0,tf,step_duration,Il,Ih,check_gs,check_gus,check_guus,check_Ih)
	conv_mem = []
	class = [[],[],[],[],[],[],[],[],[],[],[],[]]
	
	for i=1:length(p_list)
		I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p_list[i]
		if check_Ih
			Ih_ = ceil(I_list[i])+1
		else
			Ih_ = Ih
		end
		conv = classify_pattern_(t0,tf,step_duration,p_list[i],Il,Ih_)
		append!(conv_mem,conv)
		if conv != NaN
			append!(class[conv+2],i)
		end
	end
	
	
	for i=1:length(class)
		if length(class[i])>1
			X= []
			Y= []
			for j=1:length(class[i])
				if check_gus == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][9])
				end
				if check_gs == true && check_guus ==true
					append!(X,p_list[class[i]][j][10])
					append!(Y,p_list[class[i]][j][8])
				end
				if check_gs == true && check_gus ==true
					append!(X,p_list[class[i]][j][9])
					append!(Y,p_list[class[i]][j][8])
				end
			end
			scatter!(plt,X,Y,markercolor = color_classification(conv_mem[class[i][1]]),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv_mem[class[i][1]]))
			if check_gus == true && check_guus ==true
				xaxis!("gus")
				yaxis!("gss")
			end
			if check_gs == true && check_guus ==true
				xaxis!("gus")
				yaxis!("gs")
			end
			if check_gs == true && check_gus ==true
				xaxis!("gss")
				yaxis!("gs")
			end
		end
	end
	
	#sub = plot(plt,plt_gs_guus,layout=(1,2),legend=:outertopright,size=(1000,650))
	sub=plt
	
	return sub
end

# ╔═╡ 0d9d96e0-6352-4fb4-9ebb-29ec39ee8ea0
full_Ih = plot_cond_space_I_(a,b,plot(),0.0,10000.0,300.0,0.5,NaN,true,false,true,true)

# ╔═╡ 52e1ce2b-478d-40c3-8133-841f4de50acb
begin
	plot(full_Ih,legend=:outertopright,yaxis=:none)
	conv_80 = classify_pattern_(0.0,10000.0,300.0,a[80],0.5,ceil(b[80])+1.0)
	conv_66 = classify_pattern_(0.0,10000.0,300.0,a[66],0.5,ceil(b[66])+1.0)
	scatter!([a[80][10]],[a[80][8]],markercolor = color_classification(conv_80),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:log,label=translate_classification(conv_80))
	scatter!([a[66][10]],[a[66][8]],markercolor = color_classification(conv_66),markershape = :circle,markersize = 10,markeralpha = 0.9,markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot,xaxis=:log,yaxis=:none,label=translate_classification(conv_66))
	#savefig("MQIF_4D_cond_Ih_limit_ceilplus1_guus_gs_2percent.pdf")
end

# ╔═╡ 254a738a-b1f4-4890-897d-bd810027fcac
function subplot_simu_step_full(t0,tf,step_duration,p,I0,Istep,both,i_exc,mod_ci)
	t0_step = 5000
	tf_step = t0_step+step_duration
	perc_tf=0.05
	
	t_I,I,sim = simulate_step__(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc,mod_ci)
	
	
	Iuus = compute_Iuus(sim[4,:],p)
	Ius = compute_Ius(sim[3,:],p)
	It = I.+Iuus.+Ius
	
	lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:]),minimum(sim[4,:])])-10
	lim2 = maximum([maximum(sim[1,:])+20,maximum(sim[3,:])])
	
	if both==false
		tspan = (t0,tf)
		tspan_pulse = (t0_step,tf_step)
		 f_ins_full,t_sampled_f_ins = full_frequency(t_I,sim,tspan,tspan_pulse)
		
		plot_step = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)           ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)",(minimum(t_I),maximum(t_I)))	

		plot_step_z = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)           ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		plot_freq = scatter(t_sampled_f_ins,f_ins_full,label="Frequency" ,markershape = :circle,markersize = 2,markeralpha = 0.9,markercolor = RGB(0.5,0.5,1),markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot)
		yaxis!("f [Hz]")	
		xaxis!("Time (ms)",(minimum(t_I),maximum(t_I)))	
		
		plot_freq_z = scatter(t_sampled_f_ins,f_ins_full,label="Frequency",markershape = :circle,markersize = 2,markeralpha = 0.9,markercolor = RGB(0.5,0.5,1),markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot)
		yaxis!("f [Hz]")	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		plot_v = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)     ")
		yaxis!("Voltage",(lim1,lim2))	
		xaxis!("Time (ms)",(minimum(t_I),maximum(t_I))	)	
		
		plot_zoom = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)     ")
				#plot!(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)")
		yaxis!("Voltage",(lim1,lim2))	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		
		sub = plot(plot_step,plot_v,plot_freq,plot_step_z,plot_zoom,plot_freq_z,layout=@layout([a;b;c;d;e;f ]),linewidth = 1.5,legend=:outertopright)
	else
		sub = plot()
		#sub = plot(plot_step,plot_v,title_zoom,plot_step_z,plot_zoom,layout=@layout([a{0.05h};b;c{0.01h};d{0.05h};e ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ 318085c4-2b3d-4970-a65a-71300f6d6743
begin
	title_4 = plot(title = translate_classification(conv), grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=20,titlefontcolor=color_classification(conv))
	full_4 = subplot_simu_step_full(0.0,20000.0,pulse_dur,p_level4,Il,Ih,false,false,true)
	plot(title_4,full_4,
	layout = @layout([A{0.03h} ;B ]),size=(800,800))
	#savefig("MQIF_4D_sim_stable.png")
end

# ╔═╡ 717dcd6b-9058-450d-b545-6c7e280d7cd3
begin
	p_test =a[80]
	#limits_I_ = limits_I_effect(p_test,Iapp)
	conv_test = classify_pattern_(0.0,10000.0,300.0,p_test,0.5,ceil(b[80])+1.0)
	translate_classification(conv_test)

	title_test = plot(title = translate_classification(conv_test), grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=20,titlefontcolor=color_classification(conv_test))
	full_test = subplot_simu_step_full(0.0,10000.0,300.0,p_test,0.5,ceil(b[80])+1,false,false,true)
	plot(title_test,full_test,
	layout = @layout([A{0.03h} ;B ]),size=(800,800))
	#savefig("MQIF_4D_sim_accel_300_0008_000015.pdf")
end

# ╔═╡ 1b6fca49-4a8a-4a92-a9ab-7aa33fd7f0cb
function subplot_simu_step_(t0,tf,step_duration,p,I0,Istep,both,i_exc,mod_ci)
	t0_step = 5000
	tf_step = t0_step+step_duration
	perc_tf=0.05
	
	t_I,I,sim = simulate_step__(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc,mod_ci)
	
	
	Iuus = compute_Iuus(sim[4,:],p)
	Ius = compute_Ius(sim[3,:],p)
	It = I.+Iuus.+Ius
	
	lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:]),minimum(sim[4,:])])-10
	lim2 = maximum([maximum(sim[1,:])+20,maximum(sim[3,:])])
	
	if both==false
		tspan = (t0,tf)
		tspan_pulse = (t0_step,tf_step)
		 f_ins_full,t_sampled_f_ins = full_frequency(t_I,sim,tspan,tspan_pulse)
		
		plot_step = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)           ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)",(minimum(t_I),maximum(t_I)))	

		plot_step_z = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)           ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		plot_freq = scatter(t_sampled_f_ins,f_ins_full,label="Frequency" ,markershape = :circle,markersize = 2,markeralpha = 0.9,markercolor = RGB(0.5,0.5,1),markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot)
		yaxis!("f [Hz]")	
		xaxis!("Time (ms)",(minimum(t_I),maximum(t_I)))	
		
		plot_freq_z = scatter(t_sampled_f_ins,f_ins_full,label="Frequency",markershape = :circle,markersize = 2,markeralpha = 0.9,markercolor = RGB(0.5,0.5,1),markerstrokewidth = 0,markerstrokealpha = 0., markerstrokecolor = :black,markerstrokestyle = :dot)
		yaxis!("f [Hz]")	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		#sub = plot(plot_step,plot_freq,plot_step_z,plot_freq_z,layout=@layout([a;b;c;d ]),linewidth = 1.5,legend=:outertopright)
		sub = plot(plot_step,plot_freq,plot_step_z,plot_freq_z,layout=@layout([a;b;c;d ]),linewidth = 1.5,legend=:outertopright)
	else
		plot_step = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)      ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)")	

		plot_step_z = plot(t_I,I,linecolor=RGB(0.7,0,0.1),linealpha=0.8,label="I(t)      ")
		plot!(t_I,It,linecolor=RGB(1,0.5,0.7),linewidth=1,label="It(t)     ")
		plot!(t_I,find_local_bif(p).*ones(size(t_I)),linecolor=RGB(1,0.6,0.4), linewidth=1,linestyle=:dot,label="2D bif  ")
		yaxis!("I")	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		plot_v = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Voltage",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_zoom = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
				#plot!(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)")
		yaxis!("Voltage",(lim1,lim2))	
		xaxis!("Time (ms)",(tf_step-tf*perc_tf,tf_step +tf*perc_tf))	
		
		title_zoom = plot(title = "Zoom", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:right,titlefontsize=14,titlefontcolor=RGB(0,0.4,0.95))

		sub = plot(plot_step,plot_v,plot_step_z,plot_zoom,layout=@layout([a{0.1h};b;c{0.1h};d ]),linewidth = 1.5,legend=:outertopright)
		#sub = plot(plot_step,plot_v,title_zoom,plot_step_z,plot_zoom,layout=@layout([a{0.05h};b;c{0.01h};d{0.05h};e ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ d6e677ab-eb7d-4097-be0c-bd25895d71be
function subplot_simu_step__(t0,tf,step_duration,p,I0,Istep,both,i_exc,mod_ci)
	subplot_simu_step_(t0,tf,step_duration,p,I0,Istep,both,i_exc,mod_ci)
	t0_step = 5000
	tf_step = t0_step+step_duration
	perc_tf=0.05
	
	t_I,I,sim = simulate_step_(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc,mod_ci)
	
	
	Iuus = compute_Iuus(sim[4,:],p)
	Ius = compute_Ius(sim[3,:],p)
	It = I.+Iuus.+Ius
	
	lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:]),minimum(sim[4,:])])-10
	lim2 = maximum([maximum(sim[1,:])+20,maximum(sim[3,:])])
	
	if both==false
		tspan = (t0,tf)
		tspan_pulse = (t0_step,tf_step)
		f_ins_full,t_sampled_f_ins = full_frequency(t_I,sim,tspan,tspan_pulse)
		
		return f_ins_full,t_sampled_f_ins
	end
end

# ╔═╡ 77e13ad3-7e1f-41e7-b1b2-d3b10778bb6f
function subplot_simu_step(t0,tf,step_duration,p,I0,Istep,both,i_exc)
	t0_step = round((tf-t0)/2)-round(step_duration/2)
	tf_step = round((tf-t0)/2)+round(step_duration/2)
	
	t_I,I,sim = simulate_step(t0,tf,t0_step,tf_step,p,I0,Istep,i_exc)
	
	plot_step = plot(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)      ")
	yaxis!("I")	
	xaxis!("Time (ms)")	
	
	plot_step_z = plot(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)      ")
	yaxis!("I")	
	xaxis!("Time (ms)",(t0_step -tf*0.01,tf_step +tf*0.01))	
	
	if both==false
		plot_v = plot(t_I,sim[1,:],label="V(t)")
		lim1 = minimum([minimum(sim[1,:]),minimum(sim[3,:]),minimum(sim[4,:])])
		lim2 = maximum([maximum(sim[1,:]),maximum(sim[4,:])])
		yaxis!("V",(lim1,lim2))	
		xaxis!("Time (ms)")	

		plot_vs = plot(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
		yaxis!("Vs",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_vus = plot(t_I,sim[3,:],linecolor=RGB(0,0.7,0.1),label="Vus(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	
		
		plot_vuus = plot(t_I,sim[4,:],linecolor=RGB(0,0.7,0.1),label="Vuus(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_step,plot_v,plot_vs,plot_vus,plot_vuus,layout=@layout([a{0.2h};b ;c;d;e]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	
		
		plot_zoom = plot(t_I,sim[1,:],label="V(t)")
				plot!(t_I,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t_I,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t_I,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
				#plot!(t_I,I,linecolor=RGB(0.7,0,0.1),label="I(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)",(t0_step -tf*0.01,tf_step +tf*0.01))	
		
		title_zoom = plot(title = "Zoom", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:right,titlefontsize=14,titlefontcolor=RGB(0,0.4,0.95))

		sub = plot(plot_step,plot_v,plot_step_z,plot_zoom,layout=@layout([a{0.1h};b;c{0.1h};d ]),linewidth = 1.5,legend=:outertopright)
		#sub = plot(plot_step,plot_v,title_zoom,plot_step_z,plot_zoom,layout=@layout([a{0.05h};b;c{0.01h};d{0.05h};e ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ 4106e96b-9a53-4a44-a063-97837d57bfc8
function simu_pulse(t0,tf,t_spike,n_spike,p,I0,Ispike,spike_duration,i_exc)
	tspan0=(t0,t_spike)
	
	if i_exc == false
		u00 =[-40,-40,-40,-40]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-6,abstol=1e-6)
	u0s=[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	
	v = []
	vs = []
	vus = []
	vuus = []
	t = []
	I = []
	append!(v,sol0[1,:])
	append!(vs,sol0[2,:])
	append!(vus,sol0[3,:])
	append!(vuus,sol0[4,:])
	append!(t,sol0.t)
	append!(I,I0.*ones(length(sol0.t)))
	
	for i=1:n_spike
		tspan_spike=(i*t_spike,i*t_spike+spike_duration)
		p_spike = change_I_(Ispike,p)
		
		probs = ODEProblem(MQIF_4D!,u0s,tspan_spike,p_spike,callback=cb)
		sols = solve(probs,DP5(),reltol=1e-6,abstol=1e-6)
		u0r = [sols[1,end],sols[2,end],sols[3,end],sols[4,end]]
		
		append!(v,sols[1,:])
		append!(vs,sols[2,:])
		append!(vus,sols[3,:])
		append!(vuus,sols[4,:])
		append!(t,sols.t)
		append!(I,Ispike.*ones(length(sols.t)))
		
		if (i+1)*t_spike<=tf
			tspan_rest=(i*t_spike+spike_duration,(i+1)*t_spike)
		else
			tspan_rest=(i*t_spike+spike_duration,tf)
		end
		
		probr = ODEProblem(MQIF_4D!,u0r,tspan_rest,p0,callback=cb)
		solr = solve(probr,DP5(),reltol=1e-6,abstol=1e-6)
		
		u0s =[solr[1,end],solr[2,end],solr[3,end],solr[4,end]]
		
		append!(v,solr[1,:])
		append!(vs,solr[2,:])
		append!(vus,solr[3,:])
		append!(vuus,solr[4,:])
		append!(t,solr.t)
		append!(I,I0.*ones(length(solr.t)))
	end
	
	sim=zeros(4,length(t))
	sim[1,:] = v
	sim[2,:] = vs
	sim[3,:] = vus
	sim[4,:] = vuus
	
	return t,I,sim
end

# ╔═╡ 043c2115-ea95-4b53-877a-6e00ff92f972
function subplot_simu_spike(t0,tf,t_spike,p,I0,Ispike,spike_duration,both,i_exc)
	n_spike = floor((tf-t0)/t_spike)-1
	
	t,I,sim = simu_pulse(t0,tf,t_spike,n_spike,p,I0,Ispike,spike_duration,i_exc)
	
	plot_spike = plot(t,I,linecolor=RGB(0.7,0,0.1),label="I(t)      ")
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
		
		plot_vuus = plot(t,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_spike,plot_v,plot_vs,plot_vus,plot_vuus,layout=@layout([a{0.1h};b ;c;d;e]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t,sim[1,:],label="V(t)")
				plot!(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	

		sub = plot(plot_spike,plot_v,layout=@layout([a{0.1h};b ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ 9a565434-315f-45d0-b64b-d8f21bf76288
function simu_ud_(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,i_exc)
	tspan0=(t0,t_u[1])

	if i_exc == false
		u00 =[-40,-40,-40,-40]
	else
		u00 =[-10,-20,-30,-50]
	end
	
	p0=change_I_(I0,p)
	
	prob0 = ODEProblem(MQIF_4D!,u00,tspan0,p0,callback=cb)
	sol0 = solve(prob0,DP5(),reltol=1e-6,abstol=1e-6)
	u0u=[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	
	v = []
	vs = []
	vus = []
	vuus = []
	t = []
	I = []
	append!(v,sol0[1,:])
	append!(vs,sol0[2,:])
	append!(vus,sol0[3,:])
	append!(vuus,sol0[3,:])
	append!(t,sol0.t)
	append!(I,I0.*ones(length(sol0.t)))
	
	for i=1:length(t_u)
		## up
		tspan_u=(t_u[i],t_u[i]+spike_duration)
		p_u = change_I_(Ispike,p)
		
		probu = ODEProblem(MQIF_4D!,u0u,tspan_u,p_u,callback=cb)
		solu = solve(probu,DP5(),reltol=1e-5,abstol=1e-5)
		u0r = [solu[1,end],solu[2,end],solu[3,end],solu[4,end]]
		
		append!(v,solu[1,:])
		append!(vs,solu[2,:])
		append!(vus,solu[3,:])
		append!(vuus,solu[4,:])
		append!(t,solu.t)
		append!(I,Ispike.*ones(length(solu.t)))
		
		## rest
		tspan_rest=(t_u[i]+spike_duration,t_d[i])
		
		probr = ODEProblem(MQIF_4D!,u0r,tspan_rest,p0,callback=cb)
		solr = solve(probr,DP5(),reltol=1e-6,abstol=1e-6)
		
		u0d =[solr[1,end],solr[2,end],solr[3,end],solr[4,end]]
		
		append!(v,solr[1,:])
		append!(vs,solr[2,:])
		append!(vus,solr[3,:])
		append!(vuus,solr[4,:])
		append!(t,solr.t)
		append!(I,I0.*ones(length(solr.t)))
		
		## down
		tspan_d =(t_d[i],t_d[i]+spike_duration)
		p_d = change_I_(-Ispike,p)
		
		probd = ODEProblem(MQIF_4D!,u0d,tspan_d,p_d,callback=cb)
		sold = solve(probd,DP5(),reltol=1e-6,abstol=1e-6)
		u0r = [sold[1,end],sold[2,end],sold[3,end],sold[4,end]]
		
		append!(v,sold[1,:])
		append!(vs,sold[2,:])
		append!(vus,sold[3,:])
		append!(vuus,sold[4,:])
		append!(t,sold.t)
		append!(I,(-Ispike).*ones(length(sold.t)))
		
		
		if t_d[i]+spike_duration < tf
			if i<length(t_u)
				tspan_rest=(t_d[i]+spike_duration,t_u[i+1])
			else
				tspan_rest=(t_d[i]+spike_duration,tf)
			end

			probr = ODEProblem(MQIF_4D!,u0r,tspan_rest,p0,callback=cb)
			solr = solve(probr,DP5(),reltol=1e-6,abstol=1e-6)

			u0u =[solr[1,end],solr[2,end],solr[3,end],solr[4,end]]

			append!(v,solr[1,:])
			append!(vs,solr[2,:])
			append!(vus,solr[3,:])
			append!(vuus,solr[4,:])
			append!(t,solr.t)
			append!(I,I0.*ones(length(solr.t)))
		end
	end
	
	sim=zeros(4,length(t))
	sim[1,:] = v
	sim[2,:] = vs
	sim[3,:] = vus
	sim[4,:] = vuus
	
	return t,I,sim
end

# ╔═╡ 00a3e093-719d-4205-a439-3d99c98fc52b
function subplot_simu_ud(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,both,i_exc)

	t,I,sim = simu_ud_(t0,tf,t_u,t_d,p,I0,Ispike,spike_duration,i_exc)
	
	plot_ud = plot(t,I,linecolor=RGB(0.7,0,0.1),label="I(t)      ")
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
		
		plot_vuus = plot(t,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Vus",(lim1,lim2))	
		xaxis!("Time (ms)")	

		sub = plot(plot_ud,plot_v,plot_vs,plot_vus,plot_vuus,layout=@layout([a{0.2h};b ;c;d;e]),linewidth = 1.5,legend=:outertopright)
	else
		plot_v = plot(t,sim[1,:],label="V(t)")
				plot!(t,sim[2,:],linecolor=RGB(0,0.7,0.1),label="Vs(t)")
				plot!(t,sim[3,:],linecolor=RGB(1,0.3,0.1),label="Vus(t)")
				plot!(t,sim[4,:],linecolor=RGB(0.58,0.34,0.69),label="Vuus(t)")
		yaxis!("Voltage")	
		xaxis!("Time (ms)")	

		sub = plot(plot_ud,plot_v,layout=@layout([a{0.2h};b ]),linewidth = 1.5,legend=:outertopright)
	end
	
	return sub
end

# ╔═╡ Cell order:
# ╟─eb93ed6e-a1f2-11eb-00a6-d5bbd921a4be
# ╠═eb47cf0d-f73a-4719-afda-3188b2998e52
# ╟─b9b842e4-7861-4bc5-9956-2db1d09c1c7d
# ╠═d49d2169-a568-4319-9287-d13044c30bd4
# ╟─78d2b9c6-9953-49dd-81a2-e1f22c887d30
# ╟─97629afb-5289-422a-bb17-5c1a0afebc00
# ╟─94979c24-7cc6-4c67-a0fa-13e0ea2c51f7
# ╟─3ae16e48-9afb-4d31-9123-cf6edb092b51
# ╟─9a5b6a2b-4080-4128-aabb-c1846eab7a40
# ╟─c7930bda-a02c-4ec1-a975-83a61de2592a
# ╟─042c10bf-0755-48f0-be92-f7451765a816
# ╟─d21074d8-6462-4ec1-b2a2-dd46028f675c
# ╟─139362f6-4a22-48a0-8d35-219deb803372
# ╟─58cbec76-b9bd-4e3c-bd7c-295bfd036f5b
# ╟─25553d60-2155-45ad-a6cd-d7bfb247e7c8
# ╟─be4c0e03-30c1-4b24-b491-4c8217ee27c8
# ╟─c5bc2248-002b-4da9-9640-42db8bca0ee6
# ╟─774f5bdf-863d-48f7-8115-625a06c171a3
# ╠═1b9fd190-0614-4730-b9f4-b48619abd92d
# ╟─4d024409-a4f4-4edd-9039-988e3c37c09d
# ╟─7a533a24-1a7d-4502-bfc9-075ed377866d
# ╟─455c242e-c3fc-4b2f-ae99-a15897ff4613
# ╟─9d5d2669-dba9-4c00-9b5c-5e27428d05b1
# ╟─acec6140-f6b1-446f-8401-777d78c648e5
# ╟─507962d3-b77c-478a-b39b-ebb75c0e98c5
# ╟─4e0b1fd0-107d-42c1-9e9c-c0fb0e7de3df
# ╟─cbae9324-dd13-4769-a381-8089c1114a98
# ╟─50a836a6-c023-49ab-8aad-79983dd7c2b8
# ╟─1cfc2491-8e71-42ac-a8c3-5001adbe3597
# ╟─5f609d44-ea73-4c26-95e8-1158e0d55e44
# ╟─d1c422a3-82ab-42b9-8a17-be40b3f8664b
# ╟─1c92ffd6-3284-444e-bbe6-22103df0dbff
# ╟─9fa7fb2b-2d9c-4775-9a02-f98bb27133d9
# ╟─8d571a7a-78ac-4841-8ae1-6c4a6adcc650
# ╟─fe3f6a22-5450-4983-bfa2-d9440f72e4f3
# ╟─8172c9ee-117c-468a-9c37-71f1de0708dd
# ╟─07633dbf-9633-45ef-84a9-3f839122386f
# ╟─120b0cc1-d150-4720-97dd-9e473117bd38
# ╟─02a4331e-ad81-49fa-bcc5-b2131af42691
# ╟─f1226b78-526e-46fb-96d4-a7da54c6126e
# ╟─41fcfdfe-0f90-4362-9c95-c8ff6201f643
# ╟─9d241efb-d465-4ba1-987a-f509ba0b81aa
# ╟─cb1e1b89-f2c3-4f3a-ba56-dd1c4f6d2c98
# ╟─7530943d-039a-4467-9c7d-8f565c5e5298
# ╟─0d58c293-3a2b-4518-ac0f-f26620ca4e44
# ╟─36f2a546-2716-494c-ad21-09c21f78cb8a
# ╟─411e46b8-8ebc-4ef7-bb0c-688020416e5b
# ╟─44d2c44b-44a3-490d-a933-12f74d0ef5b1
# ╟─47bccea1-4cb9-40ac-a5dc-62322512ce1a
# ╠═ae088c53-ab45-4dd8-8e01-384940acb491
# ╟─be306caa-d310-4a27-bfcf-16b5cff5590b
# ╟─006de583-9601-432c-95a7-a3f786b0c937
# ╟─61645246-ac2d-4d11-9252-ab4929841135
# ╟─20c01555-448a-48b8-ad43-345ad4c82c34
# ╟─a7234266-83ee-4f54-b3a8-3d2ada1c2a00
# ╟─7ab02730-0c95-445f-a1f3-f3d645776f9b
# ╟─3e13d7ce-1056-4562-86ad-412598a9eadf
# ╟─0ef6fa1a-efe7-441c-8360-9afba715e066
# ╟─db5df20a-fea3-4f51-827f-5ff657bd21a4
# ╟─17856e2b-0541-42fb-b393-9858be2411be
# ╟─35a023dc-16c4-4054-903a-1659474ac2a5
# ╟─e784cdff-fab4-4344-883f-6fcc17af8ad0
# ╟─80acbb79-335f-4f40-90a5-fe3a03c937aa
# ╟─14e0a258-d957-45cc-8543-9203fea533be
# ╟─67122038-ec98-4c9f-9d18-51f02e28a183
# ╟─33117e03-5d89-47fa-90a5-e89dc2071512
# ╟─4f69ac27-61e1-4ec2-92ce-4bd3b9d49bf8
# ╟─0d082a4b-9fbf-4045-b710-1c2930c1c530
# ╟─254f1c1d-5f79-4078-b686-51f2d1f6f1a3
# ╟─316a370c-e8bf-4da6-8782-381fba2e0d61
# ╟─4f6c0eb5-819d-4416-9449-e74f95c373f0
# ╟─2a6eb3be-ce0f-4e61-ae4f-d3e4abfa670e
# ╟─7f690094-96ac-439a-bf4a-8023cfdff5a0
# ╟─de46282f-4ce1-44f8-b645-c23ee550be02
# ╟─20bfa55e-0ce0-46ae-8e7b-dc28004d5a38
# ╟─a2b418be-4fb6-4abf-95aa-f914d29ee325
# ╟─04ed6d51-354a-4432-904d-76a72d7fa322
# ╟─57202474-8e51-43b3-8a71-17cc786aead6
# ╟─61bff036-322e-4e53-9114-c99c6457e115
# ╟─22a4c6d3-ab80-4f99-9b04-8f467edaa965
# ╟─7d63576b-39f7-41e0-a82b-ad570d31f8e5
# ╟─6146153e-49d6-4b51-9b03-00c56c5c08c8
# ╟─bcd46d0c-4738-4885-b6d4-ffe13143fce1
# ╟─26dd3baf-c6ed-47c1-8337-69a0589a237a
# ╟─b77a655e-4d54-4e3c-8bdb-a557141aa774
# ╟─f1789d85-45bc-4544-a175-197225e8ed6f
# ╟─481602fa-f9e9-450f-a831-f0caf68eca54
# ╟─621533e4-1c67-4436-baf8-8d55598c1184
# ╟─4c197fd0-9ae6-447f-90d4-4e6f395613fd
# ╟─7307d524-5bd4-4a24-9fd0-99313af51db9
# ╟─90318732-189c-49e0-b24a-cc499b19f152
# ╟─e02a4195-9b39-4388-84e6-fad12586ac5e
# ╟─f2a9315d-b76d-4017-9bc2-5280f1dc8a97
# ╟─eae50c5b-81bd-41f4-995c-4df44d541c05
# ╟─8dfb9cd4-d16b-49a7-81e1-5ba4085c614a
# ╠═0d3b4e72-8498-4c99-b380-cc269dd74b1c
# ╟─8d5b0ab0-0825-446e-a6a3-132398e2ec04
# ╠═33915b18-1f42-4b0e-8f15-2d35dea8451f
# ╠═1fa67863-e6c0-48a3-a937-a8c42cd49c29
# ╠═672f4668-856f-4da7-977c-a1f2f28bd68e
# ╠═e7e0b571-c282-4db7-8a0c-842b3fe04f29
# ╠═497f8d3c-afcb-42ce-a50f-00b5bc6ecbb1
# ╠═5c8ec9ad-23c0-4cb7-a5cc-fa95d4a1368d
# ╠═cef9fb8f-8824-452a-b335-2e9f31471c92
# ╠═318085c4-2b3d-4970-a65a-71300f6d6743
# ╟─12043a36-7ae2-4f1d-a5a6-21827356335e
# ╠═f75aee30-303c-48f4-a554-9d144ff3611d
# ╟─5606e351-b44c-46aa-b5d9-b7d5b43928da
# ╟─048a13f5-b5e2-46b6-b370-6c388dcf2d4d
# ╟─62899eda-55e0-4ebf-bfdf-86835869e528
# ╟─047609ef-d362-431e-a16c-50570f2eaf92
# ╟─9d4c91d7-2e51-4cb2-9ed6-0f1f67c694c7
# ╟─b97e4725-8b6d-46d5-ad95-f2ae8180079c
# ╠═0f7db3d9-4037-4aea-9314-2d83c7e32f14
# ╠═0d9d96e0-6352-4fb4-9ebb-29ec39ee8ea0
# ╠═52e1ce2b-478d-40c3-8133-841f4de50acb
# ╠═bdb2bf03-a44e-4c1a-b303-8533ef3907af
# ╠═78ded07f-ceac-4ac4-82bd-991da87d8d78
# ╠═fdc42466-d089-4f44-8fd1-0d71de8e9514
# ╠═6b0d2310-8252-4790-8fea-4a542ab0cf78
# ╠═3f7ff72d-971a-492c-8d1c-80359c573c15
# ╟─1a665a69-8122-4547-ba84-223bc4aba2ab
# ╟─ee8900b4-b0a1-4f2d-9f71-53e9ce6edf72
# ╟─932f0f93-1cb5-44ef-97be-cd36f352805d
# ╠═fc7cb2be-9810-4dbc-bce3-52f1aa45d4ce
# ╟─d5f0e41d-9b82-4c5a-af89-3560464895d0
# ╟─46a50de2-2aa3-462e-9fc4-9dd96b3b52bc
# ╟─a6632203-f59c-4793-ad46-4b508214d08d
# ╟─d2e58edd-5759-42f3-920c-0ab640c063fd
# ╟─990efd90-eb3c-4f1d-8f2b-7de99916cbf4
# ╠═d99f7e3c-3de6-4a6d-86d6-2d63755e8e28
# ╟─3a8f030f-4da3-4d3b-8fb2-000dd162823b
# ╟─2918e0cb-8d0a-403d-a04d-679df9660b5c
# ╠═cbe7686f-4d76-4bd0-8983-c0a5633c22d6
# ╠═cee6d373-c5f9-441f-ba98-b7c1cd2ff394
# ╠═3b5c89f8-f251-4016-943b-7cdb4b9cddb5
# ╟─34ed2f38-0376-454f-877f-2e6ac7fb34be
# ╠═20f0f124-055c-44af-ae46-26507b04e3a7
# ╠═525bfa60-4254-4128-a090-ac63da977919
# ╠═7dd65341-de1f-49a0-bd39-d374e287f23a
# ╟─549bfbf4-db61-4887-9334-aa163a0153da
# ╠═6e99ef78-e4cf-4c91-9e26-09a6cf31b165
# ╠═915b6ea6-8ec6-4a34-93d0-5c8d2596980e
# ╠═4bf33fbd-e7db-4dd1-a516-a35657f4055a
# ╠═bf971085-7df9-48e0-a126-f0fd527248c0
# ╠═833e6454-c96c-4d40-838a-8b631f7d171b
# ╠═adb918d6-6614-4749-8f50-222138ac2518
# ╠═dce92cf6-b5f0-40e7-b0b9-79eed0a74ffe
# ╠═717dcd6b-9058-450d-b545-6c7e280d7cd3
# ╠═2cedbac4-9944-4299-a6ec-331f7cbf1e2c
# ╠═bf8883c8-1387-4b4c-aa0a-dfe184d48816
# ╟─3ee5ac17-abd8-4dd7-8ec5-fa94dd93f013
# ╠═212a5853-8072-4614-bf6b-15c3029e26a9
# ╠═00505f0c-12a6-4144-b4a8-46c3d9ead091
# ╟─89328195-8ed9-47ab-912c-16fbb12bf5eb
# ╠═b459c6fe-696a-4af2-80b4-a7c81f7136ea
# ╠═4d7c6062-b61e-44bc-a90d-7fd8f62b178e
# ╟─a35ef1ed-b9b5-4702-aa91-555a28911578
# ╟─842a06d5-72ee-44c6-b38c-f549e439b95c
# ╟─18f4d337-2549-4dc9-97f0-abb3130fef8f
# ╟─953b5d8c-5fbd-4ed0-8b46-21e4bee0e565
# ╟─9688f176-b51e-4dcd-bf21-fbe5af0c366c
# ╟─3d9b0adc-389b-4a58-869c-b1a5c1af22a6
# ╟─a7e59710-2ee0-47fa-a7ba-63054e8cc708
# ╟─1b4a0e49-b64c-4dd7-a03c-94e2611613f1
# ╟─4471dd69-4a70-4f8c-aacf-589d26a53202
# ╟─97d078a4-2005-443b-ac35-3802c85dfa21
# ╟─13b1f81a-cec8-4ba3-b2d7-4beaf3103c39
# ╟─b66d7b10-7426-4f2b-a4bf-b55772244507
# ╟─5bc8263f-184c-4549-8f2b-681bdf697ffd
# ╟─d6e677ab-eb7d-4097-be0c-bd25895d71be
# ╠═254a738a-b1f4-4890-897d-bd810027fcac
# ╟─1b6fca49-4a8a-4a92-a9ab-7aa33fd7f0cb
# ╟─77e13ad3-7e1f-41e7-b1b2-d3b10778bb6f
# ╟─4106e96b-9a53-4a44-a063-97837d57bfc8
# ╟─043c2115-ea95-4b53-877a-6e00ff92f972
# ╟─9a565434-315f-45d0-b64b-d8f21bf76288
# ╟─00a3e093-719d-4205-a439-3d99c98fc52b
