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
end

# ╔═╡ eb93ed6e-a1f2-11eb-00a6-d5bbd921a4be
md" #### Packages"

# ╔═╡ b9b842e4-7861-4bc5-9956-2db1d09c1c7d
md" #### Problem parameters"

# ╔═╡ d49d2169-a568-4319-9287-d13044c30bd4
p=(20,-40.0,-38.4,-10.0,-50,1.0,1.0,0.5,0.015,0.0015,10.0,100.0,1000.0)

# ╔═╡ c3cd5d9c-3192-429e-bcd5-86ac8bf38460
#p=(70,-40.0,-38.4,-19.0,-50,1.0,1.0,0.5,0.1,0.01,10.0,100.0,1000.0)

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

# ╔═╡ 948e5dfa-573f-48de-8899-900405f28bf7
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

# ╔═╡ 566da181-339d-4a3d-abb1-eff5828b4db1
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

# ╔═╡ 90916e6d-7004-4f78-b1ac-4a9f2ea97561
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

# ╔═╡ 774f5bdf-863d-48f7-8115-625a06c171a3
function find_global_bif(p)
	I,v0,vs0,vus0,vuus0,C,gf,gs,gus,guus,ts,tus,tuus = p
	
	numerator = gf*gs*(v0-vs0)^2 +gf*gus*(v0-vus0)^2 +gf*guus*(v0-vuus0)^2 -gs*gus*(vs0-vus0)^2 -gs*guus*(vs0-vuus0)^2 -gus*guus*(vus0-vuus0)^2
	Ibif = numerator /(gf-gs-gus-guus)
end

# ╔═╡ 78d2b9c6-9953-49dd-81a2-e1f22c887d30
Ibif = find_global_bif(change_guus_(0.00015,p))

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

# ╔═╡ fd2dd460-60a1-40de-b05e-2b302b8ad40d
gr()

# ╔═╡ cd5f8df3-1b3a-4343-ab75-172b3fb702de
begin
	ref_time = plot(sol.t,sol[1,:],linecolor="blue",label="V(t)")
	plot!(sol.t,sol[2,:],linecolor="orange",label="Vs(t)")
	plot!(sol.t,sol[3,:],linecolor="pink",label="Vss(t)")
	plot!(sol.t,sol[4,:],linecolor="black",label="Vus(t)")
	#plot!(sol.t,p[3]ones(size(sol.t)),linecolor="red",label="Vs0")
	#plot!(sol.t,p[4]ones(size(sol.t)),linecolor="magenta",label="Vus0")
	#plot!(sol.t,p[5]ones(size(sol.t)),linecolor="green",label="Vuus0")
	#scatter!(sol.t[ind_I_0],v[ind_I_0])
	#scatter!([sol.t[ind1],sol.t[ind2]],[v[ind1],v[ind2]])
	yaxis!("Voltage")
	xaxis!("Time [ms]")
end

# ╔═╡ 9d4d383d-3b08-4257-b59f-5cd7191213d9
begin
	plot(sol.t,sol[3,:],linecolor="pink",label="Vus(t)")
	plot!(sol.t,p[4]ones(size(sol.t)),linecolor="magenta",label="Vus0",size=(300,200))
	plot!(sol.t,sol[4,:],linecolor="orange",label="Vuus(t)")
	plot!(sol.t,p[5]ones(size(sol.t)),linecolor="red",label="Vuus0",size=(300,200))
end

# ╔═╡ 35f9c65b-79a7-4c54-b9f2-436dfed23fde
plotly()

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

# ╔═╡ a9db6e66-9353-44dc-9e68-8ead6251633e
begin
	ref_I = plot(sol.t,Iuus.+p[1].+Ius,label="It",linewidth=1.5)
	#plot!(sol.t,Iuus,label="Iuus")
	#plot!(sol.t,Ius.+p[1],label="I + Ius")
	
	#scatter!(sol.t[ind_I_0],It[ind_I_0])
	yaxis!("Current")
	xaxis!("Time [ms]")
end

# ╔═╡ 3c1c17a8-5e00-4f19-aed1-658a4a791e6a
begin
	plot(ref_time,ref_I,layout=(2,1),size=(800,800),legend=:outertopright)
	#savefig("MQIF_4D_const_current.pdf")
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

# ╔═╡ ff43b14d-da17-4b81-acc6-88ce5d772a90
ind1,ind2 = find_ind_period(sol[1,:],sol[4,:],sol.t)

# ╔═╡ 0b575c04-75e2-4ffd-8b61-8b67103ee35e
plotly()

# ╔═╡ 58f8d821-7868-4ca9-8ee5-eb527ccc90b0
begin
	plot_period = plot(sol.t[ind1:ind2],sol[1,ind1:ind2],linecolor="blue",label="V(t)")
	plot!(sol.t[ind1:ind2],sol[2,ind1:ind2],linecolor="orange",label="Vs(t)")
	plot!(sol.t[ind1:ind2],sol[3,ind1:ind2],linecolor="pink",label="Vss(t)")
	plot!(sol.t[ind1:ind2],sol[4,ind1:ind2],linecolor="black",label="Vus(t)")
	plot!(sol.t[ind1:ind2],p[3]ones(size(sol.t[ind1:ind2])),linecolor="red",label="Vs0",linestyle=:dash)
	plot!(sol.t[ind1:ind2],p[4]ones(size(sol.t[ind1:ind2])),linecolor="magenta",label="Vss0",linestyle=:dash)
	plot!(sol.t[ind1:ind2],p[5]ones(size(sol.t[ind1:ind2])),linecolor="green",label="Vus0",legend=:outertopright,linestyle=:dash)
	xaxis!("Time [ms]")
	yaxis!("Voltage")
end

# ╔═╡ 5e8b79ec-c6fa-4228-a504-15c66c437816
begin
	plot_period_I = plot(sol.t[ind1:ind2],(Iuus.+p[1].+Ius)[ind1:ind2],label="It",legend=:outertopright)
	plot!(sol.t[ind1:ind2],(Ius)[ind1:ind2],label="Iss")
	plot!(sol.t[ind1:ind2],Iuus[ind1:ind2],label="Ius")
	xaxis!("Time [ms]")
	yaxis!("Current")
	
	#scatter!(sol.t[ind1:ind2],ones(size(sol.t[ind1:ind2])))
end

# ╔═╡ b13acccd-827e-48c0-bcc1-21b36939abf5
begin
	plot(plot_period, plot_period_I,layout=(2,1),size=(800,800))
	#savefig("MQIF_4D_const_current_zoom.pdf")
end

# ╔═╡ 0fee90a2-c022-477f-87d1-03c09e5a8f54
plot(sol.t[ind1:ind2],Iuus[ind1:ind2],label="Iuus",size=(400,200))

# ╔═╡ 16b1fac2-35da-4997-bf27-5a1a8ce9a198
plot(sol.t[ind1:ind2],Ius[ind1:ind2],label="Ius",size=(400,200))

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
		
	pp1 =plot(v,vuus,V1_null_,st=:surface,c=:orange,colorbar_entry=false,xlabel="V", ylabel="Vuus",zlabel="Vs",camera=(60,20),size=(500,400)) #yellow-orange
	plot!(v,vuus,V2_null_,st=:surface,c=:orange,colorbar_entry=false,label="dV/dt = 0") #yellow-orange
	plot!([inter_V_Vs1,inter_V_Vs2],[vuus,vuus],[inter_V_Vs1,inter_V_Vs2],linecolor=RGB(1,0.5,0.7),linewidth=3,label="dVs/dt = 0") #red 
	plot!([v,v],[v,v],[inter_V_Vuus1,inter_V_Vuus2],linecolor=RGB(0.31,0.66,1),linewidth=3,label="dVus/dt = 0") #blue
	plot!([vus_level.*ones(size(vuus)),vus_level.*ones(size(vuus))],[vuus,vuus],[inter_V_Vus1,inter_V_Vus2],linecolor=:green,linewidth=3,label="dVss/dt = 0") #blue
	#title!("Vss = $(vus_level)")
	
	pp1_side = plot(pp1,camera=(0,0),legend=:outertopright)
	pp1_oside = plot(pp1,camera=(30,20),legend=:outertopright,colorbar=true)
	#pp1_supp = plot()
	
	#plot(pp1,pp1_side,pp1_oside,pp1_supp,layout=(2,2))
	
	vs_lim1 = maximum(inter_V_Vuus1)
	vs_lim2 = minimum(inter_V_Vuus2)
	return pp1
end

# ╔═╡ 0e22852c-07a5-4bcc-97cf-23fc10b006bd
gr()

# ╔═╡ aaf2d6db-ded3-4aca-8b9a-855c97591fb9
begin 
	v=collect(range(-60,stop=0,length=150))
	vuus=collect(range(-60,stop=30,length=50))
end

# ╔═╡ f94c1f5d-1d58-4ccf-bb35-c2f1570ecd25
@bind vus_level1 html"<input type=range min=-60 max=0 step=10>"

# ╔═╡ e21faa25-c167-4488-8812-7ddd2a46fd08
@bind vus_level2 html"<input type=range min=-60 max=0 step=10>"

# ╔═╡ 0d49fdd3-2224-40c9-b258-b20676798e0a
vus_level1,vus_level2

# ╔═╡ c694d484-2013-465f-8cee-747fd6fcc308
md"""Vus = $vus_level1"""

# ╔═╡ 45fa3d0b-7fef-4261-8842-a0080101f721
begin
	pp1,vslim11,vslim21 = plot_3D_pp(v,vus_level1,vuus,p)
	plot!(pp1,size=(500,400))
end

# ╔═╡ 5c5831f1-96aa-436d-9596-70e07f45092f
md""" Vus = $vus_level2"""

# ╔═╡ cba99b78-9f68-4c65-847a-8bc027d5e43e
begin
	pp2,vslim12,vslim22 = plot_3D_pp(v,vus_level2,vuus,p)
	plot!(pp2,size=(500,400))
end

# ╔═╡ b09cdea7-75a3-4d44-a233-247ee7665975
function plot_potentials(ind,v_p,vs_p,vus_p,vuus_p,t_p)
	plt = plot(t_p[1:ind],vus_p[1:ind],linecolor="pink",label="Vss(t)")
	xaxis!("Time [ms]")
	yaxis!("Voltage")
	return plt
end

# ╔═╡ 098b56bb-0cc5-4623-ba20-66267734b2d0
function make_gif_pp_3D(v,v_p,vs_p,vus_p,vuus_p,vuus,t_p,p)
	plt_3D_pp= []
	plt_pot=[]
	ind=[]
	t = t_p .-t_p[1]
	t_prev=[0]
	for i=1:length(vus_p)
		if t[i]>t_prev[1]+20
			push!(plt_3D_pp,plot_3D_pp_gif(v,vus_p[i],vuus,p))
			t_prev=[t[i]]
			push!(ind,i)
			push!(plt_pot,plot_potentials(i,v_p,vs_p,vus_p,vuus_p,t_p))
		end
	end
	return plt_3D_pp,plt_pot
end

# ╔═╡ 1a1c3de2-8d82-404f-804c-65e43a90db33
bob=[]

# ╔═╡ d670a124-b520-462b-82b0-bf720a2e4d05
begin
	push!(bob,1)
end

# ╔═╡ d1a284ab-9206-46f5-a16e-3b8831a157a9
begin
	plot_3D_pp_gif(v,-40,vuus,p)
end

# ╔═╡ f6e3cb8e-868a-4206-b546-2d69d2ef7b16
md""" Vus = $vus_level1"""

# ╔═╡ 835938d0-75b8-4e5c-85c3-daa4c4748acf
plot(pp1,camera=(90,0),legend=:outertopright)

# ╔═╡ 060f260c-7506-4e99-be95-98a9f1d1568b
md""" Vus = $vus_level2"""

# ╔═╡ c8ce3b6d-607d-44e6-937c-06c848fe4967
begin
	plot(pp2,camera=(90,0),legend=:outertopright)
	md"""view in Vuus;Vs from 3D"""
end

# ╔═╡ 7b9ce5de-ca70-4dd6-bbac-1e087b58b7c5
begin
	vs_max1 = zeros(length(vuus))
	vs_min1 = zeros(length(vuus))
	vs_max2 = zeros(length(vuus))
	vs_min2 = zeros(length(vuus))
	for i=1:length(vuus)
		vs_max1[i],vs_min1[i] = compute_vs_extremas(vus_level1,vuus[i],p)
		vs_max2[i],vs_min2[i] = compute_vs_extremas(vus_level2,vuus[i],p)
	end
	plot(vuus,vs_max1,linecolor=:red,label = "Extrema ordinate, vus = $vus_level1")
	plot!(vuus,vs_min1,linecolor=:red,label="Extrema ordinate, vus = $vus_level1")
	plot!(vuus,vs_max1,fill=(vs_min1, 0.3, :orange),linealpha=0,label="Excitability well, vus = $vus_level1")
	plot!(vuus,vs_max2,linecolor=RGB(1,0.3,0),label = "Extrema ordinate, vus = $vus_level2")
	plot!(vuus,vs_min2,linecolor=RGB(1,0.3,0),label="Extrema ordinate, vus = $vus_level2")
	plot!(vuus,vs_max2,fill=(vs_min2, 0.3, :yellow),linealpha=0,label="Excitability well, vus = $vus_level2",size=(500,400))
	yaxis!("Vs",(vslim21,vslim11))
	xaxis!("Vuus")
end

# ╔═╡ 3c9cbb99-1fd5-4dd2-a77b-711080fecf43
plot([Shape([1,2,3], [1,1,i]) for i=1:5], fill_z=[1,2,3,4,5],size=(200,100))

# ╔═╡ 17856e2b-0541-42fb-b393-9858be2411be
md" #### Excitability well Plane"

# ╔═╡ 35a023dc-16c4-4054-903a-1659474ac2a5
md"""Study of one period of the reference simulation"""

# ╔═╡ 904855d1-66fc-42f9-95e1-d8c6b43d8b9d
begin 
	v_p = sol[1,ind1:ind2]
	vs_p = sol[2,ind1:ind2]
	vus_p = sol[3,ind1:ind2]
	vuus_p = sol[4,ind1:ind2]
	t_p = sol.t[ind1:ind2]
	It_f = Iuus.+p[1].+Ius
	It_p = It_f[ind1:ind2]
end

# ╔═╡ 701ee8d1-002f-40ba-9526-a368faf81fe8
begin
	md"""Gif of 3D pp"""
	#plt_3D_pp = [plot_3D_pp_gif(v,vus_p[i],vuus,p) for i=1:length(vus_p)]
	plt_3D_pp = make_gif_pp_3D(v,v_p,vs_p,vus_p,vuus_p,vuus,t_p,p)
end

# ╔═╡ bb5563b2-bdb8-4e3e-9c89-86ba126587bc
length(plt_3D_pp[1])

# ╔═╡ 0056a35e-1ec4-49c8-b2b0-98aea93c2bae
begin
	for i=1:length(plt_3D_pp[1])
		plot(plt_3D_pp[1][i],plt_3D_pp[2][i],layout=@layout([a{0.7w} b]),size=(1000,850),legend=:outertopright)
		savefig("gif_4D_pp_$i.png")
	end
end

# ╔═╡ df1a06db-dbac-4ca3-b302-25c0dce920d6
length(plt_3D_pp)

# ╔═╡ f2022bbf-d45b-43f8-b076-b7c9babe29e2
t = t_p .-t_p[1]

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
	plt = plot(vuus_,vs_max,linecolor=:red,label = "Extrema ordinate, vss = $(Int(round(vus_level)))")
	plot!(vuus_,vs_min,linecolor=:red,label="Extrema ordinate, vss = $(Int(round(vus_level)))")
	plot!(collect(range(minimum(vuus),stop=maximum(vuus),length=length(vuus_))),vs_max,fill=(vs_min, 0.5, :orange),linealpha=0,label="Excitability well, vss = $(Int(round(vus_level)))",size=(500,400))
	plot!(vuus[1:i],vs[1:i],line_z=t,c=:cyclic_mygbm_30_95_c78_n256,label="Trajectory")	
	yaxis!("Vs",(-80,0))
	xaxis!("Vus")
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
	plot!(t[1:i],vus[1:i],linecolor="pink",label="Vss(t)")
	plot!(t[1:i],vuus[1:i],linecolor="black",label="Vus(t)")
	plot!(t[1:i],p[3]ones(size(t[1:i])),linecolor="red",linestyle=:dash,label="Vs0")
	plot!(t[1:i],p[4]ones(size(t[1:i])),linecolor="magenta",linestyle=:dash,label="Vss0")
	plot!(t[1:i],p[5]ones(size(t[1:i])),linecolor="green",label="Vus0",linestyle=:dash,legend=:outertopright)
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
	
	p2D = p_4D_to_2D(p4D)
	
	if i<length(sample_index)
		plot_test = plot(t[1:sample_index[i+1]],v[1:sample_index[i+1]],linecolor="blue",label="V(t)")
	plot!(t[1:sample_index[i+1]],vs[1:sample_index[i+1]],linecolor="orange",label="Vs(t)")
	plot!(t[1:sample_index[i+1]],vus[1:sample_index[i+1]],linecolor="pink",label="Vss(t)")
	plot!(t[1:sample_index[i+1]],vuus[1:sample_index[i+1]],linecolor="black",label="Vus(t)")
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
		plot!(t_,vus_,linecolor="pink",label="Vss(t)")
		plot!(t_,vuus_,linecolor="black",label="Vus(t)")
	end
	
	xaxis!("Time")
	yaxis!("Voltage",(-80,0))
	title!("  B. ",titlelocation=:left,titlefontsize=18)
	return plot_test
end

# ╔═╡ 9bba65a2-b723-4f85-9301-673637478f51
gr()

# ╔═╡ e5a9ae32-8fa4-49bb-b3da-3cbbe31ef3c1
begin
	recon_It,sample_index = sample_It(30,100,70,t_p,v_p,vs_p,vus_p,vuus_p,It_p)
	plot(t_p,It_p)
	plot!(t_p,recon_It)
	scatter!([t_p[last_max_It(v_p,vs_p,vus_p,vuus_p,p)[1]]],[It_p[last_max_It(v_p,vs_p,vus_p,vuus_p,p)[1]]])
end

# ╔═╡ f3355244-baf2-4203-9073-7c524eed3932
begin
	v_potential = collect(range(-80,stop=0,length=150))
	vs_null_2D= Vsnullcline_2D(v_potential,p_4D_to_2D(p))
	v_null_1_list,v_null_2_list,traj_v_list,traj_vs_list,fp1_list,fp2_list, stab_list=give_me_phase_planes_l(v_potential,p,recon_It,v_p,vs_p,sample_index) 
	plot_list = [make_my_phase_plane(v_potential,v_null_1_list,v_null_2_list,vs_null_2D,traj_v_list,traj_vs_list,fp1_list,fp2_list,stab_list,recon_It,t_p,sample_index,p,i) for i=1:length(traj_v_list)]
	plot_t = [make_my_time_plot(v_potential,v_null_1_list,v_null_2_list,vs_null_2D,v_p,vs_p,vus_p,vuus_p,fp1_list,fp2_list,stab_list,recon_It,t_p,sample_index,p,i) for i=1:length(traj_v_list)]
	

	#savefig("simu2D_for_spikelat_3D_phase.pdf")
end

# ╔═╡ 3d98e1c7-2f93-46d1-822d-695035bcc1f8
begin
	plt_ = [gif_well(v_p,vs_p,vus_p,vuus_p,t_p,p,sample_index[i]) for i=1:length(sample_index)]
	plt_t = [gif_time(v_p,vs_p,vus_p,vuus_p,t_p,p,sample_index[i]) for i=1:length(sample_index)]
	plt_I = [gif_It_time(It_p,t_p,p,sample_index[i]) for i=1:length(sample_index)]
end

# ╔═╡ f63f227d-6370-42f2-8585-34ef60cadd51
plt_[100]

# ╔═╡ 6d664df4-946b-4a8d-a5dd-acf3f15719da
begin
	for i=1:length(plot_list)
		z=i
		pp_t = plot(plot_list[z],plot_t[z],layout=(1,2),size=(800,500))
		well_I = plot(plt_[z],plt_I[z],layout=(1,2),size=(800,500))
		plot(pp_t,well_I,layout=(2,1),size=(1000,750))
		#savefig("gif_4D_$z.png")
	end
end

# ╔═╡ 7c8232d3-535d-4909-ac65-3036668cd0cd
bob[200]

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
	f_intra[ind_t_period_s[end]:end].=f_inter_ref
	
end

# ╔═╡ b8fc2a68-9401-4915-880e-4efa4f5c2a3f
f_intra_full

# ╔═╡ f508c845-1dc1-4999-9d79-952a9c577301
begin
	plt_f_intra = plot(t_period,f_intra,label="Burst instant frequency",size=(400,200))
	xaxis!("Time (ms)")
	yaxis!("Frequency [Hz]")
	
	ind_end_burst = findfirst(x -> x>=t_burst+t_p[1],t_p)
	plt_f_intra_It = plot(t_p[1:ind_end_burst].-t_p[1],It_p[1:ind_end_burst],linecolor="green",label="It")
	xaxis!("Time (ms)")
	yaxis!("Current")
	
	plot(plt_f_intra,plt_f_intra_It,layout=(2,1))
	#savefig("MQIF_4D_intraburst_f_.pdf")
end

# ╔═╡ 86acc7c8-971e-4867-8c1d-b10bd7a11d7f
gr()

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

# ╔═╡ 9d3630b9-15a0-4761-bfba-04cb4e90ad5c
begin
	plot(I_global_bif,f_intra_av,label="Intraburst average frequency")
	plot!(I_global_bif,f_intra_max,fill = (f_intra_min, 0.5, :pink),linecolor=:purple,linealpha=0.1,label="Min max interburst frequency",size=(450,350),legend=:topleft)
	#plot!(I_global_bif,f_intra_av_r,label="Intraburst average frequency")
	#plot!(I_global_bif,f_intra_max_r,fill = (f_intra_min_r, 0.5, :orange),linecolor=:orange,linealpha=0.1,label="Min max interburst frequency",size=(450,350),legend=:bottomright)
	yaxis!("Frequency [Hz]")
	xaxis!("Current",(10,50))
	#savefig("MQIF_4D_intraburst_I_f.pdf")
end

# ╔═╡ 42eaf693-f686-4dd9-b1ad-70763df7f35c
begin 
	plot_fr_inter = plot(I_global_bif,f_inter,linealpha=0.5,label="Interburst frequency     ",size=(450,350))
	yaxis!("Frequency [Hz]",(0,5))
	plot_n_fr = scatter(I_global_bif,n_spike,markersize=1,size=(450,350),label="Intraburst spike number")
	yaxis!("Number of spike",(20,40))
	plot(plot_fr_inter,plot_n_fr,layout=(2,1),legend=:outertopright)
	xaxis!("Current",(10,50))
	#savefig("MQIF_4D_interburst_I_f_.pdf")
end

# ╔═╡ b6f58bd2-9695-440e-860d-f7a778da90cd
f_inter

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

# ╔═╡ 09a54577-14de-4154-905d-d82a2e53f5f3
begin
			ind_spikeV_test = find_Vspikes(sol_test.t,sol_test[1,:],maximum(tspan)) 
			ind_spikeVuus_test = find_Vuusspikes(sol_test.t,sol_test[4,:],maximum(tspan))
			
		delta_t_Vuus_test = sol_test.t[ind_spikeVuus_test[2:end]] - sol_test.t[ind_spikeVuus_test[1:end-1]]
		delta_i_Vuus_test = ind_spikeVuus_test[end]-ind_spikeVuus_test[end-1]
			t_period_inter_test= mean(delta_t_Vuus_test)
			
			
			delta_t_V_test = sol_test.t[ind_spikeV_test[2:end]] - sol_test.t[ind_spikeV_test[1:end-1]]
			i_p_n_1_test =findfirst(x-> x>=ind_spikeVuus_test[end-1],ind_spikeV_test)
			i_p_n_2_test =length(ind_spikeV_test) - (findfirst(x-> x<=ind_spikeVuus_test[end],reverse(ind_spikeV_test))-1)
			
			n_spike_test= i_p_n_2_test-i_p_n_1_test+1

			if n_spike_test >1
				t_period_intra1_test=maximum(delta_t_V_test[i_p_n_1_test:i_p_n_2_test-1])
				t_period_intra2_test=minimum(delta_t_V_test[i_p_n_1_test:i_p_n_2_test-1])
				
				ar = delta_t_V_test[i_p_n_1_test:i_p_n_2_test-1]
				t_period_intra_av_test = mean(ar) 
			else
				t_period_intra1_test=NaN
				t_period_intra2_test=NaN
				t_period_intra_av_test = NaN
			end
	plot(sol_test.t,sol_test[1,:])
	plot!(sol_test.t,sol_test[4,:])
	scatter!(sol_test.t[ind_spikeVuus_test],sol_test[4,ind_spikeVuus_test])
	scatter!(sol_test.t[ind_spikeV_test],sol_test[1,ind_spikeV_test])
	scatter!([sol_test.t[ind_spikeV_test[i_p_n_2_test]]],[sol_test[1,ind_spikeV_test[i_p_n_2_test]]])
	scatter!([sol_test.t[ind_spikeV_test[i_p_n_1_test]]],[sol_test[1,ind_spikeV_test[i_p_n_1_test]]],size=(200,200))

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

# ╔═╡ 86ed9930-0f8d-4fa3-97c2-493f61dc41d3
begin
	plot(list_I_stable_global_bif,list_V_stable_global_bif,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(list_I_saddle_global_bif,list_V_saddle_global_bif,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	#plot!(list_I_unstable_global_bif,list_V_unstable_global_bif,linecolor="red",linewidth=1.5,label="Unstable node")
	scatter!([Ibif],[Vbif],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	#plot!(I_global_bif_cycle,V_cycle_global_bif[:,1])
	plot!(I_global_bif,V_cycle_global_bif[:,2],label="Limit cycle minimum")
	#yaxis!((-50,-30))
	xaxis!("I")
	yaxis!("V")
	#savefig("MQIF_4D_bifI_.pdf")
end

# ╔═╡ 74e4e531-3f7b-4d89-b9d4-14700c12fae2
md"""Todo : 

-study the impact of Vus0, gus, 

-update the results pdf"""

# ╔═╡ 4471dd69-4a70-4f8c-aacf-589d26a53202
md"#### Pattern simulations" 

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
	sol0 = solve(prob0,DP5(),reltol=1e-6,abstol=1e-6)
	
	tspan1=(t0_step,tf_step)
	u01 =[sol0[1,end],sol0[2,end],sol0[3,end],sol0[4,end]]
	p1=change_I_(Istep,p)
	
	prob1 = ODEProblem(MQIF_4D!,u01,tspan1,p1,callback=cb)
	sol1 = solve(prob1,DP5(),reltol=1e-6,abstol=1e-6)
	
	tspan2=(tf_step,tf)
	u02 =[sol1[1,end],sol1[2,end],sol1[3,end],sol1[4,end]]
	p2=change_I_(I0,p)
	
	prob2 = ODEProblem(MQIF_4D!,u02,tspan2,p2,callback=cb)
	sol2 = solve(prob2,DP5(),reltol=1e-6,abstol=1e-6)

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

# ╔═╡ 14ae9ded-d965-44c9-bdab-7053feaee9ad
begin
	sub =subplot_simu_step(0.0,10000.0,1200.0,p,0.0,20.0,true,false)	
	#savefig("MQIF_4D-simustep.pdf")
end

# ╔═╡ 0ee33d9f-2726-467b-b2b8-7daacfe9b9e1
begin
	subs = subplot_simu_spike(0.0,10000.0,1500.0,p,1.0,20.0,100,true,false)
	#savefig("MQIF_4D-simupulse.pdf")
end

# ╔═╡ 2192bcbc-4cba-42b3-a2b5-15918a33c512
begin
	t_u=[250,2000,7000]
	delta_ud = [600,1000,2500]
	t_d=t_u+delta_ud
	sub_ud =subplot_simu_ud(0.0,10000.0,t_u,t_d,p,0.4,20.0,100.0,true,false)
	#savefig("MQIF_4D-simuud.pdf")
end

# ╔═╡ 219d865e-0f6d-4c44-bcfa-ce36c2ac64b7
md"#### Dynamic input conductances"

# ╔═╡ 3bb8bfd4-cd13-471f-a80c-9beb483ca5aa
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
	if norm==true
		n_dic_gf = zeros(size(dic_gf))
		n_dic_gs = zeros(size(dic_gs))
		n_dic_gus = zeros(size(dic_gus))
		n_dic_guus = zeros(size(dic_guus))
		
		max_dic_gf=maximum(dic_gf)
		min_dic_gf=minimum(dic_gf)
		n_dic_gf .= (dic_gf.-min_dic_gf)./((max_dic_gf.-min_dic_gf)./2).-1
		
		max_dic_gs=maximum(dic_gs)
		min_dic_gs=minimum(dic_gs)
		n_dic_gs .= (dic_gs.-min_dic_gs)./((max_dic_gs.-min_dic_gs)./2).-1
		
		max_dic_gus=maximum(dic_gus)
		min_dic_gus=minimum(dic_gus)
		n_dic_gus .= (dic_gus.-min_dic_gus)./((max_dic_gus.-min_dic_gus)./2).-1
	end
	dic_plot = plot(v,dic_gf,lincolor="blue",label="gf")
	plot!(v,dic_gs,lincolor="red",label="gs")
	plot!(v,dic_gus,lincolor="green",label="gus")
	plot!(v,dic_guus,lincolor="purple",label="guus")
	scatter!([v0],[jacobian(v0,p)[1,1]],markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = "blue",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
	scatter!([vs0],[jacobian(vs0,p)[1,2]],markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = "red",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
	scatter!([vus0],[jacobian(vus0,p)[1,3]],markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = "green",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
	scatter!([vuus0],[jacobian(vuus0,p)[1,4]],markershape = :circle,markersize = 2,markeralpha = 0.6,markercolor = "purple",markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
	xaxis!("V")
	yaxis!("DIC",(-50,50))
	
	n_plot = plot(v,n_dic_gf,lincolor="blue",label="gf")
	plot!(v,n_dic_gs,lincolor="red",label="gs")
	plot!(v,n_dic_gus,lincolor="green",label="gus")
	plot!(v,n_dic_guus,lincolor="purple",label="guus")
	
	return dic_plot
end

# ╔═╡ f5562329-04fa-4a9e-b20c-f8fbe23f59ec
dic(p,true)

# ╔═╡ 35cc5007-235e-4671-b29c-c6016ca25259
begin
	p_tonic=change_gus_(0.0015,change_guus_(0.0028,p))
	dic(p_tonic,true)
end

# ╔═╡ 884e475e-b8a4-43b3-a596-803fbd0624e6
p_tonic

# ╔═╡ Cell order:
# ╠═eb93ed6e-a1f2-11eb-00a6-d5bbd921a4be
# ╠═eb47cf0d-f73a-4719-afda-3188b2998e52
# ╟─b9b842e4-7861-4bc5-9956-2db1d09c1c7d
# ╠═d49d2169-a568-4319-9287-d13044c30bd4
# ╠═c3cd5d9c-3192-429e-bcd5-86ac8bf38460
# ╠═78d2b9c6-9953-49dd-81a2-e1f22c887d30
# ╟─97629afb-5289-422a-bb17-5c1a0afebc00
# ╟─94979c24-7cc6-4c67-a0fa-13e0ea2c51f7
# ╟─3ae16e48-9afb-4d31-9123-cf6edb092b51
# ╟─9a5b6a2b-4080-4128-aabb-c1846eab7a40
# ╟─c7930bda-a02c-4ec1-a975-83a61de2592a
# ╟─042c10bf-0755-48f0-be92-f7451765a816
# ╟─948e5dfa-573f-48de-8899-900405f28bf7
# ╟─566da181-339d-4a3d-abb1-eff5828b4db1
# ╟─90916e6d-7004-4f78-b1ac-4a9f2ea97561
# ╟─774f5bdf-863d-48f7-8115-625a06c171a3
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
# ╠═7530943d-039a-4467-9c7d-8f565c5e5298
# ╟─0d58c293-3a2b-4518-ac0f-f26620ca4e44
# ╟─36f2a546-2716-494c-ad21-09c21f78cb8a
# ╟─411e46b8-8ebc-4ef7-bb0c-688020416e5b
# ╟─44d2c44b-44a3-490d-a933-12f74d0ef5b1
# ╠═47bccea1-4cb9-40ac-a5dc-62322512ce1a
# ╠═ae088c53-ab45-4dd8-8e01-384940acb491
# ╟─be306caa-d310-4a27-bfcf-16b5cff5590b
# ╠═006de583-9601-432c-95a7-a3f786b0c937
# ╠═fd2dd460-60a1-40de-b05e-2b302b8ad40d
# ╠═cd5f8df3-1b3a-4343-ab75-172b3fb702de
# ╠═9d4d383d-3b08-4257-b59f-5cd7191213d9
# ╠═35f9c65b-79a7-4c54-b9f2-436dfed23fde
# ╠═a9db6e66-9353-44dc-9e68-8ead6251633e
# ╠═3c1c17a8-5e00-4f19-aed1-658a4a791e6a
# ╟─61645246-ac2d-4d11-9252-ab4929841135
# ╠═20c01555-448a-48b8-ad43-345ad4c82c34
# ╠═a7234266-83ee-4f54-b3a8-3d2ada1c2a00
# ╠═ff43b14d-da17-4b81-acc6-88ce5d772a90
# ╠═0b575c04-75e2-4ffd-8b61-8b67103ee35e
# ╠═58f8d821-7868-4ca9-8ee5-eb527ccc90b0
# ╠═5e8b79ec-c6fa-4228-a504-15c66c437816
# ╠═b13acccd-827e-48c0-bcc1-21b36939abf5
# ╟─0fee90a2-c022-477f-87d1-03c09e5a8f54
# ╟─16b1fac2-35da-4997-bf27-5a1a8ce9a198
# ╟─7ab02730-0c95-445f-a1f3-f3d645776f9b
# ╟─3e13d7ce-1056-4562-86ad-412598a9eadf
# ╠═0ef6fa1a-efe7-441c-8360-9afba715e066
# ╠═db5df20a-fea3-4f51-827f-5ff657bd21a4
# ╠═0e22852c-07a5-4bcc-97cf-23fc10b006bd
# ╠═aaf2d6db-ded3-4aca-8b9a-855c97591fb9
# ╟─f94c1f5d-1d58-4ccf-bb35-c2f1570ecd25
# ╟─e21faa25-c167-4488-8812-7ddd2a46fd08
# ╟─0d49fdd3-2224-40c9-b258-b20676798e0a
# ╟─c694d484-2013-465f-8cee-747fd6fcc308
# ╟─45fa3d0b-7fef-4261-8842-a0080101f721
# ╟─5c5831f1-96aa-436d-9596-70e07f45092f
# ╟─cba99b78-9f68-4c65-847a-8bc027d5e43e
# ╠═701ee8d1-002f-40ba-9526-a368faf81fe8
# ╠═098b56bb-0cc5-4623-ba20-66267734b2d0
# ╠═b09cdea7-75a3-4d44-a233-247ee7665975
# ╠═f2022bbf-d45b-43f8-b076-b7c9babe29e2
# ╠═1a1c3de2-8d82-404f-804c-65e43a90db33
# ╠═d670a124-b520-462b-82b0-bf720a2e4d05
# ╠═d1a284ab-9206-46f5-a16e-3b8831a157a9
# ╠═bb5563b2-bdb8-4e3e-9c89-86ba126587bc
# ╠═0056a35e-1ec4-49c8-b2b0-98aea93c2bae
# ╠═df1a06db-dbac-4ca3-b302-25c0dce920d6
# ╟─f6e3cb8e-868a-4206-b546-2d69d2ef7b16
# ╟─835938d0-75b8-4e5c-85c3-daa4c4748acf
# ╟─060f260c-7506-4e99-be95-98a9f1d1568b
# ╟─c8ce3b6d-607d-44e6-937c-06c848fe4967
# ╟─7b9ce5de-ca70-4dd6-bbac-1e087b58b7c5
# ╟─3c9cbb99-1fd5-4dd2-a77b-711080fecf43
# ╟─17856e2b-0541-42fb-b393-9858be2411be
# ╟─35a023dc-16c4-4054-903a-1659474ac2a5
# ╠═904855d1-66fc-42f9-95e1-d8c6b43d8b9d
# ╠═e784cdff-fab4-4344-883f-6fcc17af8ad0
# ╠═80acbb79-335f-4f40-90a5-fe3a03c937aa
# ╟─14e0a258-d957-45cc-8543-9203fea533be
# ╠═67122038-ec98-4c9f-9d18-51f02e28a183
# ╠═33117e03-5d89-47fa-90a5-e89dc2071512
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
# ╠═04ed6d51-354a-4432-904d-76a72d7fa322
# ╠═9bba65a2-b723-4f85-9301-673637478f51
# ╠═e5a9ae32-8fa4-49bb-b3da-3cbbe31ef3c1
# ╠═f3355244-baf2-4203-9073-7c524eed3932
# ╠═f63f227d-6370-42f2-8585-34ef60cadd51
# ╠═3d98e1c7-2f93-46d1-822d-695035bcc1f8
# ╠═6d664df4-946b-4a8d-a5dd-acf3f15719da
# ╠═7c8232d3-535d-4909-ac65-3036668cd0cd
# ╟─57202474-8e51-43b3-8a71-17cc786aead6
# ╟─61bff036-322e-4e53-9114-c99c6457e115
# ╟─22a4c6d3-ab80-4f99-9b04-8f467edaa965
# ╟─6146153e-49d6-4b51-9b03-00c56c5c08c8
# ╟─bcd46d0c-4738-4885-b6d4-ffe13143fce1
# ╟─26dd3baf-c6ed-47c1-8337-69a0589a237a
# ╠═b77a655e-4d54-4e3c-8bdb-a557141aa774
# ╠═b8fc2a68-9401-4915-880e-4efa4f5c2a3f
# ╠═f508c845-1dc1-4999-9d79-952a9c577301
# ╠═86acc7c8-971e-4867-8c1d-b10bd7a11d7f
# ╟─f1789d85-45bc-4544-a175-197225e8ed6f
# ╠═9d3630b9-15a0-4761-bfba-04cb4e90ad5c
# ╠═42eaf693-f686-4dd9-b1ad-70763df7f35c
# ╠═b6f58bd2-9695-440e-860d-f7a778da90cd
# ╟─481602fa-f9e9-450f-a831-f0caf68eca54
# ╟─09a54577-14de-4154-905d-d82a2e53f5f3
# ╟─621533e4-1c67-4436-baf8-8d55598c1184
# ╠═4c197fd0-9ae6-447f-90d4-4e6f395613fd
# ╟─7307d524-5bd4-4a24-9fd0-99313af51db9
# ╟─90318732-189c-49e0-b24a-cc499b19f152
# ╟─e02a4195-9b39-4388-84e6-fad12586ac5e
# ╠═86ed9930-0f8d-4fa3-97c2-493f61dc41d3
# ╟─74e4e531-3f7b-4d89-b9d4-14700c12fae2
# ╟─4471dd69-4a70-4f8c-aacf-589d26a53202
# ╟─13b1f81a-cec8-4ba3-b2d7-4beaf3103c39
# ╟─77e13ad3-7e1f-41e7-b1b2-d3b10778bb6f
# ╟─4106e96b-9a53-4a44-a063-97837d57bfc8
# ╟─043c2115-ea95-4b53-877a-6e00ff92f972
# ╟─9a565434-315f-45d0-b64b-d8f21bf76288
# ╟─00a3e093-719d-4205-a439-3d99c98fc52b
# ╠═14ae9ded-d965-44c9-bdab-7053feaee9ad
# ╠═0ee33d9f-2726-467b-b2b8-7daacfe9b9e1
# ╠═2192bcbc-4cba-42b3-a2b5-15918a33c512
# ╟─219d865e-0f6d-4c44-bcfa-ce36c2ac64b7
# ╟─3bb8bfd4-cd13-471f-a80c-9beb483ca5aa
# ╠═f5562329-04fa-4a9e-b20c-f8fbe23f59ec
# ╠═35cc5007-235e-4671-b29c-c6016ca25259
# ╠═884e475e-b8a4-43b3-a596-803fbd0624e6
