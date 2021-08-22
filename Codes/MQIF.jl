### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ b9abc652-6491-11eb-0e92-7b7877bf9faf
import Pkg

# ╔═╡ 50db35c0-6af4-11eb-3b81-b50d69d7a085
using LinearAlgebra

# ╔═╡ d47d0250-6491-11eb-0850-3f015c450720
using DifferentialEquations

# ╔═╡ 9b413122-8e79-11eb-39d0-739f756ee127
using Statistics

# ╔═╡ 5af7bdb0-64ac-11eb-3b4e-f7d3efeead19
using Plots 

# ╔═╡ f0a97dfe-6490-11eb-0684-c50efecf2bfd
md"# MQIF model"

# ╔═╡ 6a3dc280-6491-11eb-0ce0-6b468892be37
#Pkg.add("DifferentialEquations")

# ╔═╡ 6bb6b2a0-64ac-11eb-2b8e-b71abbdf3043
#Pkg.add("Plots")

# ╔═╡ b654ccb0-654d-11eb-2b8c-1700940eb04b
md"#### Nullclines "

# ╔═╡ 6bbef720-6548-11eb-220d-094c9f1f4504
function Vnullcline(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	
	vs1 = vs0 + sqrt((gf*(v-v0)^2 + I)/gs)
	vs2 = vs0 - sqrt((gf*(v-v0)^2 + I)/gs)
	return vs1,vs2
end

# ╔═╡ a447df80-6a25-11eb-2477-aba13d9a9d7e
function Vnullcline_(vs,p)
	I,v0,vs0,C,gf,gs,ts = p
	
	v1 = v0 + sqrt((gs*(vs-vs0)^2 - I)/gf)
	v2 = v0 - sqrt((gs*(vs-vs0)^2 - I)/gf)
	return v1,v2
end

# ╔═╡ cca2ee5e-6549-11eb-3b9e-0934cddd0aa8
function Vsnullcline(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	vs = v
	return vs
end

# ╔═╡ 779bca42-6bc5-11eb-14b5-b38a2a99db3e
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

# ╔═╡ d0569450-8b58-11eb-16b6-8762d560dafb
function change_v0(v0,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==2
			p2[i]=v0
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ c6476770-7f58-11eb-2fa2-5b79fa30d78e
function change_gf(gf,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==5
			p2[i]=gf
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 00178390-7f59-11eb-1c6c-75a2bcbedfc4
function change_gs(gs,p)
	p2=zeros(length(p))
	for i=1:length(p)
		if i==6
			p2[i]=gs
		else
			p2[i]=p[i]
		end
	end
	return p2
end

# ╔═╡ 3dd5be70-7f88-11eb-01b7-a12dc52591ae
function fixedpoints_only(p)
	I,v0,vs0,C,gf,gs,ts = p
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
	end
	
	fp1=ones(1,2).*vfp1
	fp2=ones(1,2).*vfp2	
	
	return fp1,fp2
end

# ╔═╡ 1e623420-6afb-11eb-0e56-45472f98e9fe
function jacobian(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	J = zeros(2,2)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[2,1] = 1
	J[2,2] = -1
	
	return J
end

# ╔═╡ b75c43ce-6aee-11eb-385c-3f46179d7089
function fixedpoints(p)
	I,v0,vs0,C,gf,gs,ts = p
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
		for i=1:length(fp)
			J = jacobian(fp[i],p)
			lambda1,lambda2 = eigvals(J)
			if real(lambda1)>0 && real(lambda2)>0
				stability[i]=1.0
			end
			if (real(lambda1)>0 && real(lambda2)<0) || (real(lambda1)<0 && 	
				real(lambda2)>0)
				
				stability[i]=0.0
			end
			if real(lambda1)<0 && real(lambda2)<0
				stability[i]=-1.0
			end
		end
	end
	
	fp1=ones(1,2).*vfp1
	fp2=ones(1,2).*vfp2	
	
	return fp1,fp2,stability
end

# ╔═╡ 8e67dd20-7225-11eb-1489-6f6dbe8cbd6c
function find_bifurcation(p)
	I,v0,vs0,C,gf,gs,ts = p
	vs0_bif_l = v0 - sqrt(I*(gf-gs)/(gf*gs))
	vs0_bif_h = v0 + sqrt(I*(gf-gs)/(gf*gs))
	
	return vs0_bif_l,vs0_bif_h
end

# ╔═╡ 200e05e0-6549-11eb-2151-51f90adb859f
cell_potential = collect(range(-70,stop=30,length=100))

# ╔═╡ 9f24b0c0-6654-11eb-2511-09afcd5a87e5
md"####  Main parameters"

# ╔═╡ 7af74d3e-6562-11eb-16a4-5195708272f3
begin 
	Vmax = -10 #30
	Vr = -30 #-40
	Vsr = -10 #-20	
	
	md"""Voltages used for reset"""
end

# ╔═╡ 137cf390-6643-11eb-1eb2-9f87c024400f
begin 
	Il = 0.352
	Ih = 30
	
	#for phase plane section 
	Ip = 5.0
	In = -5.0
	Inh = -30.0
	md""" Currents used for phase planes with varying current"""
end

# ╔═╡ 9bad8b60-6654-11eb-1c88-eb5853740148
#p=(Il,-40.0,-35.0,1.0,1.0,0.5,10.0)

# ╔═╡ 28be7ed2-8e21-11eb-15c3-715f4d558d6c
p=(0.0,-40.0,-35.0,1.0,1.0,0.5,10.0)

# ╔═╡ 4fa199b0-6bb7-11eb-32a6-0ff14be97216
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

# ╔═╡ 35ce42a0-6549-11eb-0201-d1a1b2f5ede1
begin 
	V1 = zeros(size(cell_potential))
	V2 = zeros(size(cell_potential))
	V3 = zeros(size(cell_potential))
	
	for i=1:length(cell_potential)
		V1[i],V2[i]=Vnullcline(cell_potential[i],p)
		V3[i]=Vsnullcline(cell_potential[i],p)
	end
	md"""Low current nullclines"""
end

# ╔═╡ c7d4eec0-654d-11eb-0d20-c5657cc52a2e
md" #### ODE problem functions"

# ╔═╡ 42327b20-64a8-11eb-1c0b-47200e371b58
function MQIF!(du,u,p,t)	
	I,v0,vs0,C,gf,gs,ts = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
end

# ╔═╡ 4f8838e0-6562-11eb-1d9b-b9cdffcece31
function spike(x)  # spikes when spike(x) goes from negative to positive
    (x[1] - Vmax)
end

# ╔═╡ 98fcb230-6562-11eb-1194-8d99fc6dad03
function reset!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
end

# ╔═╡ e2dac2be-6562-11eb-0607-a3ec86eccd70
# event when event_f(u,t) == 0
function condition(x,t,integrator) # 
    spike(x)
end

# ╔═╡ fe5732e0-6562-11eb-0a9b-a15a4acdb801
# when condition == 0 and upcrossing (from negative to positive) 
function affect!(integrator)      
    reset!(integrator.u)
end

# ╔═╡ f841cd00-6654-11eb-0a0d-f9c172a59637
md"#### Solvers & plots"

# ╔═╡ be6f1530-66ca-11eb-2c39-6b1a821a91e1
md"###### Low currrent"

# ╔═╡ c733ad50-64ab-11eb-087d-e79f7293c70b
begin
	u0=[-60.0,-40.0]
	tspan = (0.0,200.0)
	cb   = ContinuousCallback(condition,affect!,nothing)
	prob = ODEProblem(MQIF!,u0,tspan,p,callback=cb)
end

# ╔═╡ 69e526f0-64ac-11eb-21c0-ff1dcf1ef711
sol = solve(prob,dense=false,dtmax=0.01)

# ╔═╡ 333b3c00-6657-11eb-2876-536873d6a68a
md""" Function used for markers so that we keep the vectors x and y with N_wanted points in for a good display. It supposes that the steps in x are constants so that we keep N-wanted points with a constant step (which helps to see the evolution of y across x."""

# ╔═╡ 3225fb70-6657-11eb-3d47-871fb7a258b1
function split_for_markers(x,y,N_wanted)
	idx_step = Int(floor(length(x)/N_wanted))
	if idx_step == 0
		return NaN
	end
	split_x = zeros(N_wanted)
	split_y = zeros(N_wanted)
	j=1
	for i=1:N_wanted
		split_x[i]=x[j]
		split_y[i]=y[j]
		j=j+idx_step
	end
	if x[1]==0
		x_step = x[1+idx_step]
	else
		x_step = NaN
	end
	return split_x,split_y,x_step
end

# ╔═╡ d6475be0-64ae-11eb-0b97-a7455f721734
plotly()

# ╔═╡ 16c03280-6655-11eb-2a36-05188dd54bc4
md""" System trajectory for low current"""

# ╔═╡ e73b5fd0-66f9-11eb-137f-872a04c6a9cb
gr()

# ╔═╡ b9796ca0-654a-11eb-17f7-6fb798368c65
begin
	plot(cell_potential,V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines, low I ")
	plot!(cell_potential,V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines, low I ")
	plot!(cell_potential,V3,linewidth = 1.5,label="Vs-nullcline")
	plot!(sol,vars=(1,2),linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="Trajectory")
	
	sc_t,sc_V,dt_sc_v_ = split_for_markers(sol.t,sol[1,:],15)
	sc_t,sc_Vs = split_for_markers(sol.t,sol[2,:],15)
	dt_sc_v=round(dt_sc_v_*10)/10
	
	scatter!(sc_V,sc_Vs,markershape = :diamond,markersize = 2,markeralpha = 0.6,markercolor = RGB(0,0.6,0.8),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Dt = $dt_sc_v ms")
	
	scatter!([sol[1,1],sol[1,1]],[sol[2,1],sol[2,1]],markershape = :diamond,markersize = 2,markeralpha = 1,markercolor = RGB(1,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Initial")
	
	f1 = scatter!([sol[1,length(sol)],sol[1,length(sol)]],[sol[2,length(sol)],sol[2,length(sol)]],markershape = :diamond,markersize = 2,markeralpha = 1,markercolor = RGB(0.7,0,0.8),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Final",legend=:outertopright)
	
	yaxis!("Vs",(-60,0))
	xaxis!("V",(-70,0))
end

# ╔═╡ cd314fc0-66ca-11eb-1e8e-955995eaf311
md"###### High currrent"

# ╔═╡ 9cc8c8c0-663b-11eb-1dc9-fdb3997780e9
begin
	ph=zeros(length(p))
	for i=1:length(p)
		if i > 1
			ph[i]=p[i]
		else
			ph[i]=Ih
		end
	end
	md"""Change current value in the parameters vector"""
end

# ╔═╡ f0a870ae-654d-11eb-2f67-7d8909eb2993
begin	
	tspanh=(0.0,150.0)
	u0h =sol[length(sol)]
	probh = ODEProblem(MQIF!,u0h,tspanh,ph,callback=cb)
	solh = solve(probh,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
end

# ╔═╡ cbde61fe-665e-11eb-13d8-07c58a3ef63b
function disp_traj_resets(v,vs)
	idx_reset_ = findall(x->x>=Vmax-1,v)
	
	#find 1 index for each peak
	idx_reset=[]
	j=0
	for i=1:length(idx_reset_)
		if i==1
			idx_reset = push!(idx_reset,idx_reset_[i])
			j=j+1
		else
			if abs(idx_reset_[i-1]-idx_reset_[i])<=10
				idx_reset[j]=idx_reset_[i]
			else
				idx_reset = push!(idx_reset,idx_reset_[i])
				j=j+1
			end
		end
		if maximum(idx_reset)==maximum(idx_reset_)
			break
		end
	end
	#add the length of v (=length(vs)) so that the last part of the curve is added to the list 
	idx_reset = push!(idx_reset,length(v))
	
	#create a list with the voltages time evolution between 0 and the firs peak and between each other peaks
	
	id=1
	id_prev=1
	parts_v = []
	parts_vs = []
	for i=1:length(idx_reset)
		id = idx_reset[i]
		
		parts_v = push!(parts_v,v[id_prev:id])
		parts_vs = push!(parts_vs,vs[id_prev:id],)
		
		id_prev=id+1
	end
	
	return parts_v,parts_vs
end

# ╔═╡ efb781e0-66ca-11eb-2ca5-6d9d9b3af2ba
md""" System trajectory for high current"""

# ╔═╡ 45add010-6575-11eb-0ff1-c3f6c939f522
begin
	V1h = zeros(size(cell_potential))
	V2h = zeros(size(cell_potential))
	
	for i=1:length(cell_potential)
		V1h[i],V2h[i]=Vnullcline(cell_potential[i],ph)
	end

	plot(cell_potential,V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines, low I ")
	plot!(cell_potential,V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="Low I V-nullclines, low I ")
	plot!(cell_potential,V1h,linecolor=RGB(1,0.5,0),linewidth = 1.5,label="V-nullclines, high I ")
	plot!(cell_potential,V2h,linecolor=RGB(1,0.5,0),linewidth = 1.5,label="V-nullclines, high I ")
	plot!(cell_potential,V3,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")#does not depend on the current 
	plot!(solh[1,:],solh[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory")
	yaxis!((-65,-10))
	xaxis!((-65,0))
end

# ╔═╡ 468a3b70-6577-11eb-26c4-376c402f17f1
begin
	p_v,p_vs = disp_traj_resets(solh[1,:],solh[2,:])
	
	plot(cell_potential,V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,linestyle=:dot,label="V-nullclines, low I ")
	plot!(cell_potential,V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,linestyle=:dot,label="V-nullclines, low I ")
	plot!(cell_potential,V1h,linecolor=RGB(1,0.5,0),linewidth = 1.5,label="V-nullclines, high I ")
	plot!(cell_potential,V2h,linecolor=RGB(1,0.5,0),linewidth = 1.5,label="V-nullclines, high I ")
	plot!(cell_potential,V3,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	for i=1:length(p_v)
		plot!(p_v[i],p_vs[i],linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="Trajectory",legend=:outertopright)
	end
	
	scatter!([p_v[1][length(p_v[1])]],
		[p_vs[1][length(p_vs[1])]],
		marker=(:circle, 70, 1, RGB(1,0,0.5), stroke(0.1, 0, :black, :dot)),label="Reset")
	
	scatter!([p_v[2][length(p_v[2])]],
		[p_vs[2][length(p_vs[2])]],
		marker=(:circle, 70, 1, RGB(1,0,0.5), stroke(0.1, 0, :black, :dot)),label="Reset")
	
	scatter!([p_v[length(p_v)][1]],
	[p_vs[length(p_v)][1]],
	marker=(:circle, 70, 1, RGB(1,0,0.5), stroke(0.1, 0, :black, :dot)),label="Reset")
	
	scatter!([solh[1,1]],[solh[2,1]],markershape = :diamond,markersize = 2,markeralpha = 1,markercolor = RGB(1,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0.2,markerstrokecolor = :black,markerstrokestyle = :dot,label="Initial")
	
	f2 = yaxis!("Vs",(-60,0))
	xaxis!("V",(-70,0))
	
	
	
end

# ╔═╡ 8bcd0190-66ee-11eb-0dc4-f5d028d63ae2
md"###### Low current"

# ╔═╡ 8456b260-66ff-11eb-37ae-b9ca83e01afc
begin
	pl=zeros(length(p))
	for i=1:length(p)
		if i > 1
			pl[i]=p[i]
		else
			pl[i]=Il
		end
	end
	md"""Change current value in the parameters vector"""
end

# ╔═╡ c9180910-66ed-11eb-0f7d-251c6f0944c6
begin
	u0l=solh[length(solh)]
	tspanl=(0.0,150.)
	probl = ODEProblem(MQIF!,u0l,tspanl,pl,callback=cb)
	soll = solve(probl,dense=false,dtmax=0.01)
end

# ╔═╡ a9fdbc40-66ee-11eb-2110-add301449188
begin
	plot(cell_potential,V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines, low I ")
	plot!(cell_potential,V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="Low I V-nullclines, low I ")
	plot!(cell_potential,V1h,linecolor=RGB(1,0.5,0),linewidth = 1.5,linestyle=:dot,label="V-nullclines, high I ")
	plot!(cell_potential,V2h,linecolor=RGB(1,0.5,0),linewidth = 1.5,linestyle=:dot,label="V-nullclines, high I ")
	plot!(cell_potential,V3,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	plot!(soll[1,:],soll[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory")
	
	sc_tl,sc_V_l,dt_sc_v_l_ = split_for_markers(soll.t,soll[1,:],15)
	sc_tl,sc_Vs_l = split_for_markers(soll.t,soll[2,:],15)
	dt_sc_v_l=round(dt_sc_v_l_*10)/10
	
	scatter!(sc_V_l,sc_Vs_l,markershape = :diamond,markersize = 2,markeralpha = 0.6,markercolor = RGB(0,0.6,0.8),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Dt = $dt_sc_v_l ms")
	
	scatter!([soll[1,1],soll[1,1]],[soll[2,1],soll[2,1]],markershape = :diamond,markersize = 2,markeralpha = 1,markercolor = RGB(1,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Initial")
	
	f3 = scatter!([soll[1,length(soll)],soll[1,length(soll)]],[soll[2,length(soll)],soll[2,length(soll)]],markershape = :diamond,markersize = 2,markeralpha = 1,markercolor = RGB(0.7,0,0.8),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Final")
	
	yaxis!("Vs",(-60,0))
	xaxis!("V",(-70,0))
end

# ╔═╡ 50a384de-66cb-11eb-1169-059806818ba0
md""" Join low and high current to show the final result of the simulation"""

# ╔═╡ 857d7b00-6642-11eb-10a8-35dd3fee920a
begin
	time = zeros(length(sol.t)+length(solh.t)+length(soll.t))
	
	current = ones(length(sol.t)+length(solh.t)+length(soll.t)) 
	current[1:length(sol.t)].=Il
	current[(length(sol.t)+1):(length(sol.t)+length(solh.t))].=Ih
	current[(length(sol.t)+length(solh.t)+1):length(current)].=Il
	
	Vt = ones(length(sol.t)+length(solh.t)+length(soll.t)) 
	
	#start from resting point
	#Vt[1:length(sol.t)].=solh[1,1]
	
	Vst = ones(length(sol.t)+length(solh.t)+length(soll.t)) 
	
	#start from resting point
	#Vst[1:length(sol.t)].=solh[2,1]

	for i=1:length(time)
		if i <=length(sol.t)
			time[i]=sol.t[i]
			Vt[i]=sol[1,i]
			Vst[i]=sol[2,i]
		else
			if i > length(sol.t) && i <= length(solh.t)+length(sol.t)
				time[i]=solh.t[i-length(sol.t)]+time[length(sol.t)]
				
				Vt[i]=solh[1,i-length(sol.t)]
				Vst[i]=solh[2,i-length(sol.t)]
			else
				time[i]=soll.t[i-length(sol.t)-length(solh.t)]+time[length(sol.t)+length(solh.t)]
				
				Vt[i]=soll[1,i-length(sol.t)-length(solh.t)]
				Vst[i]=soll[2,i-length(sol.t)-length(solh.t)]
			end
		end
	end
end

# ╔═╡ 57dc9440-66cb-11eb-12ca-0d8366d81ffd
md"###### Final simulation"

# ╔═╡ 24f156b0-654f-11eb-11b5-71ee401f0a67
begin 
	plo = plot(time,current,linecolor=RGB(0.7,0,0.1),label="I(t)")
	
	pvh = plot(time,Vt,label="V(t)")
	
	pvsh = plot(time,Vst,linecolor=RGB(0,0.7,0.1),label="Vs(t)")
	xaxis!("Time (ms)")	
	
	mysub = plot(plo,pvh,pvsh,layout=@layout([a{0.1h};b ;c]),linewidth = 1.5,legend=:outertopright)
	#size=(1200,600)
	
end

# ╔═╡ 758691f0-66fb-11eb-0bf8-7fe8ae75ad1b
#savefig("MQIF_2D_Fig6A7A-time.pdf")

# ╔═╡ 9b9f2480-66f9-11eb-3c16-174fe3e6ac25
gr()

# ╔═╡ b6d90a30-66f5-11eb-220c-df794163eb6c
begin
	finalplot = plot(f1,f2,f3,layout=(1,3),size=(1800,600))
end

# ╔═╡ fdac2f40-66f6-11eb-21ec-c59f7fe964d4
#savefig("MQIF_2D_Fig6A7A.pdf")

# ╔═╡ 3d281e40-6a1c-11eb-32a4-219e323e9af9
md" #### Phase portraits"

# ╔═╡ 47744e20-6a29-11eb-0485-dfbff003e543
begin
	limVmin =-60.0
	limVmax =-20.0
	limVsmin=-60.0
	limVsmax=-20.0
end

# ╔═╡ 42085da0-6a29-11eb-3d97-5faa899d3c62
v_potential = collect(range(limVmin-5,stop=limVmax+5,length=100))

# ╔═╡ 31495e50-6a4d-11eb-21b8-21c20147d5c8
begin
	vs_null=zeros(size(v_potential))
	for i=1:length(v_potential)
		vs_null[i] = Vsnullcline(v_potential[i],p)
	end
end

# ╔═╡ aeab4d6e-6a1d-11eb-1e4c-b3cfe5e2d195
vs_potential = collect(range(limVsmin-5,stop=limVsmax+5,length=100))

# ╔═╡ 62d7e340-6a22-11eb-0270-f1d22b89fe1c
md""" Mesh """

# ╔═╡ 6cc0a950-6a22-11eb-0cd5-f54f55ab18f1
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# ╔═╡ 4190fd10-7f59-11eb-18ce-b9ad0b8b8971
x, y = meshgrid(limVmin:3.0:limVmax, limVsmin:3.0:limVsmax)

# ╔═╡ 3ac785c0-6b8c-11eb-12a0-ab81fa7fc6d1
xx, yy = meshgrid(limVmin:4.0:limVmax, limVsmin:4.0:limVsmax)

# ╔═╡ 42679f30-6a2a-11eb-3476-57a1d9ae121d
scale=[0.015,0.5]

# ╔═╡ 92f34f60-6a27-11eb-37cf-7773810bca67
function return_MQIF_gradient(u,p,t)
	
	grad = [0.0,0.0]
	MQIF!(grad,u,p,0.0)	
	
	return grad
end

# ╔═╡ 335a246e-6a22-11eb-0035-c564fda8d2fa
md"###### I=0"

# ╔═╡ 3ce5efb0-6a1d-11eb-0c9a-3b230f943f01
begin 
	Vl1 = zeros(size(v_potential))
	Vl2 = zeros(size(v_potential))
	
	for i=1:length(v_potential)
		Vl1[i],Vl2[i]=Vnullcline(v_potential[i],p)
	end
	md"""Nullclines for I = 0 """ #Il
end

# ╔═╡ 7fab57e0-6a22-11eb-2004-09b13178e043
begin
	dx_l = zeros(size(x))
	dy_l = zeros(size(x))
	for k=1:length(x)
		dx_l[k],dy_l[k]=return_MQIF_gradient([x[k],y[k]],p,0.0)	
	end
	dx_l = dx_l*scale[1]
	dy_l = dy_l*scale[2]
	md""" Gradient"""
end

# ╔═╡ 4375ad70-6a22-11eb-3579-c581e408c2a0
md"###### I=5"

# ╔═╡ 7e404870-6a1d-11eb-2d4d-8d972c5a0993
begin
	#Change current value in the parameters vector
	pp=zeros(length(p))
	for i=1:length(p)
		if i > 1
			pp[i]=p[i]
		else
			pp[i]=Ip
		end
	end
	
	#V-nullclines computation
	Vp1 = zeros(size(v_potential))
	Vp2 = zeros(size(v_potential))
	
	for i=1:length(v_potential)
		Vp1[i],Vp2[i]=Vnullcline(v_potential[i],pp)
	end
	Vp1
	Vp2
	md"""Nullclines for I = 5""" #Ip
end

# ╔═╡ 2a02c4c0-6a29-11eb-3c98-7135a450b2c2
begin
	dx_p = zeros(size(x))
	dy_p = zeros(size(x))
	for k=1:length(x)
		dx_p[k],dy_p[k]=return_MQIF_gradient([x[k],y[k]],pp,0.0)	
	end
	dx_p = dx_p*scale[1]
	dy_p = dy_p*scale[2]
	md""" Gradient"""
end

# ╔═╡ 4af5d7f0-6a22-11eb-11c1-07234671e997
md"###### I=30"

# ╔═╡ 7f1b3f20-6a1d-11eb-31b5-8f6eda38efb3
begin	
	#V-nullclines computation
	Vh1 = zeros(size(v_potential))
	Vh2 = zeros(size(v_potential))
	
	for i=1:length(v_potential)
		Vh1[i],Vh2[i]=Vnullcline(v_potential[i],ph)
	end
	Vh1
	Vh2
	md"""Nullclines for I = 30""" #Ih
end

# ╔═╡ a1af8f20-6a2a-11eb-099c-053cd7e54835
begin
	dx_h = zeros(size(xx))
	dy_h = zeros(size(xx))
	for k=1:length(xx)
		dx_h[k],dy_h[k]=return_MQIF_gradient([xx[k],yy[k]],ph,0.0)	
	end
	dx_h = dx_h*scale[1]
	dy_h = dy_h*scale[2]
	md""" Gradient"""
end

# ╔═╡ de713f30-6a25-11eb-1c55-d7ef1f344cd0
md"###### I=-5"

# ╔═╡ de03fe72-6a25-11eb-0d51-6743c6138f0c
begin 
	#Change current value in the parameters vector
	pn=zeros(length(p))
	for i=1:length(p)
		if i > 1
			pn[i]=p[i]
		else
			pn[i]=In
		end
	end
	
	#V-nullclines computation
	Vn1 = zeros(size(vs_potential))
	Vn2 = zeros(size(vs_potential))
	
	for i=1:length(vs_potential)
		Vn1[i],Vn2[i]=Vnullcline_(vs_potential[i],pn)
	end
	Vn1
	Vn2
	md"""Nullclines for I = -5""" #In
end

# ╔═╡ b3f3c480-6a2a-11eb-1485-877a22b19db1
begin
	dx_n = zeros(size(x))
	dy_n = zeros(size(x))
	for k=1:length(x)
		dx_n[k],dy_n[k]=return_MQIF_gradient([x[k],y[k]],pn,0.0)	
	end
	dx_n = dx_n*scale[1]
	dy_n = dy_n*scale[2]
	md""" Gradient"""
end

# ╔═╡ ddcf32d0-6a25-11eb-02ba-6f8e3097569a
md"###### I=-30"

# ╔═╡ dd713450-6a25-11eb-3c88-8fc64c499bbe
begin 
	#Change current value in the parameters vector
	pnh=zeros(length(p))
	for i=1:length(p)
		if i > 1
			pnh[i]=p[i]
		else
			pnh[i]=Inh
		end
	end
	
	#V-nullclines computation
	Vnh1 = zeros(size(vs_potential))
	Vnh2 = zeros(size(vs_potential))
	
	for i=1:length(vs_potential)
		Vnh1[i],Vnh2[i]=Vnullcline_(vs_potential[i],pnh)
	end
	Vnh1
	Vnh2
	md"""Nullclines for I = -30""" #Inh
end

# ╔═╡ c868204e-6a2a-11eb-0f6a-6d945a550bf5
begin
	dx_nh = zeros(size(x))
	dy_nh = zeros(size(x))
	for k=1:length(x)
		dx_nh[k],dy_nh[k]=return_MQIF_gradient([x[k],y[k]],pnh,0.0)	
	end
	dx_nh = dx_nh*scale[1]
	dy_nh = dy_nh*scale[2]
	md""" Gradient"""
end

# ╔═╡ 500fd7e0-6a22-11eb-3dff-27c30b5d8460
md"###### Plots"

# ╔═╡ b9f04b50-6a21-11eb-031d-fdb7c4647aed
gr()

# ╔═╡ 90cffe60-6a20-11eb-159d-71fd0b3d58ee
begin 		
	plot(Vnh1,vs_potential,linecolor=RGB(0.86,0.06,0.24),label="V-nullclines")
	plot!(Vnh2,vs_potential,linecolor=RGB(0.86,0.06,0.24),label="V-nullclines")
	plot!(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	pin30 = quiver!(x, y, quiver=(dx_nh, dy_nh),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	
	nhfp1,nhfp2,nhstab=fixedpoints(pnh)
	nhfp=[nhfp1,nhfp2]
	if length(nhstab)>1
		for i=1:length(nhstab)
			if nhstab[i] > 0
				scatter!([nhfp[i][1]],[nhfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if nhstab[i] < 0
				scatter!([nhfp[i][1]],[nhfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
			if nhstab[i] == 0
				scatter!([nhfp[i][1]],[nhfp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright)#saddle
			end
		end
	end
	
	title!("I = $Inh ")
	yaxis!("Vs",(limVsmin,limVsmax))
	xaxis!("V",(limVmin,limVmax))
	
	#savefig("MQIF_2D_phaseportraits-30.png")
	
end 

# ╔═╡ 46d34140-6a2b-11eb-053b-27fa5e6b78bd
begin
	plot(Vn1,vs_potential, linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(Vn2,vs_potential,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	pin5 = quiver!(x, y, quiver=(dx_n, dy_n),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	
	nfp1,nfp2,nstab=fixedpoints(pn)
	nfp=[nfp1,nfp2]
	if length(nstab)>1
		for i=1:length(nstab)
			if nstab[i] > 0
				scatter!([nfp[i][1]],[nfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if nstab[i] < 0
				scatter!([nfp[i][1]],[nfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
			if nstab[i] == 0
				scatter!([nfp[i][1]],[nfp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright)#saddle
			end
		end
	end
	
	title!("I = $In ")
	yaxis!("Vs",(limVsmin,limVsmax))
	xaxis!("V",(limVmin,limVmax))
	
	#savefig("MQIF_2D_phaseportraits-5.png")
end

# ╔═╡ 8ae29620-6a52-11eb-1288-05ff0a31d717
begin
	plot(v_potential,Vl1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,Vl2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V3,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	pi0 = quiver!(x, y, quiver=(dx_l, dy_l),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	
	lfp1,lfp2,lstab=fixedpoints(pl)
	lfp=[lfp1,lfp2]
	if length(lstab)>1
		for i=1:length(lstab)
			if lstab[i] > 0
				scatter!([lfp[i][1]],[lfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if lstab[i] < 0
				scatter!([lfp[i][1]],[lfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
			if lstab[i] == 0
				scatter!([lfp[i][1]],[lfp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright)#saddle
			end
		end
	end
	
	tspan0=(0.0,300.0)
	u00 =[-20,-40]
	prob0 = ODEProblem(MQIF!,u00,tspan0,pl,callback=cb)
	sol0 = solve(prob0,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
	plot!(sol0[1,:],sol0[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory",linecolor=RGB(0,0,0))
	
	title!("I = $Il ")
	yaxis!("Vs",(limVsmin,limVsmax))
	xaxis!("V",(limVmin,limVmax))
	
	#savefig("MQIF_2D_phaseportraits0.png")
end

# ╔═╡ 99a2fbf0-6a52-11eb-3dce-57131449d473
begin
	plot(v_potential,Vp1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,Vp2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	pi5 = quiver!(x, y, quiver=(dx_p, dy_p),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	
	pfp1,pfp2,pstab=fixedpoints(pp)
	pfp=[pfp1,pfp2]
	if length(pstab)>1
		for i=1:length(pstab)
			if pstab[i] > 0
				scatter!([pfp[i][1]],[pfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if pstab[i] < 0
				scatter!([pfp[i][1]],[pfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
			if pstab[i] == 0
				scatter!([pfp[i][1]],[pfp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
			end
		end
	end
	
	title!("I = $Ip ")
	yaxis!("Vs",(limVsmin,limVsmax))
	xaxis!("V",(limVmin,limVmax))
	
	#savefig("MQIF_2D_phaseportraits5.png")
end

# ╔═╡ acbe7e80-6a52-11eb-0383-71d5ce18e3ac
begin
	plot(v_potential,Vh1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,Vh2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	quiver!(xx, yy, quiver=(dx_h, dy_h),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	
	hfp1,hfp2,hstab=fixedpoints(ph)
	hfp=[hfp1,hfp2]
	if length(hstab)>1
		for i=1:length(hstab)
			if Int(hstab[i]) > 0
				scatter!([hfp[i][1]],[hfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if Int(hstab[i]) < 0
				scatter!([hfp[i][1]],[hfp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
			if Int(hstab[i]) == 0
				scatter!([hfp[i][1]],[hfp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright)#saddle
			end
		end
	end
	
	
	title!("I = $Ih ")
	yaxis!("Vs",(limVsmin,limVsmax))
	xaxis!("V",(limVmin,limVmax))
	
	#savefig("MQIF_2D_phaseportraits30.png")
end

# ╔═╡ c87b36fe-6b86-11eb-0a00-63d6799af669
md" #### Bifurcation diagram with I"

# ╔═╡ 15d6004e-6b89-11eb-14f2-216a730247b3
md""" There are 3 parts in the bifurcation diagram

               I <  0.355 
Stable resting state  

      0.355 <= I <= 25    
Stable resting state, stable cycle
   
         25 <  I          
Stable cycle"""

# ╔═╡ ce724fb0-6bb6-11eb-2671-cf3325b6d209
Ibis = 0.355

# ╔═╡ 602d37b2-6bb9-11eb-3b43-4ffcd0c06948
Ibif = 25

# ╔═╡ f7f54c50-6b86-11eb-0ae2-c9a0c1cbe53c
md" ###### I<0.355"

# ╔═╡ ca3bd1f2-6bb6-11eb-1670-6d6c3e2a68c8
Ipart1 = collect(range(-50,stop=Ibis-0.005,length=100))

# ╔═╡ 0b7de7c0-6bb7-11eb-1707-8da1e13d9992
begin
	list_stable_1 = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	list_saddle_1 = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	list_unstable_1 = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	for i=1:length(Ipart1)
		pI_1=change_I(Ipart1[i])
		#pI_1=change_vs0(-45,pI_1)
		fp1_1,fp2_1,stability_1 = fixedpoints(pI_1)
		fp_1=[fp1_1,fp2_1]
		for j=1:length(stability_1)
			if stability_1[j]<0
				list_stable_1[i,:]=fp_1[j][:]
			else
				if stability_1[j]==0
					list_saddle_1[i,:]=fp_1[j][:]
				else
					list_unstable_1[i,:]=fp_1[j][:]
				end
			end
		end
	end
	md""" Calculus of stable points and saddle points for each I, regenerative vs0"""
end

# ╔═╡ 7f4e4fb0-7dea-11eb-2246-5dff4477acb7
begin
	list_stable_1_resto = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	list_saddle_1_resto = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	list_unstable_1_resto = zeros(length(Ipart1),2) #column 1 : V , column 2 : Vs
	for i=1:length(Ipart1)
		pI_1=change_I(Ipart1[i])
		pI_1=change_vs0(-45,pI_1)
		fp1_1,fp2_1,stability_1 = fixedpoints(pI_1)
		fp_=[fp1_1,fp2_1]
		for j=1:length(stability_1)
			if stability_1[j]<0
				list_stable_1_resto[i,:]=fp_[j][:]
				list_unstable_1_resto[i,:]=[NaN,NaN]
			else
				if stability_1[j]==0
					list_saddle_1_resto[i,:]=fp_[j][:]
				else
					list_unstable_1_resto[i,:]=fp_[j][:]
					list_stable_1_resto[i,:]=[NaN,NaN]
				end
			end
		end
	end
	md""" Calculus of stable points and saddle points for each I, restorative vs0"""
end

# ╔═╡ 4732de90-6bb9-11eb-21bc-d9f23a7e5b02
md" ###### 0.355 <= I <= 25   "

# ╔═╡ 5690e3f0-6bb9-11eb-2ef7-a9fc4505ebfa
Ipart2 = collect(range(Ibis,stop=Ibif-0.005,length=100))

# ╔═╡ 825c4920-6bb9-11eb-1867-6b8e08310024
begin
	list_stable_2 = zeros(length(Ipart2),2) #column 1 : V , column 2 : Vs
	list_cycleV_2 = zeros(length(Ipart2),2) #column 1 : max, column 2 : min
	list_cycleVs_2 = zeros(length(Ipart2),2) #column 1 : max, column 2 : min
	cycle_freq_2 = zeros(length(Ipart2))
	cycle_m_freq_2 = zeros(length(Ipart2))
	d_reset_min_2 = zeros(length(Ipart2))
	
	list_saddle_2 = zeros(length(Ipart2),2) #column 1 : V , column 2 : Vs
	
	tspan_2=(0.0,200.0)
	u0_2 =[-20,-40]
	
	for i=1:length(Ipart2)
		pI_2=change_I(Ipart2[i])
		fp1_2,fp2_2,stability_2 = fixedpoints(pI_2)
		fp_2=[fp1_2,fp2_2]
		for j=1:length(stability_2)
			if stability_2[j]<0
				list_stable_2[i,:]=fp_2[j][:]
			else
				list_saddle_2[i,:]=fp_2[j][:]
			end
		end
		
		prob_2 = ODEProblem(MQIF!,u0_2,tspan_2,pI_2,callback=cb)
		sol_2 = solve(prob_2,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_2 = findfirst(x -> x>=maximum(tspan_2)/10, sol_2.t)
		list_cycleV_2[i,1] = maximum(sol_2[1,ind_tmin_2:length(sol_2)])
		list_cycleV_2[i,2] = minimum(sol_2[1,ind_tmin_2:length(sol_2)])
		list_cycleVs_2[i,1] = maximum(sol_2[2,ind_tmin_2:length(sol_2)])
		list_cycleVs_2[i,2] = minimum(sol_2[2,ind_tmin_2:length(sol_2)])
		ind_spike_2 = findall(x-> x==Vr,sol_2[1,ind_tmin_2:length(sol_2)]).+ind_tmin_2
		t_spike_2 = sol_2.t[ind_spike_2[2:end]]-sol_2.t[ind_spike_2[1:end-1]]
		cycle_freq_2[i]=1/mean(t_spike_2)
		dt_1_end_spikes_2 = sol_2.t[ind_spike_2[end]]-sol_2.t[ind_spike_2[1]]
		cycle_m_freq_2[i]= length(ind_spike_2)/dt_1_end_spikes_2
		
		d_reset_min_2[i] = sqrt((pI_2[2]-Vr)^2 + (pI_2[3]+sqrt(pI_2[1]/pI_2[6])-Vsr)^2)
	end	
	md""" Calculus of stable points, saddle points, frequency, maximum and minimum values of the limit cycle for each I, regenerative vs0"""
end

# ╔═╡ c676b620-7dea-11eb-2d3f-0dd5a7da39db
begin
	list_stable_2_resto = zeros(length(Ipart2),2) #column 1 : V , column 2 : Vs
	list_cycleV_2_resto = zeros(length(Ipart2),2) #column 1 : max, column 2 : min
	list_cycleVs_2_resto = zeros(length(Ipart2),2) #column 1 : max, column 2 : min
	
	list_saddle_2_resto = zeros(length(Ipart2),2) #column 1 : V , column 2 : Vs
	list_unstable_2_resto = zeros(length(Ipart2),2) #column 1 : V , column 2 : Vs
	cycle_freq_2_resto = zeros(length(Ipart2))
	d_reset_min_2_resto = zeros(length(Ipart2))
	
	for i=1:length(Ipart2)
		pI_2=change_I(Ipart2[i])
		pI_2=change_vs0(-45,pI_2)
		fp1_2,fp2_2,stability_2 = fixedpoints(pI_2)
		fp_2=[fp1_2,fp2_2]
		for j=1:length(stability_2)
			if stability_2[j]<0
				list_stable_2_resto[i,:]=fp_2[j][:]
				list_unstable_2_resto[i,:]=[NaN,NaN]
			else
				if stability_2[j]==0
					list_saddle_2_resto[i,:]=fp_2[j][:]
				else
					list_unstable_2_resto[i,:]=fp_2[j][:]
					list_stable_2_resto[i,:]=[NaN,NaN]
				end
			end
		end
		
		prob_2 = ODEProblem(MQIF!,u0_2,tspan_2,pI_2,callback=cb)
		sol_2 = solve(prob_2,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_2 = findfirst(x -> x>=3*maximum(tspan_2)/10, sol_2.t)
		list_cycleV_2_resto[i,1] = maximum(sol_2[1,ind_tmin_2:length(sol_2)])
		list_cycleV_2_resto[i,2] = minimum(sol_2[1,ind_tmin_2:length(sol_2)])
		list_cycleVs_2_resto[i,1] = maximum(sol_2[2,ind_tmin_2:length(sol_2)])
		list_cycleVs_2_resto[i,2] = minimum(sol_2[2,ind_tmin_2:length(sol_2)])
		if list_cycleV_2_resto[i,1]-list_cycleV_2_resto[i,2]<= 30
			list_cycleV_2_resto[i,:]=[NaN,NaN] 
			list_cycleVs_2_resto[i,:]=[NaN,NaN] 
			cycle_freq_2_resto[i]=NaN
		else
			ind_spike_2_resto = findall(x-> x==Vr,sol_2[1,ind_tmin_2:length(sol_2)]).+ind_tmin_2
			t_spike_2_resto = sol_2.t[ind_spike_2_resto[2:end]]-sol_2.t[ind_spike_2_resto[1:end-1]]
			cycle_freq_2_resto[i]=1/mean(t_spike_2_resto)
		end
		
		d_reset_min_2_resto[i] = sqrt((pI_2[2]-Vr)^2 + (pI_2[3]+sqrt(pI_2[1]/pI_2[6])-Vsr)^2)
	end	
	md""" Calculus of stable points, saddle points, frequency, maximum and minimum values of the limit cycle for each I, restorative vs0"""
end

# ╔═╡ 53576700-6bbd-11eb-3f52-a7e87ea6224d
md" ######  25 <  I    "

# ╔═╡ 640ebad0-6bbd-11eb-2603-b997744c46a4
Ipart3 = collect(range(Ibif+0.005,stop=200,length=70))

# ╔═╡ 82f18770-6bbd-11eb-0f2a-9db245be295d
begin
	list_cycleV_3 = zeros(length(Ipart3),2) #column 1 : max, column 2 : min
	list_cycleVs_3 = zeros(length(Ipart3),2) #column 1 : max, column 2 : min
	cycle_freq_3 = zeros(length(Ipart3))
	cycle_m_freq_3 = zeros(length(Ipart3))
	d_reset_min_3 = zeros(length(Ipart3))
	
	tspan_3=(0.0,150.0)
	u0_3 =[-20,-40]
	
	for i=1:length(Ipart3)
		pI_3=change_I(Ipart3[i])
		prob_3 = ODEProblem(MQIF!,u0_3,tspan_3,pI_3,callback=cb)
		sol_3 = solve(prob_3,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_3 = findfirst(x -> x>=maximum(tspan_3)/10, sol_3.t)
		list_cycleV_3[i,1] = maximum(sol_3[1,ind_tmin_3:length(sol_3)])
		list_cycleV_3[i,2] = minimum(sol_3[1,ind_tmin_3:length(sol_3)])
		list_cycleVs_3[i,1] = maximum(sol_3[2,ind_tmin_3:length(sol_3)])
		list_cycleVs_3[i,2] = minimum(sol_3[2,ind_tmin_3:length(sol_3)])
		ind_spike_3 = findall(x-> x==Vr,sol_3[1,ind_tmin_3:length(sol_3)]).+ind_tmin_3
		t_spike_3 = sol_3.t[ind_spike_3[2:end]]-sol_3.t[ind_spike_3[1:end-1]]
		cycle_freq_3[i]=1/mean(t_spike_3)
		dt_1_end_spikes_3 = sol_3.t[ind_spike_3[end]]-sol_3.t[ind_spike_3[1]]
		cycle_m_freq_3[i]= length(ind_spike_3)/dt_1_end_spikes_3
		
		d_reset_min_3[i] = sqrt((pI_3[2]-Vr)^2 + (pI_3[3]+sqrt(pI_3[1]/pI_3[6])-Vsr)^2)
	end	
	md""" Calculus of frequency, maximum and minimum values of the limit cycle for each I, regenerative vs0"""
end

# ╔═╡ 4f09f010-7deb-11eb-0aca-fd4072ea5648
begin
	list_cycleV_3_resto = zeros(length(Ipart3),2) #column 1 : max, column 2 : min
	list_cycleVs_3_resto = zeros(length(Ipart3),2) #column 1 : max, column 2 : min
	cycle_freq_3_resto = zeros(length(Ipart3))
	d_reset_min_3_resto = zeros(length(Ipart3))
	
	for i=1:length(Ipart3)
		pI_3=change_I(Ipart3[i])
		pI_3=change_vs0(-45,pI_3)
		prob_3 = ODEProblem(MQIF!,u0_3,tspan_3,pI_3,callback=cb)
		sol_3 = solve(prob_3,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_3 = findfirst(x -> x>=maximum(tspan_3)/10, sol_3.t)
		list_cycleV_3_resto[i,1] = maximum(sol_3[1,ind_tmin_3:length(sol_3)])
		list_cycleV_3_resto[i,2] = minimum(sol_3[1,ind_tmin_3:length(sol_3)])
		list_cycleVs_3_resto[i,1] = maximum(sol_3[2,ind_tmin_3:length(sol_3)])
		list_cycleVs_3_resto[i,2] = minimum(sol_3[2,ind_tmin_3:length(sol_3)])
		ind_spike_3_resto = findall(x-> x==Vr,sol_3[1,ind_tmin_3:length(sol_3)]).+ind_tmin_3
		t_spike_3_resto = sol_3.t[ind_spike_3_resto[2:end]]-sol_3.t[ind_spike_3_resto[1:end-1]]
		cycle_freq_3_resto[i]=1/mean(t_spike_3_resto)
		
		d_reset_min_3_resto[i] = sqrt((pI_3[2]-Vr)^2 + (pI_3[3]+sqrt(pI_3[1]/pI_3[6])-Vsr)^2)
	end	
	md""" Calculus of frequency, maximum and minimum values of the limit cycle for each I, restorative vs0"""
end

# ╔═╡ 6a811010-6bc3-11eb-202f-23e520ba250f
gr()

# ╔═╡ 96c46370-7e92-11eb-0e5c-8fb1ebe10c89
md" ######  Regenerative vs0 (vs0>v0)    "

# ╔═╡ fb4ddda0-6bbc-11eb-0acf-db128b83a622
begin
	Ipart12 = zeros(1,length(Ipart1)+length(Ipart2))
	list_stable_12 = zeros(length(Ipart1)+length(Ipart2),2)
	list_saddle_12 = zeros(length(Ipart1)+length(Ipart2),2)
	for i=1:length(Ipart1)+length(Ipart2)
		if i<=length(Ipart1)
			Ipart12[i]=Ipart1[i]
			list_stable_12[i,:]=list_stable_1[i,:]
			list_saddle_12[i,:]=list_saddle_1[i,:]
		else
			Ipart12[i]=Ipart2[i-length(Ipart1)]
			list_stable_12[i,:]=list_stable_2[i-length(Ipart1),:]
			list_saddle_12[i,:]=list_saddle_2[i-length(Ipart1),:]
		end
	end
	
	Ipart23 = zeros(1,length(Ipart2)+length(Ipart3))
	list_cycleV_23 = zeros(length(Ipart2)+length(Ipart3),2)
	list_cycleVs_23 = zeros(length(Ipart2)+length(Ipart3),2)
	cycle_freq_23 = zeros(length(Ipart2)+length(Ipart3))
	d_reset_min_23 = zeros(length(Ipart2)+length(Ipart3))
	for i=1:length(Ipart2)+length(Ipart3)
		if i<=length(Ipart2)
			Ipart23[i]=Ipart2[i]
			list_cycleV_23[i,:]=list_cycleV_2[i,:]
			list_cycleVs_23[i,:]=list_cycleVs_2[i,:]
			cycle_freq_23[i]=cycle_freq_2[i]
			d_reset_min_23[i]=d_reset_min_2[i]
		else
			Ipart23[i]=Ipart3[i-length(Ipart2)]
			list_cycleV_23[i,:]=list_cycleV_3[i-length(Ipart2),:]
			list_cycleVs_23[i,:]=list_cycleVs_3[i-length(Ipart2),:]
			cycle_freq_23[i]=cycle_freq_3[i-length(Ipart2)]
			d_reset_min_23[i]=d_reset_min_3[i-length(Ipart2)]
		end
	end
	
	bifIV =plot(Ipart12', list_stable_12[:,1],label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)
	plot!(Ipart12', list_saddle_12[:,1],label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)
	scatter!([Ibif],[-45],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	
	plot!(Ipart23', list_cycleV_23[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	plot!(Ipart23', list_cycleV_23[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,linestyle=:dot,legend=:best)
	yaxis!("V",(-65,-10))
	
	bifIVs = plot(Ipart12', list_stable_12[:,2],label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)
	plot!(Ipart12', list_saddle_12[:,2],label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)
	scatter!([Ibif],[-45],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	
	plot!(Ipart23', list_cycleVs_23[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0,0.4,1),linewidth=1.5,linestyle=:dot)
	plot!(Ipart23', list_cycleVs_23[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0,0.2,0.7),linewidth=1.5,linestyle=:dot,legend=:bottomright)
	yaxis!("Vs",(-65,-10))
	
	plot(bifIV,bifIVs,size=(800,400))
	xaxis!("I")
	
	#savefig("bifIV_regenerative_vs0.pdf")
	
end 

# ╔═╡ a71ea31e-7e92-11eb-0b09-11bca509a568
md" ######  Restorative vs0 (vs0<v0)    "

# ╔═╡ c445afe0-7deb-11eb-3726-0b309a0c35ec
begin
	list_stable_12_resto = zeros(length(Ipart1)+length(Ipart2),2)
	list_unstable_12_resto = zeros(length(Ipart1)+length(Ipart2),2)
	list_saddle_12_resto = zeros(length(Ipart1)+length(Ipart2),2)
	for i=1:length(Ipart1)+length(Ipart2)
		if i<=length(Ipart1)
			list_stable_12_resto[i,:]=list_stable_1_resto[i,:]
			list_unstable_12_resto[i,:]=list_unstable_1_resto[i,:]
			list_saddle_12_resto[i,:]=list_saddle_1_resto[i,:]
		else
			list_stable_12_resto[i,:]=list_stable_2_resto[i-length(Ipart1),:]
			list_unstable_12_resto[i,:]=list_unstable_2_resto[i-length(Ipart1),:]
			list_saddle_12_resto[i,:]=list_saddle_2_resto[i-length(Ipart1),:]
		end
	end
	
	list_cycleV_23_resto = zeros(length(Ipart2)+length(Ipart3),2)
	list_cycleVs_23_resto = zeros(length(Ipart2)+length(Ipart3),2)
	cycle_freq_23_resto = zeros(length(Ipart2)+length(Ipart3))
	d_reset_min_23_resto = zeros(length(Ipart2)+length(Ipart3))
	for i=1:length(Ipart2)+length(Ipart3)
		if i<=length(Ipart2)
			list_cycleV_23_resto[i,:]=list_cycleV_2_resto[i,:]
			list_cycleVs_23_resto[i,:]=list_cycleVs_2_resto[i,:]
			cycle_freq_23_resto[i] = cycle_freq_2_resto[i]
			d_reset_min_23_resto[i]=d_reset_min_2_resto[i]
		else
			list_cycleV_23_resto[i,:]=list_cycleV_3_resto[i-length(Ipart2),:]
			list_cycleVs_23_resto[i,:]=list_cycleVs_3_resto[i-length(Ipart2),:]
			cycle_freq_23_resto[i]=cycle_freq_3_resto[i-length(Ipart2)]
			d_reset_min_23_resto[i]=d_reset_min_3_resto[i-length(Ipart2)]
		end
	end
	
	bifIV_resto =plot(Ipart12', list_stable_12_resto[:,1],label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)
	plot!(Ipart12', list_unstable_12_resto[:,1],label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	plot!(Ipart12', list_saddle_12_resto[:,1],label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)
	scatter!([Ibif],[-35],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	
	plot!(Ipart23', list_cycleV_23_resto[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0.8,0,0.7),linewidth=1.5,linestyle=:dot)
	plot!(Ipart23', list_cycleV_23_resto[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0.5,0,0.4),linewidth=1.5,linestyle=:dot,legend=:best)
	yaxis!("V",(-65,-10))
	
	bifIVs_resto = plot(Ipart12', list_stable_12_resto[:,2],label="Stable state",linecolor=RGB(0,0.5,0),linewidth=1.5)
	plot!(Ipart12', list_unstable_12_resto[:,2],label="Unstable state",linecolor=RGB(0.7,0,0),linewidth=1.5)
	plot!(Ipart12', list_saddle_12_resto[:,2],label="Saddle",linecolor=RGB(0.8,0.6,0),linewidth=1.5,linestyle=:dashdot)
	scatter!([Ibif],[-35],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="SN bifurcation")
	
	plot!(Ipart23', list_cycleVs_23_resto[:,1],label="Stable limit cycle : maximum",linecolor=RGB(0,0.4,1),linewidth=1.5,linestyle=:dot)
	plot!(Ipart23', list_cycleVs_23_resto[:,2],label="Stable limit cycle : minimum",linecolor=RGB(0,0.2,0.7),linewidth=1.5,linestyle=:dot,legend=:bottomright)
	yaxis!("Vs",(-65,-10))
	
	bifI_resto = plot(bifIV_resto,bifIVs_resto,size=(800,400))
	plot(bifI_resto)
	xaxis!("I")
	
	#savefig("bifIV_restorative_vs0.pdf")
	
end 

# ╔═╡ 14386390-8e80-11eb-36ed-812d431e537b
md" ######  Limit cycle Frequency evolution with I    "

# ╔═╡ 0f30b320-8e80-11eb-184c-a58dfa5bfdf2
begin
	plot(Ipart23',cycle_freq_23.*1000,label="Regenerative",linewidth=2)
	plot!(Ipart23',cycle_freq_23_resto.*1000,label="Restorative",linewidth=2,legend=:outertopright)	
	xaxis!("I")
	yaxis!("Frequency [Hz]")
	
	#savefig("MQIF-frequency_variation_with_I.pdf")
end

# ╔═╡ 6f36d9f0-8ff4-11eb-2186-0938b0a1a009
begin
	d_reset_min_23_resto_ = zeros(length(d_reset_min_23_resto))
	d_reset_min_23_ = zeros(length(d_reset_min_23))
	for i=1:length(d_reset_min_23_resto)
		k=1
		d_reset_min_23_resto_[i] = k/ d_reset_min_23_resto[i]
		d_reset_min_23_[i] = k/ d_reset_min_23[i]
	end
	vd_reset_min_23 = zeros(length(d_reset_min_23))
	vd_reset_min_23_resto = zeros(length(d_reset_min_23_resto))
	for i=1:length(d_reset_min_23_resto)
		vd_reset_min_23[i]=sqrt(d_reset_min_23[i]^2+10^2)
		vd_reset_min_23_resto[i]=sqrt(d_reset_min_23_resto[i]^2+10^2)
	end
	
	vd_reset_min_23_ = zeros(length(vd_reset_min_23))
	vd_reset_min_23_resto_ = zeros(length(vd_reset_min_23_resto))
	for i=1:length(d_reset_min_23_resto)
		k=1
		vd_reset_min_23_resto_[i] = k/ vd_reset_min_23_resto[i]
		vd_reset_min_23_[i] = k/ vd_reset_min_23[i]
	end
	
	plot(Ipart23',d_reset_min_23_,label="Dist Regenerative",linewidth=2,legend=:outertopright)
	plot!(Ipart23',d_reset_min_23_resto_,label="Dist Restorative",linewidth=2,legend=:outertopright)
	plot!(Ipart23',vd_reset_min_23_,label="VDist Regenerative",linewidth=2,legend=:outertopright)
	plot!(Ipart23',vd_reset_min_23_resto_,label="VDist Restorative",linewidth=2,legend=:outertopright)
	
	xaxis!("I")
	yaxis!("Inverse of distance nullcline minimum to reset")
end

# ╔═╡ 178a21b0-6bc5-11eb-37dd-1f23c2e9a31a
md" #### Bifurcation with Vs0"

# ╔═╡ a9f755ce-6bdf-11eb-114e-6b6938c70ad6
x_vs0, y_vs0 = meshgrid(-70:5.0:0, -70:5.0:0)

# ╔═╡ 665b7050-6bc5-11eb-393d-990a0aafbb60
begin 
	#initial = -35.0
	pl_1=change_I(1)
	pvs0_v0 = change_vs0(-40,pl_1)
	vs0_v0_V1 = zeros(size(cell_potential))
	vs0_v0_V2 = zeros(size(cell_potential))
	
	for i=1:length(cell_potential)
		vs0_v0_V1[i],vs0_v0_V2[i]=Vnullcline(cell_potential[i],pvs0_v0)
	end
	md"""Nullclines for vs0 = -40 (= v0)""" 
end

# ╔═╡ 3f1e6e10-6bdf-11eb-3518-915c7e60782f
begin
	dx_vs0_v0 = zeros(size(x_vs0))
	dy_vs0_v0 = zeros(size(x_vs0))
	for k=1:length(x_vs0)
		dx_vs0_v0[k],dy_vs0_v0[k]=return_MQIF_gradient([x_vs0[k],y_vs0[k]],pvs0_v0,0.0)	
	end
	dx_vs0_v0 = dx_vs0_v0*scale[1]
	dy_vs0_v0 = dy_vs0_v0*scale[2]
	md""" Gradient"""
end

# ╔═╡ 033caa5e-6edc-11eb-0cf7-017566eccf30
begin
	u0_vs0_1 =[-20,-40]
	u0_vs0_2 =[-60,-40]
	tspan_vs0=(0.0,300.0)
	prob_vs0_v0_1 = ODEProblem(MQIF!,u0_vs0_1,tspan_vs0,pvs0_v0,callback=cb)
	sol_vs0_v0_1 = solve(prob_vs0_v0_1,DP5(),reltol=1e-6,abstol=1e-6)
	prob_vs0_v0_2 = ODEProblem(MQIF!,u0_vs0_2,tspan_vs0,pvs0_v0,callback=cb)
	sol_vs0_v0_2 = solve(prob_vs0_v0_2,DP5(),reltol=1e-6,abstol=1e-6)
	md"""Trajectories"""
end

# ╔═╡ 51709810-6bdd-11eb-3f50-41c504469836
begin 
	#initial = -35.0
	pvs0_h = change_vs0(-35,pl_1)
	vs0_h_V1 = zeros(size(cell_potential))
	vs0_h_V2 = zeros(size(cell_potential))
	
	for i=1:length(cell_potential)
		vs0_h_V1[i],vs0_h_V2[i]=Vnullcline(cell_potential[i],pvs0_h)
	end
	md"""Nullclines for vs0 = -35""" 
end

# ╔═╡ 477972c0-6be0-11eb-1ca3-8157c59417fa
begin
	dx_vs0_h = zeros(size(x_vs0))
	dy_vs0_h = zeros(size(x_vs0))
	for k=1:length(x_vs0)
		dx_vs0_h[k],dy_vs0_h[k]=return_MQIF_gradient([x_vs0[k],y_vs0[k]],pvs0_h,0.0)	
	end
	dx_vs0_h = dx_vs0_h*scale[1]
	dy_vs0_h = dy_vs0_h*scale[2]
	md""" Gradient"""
end

# ╔═╡ 4f91f1f0-6edb-11eb-3f31-7927596e94ea
begin
	prob_vs0_h_1 = ODEProblem(MQIF!,u0_vs0_1,tspan_vs0,pvs0_h,callback=cb)
	sol_vs0_h_1 = solve(prob_vs0_h_1,DP5(),reltol=1e-6,abstol=1e-6)
	prob_vs0_h_2 = ODEProblem(MQIF!,u0_vs0_2,tspan_vs0,pvs0_h,callback=cb)
	sol_vs0_h_2 = solve(prob_vs0_h_2,DP5(),reltol=1e-6,abstol=1e-6)
	md"""Trajectories"""
end

# ╔═╡ 032b1940-6bde-11eb-0918-17fed0abd003
begin 
	#initial = -35.0
	pvs0_l = change_vs0(-45,pl_1)
	vs0_l_V1 = zeros(size(cell_potential))
	vs0_l_V2 = zeros(size(cell_potential))
	
	for i=1:length(cell_potential)
		vs0_l_V1[i],vs0_l_V2[i]=Vnullcline(cell_potential[i],pvs0_l)
	end
	md"""Nullclines for vs0 = -45""" 
end

# ╔═╡ 37137dd0-6be1-11eb-1757-ab7ac23e4d7c
begin
	dx_vs0_l = zeros(size(x_vs0))
	dy_vs0_l = zeros(size(x_vs0))
	for k=1:length(x_vs0)
		dx_vs0_l[k],dy_vs0_l[k]=return_MQIF_gradient([x_vs0[k],y_vs0[k]],pvs0_l,0.0)	
	end
	dx_vs0_l = dx_vs0_l*scale[1]
	dy_vs0_l = dy_vs0_l*scale[2]
	md""" Gradient"""
end

# ╔═╡ 84d16cf0-6edc-11eb-1078-8b5c6781e8a4
begin
	prob_vs0_l_1 = ODEProblem(MQIF!,u0_vs0_1,tspan_vs0,pvs0_l,callback=cb)
	sol_vs0_l_1 = solve(prob_vs0_l_1,DP5(),reltol=1e-6,abstol=1e-6)
	prob_vs0_l_2 = ODEProblem(MQIF!,u0_vs0_2,tspan_vs0,pvs0_l,callback=cb)
	sol_vs0_l_2 = solve(prob_vs0_l_2,DP5(),reltol=1e-6,abstol=1e-6)
	md"""Trajectories"""
end

# ╔═╡ 8279a940-6bc6-11eb-1151-c1bc4a7882e1
gr()

# ╔═╡ 327db2c0-6bc5-11eb-3258-7314d88c9039
begin
	plotvs0_h = plot(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_h_V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1,label="V-nullclines Vs0>V0")
	plot!(cell_potential,vs0_h_V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1,label="V-nullclines Vs0>V0")
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))
	
	plotvs0_v0 = plot(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_v0_V1,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="V-nullclines Vs0=V0")
	plot!(cell_potential,vs0_v0_V2,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="V-nullclines Vs0=V0")
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))
	
	plotvs0_l = plot(v_potential,vs_null,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_l_V1,linecolor=RGB(1,0.6,0.1),linewidth = 1.5,label="V-nullclines Vs0<V0")
	plot!(cell_potential,vs0_l_V2,linecolor=RGB(1,0.6,0.1),linewidth = 1.5,label="V-nullclines Vs0<V0")
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))

	plot(plotvs0_h,plotvs0_v0,plotvs0_l,layout=(1,3),legend=:outertopright,size=(1200,400))
	title!("I=1")
	
end

# ╔═╡ 31f6a870-6bc5-11eb-3afa-230cb80e952b
begin
	plot(v_potential,vs_null,linecolor=RGB(0,0.7,0),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_h_V1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1,label="V-nullclines Vs0>V0")
	plot!(cell_potential,vs0_h_V2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1,label="V-nullclines Vs0>V0")
	quiver!(x_vs0, y_vs0, quiver=(dx_vs0_h, dy_vs0_h),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	plot!(sol_vs0_h_1[1,:],sol_vs0_h_1[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 1 : cycle",linecolor=RGB(0,0,0))
	plot!(sol_vs0_h_2[1,:],sol_vs0_h_2[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 2 : rest",linecolor=RGB(0,0,0))
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))
end

# ╔═╡ 32473970-6bc5-11eb-1e25-437e376ea951
begin
	plot(v_potential,vs_null,linecolor=RGB(0,0.7,0),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_v0_V1,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="V-nullclines Vs0=V0",linestyle=:dot)
	plot!(cell_potential,vs0_v0_V2,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="V-nullclines Vs0=V0",linestyle=:dot)
	quiver!(x_vs0, y_vs0, quiver=(dx_vs0_v0, dy_vs0_v0),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	plot!(sol_vs0_v0_1[1,:],sol_vs0_v0_1[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 1 : cycle",linecolor=RGB(0,0,0))
	plot!(sol_vs0_v0_2[1,:],sol_vs0_v0_2[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 2 : cycle",linecolor=RGB(0,0,0))
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))
end

# ╔═╡ 607d2540-6be1-11eb-04d7-719327df0298
begin
	plot(v_potential,vs_null,linecolor=RGB(0,0.7,0),linewidth = 1.5,label="Vs-nullclines")
	plot!(cell_potential,vs0_l_V1,linecolor=RGB(1,0.6,0.1),linewidth = 1.5,label="V-nullclines Vs0<V0")
	plot!(cell_potential,vs0_l_V2,linecolor=RGB(1,0.6,0.1),linewidth = 1.5,label="V-nullclines Vs0<V0")
	quiver!(x_vs0, y_vs0, quiver=(dx_vs0_l, dy_vs0_l),label="test",linecolor=RGB(0,0.6,0.8),linewidth=1.5,size=(500,500),legend=:outertopright)
	plot!(sol_vs0_l_1[1,:],sol_vs0_l_1[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 1 : rest",linecolor=RGB(0,0,0))
	plot!(sol_vs0_l_2[1,:],sol_vs0_l_2[2,:],linewidth = 1.5,legend=:outertopright,label="Trajectory 2 : rest",linecolor=RGB(0,0,0))
	yaxis!("Vs",(-70,0))
	xaxis!("V",(-70,0))
end

# ╔═╡ 3c0c2f40-7225-11eb-2f84-811e39d11ccb
begin
	I_vs0_ref = 0.1
	I_vs0_l = 5
	I_vs0_h =20
	I_vs0_nl = -5
	I_vs0_nh = -20
	
	p_I_vs0_ref=change_I(I_vs0_ref)
	p_I_vs0_l=change_I(I_vs0_l)
	p_I_vs0_h=change_I(I_vs0_h)
	p_I_vs0_nl=change_I(I_vs0_nl)
	p_I_vs0_nh=change_I(I_vs0_nh)
	
	vs0_bif_ref = zeros(1,2)
	vs0_bif_l = zeros(1,2)
	vs0_bif_h = zeros(1,2)
	
	vs0_bif_ref[1],vs0_bif_ref[2]=find_bifurcation(p_I_vs0_ref)
	vs0_bif_l[1],vs0_bif_l[2]=find_bifurcation(p_I_vs0_l)
	vs0_bif_h[1],vs0_bif_h[2]=find_bifurcation(p_I_vs0_h)
	
	vs0_min = -50
	vs0_max = -30
	md"""Currents used for bifurcation diagram with vs0 and their corresponding value of bifurcation with vs0"""
end

# ╔═╡ bf34e5d0-6f94-11eb-15c6-1b081518100c
md""" We will evaluate the bifurcation diagram with vs0 as a bifurcation parameter with 

               I =  0.1 

"""

# ╔═╡ c31d95a0-729d-11eb-07c7-d9c0db5dc535
begin
	# I is  very low
	lim_vs0_range_1=vs0_bif_ref[1]-abs(vs0_bif_ref[1]/1000)
	vs0_ref_range_1 = collect(range(vs0_min,stop=lim_vs0_range_1,length=50))
	
	
	lim_vs0_range_2_1=vs0_bif_ref[1]-abs(vs0_bif_ref[1]/1000)
	lim_vs0_range_2_2=vs0_bif_ref[2]+abs(vs0_bif_ref[2]/1000)
	vs0_ref_range_2 = collect(range(lim_vs0_range_2_1,stop=lim_vs0_range_2_2,length=30))
	
	lim_vs0_range_3=vs0_bif_ref[2]+abs(vs0_bif_ref[2]/1000)
	vs0_ref_range_3 = collect(range(lim_vs0_range_3,stop=vs0_max,length=50))
	
						vs0_ref_range=zeros(length(vs0_ref_range_1)+length(vs0_ref_range_2)+length(vs0_ref_range_3),1)
	
	for i=1:length(vs0_ref_range)
		if i<=length(vs0_ref_range_1)
			vs0_ref_range[i] = vs0_ref_range_1[i]
		else
			if i<=length(vs0_ref_range_1)+length(vs0_ref_range_2)
				vs0_ref_range[i] = vs0_ref_range_2[i-length(vs0_ref_range_1)]
			else
				vs0_ref_range[i] = vs0_ref_range_3[i-length(vs0_ref_range_1)-length(vs0_ref_range_2)]
			end
		end
	end
	
	md"""Range of vs0 values separated in 3 parts according to the bifurcations"""
end

# ╔═╡ 49d451f0-729f-11eb-3996-054cabdfbf5a
begin
	list_stable_vs0_ref = zeros(size(vs0_ref_range))
	list_saddle_vs0_ref = zeros(size(vs0_ref_range))
	list_unstable_vs0_ref = zeros(size(vs0_ref_range))
	list_cycleV_max_vs0_ref = zeros(size(vs0_ref_range))
	list_cycleV_min_vs0_ref = zeros(size(vs0_ref_range))
	list_cycleVs_vs0_ref = zeros(length(vs0_ref_range),2)
	p_used_vs0_ref = zeros(length(vs0_ref_range),length(p))
	
	tspan_vs0_ref=(0.0,150.0)
	u0_vs0_ref =[-20,-40]
	
	for i=1:length(vs0_ref_range)
		vs0_ref_ = vs0_ref_range[i]
		p_I_vs0_ref_= change_vs0(vs0_ref_,p_I_vs0_ref)
		p_used_vs0_ref[i,:]=p_I_vs0_ref_
		vs0_ref_fp1,vs0_ref_fp2,vs0_ref_stability = fixedpoints(p_I_vs0_ref_)
		
		if length(vs0_ref_stability) >1
			vs0_ref_fp=[vs0_ref_fp1,vs0_ref_fp2]
			for k=1:length(vs0_ref_stability)
				if vs0_ref_stability[k] < 0
					list_stable_vs0_ref[i] = vs0_ref_fp[k][1]
					list_unstable_vs0_ref[i] = NaN
				end
				if vs0_ref_stability[k] == 0
					list_saddle_vs0_ref[i] = vs0_ref_fp[k][1]
				end
				if vs0_ref_stability[k] > 0
					list_unstable_vs0_ref[i] = vs0_ref_fp[k][1]
					list_stable_vs0_ref[i] = NaN
				end
			end
		else 
			list_stable_vs0_ref[i] = NaN
			list_saddle_vs0_ref[i] = NaN
			list_unstable_vs0_ref[i] = NaN
		end
		
			
		prob_vs0_ref = ODEProblem(MQIF!,u0_vs0_ref,tspan_vs0_ref,p_I_vs0_ref_,callback=cb)
			
		sol_vs0_ref = solve(prob_vs0_ref,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_vs0_ref = findfirst(x -> x>=maximum(tspan_vs0_ref)/2, sol_vs0_ref.t)
			
		sol_vs0_ref_max_V = maximum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_V = minimum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_max_Vs = maximum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_Vs = minimum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
			
		if sol_vs0_ref_max_V-sol_vs0_ref_min_V>20
			list_cycleV_max_vs0_ref[i]=sol_vs0_ref_max_V
			list_cycleV_min_vs0_ref[i]=sol_vs0_ref_min_V
			list_cycleVs_vs0_ref[i,1]=sol_vs0_ref_max_Vs
			list_cycleVs_vs0_ref[i,2]=sol_vs0_ref_min_Vs
		else 
			list_cycleV_max_vs0_ref[i]=NaN
			list_cycleV_min_vs0_ref[i]=NaN
			list_cycleVs_vs0_ref[i,1]=NaN
			list_cycleVs_vs0_ref[i,2]=NaN
		end
	end
	
	
	md""" Calculus of stable points, saddle points, maximum and minimum values of the limit cycle for each I"""
end

# ╔═╡ bfdab820-72c8-11eb-0d57-0d8988ef5aa7
begin	
	bif_vs0_ref = plot(vs0_ref_range,list_stable_vs0_ref,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(vs0_ref_range,list_saddle_vs0_ref,linecolor="orange",linewidth=1.5,label="Saddle node")
	plot!(vs0_ref_range,list_cycleV_max_vs0_ref,linecolor="pink",linewidth=2.5,label="Maximum of the cycle",linestyle=:dashdot)
	plot!(vs0_ref_range,list_cycleV_min_vs0_ref,linecolor="pink",linewidth=2.5,label="Minimum of the cycle",linestyle=:dashdot)
	xaxis!("Vs0")
	yaxis!("V")
	title!("Bifurcation diagram, I = $I_vs0_ref")
end

# ╔═╡ a01bb0a0-745a-11eb-1554-2b2749deee03
begin
	ind_vs0_v0_approx_ref = findfirst(x ->abs(x-p[2])<0.02, vs0_ref_range)
	ind_vs0_h_approx_ref = findfirst(x ->abs(x-(-35))<0.5, vs0_ref_range)
	ind_vs0_l_approx_ref = findfirst(x ->abs(x-(-45))<0.5, vs0_ref_range)
	
	vs0_v0_approx_ref=vs0_ref_range[ind_vs0_v0_approx_ref[1]]
	vs0_h_approx_ref=vs0_ref_range[ind_vs0_h_approx_ref[1]]
	vs0_l_approx_ref=vs0_ref_range[ind_vs0_l_approx_ref[1]]
	
	vnull_vs0_v0_ref1=zeros(size(v_potential))
	vnull_vs0_v0_ref2=zeros(size(v_potential))
	vnull_vs0_h_ref1=zeros(size(v_potential))
	vnull_vs0_h_ref2=zeros(size(v_potential))
	vnull_vs0_l_ref1=zeros(size(v_potential))
	vnull_vs0_l_ref2=zeros(size(v_potential))
	
	for i=1:length(v_potential)
		vnull_vs0_v0_ref1[i],vnull_vs0_v0_ref2[i] = Vnullcline(v_potential[i],p_used_vs0_ref[ind_vs0_v0_approx_ref[1],:])
		vnull_vs0_h_ref1[i],vnull_vs0_h_ref2[i] = Vnullcline(v_potential[i],p_used_vs0_ref[ind_vs0_h_approx_ref[1],:])
		vnull_vs0_l_ref1[i],vnull_vs0_l_ref2[i] = Vnullcline(v_potential[i],p_used_vs0_ref[ind_vs0_l_approx_ref[1],:])
	end
	
	phase_vs0_v0_ref = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_v0_ref1,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_v0_ref2,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_v0_approx_ref))")
	
	phase_vs0_h_ref = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_h_ref1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_h_ref2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_ref[ind_vs0_h_approx_ref]],[list_stable_vs0_ref[ind_vs0_h_approx_ref]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_ref[ind_vs0_h_approx_ref]],[list_saddle_vs0_ref[ind_vs0_h_approx_ref]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_h_approx_ref))")
	
	phase_vs0_l_ref = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_l_ref1,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_l_ref2,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_ref[ind_vs0_l_approx_ref]],[list_stable_vs0_ref[ind_vs0_l_approx_ref]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_ref[ind_vs0_l_approx_ref]],[list_saddle_vs0_ref[ind_vs0_l_approx_ref]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_l_approx_ref))")
	
	phase_vs0_ref = plot(phase_vs0_l_ref,phase_vs0_v0_ref,phase_vs0_h_ref,layout=(3,1),linewidth = 1.5)
	
	
	plot(bif_vs0_ref,phase_vs0_ref,layout=(1,2))
	#savefig("MQIF_2D_bif_vs0_lowI.pdf")
end

# ╔═╡ ef10370e-74fc-11eb-0893-23b71d91939d
#savefig("MQIF_2D_bif_vs0_lowI.pdf")

# ╔═╡ c3271fe0-6f9a-11eb-310f-5b846b6f646f
md"""

      		   I = 5    
 """

# ╔═╡ a7b2f5a0-72e6-11eb-2f6b-5f6177cce876
begin
	# I is low
	lim_vs0_l_range_1=vs0_bif_l[1]-abs(vs0_bif_l[1]/1000)
	vs0_l_range_1 = collect(range(vs0_min,stop=lim_vs0_l_range_1,length=100))
	
	
	lim_vs0_l_range_2_1=vs0_bif_l[1]-abs(vs0_bif_l[1]/1000)
	lim_vs0_l_range_2_2=vs0_bif_l[2]+abs(vs0_bif_l[2]/1000)
	vs0_l_range_2 = collect(range(lim_vs0_l_range_2_1,stop=lim_vs0_l_range_2_2,length=30))
	
	lim_vs0_l_range_3=vs0_bif_l[2]+abs(vs0_bif_l[2]/1000)
	vs0_l_range_3 = collect(range(lim_vs0_l_range_3,stop=vs0_max,length=50))
	
	vs0_l_range=zeros(length(vs0_l_range_1)+length(vs0_l_range_2)+length(vs0_l_range_3),1)
	
	for i=1:length(vs0_l_range)
		if i<=length(vs0_l_range_1)
			vs0_l_range[i] = vs0_l_range_1[i]
		else
			if i<=length(vs0_l_range_1)+length(vs0_l_range_2)
				vs0_l_range[i] = vs0_l_range_2[i-length(vs0_l_range_1)]
			else
				vs0_l_range[i] = vs0_l_range_3[i-length(vs0_l_range_1)-length(vs0_l_range_2)]
			end
		end
	end
	md"""Range of vs0 values separated in 3 parts according to the bifurcations"""
end

# ╔═╡ 287ffb62-72e7-11eb-3460-6bc59ea30e12
begin
	list_stable_vs0_l = zeros(size(vs0_l_range))
	list_saddle_vs0_l = zeros(size(vs0_l_range))
	list_unstable_vs0_l = zeros(size(vs0_l_range))
	list_cycleV_max_vs0_l = zeros(size(vs0_l_range))
	list_cycleV_min_vs0_l = zeros(size(vs0_l_range))
	list_cycleVs_vs0_l = zeros(length(vs0_l_range),2)
	p_used_vs0_l = zeros(length(vs0_l_range),length(p))
	
	tspan_vs0_l=(0.0,150.0)
	u0_vs0_l =[-20,-40]
	
	for i=1:length(vs0_l_range)
		vs0_ref_ = vs0_l_range[i]
		p_I_vs0_ref_= change_vs0(vs0_ref_,p_I_vs0_l)
		p_used_vs0_l[i,:]=p_I_vs0_ref_
		vs0_ref_fp1,vs0_ref_fp2,vs0_ref_stability = fixedpoints(p_I_vs0_ref_)
		
		if length(vs0_ref_stability) >1
			vs0_ref_fp=[vs0_ref_fp1,vs0_ref_fp2]
			for k=1:length(vs0_ref_stability)
				if vs0_ref_stability[k] < 0
					list_stable_vs0_l[i] = vs0_ref_fp[k][1]
					list_unstable_vs0_l[i] = NaN
				end
				if vs0_ref_stability[k] == 0
					list_saddle_vs0_l[i] = vs0_ref_fp[k][1]
				end
				if vs0_ref_stability[k] > 0
					list_unstable_vs0_l[i] = vs0_ref_fp[k][1]
					list_stable_vs0_l[i] = NaN
				end
			end
		else 
			list_stable_vs0_l[i] = NaN
			list_saddle_vs0_l[i] = NaN
			list_unstable_vs0_l[i] = NaN
		end
		
			
		prob_vs0_ref = ODEProblem(MQIF!,u0_vs0_l,tspan_vs0_l,p_I_vs0_ref_,callback=cb)
			
		sol_vs0_ref = solve(prob_vs0_ref,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_vs0_ref = findfirst(x -> x>=maximum(tspan_vs0_ref)/2, sol_vs0_ref.t)
			
		sol_vs0_ref_max_V = maximum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_V = minimum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_max_Vs = maximum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_Vs = minimum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
			
		if sol_vs0_ref_max_V-sol_vs0_ref_min_V>20
			list_cycleV_max_vs0_l[i]=sol_vs0_ref_max_V
			list_cycleV_min_vs0_l[i]=sol_vs0_ref_min_V
			list_cycleVs_vs0_l[i,1]=sol_vs0_ref_max_Vs
			list_cycleVs_vs0_l[i,2]=sol_vs0_ref_min_Vs
		else 
			list_cycleV_max_vs0_l[i]=NaN
			list_cycleV_min_vs0_l[i]=NaN
			list_cycleVs_vs0_l[i,1]=NaN
			list_cycleVs_vs0_l[i,2]=NaN
		end
	end
	
	
	md""" Calculus of stable points, saddle points, maximum and minimum values of the limit cycle for each I"""
end

# ╔═╡ f510c1a0-72e7-11eb-148e-e1bebb6e8c44
begin
	bif_vs0_l= plot(vs0_l_range,list_stable_vs0_l,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(vs0_l_range,list_unstable_vs0_l,linecolor="red",linewidth=1.5,label="Unstable node")
	plot!(vs0_l_range,list_saddle_vs0_l,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	plot!(vs0_l_range,list_cycleV_max_vs0_l,linecolor="pink",linewidth=1.5,label="Maximim of the cycle")
	plot!(vs0_l_range,list_cycleV_min_vs0_l,linecolor="pink",linewidth=1.5,label="Minimum of the cycle")
	xaxis!("Vs0")
	yaxis!("V")
	title!("Bifurcation diagram, I = $I_vs0_l")
end

# ╔═╡ 012613fe-74fe-11eb-2e19-ddee86dd3c7f
begin
	ind_vs0_v0_approx_l = findfirst(x ->abs(x-p[2])<0.1, vs0_l_range)
	ind_vs0_h_approx_l = findfirst(x ->abs(x-(-35))<0.5, vs0_l_range)
	ind_vs0_l_approx_l = findfirst(x ->abs(x-(-45))<0.5, vs0_l_range)
	
	vs0_v0_approx_l=vs0_l_range[ind_vs0_v0_approx_l[1]]
	vs0_h_approx_l=vs0_l_range[ind_vs0_h_approx_l[1]]
	vs0_l_approx_l=vs0_l_range[ind_vs0_l_approx_l[1]]
	
	vnull_vs0_v0_l1=zeros(size(v_potential))
	vnull_vs0_v0_l2=zeros(size(v_potential))
	vnull_vs0_h_l1=zeros(size(v_potential))
	vnull_vs0_h_l2=zeros(size(v_potential))
	vnull_vs0_l_l1=zeros(size(v_potential))
	vnull_vs0_l_l2=zeros(size(v_potential))
	
	for i=1:length(v_potential)
		vnull_vs0_v0_l1[i],vnull_vs0_v0_l2[i] = Vnullcline(v_potential[i],p_used_vs0_l[ind_vs0_v0_approx_l[1],:])
		vnull_vs0_h_l1[i],vnull_vs0_h_l2[i] = Vnullcline(v_potential[i],p_used_vs0_l[ind_vs0_h_approx_l[1],:])
		vnull_vs0_l_l1[i],vnull_vs0_l_l2[i] = Vnullcline(v_potential[i],p_used_vs0_l[ind_vs0_l_approx_l[1],:])
	end
	
	phase_vs0_v0_l = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_v0_l1,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_v0_l2,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_v0_approx_ref))")
	
	phase_vs0_h_l = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_h_l1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_h_l2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_l[ind_vs0_h_approx_l]],[list_stable_vs0_ref[ind_vs0_h_approx_ref]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_l[ind_vs0_h_approx_l]],[list_saddle_vs0_l[ind_vs0_h_approx_l]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_h_approx_ref))")
	
	phase_vs0_l_l = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_l_l1,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_l_l2,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_l[ind_vs0_l_approx_l]],[list_stable_vs0_l[ind_vs0_l_approx_l]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_l[ind_vs0_l_approx_l]],[list_saddle_vs0_l[ind_vs0_l_approx_l]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_l_approx_l))")
	
	phase_vs0_l = plot(phase_vs0_l_l,phase_vs0_v0_l,phase_vs0_h_l,layout=(3,1),linewidth = 1.5)
	
	
	plot(bif_vs0_l,phase_vs0_l,layout=(1,2))
	#savefig("MQIF_2D_bif_vs0_middleI.pdf")
end

# ╔═╡ 36c49b72-8504-11eb-2652-bd34549800bd
#savefig("MQIF_2D_bif_vs0_middleI.pdf")

# ╔═╡ 6720a03e-7053-11eb-3795-9fd39a387359
md"""  
         	   I = 40       
"""

# ╔═╡ f7158060-72e9-11eb-24d6-d52b86bf3848
begin
	# I is high
	lim_vs0_h_range_1=vs0_bif_h[1]-abs(vs0_bif_h[1]/1000)
	vs0_h_range_1 = collect(range(vs0_min,stop=lim_vs0_h_range_1,length=100))
	
	
	lim_vs0_h_range_2_1=vs0_bif_h[1]-abs(vs0_bif_h[1]/1000)
	lim_vs0_h_range_2_2=vs0_bif_h[2]+abs(vs0_bif_h[2]/1000)
	vs0_h_range_2 = collect(range(lim_vs0_h_range_2_1,stop=lim_vs0_h_range_2_2,length=30))
	
	lim_vs0_h_range_3=vs0_bif_h[2]+abs(vs0_bif_h[2]/1000)
	vs0_h_range_3 = collect(range(lim_vs0_h_range_3,stop=vs0_max,length=50))
	
	vs0_h_range=zeros(length(vs0_h_range_1)+length(vs0_h_range_2)+length(vs0_h_range_3),1)
	
	for i=1:length(vs0_h_range)
		if i<=length(vs0_h_range_1)
			vs0_h_range[i] = vs0_h_range_1[i]
		else
			if i<=length(vs0_h_range_1)+length(vs0_h_range_2)
				vs0_h_range[i] = vs0_h_range_2[i-length(vs0_h_range_1)]
			else
				vs0_h_range[i] = vs0_h_range_3[i-length(vs0_h_range_1)-length(vs0_h_range_2)]
			end
		end
	end
	md"""Range of vs0 values separated in 3 parts according to the bifurcations"""
end

# ╔═╡ 5f5323d0-72ea-11eb-392e-4d095f09c02c
begin
	list_stable_vs0_h = zeros(size(vs0_h_range))
	list_saddle_vs0_h = zeros(size(vs0_h_range))
	list_unstable_vs0_h = zeros(size(vs0_h_range))
	list_cycleV_max_vs0_h = zeros(size(vs0_h_range))
	list_cycleV_min_vs0_h = zeros(size(vs0_h_range))
	list_cycleVs_vs0_h = zeros(length(vs0_h_range),2)
	p_used_vs0_h = zeros(length(vs0_h_range),length(p))
	
	tspan_vs0_h=(0.0,150.0)
	u0_vs0_h =[-20,-40]
	
	for i=1:length(vs0_h_range)
		vs0_ref_ = vs0_h_range[i]
		p_I_vs0_ref_= change_vs0(vs0_ref_,p_I_vs0_h)
		p_used_vs0_h[i,:]=p_I_vs0_ref_
		vs0_ref_fp1,vs0_ref_fp2,vs0_ref_stability = fixedpoints(p_I_vs0_ref_)
		
		if length(vs0_ref_stability) >1
			vs0_ref_fp=[vs0_ref_fp1,vs0_ref_fp2]
			for k=1:length(vs0_ref_stability)
				if vs0_ref_stability[k] < 0
					list_stable_vs0_h[i] = vs0_ref_fp[k][1]
					list_unstable_vs0_h[i] = NaN
				end
				if vs0_ref_stability[k] == 0
					list_saddle_vs0_h[i] = vs0_ref_fp[k][1]
				end
				if vs0_ref_stability[k] > 0
					list_unstable_vs0_h[i] = vs0_ref_fp[k][1]
					list_stable_vs0_h[i] = NaN
				end
			end
		else 
			list_stable_vs0_h[i] = NaN
			list_saddle_vs0_h[i] = NaN
			list_unstable_vs0_h[i] = NaN
		end
		
			
		prob_vs0_ref = ODEProblem(MQIF!,u0_vs0_h,tspan_vs0_h,p_I_vs0_ref_,callback=cb)
			
		sol_vs0_ref = solve(prob_vs0_ref,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin_vs0_ref = findfirst(x -> x>=maximum(tspan_vs0_ref)/2, sol_vs0_ref.t)
			
		sol_vs0_ref_max_V = maximum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_V = minimum(sol_vs0_ref[1,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_max_Vs = maximum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
		sol_vs0_ref_min_Vs = minimum(sol_vs0_ref[2,ind_tmin_vs0_ref:length(sol_vs0_ref)])
			
		if sol_vs0_ref_max_V-sol_vs0_ref_min_V>20
			list_cycleV_max_vs0_h[i]=sol_vs0_ref_max_V
			list_cycleV_min_vs0_h[i]=sol_vs0_ref_min_V
			list_cycleVs_vs0_h[i,1]=sol_vs0_ref_max_Vs
			list_cycleVs_vs0_h[i,2]=sol_vs0_ref_min_Vs
		else 
			list_cycleV_max_vs0_h[i]=NaN
			list_cycleV_min_vs0_h[i]=NaN
			list_cycleVs_vs0_h[i,1]=NaN
			list_cycleVs_vs0_h[i,2]=NaN
		end
	end
	
	
	md""" Calculus of stable points, saddle points, maximum and minimum values of the limit cycle for each I"""
end

# ╔═╡ f6b1c010-72ea-11eb-02ce-07d898496dcd
begin
	bif_vs0_h = plot(vs0_h_range,list_stable_vs0_h,linecolor="green",linewidth=1.5,label="Stable node")
	plot!(vs0_h_range,list_unstable_vs0_h,linecolor="red",linewidth=1.5,label="Unstable node")
	plot!(vs0_h_range,list_saddle_vs0_h,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	plot!(vs0_h_range,list_cycleV_max_vs0_h,linecolor="pink",linewidth=1.5,label="Maximum of the cycle",linestyle=:dot)
	plot!(vs0_h_range,list_cycleV_min_vs0_h,linecolor="pink",linewidth=1.5,label="Minimum of the cycle",linestyle=:dot)
	xaxis!("Vs0")
	yaxis!("V")
	title!("Bifurcation diagram, I = $I_vs0_h")
end

# ╔═╡ 4ced1360-74ff-11eb-29e7-bde7e5380971
begin
	ind_vs0_v0_approx_h = findfirst(x ->abs(x-p[2])<0.25, vs0_h_range)
	ind_vs0_h_approx_h = findfirst(x ->abs(x-(-35))<0.5, vs0_h_range)
	ind_vs0_l_approx_h = findfirst(x ->abs(x-(-45))<0.5, vs0_h_range)
	
	vs0_v0_approx_h=vs0_h_range[ind_vs0_v0_approx_h[1]]
	vs0_h_approx_h=vs0_h_range[ind_vs0_h_approx_h[1]]
	vs0_l_approx_h=vs0_h_range[ind_vs0_l_approx_h[1]]
	
	vnull_vs0_v0_h1=zeros(size(v_potential))
	vnull_vs0_v0_h2=zeros(size(v_potential))
	vnull_vs0_h_h1=zeros(size(v_potential))
	vnull_vs0_h_h2=zeros(size(v_potential))
	vnull_vs0_l_h1=zeros(size(v_potential))
	vnull_vs0_l_h2=zeros(size(v_potential))
	
	for i=1:length(v_potential)
		vnull_vs0_v0_h1[i],vnull_vs0_v0_h2[i] = Vnullcline(v_potential[i],p_used_vs0_h[ind_vs0_v0_approx_h[1],:])
		vnull_vs0_h_h1[i],vnull_vs0_h_h2[i] = Vnullcline(v_potential[i],p_used_vs0_h[ind_vs0_h_approx_h[1],:])
		vnull_vs0_l_h1[i],vnull_vs0_l_h2[i] = Vnullcline(v_potential[i],p_used_vs0_h[ind_vs0_l_approx_h[1],:])
	end
	
	phase_vs0_v0_h = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_v0_h1,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_v0_h2,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_v0_approx_h))")
	
	phase_vs0_h_h = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_h_h1,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_h_h2,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_h[ind_vs0_h_approx_h]],[list_stable_vs0_h[ind_vs0_h_approx_h]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_h[ind_vs0_h_approx_h]],[list_saddle_vs0_h[ind_vs0_h_approx_h]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_h_approx_h))")
	
	phase_vs0_l_h = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(v_potential,vnull_vs0_l_h1,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	plot!(v_potential,vnull_vs0_l_h2,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_unstable_vs0_h[ind_vs0_l_approx_h]],[list_unstable_vs0_h[ind_vs0_l_approx_h]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")
	scatter!([list_saddle_vs0_h[ind_vs0_l_approx_h]],[list_saddle_vs0_h[ind_vs0_l_approx_h]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_l_approx_h))")
	
	phase_vs0_h = plot(phase_vs0_l_h,phase_vs0_v0_h,phase_vs0_h_h,layout=(3,1),linewidth = 1.5)
	
	
	plot(bif_vs0_h,phase_vs0_h,layout=(1,2))
	#savefig("MQIF_2D_bif_vs0_highI.pdf")
end

# ╔═╡ 94822410-72eb-11eb-1b05-eb82cb443e0c
md"""  
         	   I = -5       
"""

# ╔═╡ a4d8ba40-72eb-11eb-32ac-e7fd18d0fe2e
begin
	# I is low and negative
	vs0_nl_range = collect(range(vs0_min,stop=vs0_max,length=100))
	md"""Range of vs0 values separated in 3 parts according to the bifurcations"""
end

# ╔═╡ a46d272e-72eb-11eb-13e5-4f15f3de10c3
begin
	list_stable_vs0_nl =zeros(size(vs0_nl_range))
	list_unstable_vs0_nl =zeros(size(vs0_nl_range))
	list_saddle_vs0_nl =zeros(size(vs0_nl_range))
	list_cycleV_vs0_nl =zeros(length(vs0_nl_range),2)
	list_cycleVs_vs0_nl =zeros(length(vs0_nl_range),2)
	p_used_vs0_nl = zeros(length(vs0_nl_range),length(p))
	
	tspan_vs0_nl=(0.0,150.0)
	u0_vs0_nl =[-20,-40]
	
		
		
		
		for j=1:length(vs0_nl_range)
			vs0_nl = vs0_nl_range[j]
			p_I_vs0_nl_= change_vs0(vs0_nl,p_I_vs0_nl)
			p_used_vs0_nl[j,:]=p_I_vs0_nl_
			vs0_nl_fp1,vs0_nl_fp2,vs0_nl_stability = fixedpoints(p_I_vs0_nl_)
			if length(vs0_nl_stability) >1
				vs0_nl_fp=[vs0_nl_fp1,vs0_nl_fp2]
				for k=1:length(vs0_nl_stability)
					if vs0_nl_stability[k] < 0
						list_stable_vs0_nl[j] = vs0_nl_fp[k][1]
						list_unstable_vs0_nl[j] = NaN
					end
					if vs0_nl_stability[k] == 0
						list_saddle_vs0_nl[j] = vs0_nl_fp[k][1]
					end
					if vs0_nl_stability[k] > 0
						list_unstable_vs0_nl[j] = vs0_nl_fp[k][1]
						list_stable_vs0_nl[j] = NaN
					end
				end
			else 
				list_stable_vs0_nl[j] = NaN
				list_saddle_vs0_nl[j] = NaN
				list_unstable_vs0_nl[j] = NaN
			end
			
			prob_vs0_nl = ODEProblem(MQIF!,u0_vs0_nl,tspan_vs0_nl,p_I_vs0_nl_,callback=cb)
			
			sol_vs0_nl = solve(prob_vs0_nl,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
			ind_tmin_vs0_nl = findfirst(x -> x>=maximum(tspan_vs0_nl)/2, sol_vs0_nl.t)
			
			sol_vs0_nl_max_V = maximum(sol_vs0_nl[1,ind_tmin_vs0_nl:length(sol_vs0_nl)])
			sol_vs0_nl_min_V = minimum(sol_vs0_nl[1,ind_tmin_vs0_nl:length(sol_vs0_nl)])
			sol_vs0_nl_max_Vs = maximum(sol_vs0_nl[2,ind_tmin_vs0_nl:length(sol_vs0_nl)])
			sol_vs0_nl_min_Vs = minimum(sol_vs0_nl[2,ind_tmin_vs0_nl:length(sol_vs0_nl)])
			
			if sol_vs0_nl_max_V-sol_vs0_nl_min_V>30
				list_cycleV_vs0_nl[j,1]=sol_vs0_nl_max_V
				list_cycleV_vs0_nl[j,2]=sol_vs0_nl_min_V
				list_cycleVs_vs0_nl[j,1]=sol_vs0_nl_max_Vs
				list_cycleVs_vs0_nl[j,2]=sol_vs0_nl_min_Vs
			else 
				list_cycleV_vs0_nl[j,1]=NaN
				list_cycleV_vs0_nl[j,2]=NaN
				list_cycleVs_vs0_nl[j,1]=NaN
				list_cycleVs_vs0_nl[j,2]=NaN
			end
		end
	
	
	md""" Calculus of stable points, saddle points, maximum and minimum values of the limit cycle for each I"""
end

# ╔═╡ 0dce5780-72f6-11eb-3d8d-63ad82ef0ca6
begin
	bif_vs0_nl = plot(vs0_nl_range,list_stable_vs0_nl,linecolor="green",linewidth=1.5,label="Stable node")
	real_value_unstable_vs0_nl = findall(x -> x>=0 || x<0, list_unstable_vs0_nl)
	if length(real_value_unstable_vs0_nl)>0
		plot!(vs0_nl_range,list_unstable_vs0_nl[:],linecolor="red",linewidth=1.5,label="Untable node")
	end
	real_value_cycleV_vs0_nl = findall(x -> x>=0 || x<0, list_cycleV_vs0_nl)
	if length(real_value_cycleV_vs0_nl)>0
		plot!(vs0_nl_range,list_cycleV_vs0_nl[:,1],linecolor="pink",linewidth=1.5,label="Maximum of the cycle")
		plot!(vs0_nl_range,list_cycleV_vs0_nl[:,2],linecolor="pink",linewidth=1.5,label="Minimum of the cycle")
	end
	plot!(vs0_nl_range,list_saddle_vs0_nl,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	xaxis!("Vs0")
	yaxis!("V")
	title!("Bifurcation diagram, I = $I_vs0_nl")
end

# ╔═╡ 37ca0b90-7500-11eb-0393-a7cb2ce63a08
begin
	ind_vs0_v0_approx_nl = findfirst(x ->abs(x-p[2])<0.35, vs0_nl_range)
	ind_vs0_h_approx_nl = findfirst(x ->abs(x-(-35))<0.5, vs0_nl_range)
	ind_vs0_l_approx_nl = findfirst(x ->abs(x-(-45))<0.5, vs0_nl_range)
	
	vs0_v0_approx_nl=vs0_nl_range[ind_vs0_v0_approx_nl[1]]
	vs0_h_approx_nl=vs0_nl_range[ind_vs0_h_approx_nl[1]]
	vs0_l_approx_nl=vs0_nl_range[ind_vs0_l_approx_nl[1]]
	
	vnull_vs0_v0_nl1=zeros(size(v_potential))
	vnull_vs0_v0_nl2=zeros(size(v_potential))
	vnull_vs0_h_nl1=zeros(size(v_potential))
	vnull_vs0_h_nl2=zeros(size(v_potential))
	vnull_vs0_l_nl1=zeros(size(v_potential))
	vnull_vs0_l_nl2=zeros(size(v_potential))
	
	for i=1:length(v_potential)
		vnull_vs0_v0_nl1[i],vnull_vs0_v0_nl2[i] = Vnullcline_(v_potential[i],p_used_vs0_nl[ind_vs0_v0_approx_nl[1],:])
		vnull_vs0_h_nl1[i],vnull_vs0_h_nl2[i] = Vnullcline_(v_potential[i],p_used_vs0_nl[ind_vs0_h_approx_nl[1],:])
		vnull_vs0_l_nl1[i],vnull_vs0_l_nl2[i] = Vnullcline_(v_potential[i],p_used_vs0_nl[ind_vs0_l_approx_nl[1],:])
	end
	
	phase_vs0_v0_nl = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_v0_nl1,v_potential,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_v0_nl2,v_potential,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0",legend=:outertopright,size=(750,500))
	scatter!([list_stable_vs0_nl[ind_vs0_v0_approx_nl]],[list_stable_vs0_nl[ind_vs0_v0_approx_nl]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nl[ind_vs0_v0_approx_nl]],[list_saddle_vs0_nl[ind_vs0_v0_approx_nl]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_v0_approx_nl))")
	
	phase_vs0_h_nl = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_h_nl1,v_potential,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_h_nl2,v_potential,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_nl[ind_vs0_h_approx_nl]],[list_stable_vs0_nl[ind_vs0_h_approx_nl]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nl[ind_vs0_h_approx_nl]],[list_saddle_vs0_nl[ind_vs0_h_approx_nl]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_h_approx_ref))")
	
	phase_vs0_l_nl = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_l_nl1,v_potential,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_l_nl2,v_potential,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_nl[ind_vs0_l_approx_nl]],[list_stable_vs0_nl[ind_vs0_l_approx_nl]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nl[ind_vs0_l_approx_nl]],[list_saddle_vs0_nl[ind_vs0_l_approx_nl]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_l_approx_ref))")
	
	phase_vs0_nl = plot(phase_vs0_l_nl,phase_vs0_v0_nl,phase_vs0_h_nl,layout=(3,1),linewidth = 1.5)
	
	
	plot(bif_vs0_nl,phase_vs0_nl,layout=(1,2))
end

# ╔═╡ 9ab3123e-72eb-11eb-08f4-77f036b5e088
md"""  
         	   I = -40       
"""

# ╔═╡ 8ccb3bb0-72f7-11eb-2260-857b516bae2d
begin
	# I is high and negative
	vs0_nh_range = collect(range(vs0_min,stop=vs0_max,length=100))
	md"""Range of vs0 values separated in 3 parts according to the bifurcations"""
end

# ╔═╡ a611b8fe-72f7-11eb-2386-9fb53a251d7c
begin
	list_stable_vs0_nh =zeros(size(vs0_nh_range))
	list_unstable_vs0_nh =zeros(size(vs0_nh_range))
	list_saddle_vs0_nh =zeros(size(vs0_nh_range))
	list_cycleV_vs0_nh =zeros(length(vs0_nh_range),2)
	list_cycleVs_vs0_nh =zeros(length(vs0_nh_range),2)
	p_used_vs0_nh = zeros(length(vs0_nh_range),length(p))
	
	tspan_vs0_nh=(0.0,150.0)
	u0_vs0_nh =[-20,-40]
		
		for j=1:length(vs0_nh_range)
			vs0_nh = vs0_nh_range[j]
			p_I_vs0_nh_= change_vs0(vs0_nh,p_I_vs0_nh)
			p_used_vs0_nh[j,:]=p_I_vs0_nh_
			vs0_nh_fp1,vs0_nh_fp2,vs0_nh_stability = fixedpoints(p_I_vs0_nh_)
			if length(vs0_nh_stability) >1
				vs0_nh_fp=[vs0_nh_fp1,vs0_nh_fp2]
				for k=1:length(vs0_nh_stability)
					if vs0_nh_stability[k] < 0
						list_stable_vs0_nh[j] = vs0_nh_fp[k][1]
						list_unstable_vs0_nh[j] = NaN
					end
					if vs0_nh_stability[k] == 0
						list_saddle_vs0_nh[j] = vs0_nh_fp[k][1]
					end
					if vs0_nh_stability[k] > 0
						list_unstable_vs0_nh[j] = vs0_nh_fp[k][1]
						list_stable_vs0_nh[j] = NaN
					end
				end
			else 
				list_stable_vs0_nh[j] = NaN
				list_saddle_vs0_nh[j] = NaN
				list_unstable_vs0_nh[j] = NaN
			end
			
			prob_vs0_nh = ODEProblem(MQIF!,u0_vs0_nh,tspan_vs0_nh,p_I_vs0_nh_,callback=cb)
			
			sol_vs0_nh = solve(prob_vs0_nh,dtmax=0.01,DP5(),reltol=1e-6,abstol=1e-6)
			ind_tmin_vs0_nh = findfirst(x -> x>=maximum(tspan_vs0_nh)/2, sol_vs0_nh.t)
			
			sol_vs0_nh_max_V = maximum(sol_vs0_nh[1,ind_tmin_vs0_nh:length(sol_vs0_nh)])
			sol_vs0_nh_min_V = minimum(sol_vs0_nh[1,ind_tmin_vs0_nh:length(sol_vs0_nh)])
			sol_vs0_nh_max_Vs = maximum(sol_vs0_nh[2,ind_tmin_vs0_nh:length(sol_vs0_nh)])
			sol_vs0_nh_min_Vs = minimum(sol_vs0_nh[2,ind_tmin_vs0_nh:length(sol_vs0_nh)])
			
			if sol_vs0_nh_max_V-sol_vs0_nh_min_V>30
				list_cycleV_vs0_nh[j,1]=sol_vs0_nh_max_V
				list_cycleV_vs0_nh[j,2]=sol_vs0_nh_min_V
				list_cycleVs_vs0_nh[j,1]=sol_vs0_nh_max_Vs
				list_cycleVs_vs0_nh[j,2]=sol_vs0_nh_min_Vs
			else 
				list_cycleV_vs0_nh[j,1]=NaN
				list_cycleV_vs0_nh[j,2]=NaN
				list_cycleVs_vs0_nh[j,1]=NaN
				list_cycleVs_vs0_nh[j,2]=NaN
			end
		end
	
	
	md""" Calculus of stable points, saddle points, maximum and minimum values of the limit cycle for each I"""
end

# ╔═╡ 32699c10-72f8-11eb-0590-0919e7d39a9b
begin
	bif_vs0_nh = plot(vs0_nh_range,list_stable_vs0_nh,linecolor="green",linewidth=1.5,label="Stable node")
	real_value_unstable_vs0_nh = findall(x -> x>=0 || x<0, list_unstable_vs0_nh)
	if length(real_value_unstable_vs0_nh)>0
		plot!(vs0_nh_range,list_unstable_vs0_nh[:],linecolor="red",linewidth=1.5,label="Unstable node")
	end
	real_value_cycleV_vs0_nh = findall(x -> x>=0 || x<0, list_cycleV_vs0_nh)
	if length(real_value_cycleV_vs0_nh)>0
		plot!(vs0_nh_range,list_cycleV_vs0_nh[:,1],linecolor="pink",linewidth=1.5,label="Maximum of the cycle")
		plot!(vs0_nh_range,list_cycleV_vs0_nh[:,2],linecolor="pink",linewidth=1.5,label="Minimum of the cycle")
	end
	plot!(vs0_nh_range,list_saddle_vs0_nh,linecolor="orange",linewidth=1.5,label="Saddle node",linestyle=:dashdot)
	xaxis!("Vs0")
	yaxis!("V")
	title!("Bifurcation diagram, I = $I_vs0_nh")
end

# ╔═╡ b605e9b0-7501-11eb-3cb4-19be80259632
begin
	ind_vs0_v0_approx_nh = findfirst(x ->abs(x-p[2])<0.35, vs0_nh_range)
	ind_vs0_h_approx_nh = findfirst(x ->abs(x-(-35))<0.5, vs0_nh_range)
	ind_vs0_l_approx_nh = findfirst(x ->abs(x-(-45))<0.5, vs0_nh_range)
	
	vs0_v0_approx_nh=vs0_nh_range[ind_vs0_v0_approx_nh[1]]
	vs0_h_approx_nh=vs0_nh_range[ind_vs0_h_approx_nh[1]]
	vs0_l_approx_nh=vs0_nh_range[ind_vs0_l_approx_nh[1]]
	
	vnull_vs0_v0_nh1=zeros(size(v_potential))
	vnull_vs0_v0_nh2=zeros(size(v_potential))
	vnull_vs0_h_nh1=zeros(size(v_potential))
	vnull_vs0_h_nh2=zeros(size(v_potential))
	vnull_vs0_l_nh1=zeros(size(v_potential))
	vnull_vs0_l_nh2=zeros(size(v_potential))
	
	for i=1:length(v_potential)
		vnull_vs0_v0_nh1[i],vnull_vs0_v0_nh2[i] = Vnullcline_(v_potential[i],p_used_vs0_nh[ind_vs0_v0_approx_nh[1],:])
		vnull_vs0_h_nh1[i],vnull_vs0_h_nh2[i] = Vnullcline_(v_potential[i],p_used_vs0_nh[ind_vs0_h_approx_nh[1],:])
		vnull_vs0_l_nh1[i],vnull_vs0_l_nh2[i] = Vnullcline_(v_potential[i],p_used_vs0_nh[ind_vs0_l_approx_nh[1],:])
	end
	
	phase_vs0_v0_nh = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_v0_nh1,v_potential,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_v0_nh2,v_potential,linecolor=RGB(1,0.06,1),linewidth = 1.5,label="dV/dt = 0",legend=:outertopright,size=(750,500))
	scatter!([list_stable_vs0_nh[ind_vs0_v0_approx_nh]],[list_stable_vs0_nh[ind_vs0_v0_approx_nh]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nh[ind_vs0_v0_approx_nh]],[list_saddle_vs0_nh[ind_vs0_v0_approx_nh]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_v0_approx_nl))")
	
	phase_vs0_h_nh = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_h_nh1,v_potential,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_h_nh2,v_potential,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_nh[ind_vs0_h_approx_nh]],[list_stable_vs0_nh[ind_vs0_h_approx_nh]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nh[ind_vs0_h_approx_nh]],[list_saddle_vs0_nh[ind_vs0_h_approx_nh]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_h_approx_ref))")
	
	phase_vs0_l_nh = plot(v_potential,vs_null,linecolor=RGB(0,0.6,0.8),linewidth = 1.5,label="dVs/dt = 0")
	plot!(vnull_vs0_l_nh1,v_potential,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	plot!(vnull_vs0_l_nh2,v_potential,linecolor=RGB(0.6,0.6,0.9),linewidth = 1.5,label="dV/dt = 0")
	scatter!([list_stable_vs0_nh[ind_vs0_l_approx_nh]],[list_stable_vs0_nh[ind_vs0_l_approx_nh]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")
	scatter!([list_saddle_vs0_nh[ind_vs0_l_approx_nh]],[list_saddle_vs0_nh[ind_vs0_l_approx_nh]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle",legend=:outertopright,size=(750,500))
	xaxis!("V")
	yaxis!("Vs")
	title!("Vs0 = $(round(vs0_l_approx_ref))")
	
	phase_vs0_nh = plot(phase_vs0_l_nh,phase_vs0_v0_nh,phase_vs0_h_nh,layout=(3,1),linewidth = 1.5)
	
	
	plot(bif_vs0_nh,phase_vs0_nh,layout=(1,2))
end

# ╔═╡ ae69e6c0-7de9-11eb-29ea-ab97f5a9a713
md" #### Bifurcation with conductances"

# ╔═╡ ce84ea6e-82b9-11eb-1364-b75ebf42f3f8
md"""see MQIF-conductances.jl"""

# ╔═╡ Cell order:
# ╟─f0a97dfe-6490-11eb-0684-c50efecf2bfd
# ╠═b9abc652-6491-11eb-0e92-7b7877bf9faf
# ╠═6a3dc280-6491-11eb-0ce0-6b468892be37
# ╠═6bb6b2a0-64ac-11eb-2b8e-b71abbdf3043
# ╠═50db35c0-6af4-11eb-3b81-b50d69d7a085
# ╠═d47d0250-6491-11eb-0850-3f015c450720
# ╠═9b413122-8e79-11eb-39d0-739f756ee127
# ╟─b654ccb0-654d-11eb-2b8c-1700940eb04b
# ╟─6bbef720-6548-11eb-220d-094c9f1f4504
# ╟─a447df80-6a25-11eb-2477-aba13d9a9d7e
# ╟─cca2ee5e-6549-11eb-3b9e-0934cddd0aa8
# ╟─4fa199b0-6bb7-11eb-32a6-0ff14be97216
# ╟─779bca42-6bc5-11eb-14b5-b38a2a99db3e
# ╟─d0569450-8b58-11eb-16b6-8762d560dafb
# ╟─c6476770-7f58-11eb-2fa2-5b79fa30d78e
# ╟─00178390-7f59-11eb-1c6c-75a2bcbedfc4
# ╟─b75c43ce-6aee-11eb-385c-3f46179d7089
# ╟─3dd5be70-7f88-11eb-01b7-a12dc52591ae
# ╟─1e623420-6afb-11eb-0e56-45472f98e9fe
# ╟─8e67dd20-7225-11eb-1489-6f6dbe8cbd6c
# ╠═200e05e0-6549-11eb-2151-51f90adb859f
# ╠═35ce42a0-6549-11eb-0201-d1a1b2f5ede1
# ╟─9f24b0c0-6654-11eb-2511-09afcd5a87e5
# ╠═7af74d3e-6562-11eb-16a4-5195708272f3
# ╠═137cf390-6643-11eb-1eb2-9f87c024400f
# ╠═9bad8b60-6654-11eb-1c88-eb5853740148
# ╠═28be7ed2-8e21-11eb-15c3-715f4d558d6c
# ╟─c7d4eec0-654d-11eb-0d20-c5657cc52a2e
# ╠═42327b20-64a8-11eb-1c0b-47200e371b58
# ╠═4f8838e0-6562-11eb-1d9b-b9cdffcece31
# ╠═98fcb230-6562-11eb-1194-8d99fc6dad03
# ╠═e2dac2be-6562-11eb-0607-a3ec86eccd70
# ╠═fe5732e0-6562-11eb-0a9b-a15a4acdb801
# ╟─f841cd00-6654-11eb-0a0d-f9c172a59637
# ╟─be6f1530-66ca-11eb-2c39-6b1a821a91e1
# ╠═c733ad50-64ab-11eb-087d-e79f7293c70b
# ╠═69e526f0-64ac-11eb-21c0-ff1dcf1ef711
# ╠═5af7bdb0-64ac-11eb-3b4e-f7d3efeead19
# ╟─333b3c00-6657-11eb-2876-536873d6a68a
# ╟─3225fb70-6657-11eb-3d47-871fb7a258b1
# ╟─d6475be0-64ae-11eb-0b97-a7455f721734
# ╟─16c03280-6655-11eb-2a36-05188dd54bc4
# ╠═e73b5fd0-66f9-11eb-137f-872a04c6a9cb
# ╟─b9796ca0-654a-11eb-17f7-6fb798368c65
# ╟─cd314fc0-66ca-11eb-1e8e-955995eaf311
# ╟─9cc8c8c0-663b-11eb-1dc9-fdb3997780e9
# ╠═f0a870ae-654d-11eb-2f67-7d8909eb2993
# ╟─cbde61fe-665e-11eb-13d8-07c58a3ef63b
# ╟─efb781e0-66ca-11eb-2ca5-6d9d9b3af2ba
# ╟─45add010-6575-11eb-0ff1-c3f6c939f522
# ╟─468a3b70-6577-11eb-26c4-376c402f17f1
# ╟─8bcd0190-66ee-11eb-0dc4-f5d028d63ae2
# ╟─8456b260-66ff-11eb-37ae-b9ca83e01afc
# ╟─c9180910-66ed-11eb-0f7d-251c6f0944c6
# ╟─a9fdbc40-66ee-11eb-2110-add301449188
# ╟─50a384de-66cb-11eb-1169-059806818ba0
# ╟─857d7b00-6642-11eb-10a8-35dd3fee920a
# ╟─57dc9440-66cb-11eb-12ca-0d8366d81ffd
# ╟─24f156b0-654f-11eb-11b5-71ee401f0a67
# ╠═758691f0-66fb-11eb-0bf8-7fe8ae75ad1b
# ╠═9b9f2480-66f9-11eb-3c16-174fe3e6ac25
# ╠═b6d90a30-66f5-11eb-220c-df794163eb6c
# ╠═fdac2f40-66f6-11eb-21ec-c59f7fe964d4
# ╟─3d281e40-6a1c-11eb-32a4-219e323e9af9
# ╠═47744e20-6a29-11eb-0485-dfbff003e543
# ╟─42085da0-6a29-11eb-3d97-5faa899d3c62
# ╟─31495e50-6a4d-11eb-21b8-21c20147d5c8
# ╟─aeab4d6e-6a1d-11eb-1e4c-b3cfe5e2d195
# ╟─62d7e340-6a22-11eb-0270-f1d22b89fe1c
# ╠═6cc0a950-6a22-11eb-0cd5-f54f55ab18f1
# ╟─4190fd10-7f59-11eb-18ce-b9ad0b8b8971
# ╟─3ac785c0-6b8c-11eb-12a0-ab81fa7fc6d1
# ╟─42679f30-6a2a-11eb-3476-57a1d9ae121d
# ╠═92f34f60-6a27-11eb-37cf-7773810bca67
# ╟─335a246e-6a22-11eb-0035-c564fda8d2fa
# ╟─3ce5efb0-6a1d-11eb-0c9a-3b230f943f01
# ╟─7fab57e0-6a22-11eb-2004-09b13178e043
# ╟─4375ad70-6a22-11eb-3579-c581e408c2a0
# ╟─7e404870-6a1d-11eb-2d4d-8d972c5a0993
# ╟─2a02c4c0-6a29-11eb-3c98-7135a450b2c2
# ╟─4af5d7f0-6a22-11eb-11c1-07234671e997
# ╟─7f1b3f20-6a1d-11eb-31b5-8f6eda38efb3
# ╟─a1af8f20-6a2a-11eb-099c-053cd7e54835
# ╟─de713f30-6a25-11eb-1c55-d7ef1f344cd0
# ╟─de03fe72-6a25-11eb-0d51-6743c6138f0c
# ╟─b3f3c480-6a2a-11eb-1485-877a22b19db1
# ╟─ddcf32d0-6a25-11eb-02ba-6f8e3097569a
# ╟─dd713450-6a25-11eb-3c88-8fc64c499bbe
# ╟─c868204e-6a2a-11eb-0f6a-6d945a550bf5
# ╟─500fd7e0-6a22-11eb-3dff-27c30b5d8460
# ╠═b9f04b50-6a21-11eb-031d-fdb7c4647aed
# ╟─90cffe60-6a20-11eb-159d-71fd0b3d58ee
# ╟─46d34140-6a2b-11eb-053b-27fa5e6b78bd
# ╟─8ae29620-6a52-11eb-1288-05ff0a31d717
# ╠═99a2fbf0-6a52-11eb-3dce-57131449d473
# ╟─acbe7e80-6a52-11eb-0383-71d5ce18e3ac
# ╟─c87b36fe-6b86-11eb-0a00-63d6799af669
# ╟─15d6004e-6b89-11eb-14f2-216a730247b3
# ╟─ce724fb0-6bb6-11eb-2671-cf3325b6d209
# ╟─602d37b2-6bb9-11eb-3b43-4ffcd0c06948
# ╟─f7f54c50-6b86-11eb-0ae2-c9a0c1cbe53c
# ╟─ca3bd1f2-6bb6-11eb-1670-6d6c3e2a68c8
# ╟─0b7de7c0-6bb7-11eb-1707-8da1e13d9992
# ╟─7f4e4fb0-7dea-11eb-2246-5dff4477acb7
# ╟─4732de90-6bb9-11eb-21bc-d9f23a7e5b02
# ╟─5690e3f0-6bb9-11eb-2ef7-a9fc4505ebfa
# ╟─825c4920-6bb9-11eb-1867-6b8e08310024
# ╟─c676b620-7dea-11eb-2d3f-0dd5a7da39db
# ╟─53576700-6bbd-11eb-3f52-a7e87ea6224d
# ╟─640ebad0-6bbd-11eb-2603-b997744c46a4
# ╟─82f18770-6bbd-11eb-0f2a-9db245be295d
# ╟─4f09f010-7deb-11eb-0aca-fd4072ea5648
# ╠═6a811010-6bc3-11eb-202f-23e520ba250f
# ╟─96c46370-7e92-11eb-0e5c-8fb1ebe10c89
# ╟─fb4ddda0-6bbc-11eb-0acf-db128b83a622
# ╟─a71ea31e-7e92-11eb-0b09-11bca509a568
# ╠═c445afe0-7deb-11eb-3726-0b309a0c35ec
# ╟─14386390-8e80-11eb-36ed-812d431e537b
# ╟─0f30b320-8e80-11eb-184c-a58dfa5bfdf2
# ╟─6f36d9f0-8ff4-11eb-2186-0938b0a1a009
# ╟─178a21b0-6bc5-11eb-37dd-1f23c2e9a31a
# ╟─a9f755ce-6bdf-11eb-114e-6b6938c70ad6
# ╟─665b7050-6bc5-11eb-393d-990a0aafbb60
# ╟─3f1e6e10-6bdf-11eb-3518-915c7e60782f
# ╟─033caa5e-6edc-11eb-0cf7-017566eccf30
# ╟─51709810-6bdd-11eb-3f50-41c504469836
# ╟─477972c0-6be0-11eb-1ca3-8157c59417fa
# ╟─4f91f1f0-6edb-11eb-3f31-7927596e94ea
# ╟─032b1940-6bde-11eb-0918-17fed0abd003
# ╟─37137dd0-6be1-11eb-1757-ab7ac23e4d7c
# ╟─84d16cf0-6edc-11eb-1078-8b5c6781e8a4
# ╠═8279a940-6bc6-11eb-1151-c1bc4a7882e1
# ╟─327db2c0-6bc5-11eb-3258-7314d88c9039
# ╟─31f6a870-6bc5-11eb-3afa-230cb80e952b
# ╟─32473970-6bc5-11eb-1e25-437e376ea951
# ╟─607d2540-6be1-11eb-04d7-719327df0298
# ╟─3c0c2f40-7225-11eb-2f84-811e39d11ccb
# ╟─bf34e5d0-6f94-11eb-15c6-1b081518100c
# ╟─c31d95a0-729d-11eb-07c7-d9c0db5dc535
# ╟─49d451f0-729f-11eb-3996-054cabdfbf5a
# ╟─bfdab820-72c8-11eb-0d57-0d8988ef5aa7
# ╟─a01bb0a0-745a-11eb-1554-2b2749deee03
# ╠═ef10370e-74fc-11eb-0893-23b71d91939d
# ╟─c3271fe0-6f9a-11eb-310f-5b846b6f646f
# ╟─a7b2f5a0-72e6-11eb-2f6b-5f6177cce876
# ╟─287ffb62-72e7-11eb-3460-6bc59ea30e12
# ╟─f510c1a0-72e7-11eb-148e-e1bebb6e8c44
# ╟─012613fe-74fe-11eb-2e19-ddee86dd3c7f
# ╠═36c49b72-8504-11eb-2652-bd34549800bd
# ╟─6720a03e-7053-11eb-3795-9fd39a387359
# ╟─f7158060-72e9-11eb-24d6-d52b86bf3848
# ╟─5f5323d0-72ea-11eb-392e-4d095f09c02c
# ╟─f6b1c010-72ea-11eb-02ce-07d898496dcd
# ╠═4ced1360-74ff-11eb-29e7-bde7e5380971
# ╟─94822410-72eb-11eb-1b05-eb82cb443e0c
# ╟─a4d8ba40-72eb-11eb-32ac-e7fd18d0fe2e
# ╟─a46d272e-72eb-11eb-13e5-4f15f3de10c3
# ╟─0dce5780-72f6-11eb-3d8d-63ad82ef0ca6
# ╟─37ca0b90-7500-11eb-0393-a7cb2ce63a08
# ╟─9ab3123e-72eb-11eb-08f4-77f036b5e088
# ╟─8ccb3bb0-72f7-11eb-2260-857b516bae2d
# ╟─a611b8fe-72f7-11eb-2386-9fb53a251d7c
# ╟─32699c10-72f8-11eb-0590-0919e7d39a9b
# ╟─b605e9b0-7501-11eb-3cb4-19be80259632
# ╟─ae69e6c0-7de9-11eb-29ea-ab97f5a9a713
# ╟─ce84ea6e-82b9-11eb-1364-b75ebf42f3f8
