### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c0d9fad0-828e-11eb-25b7-55870e15570f
begin
	using LinearAlgebra
	using DifferentialEquations
	using Plots
	using Statistics
end

# ╔═╡ d2da7b00-828f-11eb-351c-d5e58e150870
md" #### Packages"

# ╔═╡ c15b4350-828f-11eb-1306-b7c44b018922
md" #### Stability & nullclines functions"

# ╔═╡ 06214530-828f-11eb-2e9c-096dba628a53
function Vnullcline(v,p)
	I,v0,vs0,C,gf,gs,ts = p

	vs1 = vs0 + sqrt((gf*(v-v0)^2 + I)/gs)
	vs2 = vs0 - sqrt((gf*(v-v0)^2 + I)/gs)
	return vs1,vs2
end

# ╔═╡ 080a5082-828f-11eb-11e4-696eafae4279
function Vnullcline_(vs,p)
	I,v0,vs0,C,gf,gs,ts = p

	v1 = v0 + sqrt((gs*(vs-vs0)^2 - I)/gf)
	v2 = v0 - sqrt((gs*(vs-vs0)^2 - I)/gf)
	return v1,v2
end

# ╔═╡ 11033850-828f-11eb-2f09-4db9ddbe1d90
function Vsnullcline(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	vs = v
	return vs
end

# ╔═╡ 22f18580-828f-11eb-0d06-696d9c03d1b9
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

# ╔═╡ 2ca4bca0-828f-11eb-300e-2749156ae885
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

# ╔═╡ 35c0bcd0-828f-11eb-2f0a-65ddc147eedd
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

# ╔═╡ 49030140-828f-11eb-36ff-5d2834ec2f97
function jacobian(v,p)
	I,v0,vs0,C,gf,gs,ts = p
	J = zeros(2,2)
	J[1,1] = 2*gf*(v-v0)
	J[1,2] = -2*gs*(v-vs0)
	J[2,1] = 1
	J[2,2] = -1

	return J
end

# ╔═╡ 39090a00-828f-11eb-08aa-a74e80f4320b
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

# ╔═╡ 1e28c3f0-8ef2-11eb-3567-d9ad51b92cbd
function fixedpoints_gs_is_gf(p)
	I,v0,vs0,C,gf,gs,ts = p
	#with the vsnullcline, we know that vs = v, we then use this equality in the V-
	#nullcline and solve the system using the discriminant method
	if gs==gf
		vfp = (v0^2-vs0^2+I/gf)/(2*(v0-vs0))

		stability = zeros(1,2) #stability[1] >0 if the point (vfp1,vfp1) is unstable, 0 if it is a saddle and <0 if it is stable


		J = jacobian(vfp,p)
		lambda1,lambda2 = eigvals(J)
		if real(lambda1)>0 && real(lambda2)>0
			stability=1.0
		end
		if (real(lambda1)>0 && real(lambda2)<0) || (real(lambda1)<0 &&
			real(lambda2)>0)

			stability=0.0
		end
		if real(lambda1)<0 && real(lambda2)<0
			stability=-1.0
		end


		fp=ones(1,2).*vfp

		return fp,stability
	else
		return NaN,NaN
	end
end

# ╔═╡ 4f0899b0-828f-11eb-1585-8b45c1d5ffc0
function find_bifurcation(p)
	I,v0,vs0,C,gf,gs,ts = p
	vs0_bif_l = v0 - sqrt(I*(gf-gs)/(gf*gs))
	vs0_bif_h = v0 + sqrt(I*(gf-gs)/(gf*gs))

	return vs0_bif_l,vs0_bif_h
end

# ╔═╡ bd92a290-828f-11eb-0bd2-89e017a3729d
md" #### ODE problem functions"

# ╔═╡ b866ef60-828f-11eb-2a12-45a3edded36e
function MQIF!(du,u,p,t)
	I,v0,vs0,C,gf,gs,ts = p
 	du[1] = ( gf*(u[1]-v0)^2 - gs*(u[2]-vs0)^2 + I )/C
 	du[2] = (u[1]-u[2])/ts
end

# ╔═╡ 8ded85a0-828f-11eb-1b86-e9faa599e011
begin
	Vmax = -10 #30
	Vr = -40
	Vsr = -20
	md"""Voltages used for reset"""
end

# ╔═╡ db8039c0-828f-11eb-1f81-c9adfd61227c
function spike(x)  # spikes when spike(x) goes from negative to positive
    (x[1] - Vmax)
end

# ╔═╡ ec024cc0-828f-11eb-0043-df39ff4a0368
# event when event_f(u,t) == 0
function condition(x,t,integrator) #
    spike(x)
end

# ╔═╡ e66e13c0-828f-11eb-2710-15b8ca48500c
function reset!(x) # reset function
    x[1] = Vr
    x[2] = Vsr
end

# ╔═╡ f361f6f2-828f-11eb-24cf-db17dc893199
# when condition == 0 and upcrossing (from negative to positive)
function affect!(integrator)
    reset!(integrator.u)
end

# ╔═╡ 01d1a4ae-8290-11eb-2acc-bfa854d42c67
begin
	cb   = ContinuousCallback(condition,affect!,nothing)
	md"""Callback used for solver"""
end

# ╔═╡ ecf3951e-8290-11eb-38d5-07a7d2e1d074
md" #### Mesh & gradient functions"

# ╔═╡ 43662040-8290-11eb-0049-89099fa90249
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# ╔═╡ b327e5d0-8290-11eb-3bcd-491089f35abb
function return_MQIF_gradient(u,p,t)

	grad = [0.0,0.0]
	MQIF!(grad,u,p,0.0)

	return grad
end

# ╔═╡ 34401ad0-8290-11eb-20c9-3b683651b434
function gradient_calculation(xx,yy,p,scale)
	dx= zeros(size(xx))
	dy = zeros(size(xx))
	for k=1:length(xx)
		dx[k],dy[k]=
		return_MQIF_gradient([xx[k],yy[k]],p,0.0)
	end
	dx = dx*scale[1]
	dy = dy*scale[2]
	md""" Gradient"""
	return dx,dy
end

# ╔═╡ 70186300-8290-11eb-1c6f-75d6cb485173
begin
	limVmin = -70
	limVmax = -10
	limVsmin = limVmin
	limVsmax = limVmax
	md"""Mesh limits"""
end

# ╔═╡ 552e6ee2-8290-11eb-2e08-190e44a2318b
begin
	xx, yy = meshgrid(limVmin:4.0:limVmax, limVsmin:4.0:limVsmax)
	scale=[0.015,0.5]
	md"""Mesh used"""
end

# ╔═╡ 1acd0f40-8290-11eb-3799-a7b6fceba57d
md" #### Reference values"

# ╔═╡ 81b948f0-828f-11eb-07fa-278f0310284c
p=(0.1,-40.0,-35.0,1.0,1.0,0.5,10.0)

# ╔═╡ 16c67090-828f-11eb-15c8-d1b6cc6a1865
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

# ╔═╡ 55e642a0-828f-11eb-1583-01136f112cf0
cell_potential = collect(range(-70,stop=30,length=100))

# ╔═╡ 61d279d0-828f-11eb-001b-eb08bd724932
begin
	Vs = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		Vs[i]=Vsnullcline(cell_potential[i],p)
	end
	md"""Vs nullcline computation"""
end

# ╔═╡ 29a59b80-8291-11eb-1fc0-1d4e8fee5cc1
md" #### Convergence & limits functions for display of system behaviour with conductance change"

# ╔═╡ 03a3fe42-8291-11eb-1e64-873d92c6c4c9
function g_convergence(p,step_mesh)
	gf_mesh, gs_mesh = meshgrid(0.1:step_mesh:1, 0.1:step_mesh:1)

	stable_conv = zeros(size(gf_mesh))
	cycle_conv = zeros(size(gf_mesh))
	no_conv = zeros(size(gf_mesh))
	p_g = zeros(length(gf_mesh),length(p))
	saddle = zeros(size(gf_mesh))

	tspan=(0.0,500.0)

	for i=1:length(gf_mesh)
		p_g_= change_gf(gf_mesh[i],p)
		p_g_= change_gs(gs_mesh[i],p_g_)
		p_g[i,:]=p_g_

		max_low_nullcline = minimum(Vnullcline(p_g_[2],p_g_))
		u0=[-20,max_low_nullcline-1]

		prob = ODEProblem(MQIF!,u0,tspan,p_g_,callback=cb)

		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=3*maximum(tspan)/10, sol.t)
		ind_tperiod = findall(x -> x>=5*maximum(tspan)/10, sol.t)
		if maximum(size(ind_tperiod))>0
			test = sol[1,ind_tperiod[1]:end]
			ind_period = findall(x -> x>=-10.0-0.01,test)#is a spike
			if maximum(size(ind_period))>1
				step_period = ind_period[2:end] - ind_period[1:end-1]
				if abs(mean(step_period)-maximum(step_period))<3*maximum(step_period)/10
					cycle_conv[i]=1
				end
			end
		end
		end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
		if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5
			stable_conv[i]=1
		end

		u0 = [-60,max_low_nullcline-1]

		prob = ODEProblem(MQIF!,u0,tspan,p_g_,callback=cb)

		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=3*maximum(tspan)/10, sol.t)

		if p_g_[3]<p_g_[2]
			end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
			if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5
				stable_conv[i]=1
			else
				#if cycle_conv[i]==0
				if end_gradient[1]<-100 && sol[1,end] < -100
					no_conv[i]=1
				end
			end
		end

		if p_g_[3]>p_g_[2]
			end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
			if gf_mesh[i]==gs_mesh[i]
				#no convergence or cycle (stable node at infinity)
				is_below_max_low_nullcline = findall(x -> x<max_low_nullcline-5,sol[2,ind_tmin:end])
				if maximum(size(is_below_max_low_nullcline)) >0
					#no cycle for these IC
					#no convergence as no stable node
					no_conv[i]=1
				else
					#cycle
				end
			else
				#compute FP and their stability
				fp1,fp2,stab=fixedpoints(p_g_)
				fp=[fp1,fp2]
				if length(stab)>1
					ind_stable = findall(x -> x<0,stab)
					ind_unstable = findall(x -> x>0,stab)

					if length(ind_stable)>0
						if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5
							#stable
							stable_conv[i]=1
						end
					end
					if length(ind_unstable)>0
						if abs(end_gradient[1])>0.1
							#unstable
							no_conv[i]=1
						end
					end
				else
					is_below_max_low_nullcline = findall(x -> x<max_low_nullcline-5,sol[2,ind_tmin:end])
					if abs(end_gradient[1])>0.1 && maximum(size(is_below_max_low_nullcline)) >0
						#unstable
						no_conv[i]=1
					end
				end
			end
		end
	end

	return gf_mesh,gs_mesh,stable_conv,cycle_conv,no_conv
end

# ╔═╡ a7706122-8eef-11eb-2437-2705f6dec698
function g_convergence_st(p,step_mesh)
	gf_mesh, gs_mesh = meshgrid(step_mesh:step_mesh:1, step_mesh:step_mesh:1)

	stable_conv = zeros(size(gf_mesh))
	cycle_conv = zeros(size(gf_mesh))
	no_conv = zeros(size(gf_mesh))
	p_g = zeros(length(gf_mesh),length(p))
	saddle = zeros(size(gf_mesh))

	tspan=(0.0,500.0)

	for i=1:length(gf_mesh)
		p_g_= change_gf(gf_mesh[i],p)
		p_g_= change_gs(gs_mesh[i],p_g_)
		p_g[i,:]=p_g_

		max_low_nullcline = minimum(Vnullcline(p_g_[2],p_g_))
		u0=[-20,max_low_nullcline-1]

		prob = ODEProblem(MQIF!,u0,tspan,p_g_,callback=cb)

		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=3*maximum(tspan)/10, sol.t)
		ind_tperiod = findall(x -> x>=5*maximum(tspan)/10, sol.t)
		if maximum(size(ind_tperiod))>0
			test = sol[1,ind_tperiod[1]:end]
			ind_period = findall(x -> x>=-10.0-0.01,test)#is a spike
			if maximum(size(ind_period))>1
				step_period = ind_period[2:end] - ind_period[1:end-1]
				if abs(mean(step_period)-maximum(step_period))<3*maximum(step_period)/10
					cycle_conv[i]=1
				end
			end
		end
		end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
		if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5
			stable_conv[i]=1
		end


		u0 = [-60,max_low_nullcline-1]
		if gf_mesh[i]!=gs_mesh[i]
			fp1,fp2,stab=fixedpoints(p_g_)
			fp=[fp1,fp2]
			if length(stab)>1
				ind_saddle = findall(x -> x==0,stab)
				if length(ind_saddle)>0
					saddle = fp[ind_saddle[1]]
					if saddle[1]<p[2]
						u0 = [saddle[1]-5,saddle[2]-2]
					end
				end
			end
		else
			fp,stab = fixedpoints_gs_is_gf(p_g_)
			if stab==0
				saddle = fp
				if saddle[1]<p[2]
					u0 = [saddle[1]-5,saddle[2]-2]
				end
			end
		end


		prob = ODEProblem(MQIF!,u0,tspan,p_g_,callback=cb)

		sol = solve(prob,DP5(),reltol=1e-6,abstol=1e-6)
		ind_tmin = findfirst(x -> x>=3*maximum(tspan)/10, sol.t)

		if p_g_[3]<p_g_[2]
			end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
			end_gradient_prev = return_MQIF_gradient(sol[:,end-2],p_g_,0.0)
			if length(stab)>1
				ind_stable = findall(x -> x<0,stab)
				if length(ind_stable)>0
					ind_stable = ind_stable[1][2]
				else
					ind_stable = []
				end
			else
				ind_stable = []
			end
			
			if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5 && (length(ind_stable)>0)
				stableV = fp[ind_stable][1]
				if abs(sol[1,end]-stableV)<2
					stable_conv[i]=1
				end
			else
				if ((end_gradient[1]<-1)||(end_gradient[2]<-10)||(end_gradient[1]<end_gradient_prev[1])) && sol[1,end] < -100
					no_conv[i]=1
				end
			end
		end

		if p_g_[3]>p_g_[2]
			end_gradient = return_MQIF_gradient(sol[:,end],p_g_,0.0)
			if gf_mesh[i]==gs_mesh[i]
				#no convergence or cycle (stable node at infinity)
				is_below_max_low_nullcline = findall(x -> x<max_low_nullcline-5,sol[2,ind_tmin:end])
				if maximum(size(is_below_max_low_nullcline)) >0
					#no cycle for these IC
					#no convergence as no stable node
					no_conv[i]=1
				else
					#cycle
				end
			else
				#compute FP and their stability
				fp1,fp2,stab=fixedpoints(p_g_)
				fp=[fp1,fp2]
				if length(stab)>1
					ind_stable = findall(x -> x<0,stab)
					ind_unstable = findall(x -> x>0,stab)

					if length(ind_stable)>0
						if abs(end_gradient[1])<0.1 && abs(end_gradient[2])<0.5
							#stable
							stable_conv[i]=1
						end
					end
					if length(ind_unstable)>0
						if abs(end_gradient[1])>0.1
							#unstable
							no_conv[i]=1
						end
					end
				else
					is_below_max_low_nullcline = findall(x -> x<max_low_nullcline-5,sol[2,ind_tmin:end])
					if abs(end_gradient[1])>0.1 && maximum(size(is_below_max_low_nullcline)) >0
						#unstable
						no_conv[i]=1
					end
				end
			end
		end
	end

	return gf_mesh,gs_mesh,stable_conv,cycle_conv,no_conv
end

# ╔═╡ 15f4b6c0-8291-11eb-204a-392f432709f3
function stability_with_g(gf_mesh,gs_mesh,stable_conv,cycle_conv,no_conv)
	stable_only = zeros(size(stable_conv))
	cycle_only = zeros(size(stable_conv))
	bistable = zeros(size(stable_conv))
	no_conv_only = zeros(size(stable_conv))
	no_conv_cycle = zeros(size(stable_conv))
	no_conv_stable = zeros(size(stable_conv))

	for i=1:length(gf_mesh)
		if stable_conv[i] == 1
			if cycle_conv[i] ==1
				bistable[i]=1
			else
				if no_conv[i] ==1
					no_conv_stable[i]=1
				else
					stable_only[i]=1
				end
			end
		else
			if cycle_conv[i] ==1
				if no_conv[i]==1
					no_conv_cycle[i]=1
				else
					cycle_only[i]=1
				end
			else
				if no_conv[i]==1
					no_conv_only[i]=1
				end
			end
		end
	end

	return stable_only,cycle_only,bistable,no_conv_only,no_conv_cycle,no_conv_stable
end

# ╔═╡ ae142db0-8420-11eb-338c-7b9d5123aa59
function up_down_fill_limis(i,array,up,down,step_gf,step_gs,gf_mesh,gs_mesh)
	futur_index = collect(Int(i):step_gs:Int(i+((length(gf_mesh)/step_gs)-1)*step_gs))

	if maximum(array[futur_index])>0
		is_array = findall(x -> x>0,array[futur_index])
		if size(is_array)[1]>0
			down_=gs_mesh[futur_index[is_array[1]]]*100
			up_=gs_mesh[futur_index[is_array[end]]]*100
			if down_ > minimum(gf_mesh)*100 && down_ <=100
				down[i]=down_-step_gf*100/2
			else
				down[i]=down_
			end
			if up_ >= minimum(gf_mesh)*100 && up_ <100
				up[i]=up_+step_gf*100/2
			else
				up[i]=up_
			end
		end
	end

	return up,down
end

# ╔═╡ b68fb0d0-82af-11eb-3009-f92bad601035
function check_right_neigh_fill_limits(up,down,step_gs,gf_mesh)
	for j=1:step_gs-1
		if up[j]==down[j]
			if up[j+1]>minimum(gf_mesh)*100 && down[j+1]>minimum(gf_mesh)*100
				if j<step_gs-1
					if up[j+2]==up[j+1]
						up[j] = up[j+1]
						down[j] = up[j]
					end
					if down[j+2]==down[j+1]
						down[j] = down[j+1]
						up[j] = down[j]
					end
				end
			end
		end
	end
	return up,down
end

# ╔═╡ ac79602e-82b1-11eb-3faf-272868f80428
function check_left_neigh_fill_limits(up,down,step_gs,gf_mesh)
	for j=0:step_gs-2
		if up[step_gs-j]==down[step_gs-j]
			if up[step_gs-j-1]>minimum(gf_mesh)*100 && down[step_gs-j-1]>minimum(gf_mesh)*100
				if up[step_gs-j-1]==100
					up[step_gs-j] = up[step_gs-j-1]
					down[step_gs-j] = up[step_gs-j]
				end
			end
		end
	end
	return up,down
end

# ╔═╡ 4a982f20-82ae-11eb-0f18-d142243d24f9
function stability_fill_limits(gf_mesh,gs_mesh,stable_only,cycle_only,bistable,no_conv_only,no_conv_cycle,no_conv_stable)
	step_gs = findfirst(x -> x==1, gf_mesh)
	step_gf = gf_mesh[2]-gf_mesh[1]

	stable_up = ones(step_gs)*minimum(gf_mesh)*100
	stable_down = ones(step_gs)*minimum(gf_mesh)*100

	cycle_up = ones(step_gs)*minimum(gf_mesh)*100
	cycle_down = ones(step_gs)*minimum(gf_mesh)*100

	bistable_up = ones(step_gs)*minimum(gf_mesh)*100
	bistable_down = ones(step_gs)*minimum(gf_mesh)*100

	no_conv_only_up = ones(step_gs)*minimum(gf_mesh)*100
	no_conv_only_down = ones(step_gs)*minimum(gf_mesh)*100

	no_conv_cycle_up = ones(step_gs)*minimum(gf_mesh)*100
	no_conv_cycle_down = ones(step_gs)*minimum(gf_mesh)*100

	no_conv_stable_up = ones(step_gs)*minimum(gf_mesh)*100
	no_conv_stable_down = ones(step_gs)*minimum(gf_mesh)*100

	for i=1:step_gs #loop on gf
		stable_up,stable_down = up_down_fill_limis(i,stable_only,stable_up,stable_down,step_gf,step_gs,gf_mesh,gs_mesh)
		cycle_up,cycle_down = up_down_fill_limis(i,cycle_only,cycle_up,cycle_down,step_gf,step_gs,gf_mesh,gs_mesh)

		bistable_up,bistable_down = up_down_fill_limis(i,bistable,bistable_up,bistable_down,step_gf,step_gs,gf_mesh,gs_mesh)
		no_conv_only_up,no_conv_only_down = up_down_fill_limis(i,no_conv_only,no_conv_only_up,no_conv_only_down,step_gf,step_gs,gf_mesh,gs_mesh)
		no_conv_cycle_up,no_conv_cycle_down = up_down_fill_limis(i,no_conv_cycle,no_conv_cycle_up,no_conv_cycle_down,step_gf,step_gs,gf_mesh,gs_mesh)
		no_conv_stable_up,no_conv_stable_down = up_down_fill_limis(i,no_conv_stable,no_conv_stable_up,no_conv_stable_down,step_gf,step_gs,gf_mesh,gs_mesh)


		stable_up,stable_down = check_right_neigh_fill_limits(stable_up,stable_down,step_gs,gf_mesh)
		stable_up,stable_down = check_left_neigh_fill_limits(stable_up,stable_down,step_gs,gf_mesh)

		cycle_up,cycle_down = check_right_neigh_fill_limits(cycle_up,cycle_down,step_gs,gf_mesh)
		cycle_up,cycle_down = check_left_neigh_fill_limits(cycle_up,cycle_down,step_gs,gf_mesh)

		bistable_up,bistable_down = check_right_neigh_fill_limits(bistable_up,bistable_down,step_gs,gf_mesh)
		bistable_up,bistable_down = check_left_neigh_fill_limits(bistable_up,bistable_down,step_gs,gf_mesh)

		no_conv_only_up,no_conv_only_down = check_right_neigh_fill_limits(no_conv_only_up,no_conv_only_down,step_gs,gf_mesh)
		no_conv_only_up,no_conv_only_down = check_left_neigh_fill_limits(no_conv_only_up,no_conv_only_down,step_gs,gf_mesh)

		no_conv_cycle_up,no_conv_cycle_down = check_right_neigh_fill_limits(no_conv_cycle_up,no_conv_cycle_down,step_gs,gf_mesh)
		no_conv_cycle_up,no_conv_cycle_down = check_left_neigh_fill_limits(no_conv_cycle_up,no_conv_cycle_down,step_gs,gf_mesh)

		no_conv_stable_up,no_conv_stable_down = check_right_neigh_fill_limits(no_conv_stable_up,no_conv_stable_down,step_gs,gf_mesh)
		no_conv_stable_up,no_conv_stable_down = check_left_neigh_fill_limits(no_conv_stable_up,no_conv_stable_down,step_gs,gf_mesh)

	end

	return stable_up,stable_down,bistable_up,bistable_down,cycle_up,cycle_down,no_conv_only_up,no_conv_only_down,no_conv_cycle_up,no_conv_cycle_down,no_conv_stable_up,no_conv_stable_down,step_gs
end

# ╔═╡ 6948a750-8291-11eb-1391-17e89878f727
md" #### Bifurcation with conductances"

# ╔═╡ 1880632e-841c-11eb-1dfd-6306b2bdd5a6
begin
	pg_lowvs0_lowI=change_vs0(-45,change_I(0.1))
	mesh_step_lowvs0_lowI = 0.1
	gf_mesh_lowvs0_lowI,gs_mesh_lowvs0_lowI,stable_conv_lowvs0_lowI,cycle_conv_lowvs0_lowI,no_conv_l_l = g_convergence_st(pg_lowvs0_lowI,mesh_step_lowvs0_lowI)
	md"""Convergence for a mesh with steps of $mesh_step_lowvs0_lowI for restorative feedback (vs0 = $(pg_lowvs0_lowI[3])) and low current (I = $(pg_lowvs0_lowI[1]))"""
end

# ╔═╡ db0fb6fe-8437-11eb-0786-f911dacc34a6
begin
	pg_lowvs0_middleI=change_vs0(-45,change_I(5))
	mesh_step_lowvs0_middleI = 0.1
	gf_mesh_lowvs0_middleI,gs_mesh_lowvs0_middleI,stable_conv_lowvs0_middleI,cycle_conv_lowvs0_middleI,no_conv_l_m = g_convergence_st(pg_lowvs0_middleI,mesh_step_lowvs0_middleI)
	md"""Convergence for a mesh with steps of $mesh_step_lowvs0_middleI for restorative feedback (vs0 = $(pg_lowvs0_middleI[3])) and low current (I = $(pg_lowvs0_middleI[1]))"""
end

# ╔═╡ ce321af0-8437-11eb-249d-2b9dc595f572
begin
	pg_lowvs0_highI=change_vs0(-45,change_I(20))
	mesh_step_lowvs0_highI = 0.02
	gf_mesh_lowvs0_highI,gs_mesh_lowvs0_highI,stable_conv_lowvs0_highI,cycle_conv_lowvs0_highI,no_conv_l_h = g_convergence_st(pg_lowvs0_highI,mesh_step_lowvs0_highI)
	md"""Convergence for a mesh with steps of $mesh_step_lowvs0_highI for restorative feedback (vs0 = $(pg_lowvs0_highI[3])) and low current (I = $(pg_lowvs0_highI[1]))"""
end

# ╔═╡ 059fc6ba-6b28-4831-af18-a7f15f6f7ecf
begin
	pg_lowvs0_highI_=change_vs0(-45,change_I(20))
	mesh_step_lowvs0_highI_ = 0.1
	gf_mesh_lowvs0_highI_,gs_mesh_lowvs0_highI_,stable_conv_lowvs0_highI_,cycle_conv_lowvs0_highI_,no_conv_l_h_ = g_convergence(pg_lowvs0_highI,mesh_step_lowvs0_highI)
	md"""Convergence for a mesh with steps of $mesh_step_lowvs0_highI for restorative feedback (vs0 = $(pg_lowvs0_highI[3])) and low current (I = $(pg_lowvs0_highI[1]))"""
end

# ╔═╡ a7a811a2-8428-11eb-3b12-31a9eb1bbb4e
begin
	# pg_highvs0_lowI=change_vs0(-35,change_I(0.1))
	# mesh_step_highvs0_lowI = 0.1
	# 	gf_mesh_highvs0_lowI,gs_mesh_highvs0_lowI,stable_conv_highvs0_lowI,cycle_conv_highvs0_lowI,no_conv_h_l = g_convergence(pg_highvs0_lowI,mesh_step_highvs0_lowI)
	md"""Convergence for a mesh with steps of $mesh_step_highvs0_lowI for restorative feedback (vs0 = $(pg_highvs0_lowI[3])) and low current (I = $(pg_highvs0_lowI[1]))"""
end

# ╔═╡ cbba3d20-8437-11eb-1fb8-dbe2998b1e4c
begin
	# pg_highvs0_middleI=change_vs0(-35,change_I(5))
	# mesh_step_highvs0_middleI = 0.1
	# 	gf_mesh_highvs0_middleI,gs_mesh_highvs0_middleI,stable_conv_highvs0_middleI,cycle_conv_highvs0_middleI,no_conv_h_m = g_convergence(pg_highvs0_middleI,mesh_step_highvs0_middleI)
	md"""Convergence for a mesh with steps of $mesh_step_highvs0_middleI for restorative feedback (vs0 = $(pg_highvs0_middleI[3])) and low current (I = $(pg_highvs0_middleI[1]))"""
end

# ╔═╡ 57e45972-8438-11eb-16ec-7f914987e9f2
begin
	 pg_highvs0_highI=change_vs0(-35,change_I(20))
	# mesh_step_highvs0_highI = 0.1
	# 	gf_mesh_highvs0_highI,gs_mesh_highvs0_highI,stable_conv_highvs0_highI,cycle_conv_highvs0_highI,no_conv_h_h = g_convergence(pg_highvs0_highI,mesh_step_highvs0_highI)
	md"""Convergence for a mesh with steps of $mesh_step_highvs0_highI for restorative feedback (vs0 = $(pg_highvs0_highI[3])) and low current (I = $(pg_highvs0_highI[1]))"""
end

# ╔═╡ e788c0a2-8291-11eb-0ad7-7f271f14ae4a
begin
	stable_only_lowvs0_lowI,cycle_only_lowvs0_lowI,bistable_lowvs0_lowI,no_conv_only_lowvs0_lowI,no_conv_cycle_lowvs0_lowI,no_conv_stable_lowvs0_lowI = stability_with_g(gf_mesh_lowvs0_lowI,gs_mesh_lowvs0_lowI,stable_conv_lowvs0_lowI,cycle_conv_lowvs0_lowI,no_conv_l_l)
end

# ╔═╡ b3f663f2-8291-11eb-3ec9-27d8e55c112b
begin
	stable_only_lowvs0_middleI,cycle_only_lowvs0_middleI,bistable_lowvs0_middleI,no_conv_only_lowvs0_middleI,no_conv_cycle_lowvs0_middleI,no_conv_stable_lowvs0_middleI = stability_with_g(gf_mesh_lowvs0_middleI,gs_mesh_lowvs0_middleI,stable_conv_lowvs0_middleI,cycle_conv_lowvs0_middleI,no_conv_l_m)
end

# ╔═╡ 71be11d2-8292-11eb-15ec-7dd60814a2b1
begin
	stable_only_lowvs0_highI,cycle_only_lowvs0_highI,bistable_lowvs0_highI,no_conv_only_lowvs0_highI,no_conv_cycle_lowvs0_highI,no_conv_stable_lowvs0_highI = stability_with_g(gf_mesh_lowvs0_highI,gs_mesh_lowvs0_highI,stable_conv_lowvs0_highI,cycle_conv_lowvs0_highI,no_conv_l_h)
end

# ╔═╡ c8503610-8295-11eb-21d4-4980deab4482
begin
	# stable_only_highvs0_lowI,cycle_only_highvs0_lowI,bistable_highvs0_lowI,no_conv_only_highvs0_lowI,no_conv_cycle_highvs0_lowI,no_conv_stable_highvs0_lowI = stability_with_g(gf_mesh_highvs0_lowI,gs_mesh_highvs0_lowI,stable_conv_highvs0_lowI,cycle_conv_highvs0_lowI,no_conv_h_l)
end

# ╔═╡ a4677250-8438-11eb-20e2-11a4f00cd9b5
begin
	# stable_only_highvs0_middleI,cycle_only_highvs0_middleI,bistable_highvs0_middleI,no_conv_only_highvs0_middleI,no_conv_cycle_highvs0_middleI,no_conv_stable_highvs0_middleI = stability_with_g(gf_mesh_highvs0_middleI,gs_mesh_highvs0_middleI,stable_conv_highvs0_middleI,cycle_conv_highvs0_middleI,no_conv_h_m)
end

# ╔═╡ c9a0b680-8438-11eb-0b45-85f6673c3f11
begin
	# stable_only_highvs0_highI,cycle_only_highvs0_highI,bistable_highvs0_highI,no_conv_only_highvs0_highI,no_conv_cycle_highvs0_highI,no_conv_stable_highvs0_highI = stability_with_g(gf_mesh_highvs0_highI,gs_mesh_highvs0_highI,stable_conv_highvs0_highI,cycle_conv_highvs0_highI,no_conv_h_h)
end

# ╔═╡ fca694c0-84b3-11eb-3872-730076ea9f10
gr()

# ╔═╡ 10129dd0-83ef-11eb-1c9c-f73b40d2c04d
begin
	plot()
	for i=1:length(gf_mesh_lowvs0_highI)
		if no_conv_l_h[i]==1.0
			scatter!([gf_mesh_lowvs0_highI[i]],[gs_mesh_lowvs0_highI[i]],markershape = :diamond,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
		else
			#scatter!([gf_mesh_lowvs0_highI[i]],[gs_mesh_lowvs0_highI[i]],markershape = :diamond,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,0,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,legend=false)
		end
	end
	plot!(legend=false,size=(200,200))
	xaxis!("gf")
end

# ╔═╡ 35d6c747-c96c-48f9-925b-67f38ab71df7
begin
	plot()
	for i=1:length(gf_mesh_lowvs0_highI_)
		if cycle_conv_lowvs0_highI_[i]==1.0
			scatter!([gf_mesh_lowvs0_highI_[i]],[gs_mesh_lowvs0_highI_[i]],markershape = :diamond,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot)
		else
			scatter!([gf_mesh_lowvs0_highI_[i]],[gs_mesh_lowvs0_highI_[i]],markershape = :diamond,markersize = 5,markeralpha = 0.6,markercolor = RGB(1,0,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,legend=false)
		end
	end
	plot!(legend=false,size=(200,200))
	xaxis!("gf")
end

# ╔═╡ 90f164e0-8291-11eb-0f54-151c81ae3bae
md" ###### Tests"

# ╔═╡ 50c95000-8451-11eb-2f56-9977c8baf8e8
begin
	cond1 =[0.1,0.11] #[gf,gs]
	cond2 =[0.1,1]
	cond3 =[1,0.1]
	cond4 =[1,0.9]
	md"""Conductance values chosen for subplot simulations"""
end

# ╔═╡ 9b6f2d0e-8293-11eb-1f6f-71301a4774bb
begin
	pg_cond1_l_l = change_gs(cond1[2],change_gf(cond1[1],pg_lowvs0_lowI))
	pg_cond2_l_l = change_gs(cond2[2],change_gf(cond2[1],pg_lowvs0_lowI))
	pg_cond3_l_l = change_gs(cond3[2],change_gf(cond3[1],pg_lowvs0_lowI))
	pg_cond4_l_l = change_gs(cond4[2],change_gf(cond4[1],pg_lowvs0_lowI))

	#V-nullclines computation
	V1_cond1_l_l = zeros(size(cell_potential))
	V2_cond1_l_l = zeros(size(cell_potential))

	V1_cond2_l_l = zeros(size(cell_potential))
	V2_cond2_l_l = zeros(size(cell_potential))

	V1_cond3_l_l = zeros(size(cell_potential))
	V2_cond3_l_l = zeros(size(cell_potential))

	V1_cond4_l_l = zeros(size(cell_potential))
	V2_cond4_l_l = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_l_l[i],V2_cond1_l_l[i]=Vnullcline(cell_potential[i],pg_cond1_l_l)
		V1_cond2_l_l[i],V2_cond2_l_l[i]=Vnullcline(cell_potential[i],pg_cond2_l_l)
		V1_cond3_l_l[i],V2_cond3_l_l[i]=Vnullcline(cell_potential[i],pg_cond3_l_l)
		V1_cond4_l_l[i],V2_cond4_l_l[i]=Vnullcline(cell_potential[i],pg_cond4_l_l)
	end
end

# ╔═╡ ae391d20-914d-11eb-18ca-533973cac20c
pg_cond1_l_l

# ╔═╡ 0a819880-8444-11eb-06bc-57b6c2bec1e7
begin
	pg_cond1_l_m = change_gs(cond1[2],change_gf(cond1[1],pg_lowvs0_middleI))
	pg_cond2_l_m = change_gs(cond2[2],change_gf(cond2[1],pg_lowvs0_middleI))
	pg_cond3_l_m = change_gs(cond3[2],change_gf(cond3[1],pg_lowvs0_middleI))
	pg_cond4_l_m = change_gs(cond4[2],change_gf(cond4[1],pg_lowvs0_middleI))

	#V-nullclines computation
	V1_cond1_l_m = zeros(size(cell_potential))
	V2_cond1_l_m = zeros(size(cell_potential))

	V1_cond2_l_m = zeros(size(cell_potential))
	V2_cond2_l_m = zeros(size(cell_potential))

	V1_cond3_l_m = zeros(size(cell_potential))
	V2_cond3_l_m = zeros(size(cell_potential))

	V1_cond4_l_m = zeros(size(cell_potential))
	V2_cond4_l_m = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_l_m[i],V2_cond1_l_m[i]=Vnullcline(cell_potential[i],pg_cond1_l_m)
		V1_cond2_l_m[i],V2_cond2_l_m[i]=Vnullcline(cell_potential[i],pg_cond2_l_m)
		V1_cond3_l_m[i],V2_cond3_l_m[i]=Vnullcline(cell_potential[i],pg_cond3_l_m)
		V1_cond4_l_m[i],V2_cond4_l_m[i]=Vnullcline(cell_potential[i],pg_cond4_l_m)
	end
end

# ╔═╡ c8482202-8451-11eb-2192-cfecd9726aa1
begin
	pg_cond1_l_h = change_gs(cond1[2],change_gf(cond1[1],pg_lowvs0_highI))
	pg_cond2_l_h = change_gs(cond2[2],change_gf(cond2[1],pg_lowvs0_highI))
	pg_cond3_l_h = change_gs(cond3[2],change_gf(cond3[1],pg_lowvs0_highI))
	pg_cond4_l_h = change_gs(cond4[2],change_gf(cond4[1],pg_lowvs0_highI))

	#V-nullclines computation
	V1_cond1_l_h = zeros(size(cell_potential))
	V2_cond1_l_h = zeros(size(cell_potential))

	V1_cond2_l_h = zeros(size(cell_potential))
	V2_cond2_l_h = zeros(size(cell_potential))

	V1_cond3_l_h = zeros(size(cell_potential))
	V2_cond3_l_h = zeros(size(cell_potential))

	V1_cond4_l_h = zeros(size(cell_potential))
	V2_cond4_l_h = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_l_h[i],V2_cond1_l_h[i]=Vnullcline(cell_potential[i],pg_cond1_l_h)
		V1_cond2_l_h[i],V2_cond2_l_h[i]=Vnullcline(cell_potential[i],pg_cond2_l_h)
		V1_cond3_l_h[i],V2_cond3_l_h[i]=Vnullcline(cell_potential[i],pg_cond3_l_h)
		V1_cond4_l_h[i],V2_cond4_l_h[i]=Vnullcline(cell_potential[i],pg_cond4_l_h)
	end
end

# ╔═╡ 6b7737fe-84be-11eb-100e-bd880b7eefdd
begin
	pg_cond1_h_l = change_gs(cond1[2],change_gf(cond1[1],pg_highvs0_lowI))
	pg_cond2_h_l = change_gs(cond2[2],change_gf(cond2[1],pg_highvs0_lowI))
	pg_cond3_h_l = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_lowI))
	pg_cond4_h_l = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_lowI))
	pg_cond5_h_l = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_lowI))
	pg_cond6_h_l = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_lowI))

	#V-nullclines computation
	V1_cond1_h_l = zeros(size(cell_potential))
	V2_cond1_h_l = zeros(size(cell_potential))

	V1_cond2_h_l = zeros(size(cell_potential))
	V2_cond2_h_l = zeros(size(cell_potential))

	V1_cond3_h_l = zeros(size(cell_potential))
	V2_cond3_h_l = zeros(size(cell_potential))

	V1_cond4_h_l = zeros(size(cell_potential))
	V2_cond4_h_l = zeros(size(cell_potential))

	V1_cond5_h_l = zeros(size(cell_potential))
	V2_cond5_h_l = zeros(size(cell_potential))

	V1_cond6_h_l = zeros(size(cell_potential))
	V2_cond6_h_l = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_h_l[i],V2_cond1_h_l[i]=Vnullcline(cell_potential[i],pg_cond1_h_l)
		V1_cond2_h_l[i],V2_cond2_h_l[i]=Vnullcline(cell_potential[i],pg_cond2_h_l)
		V1_cond3_h_l[i],V2_cond3_h_l[i]=Vnullcline(cell_potential[i],pg_cond3_h_l)
		V1_cond4_h_l[i],V2_cond4_h_l[i]=Vnullcline(cell_potential[i],pg_cond4_h_l)
		V1_cond5_h_l[i],V2_cond5_h_l[i]=Vnullcline(cell_potential[i],pg_cond5_h_l)
		V1_cond6_h_l[i],V2_cond6_h_l[i]=Vnullcline(cell_potential[i],pg_cond6_h_l)
	end
	md"""Nullclines & parameters computations for conductances chosen, regenerative feedback and low current"""
end

# ╔═╡ 4669dea0-84e6-11eb-1ca9-137c44741630
begin
	pg_cond1_h_m = change_gs(cond1[2],change_gf(cond1[1],pg_highvs0_middleI))
	pg_cond2_h_m = change_gs(cond2[2],change_gf(cond2[1],pg_highvs0_middleI))
	pg_cond3_h_m = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_middleI))
	pg_cond4_h_m = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_middleI))
	pg_cond5_h_m = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_middleI))
	pg_cond6_h_m = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_middleI))

	#V-nullclines computation
	V1_cond1_h_m = zeros(size(cell_potential))
	V2_cond1_h_m = zeros(size(cell_potential))

	V1_cond2_h_m = zeros(size(cell_potential))
	V2_cond2_h_m = zeros(size(cell_potential))

	V1_cond3_h_m = zeros(size(cell_potential))
	V2_cond3_h_m = zeros(size(cell_potential))

	V1_cond4_h_m = zeros(size(cell_potential))
	V2_cond4_h_m = zeros(size(cell_potential))

	V1_cond5_h_m = zeros(size(cell_potential))
	V2_cond5_h_m = zeros(size(cell_potential))

	V1_cond6_h_m = zeros(size(cell_potential))
	V2_cond6_h_m = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_h_m[i],V2_cond1_h_m[i]=Vnullcline(cell_potential[i],pg_cond1_h_m)
		V1_cond2_h_m[i],V2_cond2_h_m[i]=Vnullcline(cell_potential[i],pg_cond2_h_m)
		V1_cond3_h_m[i],V2_cond3_h_m[i]=Vnullcline(cell_potential[i],pg_cond3_h_m)
		V1_cond4_h_m[i],V2_cond4_h_m[i]=Vnullcline(cell_potential[i],pg_cond4_h_m)
		V1_cond5_h_m[i],V2_cond5_h_m[i]=Vnullcline(cell_potential[i],pg_cond5_h_m)
		V1_cond6_h_m[i],V2_cond6_h_m[i]=Vnullcline(cell_potential[i],pg_cond6_h_m)
	end
	md"""Nullclines & parameters computations for conductances chosen, regenerative feedback and midlle current"""
end

# ╔═╡ a662a0a0-84e9-11eb-3002-d5db70f50e88
begin
	pg_cond1_h_h = change_gs(cond1[2],change_gf(cond1[1],pg_highvs0_highI))
	pg_cond2_h_h = change_gs(cond2[2],change_gf(cond2[1],pg_highvs0_highI))
	pg_cond3_h_h = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_highI))
	pg_cond4_h_h = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_highI))
	pg_cond5_h_h = change_gs(cond3[2],change_gf(cond3[1],pg_highvs0_highI))
	pg_cond6_h_h = change_gs(cond4[2],change_gf(cond4[1],pg_highvs0_highI))

	#V-nullclines computation
	V1_cond1_h_h = zeros(size(cell_potential))
	V2_cond1_h_h = zeros(size(cell_potential))

	V1_cond2_h_h = zeros(size(cell_potential))
	V2_cond2_h_h = zeros(size(cell_potential))

	V1_cond3_h_h = zeros(size(cell_potential))
	V2_cond3_h_h = zeros(size(cell_potential))

	V1_cond4_h_h = zeros(size(cell_potential))
	V2_cond4_h_h = zeros(size(cell_potential))

	V1_cond5_h_h = zeros(size(cell_potential))
	V2_cond5_h_h = zeros(size(cell_potential))

	V1_cond6_h_h = zeros(size(cell_potential))
	V2_cond6_h_h = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond1_h_h[i],V2_cond1_h_h[i]=Vnullcline(cell_potential[i],pg_cond1_h_h)
		V1_cond2_h_h[i],V2_cond2_h_h[i]=Vnullcline(cell_potential[i],pg_cond2_h_h)
		V1_cond3_h_h[i],V2_cond3_h_h[i]=Vnullcline(cell_potential[i],pg_cond3_h_h)
		V1_cond4_h_h[i],V2_cond4_h_h[i]=Vnullcline(cell_potential[i],pg_cond4_h_h)
		V1_cond5_h_h[i],V2_cond5_h_h[i]=Vnullcline(cell_potential[i],pg_cond5_h_h)
		V1_cond6_h_h[i],V2_cond6_h_h[i]=Vnullcline(cell_potential[i],pg_cond6_h_h)
	end
	md"""Nullclines & parameters computations for conductances chosen, regenerative feedback and high current"""
end

# ╔═╡ 1ec2a3c0-844e-11eb-3e3f-eb9d83bd99b7
md"""Restorative feedback with low current (vs0 = $(pg_lowvs0_lowI[3]), I = $(pg_lowvs0_lowI[1]))"""

# ╔═╡ a00e0990-8293-11eb-14f5-958a6f557dab
begin
	plot_cond1_l_l = plot(cell_potential,V1_cond1_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond1_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond1_ll = minimum(Vnullcline(pg_cond1_l_l[2],pg_cond1_l_l))
	u0_cycle_test_cond1_ll = [-20,max_low_nullcline_cond1_ll-1]
	u0_stable_test_cond1_ll = [-60,max_low_nullcline_cond1_ll-1]
	
	if pg_cond1_l_l[5]!=pg_cond1_l_l[6]
		cond1_l_l_fp1,cond1_l_l_fp2,cond1_l_l_stab=fixedpoints(pg_cond1_l_l)
		cond1_l_l_fp=[cond1_l_l_fp1,cond1_l_l_fp2]
		if length(cond1_l_l_stab)>1
			for i=1:length(cond1_l_l_stab)
				if cond1_l_l_stab[i] > 0
					scatter!([cond1_l_l_fp[i][1]],[cond1_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond1_l_l_stab[i] < 0
					scatter!([cond1_l_l_fp[i][1]],[cond1_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond1_l_l_stab[i] == 0
					scatter!([cond1_l_l_fp[i][1]],[cond1_l_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_1_l_l = cond1_l_l_fp[i]
					if saddle_1_l_l[1]<pg_cond1_l_l[2]
						u0_stable_test_cond1_ll = [saddle_1_l_l[1]-5,saddle_1_l_l[2]-2]
					end
				end
			end
		end
	end

	#Trajectory
	tspan_cond=(0.0,500.0)

	prob_cond1_l_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_ll,tspan_cond,pg_cond1_l_l,callback=cb)
	prob_cond1_l_l_st = ODEProblem(MQIF!,u0_stable_test_cond1_ll,tspan_cond,pg_cond1_l_l,callback=cb)

	sol_cond1_l_l_ct = solve(prob_cond1_l_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond1_l_l_st = solve(prob_cond1_l_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond1_l_l_st[1,:],sol_cond1_l_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond1_l_l_ct[1,:],sol_cond1_l_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V")
	yaxis!("Vs")
	title!("gf = $(cond1[1]) ; gs = $(cond1[2])")

end

# ╔═╡ 6f19b812-8f0f-11eb-12fc-d9ba60894ce9
begin
	cond5 =[1,0.669422]
	cond6 =[0.5,0.7]
	pg_cond5_l_h = change_gs(cond5[2],change_gf(cond5[1],pg_lowvs0_highI))

	#V-nullclines computation
	V1_cond5_l_h = zeros(size(cell_potential))
	V2_cond5_l_h = zeros(size(cell_potential))

	for i=1:length(cell_potential)
		V1_cond5_l_h[i],V2_cond5_l_h[i]=Vnullcline(cell_potential[i],pg_cond5_l_h)
	end

	plot(cell_potential,V1_cond5_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond5_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond5_lh = minimum(Vnullcline(pg_cond5_l_h[2],pg_cond5_l_h))
	u0_cycle_test_cond5_lh = [-20,max_low_nullcline_cond5_lh-1]
	u0_stable_test_cond5_lh = [-60,max_low_nullcline_cond5_lh-1]

	if pg_cond5_l_h[5]!=pg_cond5_l_h[6]
		cond5_l_h_fp1,cond5_l_h_fp2,cond5_l_h_stab=fixedpoints(pg_cond5_l_h)
		cond5_l_h_fp=[cond5_l_h_fp1,cond5_l_h_fp2]
		if length(cond5_l_h_stab)>1
			for i=1:length(cond5_l_h_stab)
				if cond5_l_h_stab[i] > 0
					scatter!([cond5_l_h_fp[i][1]],[cond5_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond5_l_h_stab[i] < 0
					scatter!([cond5_l_h_fp[i][1]],[cond5_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond5_l_h_stab[i] == 0
					scatter!([cond5_l_h_fp[i][1]],[cond5_l_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_5_l_h = cond5_l_h_fp[i]
					if saddle_5_l_h[1]<pg_cond5_l_h[2]
						u0_stable_test_cond5_lh = [saddle_5_l_h[1]-3,saddle_5_l_h[2]-1]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond5_l_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond5_lh,tspan_cond,pg_cond5_l_h,callback=cb)
	prob_cond5_l_h_st = ODEProblem(MQIF!,u0_stable_test_cond5_lh,tspan_cond,pg_cond5_l_h,callback=cb)

	sol_cond5_l_h_ct = solve(prob_cond5_l_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond5_l_h_st = solve(prob_cond5_l_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond5_l_h_st[1,:],sol_cond5_l_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond5_l_h_ct[1,:],sol_cond5_l_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))


	title!("gf = $(cond5[1]) ; gs = $(cond5[2]) ; I = $(pg_cond5_l_h[1])")

end

# ╔═╡ c3a876b0-8fcc-11eb-0d3e-530aeee079b4
cond5_l_h_fp

# ╔═╡ 94308d90-8442-11eb-05e8-c7dff6cc828e
begin
	plot_cond2_l_l = plot(cell_potential,V1_cond2_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond2_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond2_ll = minimum(Vnullcline(pg_cond2_l_l[2],pg_cond2_l_l))
	u0_cycle_test_cond2_ll = [-20,max_low_nullcline_cond2_ll-1]
	u0_stable_test_cond2_ll = [-60,max_low_nullcline_cond2_ll-1]
	
	if pg_cond2_l_l[5]!=pg_cond2_l_l[6]
		cond2_l_l_fp1,cond2_l_l_fp2,cond2_l_l_stab=fixedpoints(pg_cond2_l_l)
		cond2_l_l_fp=[cond2_l_l_fp1,cond2_l_l_fp2]
		if length(cond2_l_l_stab)>1
			for i=1:length(cond2_l_l_stab)
				if cond2_l_l_stab[i] > 0
					scatter!([cond2_l_l_fp[i][1]],[cond2_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond2_l_l_stab[i] < 0
					scatter!([cond2_l_l_fp[i][1]],[cond2_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond2_l_l_stab[i] == 0
					scatter!([cond2_l_l_fp[i][1]],[cond2_l_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_2_l_l = cond2_l_l_fp[i]
					if saddle_2_l_l[1]<pg_cond2_l_l[2]
						u0_stable_test_cond2_ll = [saddle_2_l_l[1]-2,saddle_2_l_l[2]]
					end
				end
			end
		end
	end

	#Trajectory


	prob_cond2_l_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_ll,tspan_cond,pg_cond2_l_l,callback=cb)
	prob_cond2_l_l_st = ODEProblem(MQIF!,u0_stable_test_cond2_ll,tspan_cond,pg_cond2_l_l,callback=cb)

	sol_cond2_l_l_ct = solve(prob_cond2_l_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond2_l_l_st = solve(prob_cond2_l_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond2_l_l_st[1,:],sol_cond2_l_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond2_l_l_ct[1,:],sol_cond2_l_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond2[1]) ; gs = $(cond2[2])")

end

# ╔═╡ 3e263a20-8443-11eb-378b-a31cd56f3611
begin
	plot_cond3_l_l = plot(cell_potential,V1_cond3_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond3_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	
	max_low_nullcline_cond3_ll = minimum(Vnullcline(pg_cond3_l_l[2],pg_cond3_l_l))
	u0_cycle_test_cond3_ll = [-20,max_low_nullcline_cond3_ll-1]
	u0_stable_test_cond3_ll = [-60,max_low_nullcline_cond3_ll-1]

	if pg_cond3_l_l[5]!=pg_cond3_l_l[6]
		cond3_l_l_fp1,cond3_l_l_fp2,cond3_l_l_stab=fixedpoints(pg_cond3_l_l)
		cond3_l_l_fp=[cond3_l_l_fp1,cond3_l_l_fp2]
		if length(cond3_l_l_stab)>1
			for i=1:length(cond3_l_l_stab)
				if cond3_l_l_stab[i] > 0
					scatter!([cond3_l_l_fp[i][1]],[cond3_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond3_l_l_stab[i] < 0
					scatter!([cond3_l_l_fp[i][1]],[cond3_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond3_l_l_stab[i] == 0
					scatter!([cond3_l_l_fp[i][1]],[cond3_l_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_3_l_l = cond3_l_l_fp[i]
					if saddle_3_l_l[1]<pg_cond3_l_l[2]
						u0_stable_test_cond3_ll = [saddle_3_l_l[1]-2,saddle_3_l_l[2]]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond3_l_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_ll,tspan_cond,pg_cond3_l_l,callback=cb)
	prob_cond3_l_l_st = ODEProblem(MQIF!,u0_stable_test_cond3_ll,tspan_cond,pg_cond3_l_l,callback=cb)

	sol_cond3_l_l_ct = solve(prob_cond3_l_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond3_l_l_st = solve(prob_cond3_l_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond3_l_l_st[1,:],sol_cond3_l_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond3_l_l_ct[1,:],sol_cond3_l_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond3[1]) ; gs = $(cond3[2])")

end

# ╔═╡ 9522e760-8443-11eb-3c70-97d2ec7e2c79
begin
	plot_cond4_l_l = plot(cell_potential,V1_cond4_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond4_l_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond4_ll = minimum(Vnullcline(pg_cond4_l_l[2],pg_cond4_l_l))
	u0_cycle_test_cond4_ll = [-20,max_low_nullcline_cond4_ll-1]
	u0_stable_test_cond4_ll = [-60,max_low_nullcline_cond4_ll-1]
	
	if pg_cond4_l_l[5]!=pg_cond4_l_l[6]
		cond4_l_l_fp1,cond4_l_l_fp2,cond4_l_l_stab=fixedpoints(pg_cond4_l_l)
		cond4_l_l_fp=[cond4_l_l_fp1,cond4_l_l_fp2]
		if length(cond4_l_l_stab)>1
			for i=1:length(cond4_l_l_stab)
				if cond4_l_l_stab[i] > 0
					scatter!([cond4_l_l_fp[i][1]],[cond4_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond4_l_l_stab[i] < 0
					scatter!([cond4_l_l_fp[i][1]],[cond4_l_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond4_l_l_stab[i] == 0
					scatter!([cond4_l_l_fp[i][1]],[cond4_l_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_4_l_l = cond4_l_l_fp[i]
					if saddle_4_l_l[1]<pg_cond4_l_l[2]
						u0_stable_test_cond4_ll = [saddle_4_l_l[1]-2,saddle_4_l_l[2]]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond4_l_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_ll,tspan_cond,pg_cond4_l_l,callback=cb)
	prob_cond4_l_l_st = ODEProblem(MQIF!,u0_stable_test_cond4_ll,tspan_cond,pg_cond4_l_l,callback=cb)

	sol_cond4_l_l_ct = solve(prob_cond4_l_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond4_l_l_st = solve(prob_cond4_l_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond4_l_l_st[1,:],sol_cond4_l_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond4_l_l_ct[1,:],sol_cond4_l_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond4[1]) ; gs = $(cond4[2])")

end

# ╔═╡ 4d034eb0-844e-11eb-058e-9319628a79be
md"""Restorative feedback with middle current (vs0 = $(pg_lowvs0_middleI[3]), I = $(pg_lowvs0_middleI[1]))"""

# ╔═╡ 4be2f61e-8444-11eb-27a8-318d0f58fec9
begin
	plot_cond1_l_m = plot(cell_potential,V1_cond1_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond1_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond1_lm = minimum(Vnullcline(pg_cond1_l_m[2],pg_cond1_l_m))
	u0_cycle_test_cond1_lm = [-20,max_low_nullcline_cond1_lm-1]
	u0_stable_test_cond1_lm = [-60,max_low_nullcline_cond1_lm-1]

	if pg_cond1_l_m[5]!=pg_cond1_l_m[6]
		cond1_l_m_fp1,cond1_l_m_fp2,cond1_l_m_stab=fixedpoints(pg_cond1_l_m)
		cond1_l_m_fp=[cond1_l_m_fp1,cond1_l_m_fp2]
		if length(cond1_l_m_stab)>1
			for i=1:length(cond1_l_m_stab)
				if cond1_l_m_stab[i] > 0
					scatter!([cond1_l_m_fp[i][1]],[cond1_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond1_l_m_stab[i] < 0
					scatter!([cond1_l_m_fp[i][1]],[cond1_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond1_l_m_stab[i] == 0
					scatter!([cond1_l_m_fp[i][1]],[cond1_l_m_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_1_l_m = cond1_l_m_fp[i]
					if saddle_1_l_m[1]<pg_cond1_l_m[2]
						u0_stable_test_cond1_lm = [saddle_1_l_m[1]-5,saddle_1_l_m[2]-2]
					end
				end
			end
		end
	end

	#Trajectory

	prob_cond1_l_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_lm,tspan_cond,pg_cond1_l_m,callback=cb)
	prob_cond1_l_m_st = ODEProblem(MQIF!,u0_stable_test_cond1_lm,tspan_cond,pg_cond1_l_m,callback=cb)

	sol_cond1_l_m_ct = solve(prob_cond1_l_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond1_l_m_st = solve(prob_cond1_l_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond1_l_m_st[1,:],sol_cond1_l_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond1_l_m_ct[1,:],sol_cond1_l_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V")
	yaxis!("Vs")
	title!("gf = $(cond1[1]) ; gs = $(cond1[2])")

end

# ╔═╡ 146538fe-8445-11eb-3810-7fab93bfe878
begin
	plot_cond2_l_m = plot(cell_potential,V1_cond2_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond2_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond2_lm = minimum(Vnullcline(pg_cond2_l_m[2],pg_cond2_l_m))
	u0_cycle_test_cond2_lm = [-20,max_low_nullcline_cond2_lm-1]
	u0_stable_test_cond2_lm = [-60,max_low_nullcline_cond2_lm-1]

	if pg_cond2_l_m[5]!=pg_cond2_l_m[6]
		cond2_l_m_fp1,cond2_l_m_fp2,cond2_l_m_stab=fixedpoints(pg_cond2_l_m)
		cond2_l_m_fp=[cond2_l_m_fp1,cond2_l_m_fp2]
		if length(cond2_l_m_stab)>1
			for i=1:length(cond2_l_m_stab)
				if cond2_l_m_stab[i] > 0
					scatter!([cond2_l_m_fp[i][1]],[cond2_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond2_l_m_stab[i] < 0
					scatter!([cond2_l_m_fp[i][1]],[cond2_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond2_l_m_stab[i] == 0
					scatter!([cond2_l_m_fp[i][1]],[cond2_l_m_fp[i][2]],marker = (:xcross, 5, 0.6, RGB(0.5,0,0.5), stroke(0.1, 0, :black, :dot)),label="Saddle")#saddle
					saddle_2_l_m = cond2_l_m_fp[i]
					if saddle_2_l_m[1]<pg_cond2_l_m[2]
						u0_stable_test_cond2_lm = [saddle_2_l_m[1]-2,saddle_2_l_m[2]]
					end
				end
			end
		end
	end

	#Trajectory

	prob_cond2_l_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_lm,tspan_cond,pg_cond2_l_m,callback=cb)
	prob_cond2_l_m_st = ODEProblem(MQIF!,u0_stable_test_cond2_lm,tspan_cond,pg_cond2_l_m,callback=cb)

	sol_cond2_l_m_ct = solve(prob_cond2_l_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond2_l_m_st = solve(prob_cond2_l_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond2_l_m_st[1,:],sol_cond2_l_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond2_l_m_ct[1,:],sol_cond2_l_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond2[1]) ; gs = $(cond2[2])")

end

# ╔═╡ 66962930-844c-11eb-0b45-316b1c7671a2
begin
	plot_cond3_l_m = plot(cell_potential,V1_cond3_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond3_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond3_lm = minimum(Vnullcline(pg_cond3_l_m[2],pg_cond3_l_m))
	u0_cycle_test_cond3_lm = [-20,max_low_nullcline_cond3_lm-1]
	u0_stable_test_cond3_lm = [-60,max_low_nullcline_cond3_lm-1]

	if pg_cond3_l_m[5]!=pg_cond3_l_m[6]
		cond3_l_m_fp1,cond3_l_m_fp2,cond3_l_m_stab=fixedpoints(pg_cond3_l_m)
		cond3_l_m_fp=[cond3_l_m_fp1,cond3_l_m_fp2]
		if length(cond3_l_m_stab)>1
			for i=1:length(cond3_l_m_stab)
				if cond3_l_m_stab[i] > 0
					scatter!([cond3_l_m_fp[i][1]],[cond3_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond3_l_m_stab[i] < 0
					scatter!([cond3_l_m_fp[i][1]],[cond3_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond3_l_m_stab[i] == 0
					scatter!([cond3_l_m_fp[i][1]],[cond3_l_m_fp[i][2]],marker = (:xcross, 5, 0.6, RGB(0.5,0,0.5), stroke(0.1, 0, :black, :dot)),label="Saddle")#saddle
					saddle_3_l_m = cond3_l_m_fp[i]
					if saddle_3_l_m[1]<pg_cond3_l_m[2]
						u0_stable_test_cond3_lm = [saddle_3_l_m[1]-2,saddle_3_l_m[2]]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond3_l_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_lm,tspan_cond,pg_cond3_l_m,callback=cb)
	prob_cond3_l_m_st = ODEProblem(MQIF!,u0_stable_test_cond3_lm,tspan_cond,pg_cond3_l_m,callback=cb)

	sol_cond3_l_m_ct = solve(prob_cond3_l_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond3_l_m_st = solve(prob_cond3_l_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond3_l_m_st[1,:],sol_cond3_l_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond3_l_m_ct[1,:],sol_cond3_l_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond3[1]) ; gs = $(cond3[2])")

end

# ╔═╡ bcd1d05e-844c-11eb-228d-1da658340995
begin
	plot_cond4_l_m = plot(cell_potential,V1_cond4_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond4_l_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond4_lm = minimum(Vnullcline(pg_cond4_l_m[2],pg_cond4_l_m))
	u0_cycle_test_cond4_lm = [-20,max_low_nullcline_cond4_lm-1]
	u0_stable_test_cond4_lm = [-60,max_low_nullcline_cond4_lm-1]

	if pg_cond4_l_m[5]!=pg_cond4_l_m[6]
		cond4_l_m_fp1,cond4_l_m_fp2,cond4_l_m_stab=fixedpoints(pg_cond4_l_m)
		cond4_l_m_fp=[cond4_l_m_fp1,cond4_l_m_fp2]
		if length(cond4_l_m_stab)>1
			for i=1:length(cond4_l_m_stab)
				if cond4_l_m_stab[i] > 0
					scatter!([cond4_l_m_fp[i][1]],[cond4_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond4_l_m_stab[i] < 0
					scatter!([cond4_l_m_fp[i][1]],[cond4_l_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond4_l_m_stab[i] == 0
					scatter!([cond4_l_m_fp[i][1]],[cond4_l_m_fp[i][2]],marker = (:xcross, 5, 0.6, RGB(0.5,0,0.5), stroke(0.1, 0, :black, :dot)),label="Saddle")#saddle
					saddle_4_l_m = cond4_l_m_fp[i]
					if saddle_4_l_m[1]<pg_cond4_l_m[2]
						u0_stable_test_cond4_lm = [saddle_4_l_m[1]-2,saddle_4_l_m[2]]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond4_l_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_lm,tspan_cond,pg_cond4_l_m,callback=cb)
	prob_cond4_l_m_st = ODEProblem(MQIF!,u0_stable_test_cond4_lm,tspan_cond,pg_cond4_l_m,callback=cb)

	sol_cond4_l_m_ct = solve(prob_cond4_l_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond4_l_m_st = solve(prob_cond4_l_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond4_l_m_st[1,:],sol_cond4_l_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond4_l_m_ct[1,:],sol_cond4_l_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V")
	yaxis!("Vs")
	title!("gf = $(cond4[1]) ; gs = $(cond4[2])")

end

# ╔═╡ 69bece80-844e-11eb-1039-77245c88f570
md"""Restorative feedback with high current (vs0 = $(pg_lowvs0_highI[3]), I = $(pg_lowvs0_highI[1]))"""

# ╔═╡ 7f4f1a40-8451-11eb-1a9c-3378bd3385df
begin
	plot_cond1_l_h = plot(cell_potential,V1_cond1_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond1_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond1_lh = minimum(Vnullcline(pg_cond1_l_h[2],pg_cond1_l_h))
	u0_cycle_test_cond1_lh = [-20,max_low_nullcline_cond1_lh-1]
	u0_stable_test_cond1_lh = [-60,max_low_nullcline_cond1_lh-1]

	if pg_cond1_l_h[5]!=pg_cond1_l_h[6]
		cond1_l_h_fp1,cond1_l_h_fp2,cond1_l_h_stab=fixedpoints(pg_cond1_l_h)
		cond1_l_h_fp=[cond1_l_h_fp1,cond1_l_h_fp2]
		if length(cond1_l_h_stab)>1
			for i=1:length(cond1_l_h_stab)
				if cond1_l_h_stab[i] > 0
					scatter!([cond1_l_h_fp[i][1]],[cond1_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond1_l_h_stab[i] < 0
					scatter!([cond1_l_h_fp[i][1]],[cond1_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond1_l_h_stab[i] == 0
					scatter!([cond1_l_h_fp[i][1]],[cond1_l_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_1_l_h = cond1_l_h_fp[i]
					if saddle_1_l_h[1]<pg_cond1_l_h[2]
						u0_stable_test_cond1_lh = [saddle_1_l_h[1]-5,saddle_1_l_h[2]-2]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond1_l_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_lh,tspan_cond,pg_cond1_l_h,callback=cb)
	prob_cond1_l_h_st = ODEProblem(MQIF!,u0_stable_test_cond1_lh,tspan_cond,pg_cond1_l_h,callback=cb)

	sol_cond1_l_h_ct = solve(prob_cond1_l_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond1_l_h_st = solve(prob_cond1_l_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond1_l_h_st[1,:],sol_cond1_l_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond1_l_h_ct[1,:],sol_cond1_l_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V")
	yaxis!("Vs")
	title!("gf = $(cond1[1]) ; gs = $(cond1[2])")

end

# ╔═╡ f6e803f0-8451-11eb-2d23-ff84295dd4a5
begin
	plot_cond2_l_h = plot(cell_potential,V1_cond2_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond2_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")

	max_low_nullcline_cond2_lh = minimum(Vnullcline(pg_cond2_l_h[2],pg_cond2_l_h))
	u0_cycle_test_cond2_lh = [-20,max_low_nullcline_cond2_lh-1]
	u0_stable_test_cond2_lh = [-60,max_low_nullcline_cond2_lh-1]

	if pg_cond2_l_h[5]!=pg_cond2_l_h[6]
		cond2_l_h_fp1,cond2_l_h_fp2,cond2_l_h_stab=fixedpoints(pg_cond2_l_h)
		cond2_l_h_fp=[cond2_l_h_fp1,cond2_l_h_fp2]
		if length(cond2_l_h_stab)>1
			for i=1:length(cond2_l_h_stab)
				if cond2_l_h_stab[i] > 0
					scatter!([cond2_l_h_fp[i][1]],[cond2_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond2_l_h_stab[i] < 0
					scatter!([cond2_l_h_fp[i][1]],[cond2_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond2_l_h_stab[i] == 0
					scatter!([cond2_l_h_fp[i][1]],[cond2_l_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
					saddle_2_l_h = cond2_l_h_fp[i]
					if saddle_2_l_h[1]<pg_cond2_l_h[2]
						u0_stable_test_cond2_lh = [saddle_2_l_h[1]-2,saddle_2_l_h[2]]
					end
				end
			end
		end
	end

	#Trajectory
	prob_cond2_l_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_lh,tspan_cond,pg_cond2_l_h,callback=cb)
	prob_cond2_l_h_st = ODEProblem(MQIF!,u0_stable_test_cond2_lh,tspan_cond,pg_cond2_l_h,callback=cb)

	sol_cond2_l_h_ct = solve(prob_cond2_l_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond2_l_h_st = solve(prob_cond2_l_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond2_l_h_st[1,:],sol_cond2_l_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond2_l_h_ct[1,:],sol_cond2_l_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond2[1]) ; gs = $(cond2[2])")

end

# ╔═╡ 4530c1a0-8452-11eb-3c24-1905fb6b0fd9
begin
	plot_cond3_l_h = plot(cell_potential,V1_cond3_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond3_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")


	if pg_cond3_l_h[5]!=pg_cond3_l_h[6]
		cond3_l_h_fp1,cond3_l_h_fp2,cond3_l_h_stab=fixedpoints(pg_cond3_l_h)
		cond3_l_h_fp=[cond3_l_h_fp1,cond3_l_h_fp2]
		if length(cond3_l_h_stab)>1
			for i=1:length(cond3_l_h_stab)
				if cond3_l_h_stab[i] > 0
					scatter!([cond3_l_h_fp[i][1]],[cond3_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond3_l_h_stab[i] < 0
					scatter!([cond3_l_h_fp[i][1]],[cond3_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond3_l_h_stab[i] == 0
					scatter!([cond3_l_h_fp[i][1]],[cond3_l_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
				end
			end
		end
	end

	#Trajectory
	max_low_nullcline_cond3_lh = minimum(Vnullcline(pg_cond3_l_h[2],pg_cond3_l_h))
	u0_cycle_test_cond3_lh = [-20,-40]
	u0_stable_test_cond3_lh = [-60,max_low_nullcline_cond3_lh-1]

	prob_cond3_l_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_lh,tspan_cond,pg_cond3_l_h,callback=cb)
	prob_cond3_l_h_st = ODEProblem(MQIF!,u0_stable_test_cond3_lh,tspan_cond,pg_cond3_l_h,callback=cb)

	sol_cond3_l_h_ct = solve(prob_cond3_l_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond3_l_h_st = solve(prob_cond3_l_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond3_l_h_st[1,:],sol_cond3_l_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond3_l_h_ct[1,:],sol_cond3_l_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	xaxis!("V",(-65,-15))
	yaxis!("Vs",(-65,-15))
	title!("gf = $(cond3[1]) ; gs = $(cond3[2])")

end

# ╔═╡ e1214210-8452-11eb-2d21-c98ad1ab60be
begin
	plot_cond4_l_h = plot(cell_potential,V1_cond4_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,V2_cond4_l_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")


	if pg_cond4_l_h[5]!=pg_cond4_l_h[6]
		cond4_l_h_fp1,cond4_l_h_fp2,cond4_l_h_stab=fixedpoints(pg_cond4_l_h)
		cond4_l_h_fp=[cond4_l_h_fp1,cond4_l_h_fp2]
		if length(cond4_l_h_stab)>1
			for i=1:length(cond4_l_h_stab)
				if cond4_l_h_stab[i] > 0
					scatter!([cond4_l_h_fp[i][1]],[cond4_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
				end
				if cond4_l_h_stab[i] < 0
					scatter!([cond4_l_h_fp[i][1]],[cond4_l_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
				end
				if cond4_l_h_stab[i] == 0
					scatter!([cond4_l_h_fp[i][1]],[cond4_l_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
				end
			end
		end
	end

	#Trajectory
	max_low_nullcline_cond4_lh = minimum(Vnullcline(pg_cond4_l_h[2],pg_cond4_l_h))
	u0_cycle_test_cond4_lh = [-20,max_low_nullcline_cond4_lh-1]
	u0_stable_test_cond4_lh = [-60,max_low_nullcline_cond4_lh-1]

	prob_cond4_l_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_lh,tspan_cond,pg_cond4_l_h,callback=cb)
	prob_cond4_l_h_st = ODEProblem(MQIF!,u0_stable_test_cond4_lh,tspan_cond,pg_cond4_l_h,callback=cb)

	sol_cond4_l_h_ct = solve(prob_cond4_l_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	sol_cond4_l_h_st = solve(prob_cond4_l_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	plot!(sol_cond4_l_h_st[1,:],sol_cond4_l_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	plot!(sol_cond4_l_h_ct[1,:],sol_cond4_l_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))

	title!("gf = $(cond4[1]) ; gs = $(cond4[2])")

end

# ╔═╡ f0ffa2f0-84be-11eb-140e-c3afab500eb8
md"""Regenerative feedback with low current (vs0 = $(pg_highvs0_lowI[3]), I = $(pg_highvs0_lowI[1]))"""

# ╔═╡ e73ab890-84be-11eb-2e91-71c8509b17f6
begin
	# plot_cond1_h_l = plot(cell_potential,V1_cond1_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond1_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond1_h_l[5]!=pg_cond1_h_l[6]
	# 	cond1_h_l_fp1,cond1_h_l_fp2,cond1_h_l_stab=fixedpoints(pg_cond1_h_l)
	# 	cond1_h_l_fp=[cond1_h_l_fp1,cond1_h_l_fp2]
	# 	if length(cond1_h_l_stab)>1
	# 		for i=1:length(cond1_h_l_stab)
	# 			if cond1_h_l_stab[i] > 0
	# 				scatter!([cond1_h_l_fp[i][1]],[cond1_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond1_h_l_stab[i] < 0
	# 				scatter!([cond1_h_l_fp[i][1]],[cond1_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond1_h_l_stab[i] == 0
	# 				scatter!([cond1_h_l_fp[i][1]],[cond1_h_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond1_hl = minimum(Vnullcline(pg_cond1_h_l[2],pg_cond1_h_l))
	# u0_cycle_test_cond1_hl = [-20,max_low_nullcline_cond1_hl-1]
	# u0_stable_test_cond1_hl = [-60,max_low_nullcline_cond1_hl-1]
	#
	# prob_cond1_h_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_hl,tspan_cond,pg_cond1_h_l,callback=cb)
	# prob_cond1_h_l_st = ODEProblem(MQIF!,u0_stable_test_cond1_hl,tspan_cond,pg_cond1_h_l,callback=cb)
	#
	# sol_cond1_h_l_ct = solve(prob_cond1_h_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond1_h_l_st = solve(prob_cond1_h_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond1_h_l_st[1,:],sol_cond1_h_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond1_h_l_ct[1,:],sol_cond1_h_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond1[1]) ; gs = $(cond1[2])")
	md"""Phase plane for (gf,gs)1 regenerative, low current"""
end

# ╔═╡ 57d2e0ee-84bf-11eb-3e89-49952681af5c
begin
	# plot_cond2_h_l = plot(cell_potential,V1_cond2_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond2_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond2_h_l[5]!=pg_cond2_h_l[6]
	# 	cond2_h_l_fp1,cond2_h_l_fp2,cond2_h_l_stab=fixedpoints(pg_cond2_h_l)
	# 	cond2_h_l_fp=[cond2_h_l_fp1,cond2_h_l_fp2]
	# 	if length(cond2_h_l_stab)>1
	# 		for i=1:length(cond2_h_l_stab)
	# 			if cond2_h_l_stab[i] > 0
	# 				scatter!([cond2_h_l_fp[i][1]],[cond2_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond2_h_l_stab[i] < 0
	# 				scatter!([cond2_h_l_fp[i][1]],[cond2_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond2_h_l_stab[i] == 0
	# 				scatter!([cond2_h_l_fp[i][1]],[cond2_h_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond2_hl = minimum(Vnullcline(pg_cond2_h_l[2],pg_cond2_h_l))
	# u0_cycle_test_cond2_hl = [-20,max_low_nullcline_cond2_hl-1]
	# u0_stable_test_cond2_hl = [-60,max_low_nullcline_cond2_hl-1]
	#
	# prob_cond2_h_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_hl,tspan_cond,pg_cond2_h_l,callback=cb)
	# prob_cond2_h_l_st = ODEProblem(MQIF!,u0_stable_test_cond2_hl,tspan_cond,pg_cond2_h_l,callback=cb)
	#
	# sol_cond2_h_l_ct = solve(prob_cond2_h_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond2_h_l_st = solve(prob_cond2_h_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond2_h_l_st[1,:],sol_cond2_h_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond2_h_l_ct[1,:],sol_cond2_h_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond2[1]) ; gs = $(cond2[2])")
	md"""Phase plane for (gf,gs)2 regenerative, low current"""
end

# ╔═╡ a34b8eb0-84bf-11eb-344c-0dd2bf7f836e
begin
	# plot_cond3_h_l = plot(cell_potential,V1_cond3_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond3_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond3_h_l[5]!=pg_cond3_h_l[6]
	# 	cond3_h_l_fp1,cond3_h_l_fp2,cond3_h_l_stab=fixedpoints(pg_cond3_h_l)
	# 	cond3_h_l_fp=[cond3_h_l_fp1,cond3_h_l_fp2]
	# 	if length(cond3_h_l_stab)>1
	# 		for i=1:length(cond3_h_l_stab)
	# 			if cond3_h_l_stab[i] > 0
	# 				scatter!([cond3_h_l_fp[i][1]],[cond3_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond3_h_l_stab[i] < 0
	# 				scatter!([cond3_h_l_fp[i][1]],[cond3_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond3_h_l_stab[i] == 0
	# 				scatter!([cond3_h_l_fp[i][1]],[cond3_h_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond3_hl = minimum(Vnullcline(pg_cond3_h_l[2],pg_cond3_h_l))
	# u0_cycle_test_cond3_hl = [-20,max_low_nullcline_cond3_hl-1]
	# u0_stable_test_cond3_hl = [-60,max_low_nullcline_cond3_hl-1]
	#
	# prob_cond3_h_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_hl,tspan_cond,pg_cond3_h_l,callback=cb)
	# prob_cond3_h_l_st = ODEProblem(MQIF!,u0_stable_test_cond3_hl,tspan_cond,pg_cond3_h_l,callback=cb)
	#
	# sol_cond3_h_l_ct = solve(prob_cond3_h_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond3_h_l_st = solve(prob_cond3_h_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond3_h_l_st[1,:],sol_cond3_h_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond3_h_l_ct[1,:],sol_cond3_h_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond3[1]) ; gs = $(cond3[2])")
	md"""Phase plane for (gf,gs)3 regenerative, low current"""
end

# ╔═╡ 1111ab50-84c0-11eb-192d-9f2bf6492011
begin
	# plot_cond4_h_l = plot(cell_potential,V1_cond4_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond4_h_l,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond4_h_l[5]!=pg_cond4_h_l[6]
	# 	cond4_h_l_fp1,cond4_h_l_fp2,cond4_h_l_stab=fixedpoints(pg_cond4_h_l)
	# 	cond4_h_l_fp=[cond4_h_l_fp1,cond4_h_l_fp2]
	# 	if length(cond4_h_l_stab)>1
	# 		for i=1:length(cond4_h_l_stab)
	# 			if cond4_h_l_stab[i] > 0
	# 				scatter!([cond4_h_l_fp[i][1]],[cond4_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond4_h_l_stab[i] < 0
	# 				scatter!([cond4_h_l_fp[i][1]],[cond4_h_l_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond4_h_l_stab[i] == 0
	# 				scatter!([cond4_h_l_fp[i][1]],[cond4_h_l_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond4_hl = minimum(Vnullcline(pg_cond4_h_l[2],pg_cond4_h_l))
	# u0_cycle_test_cond4_hl = [-20,max_low_nullcline_cond4_hl-1]
	# u0_stable_test_cond4_hl = [-60,max_low_nullcline_cond4_hl-1]
	#
	# prob_cond4_h_l_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_hl,tspan_cond,pg_cond4_h_l,callback=cb)
	# prob_cond4_h_l_st = ODEProblem(MQIF!,u0_stable_test_cond4_hl,tspan_cond,pg_cond4_h_l,callback=cb)
	#
	# sol_cond4_h_l_ct = solve(prob_cond4_h_l_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond4_h_l_st = solve(prob_cond4_h_l_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond4_h_l_st[1,:],sol_cond4_h_l_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond4_h_l_ct[1,:],sol_cond4_h_l_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond4[1]) ; gs = $(cond4[2])")
	md"""Phase plane for (gf,gs)4 regenerative, low current"""
end

# ╔═╡ fbad5f7e-84e6-11eb-2195-5f36dad7ead0
md"""Regenerative feedback with middle current (vs0 = $(pg_highvs0_middleI[3]), I = $(pg_highvs0_middleI[1]))"""

# ╔═╡ 0d119520-84e7-11eb-1f47-93c154733aa7
begin
	# plot_cond1_h_m = plot(cell_potential,V1_cond1_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond1_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond1_h_m[5]!=pg_cond1_h_m[6]
	# 	cond1_h_m_fp1,cond1_h_m_fp2,cond1_h_m_stab=fixedpoints(pg_cond1_h_m)
	# 	cond1_h_m_fp=[cond1_h_m_fp1,cond1_h_m_fp2]
	# 	if length(cond1_h_m_stab)>1
	# 		for i=1:length(cond1_h_m_stab)
	# 			if cond1_h_m_stab[i] > 0
	# 				scatter!([cond1_h_m_fp[i][1]],[cond1_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond1_h_m_stab[i] < 0
	# 				scatter!([cond1_h_m_fp[i][1]],[cond1_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond1_h_m_stab[i] == 0
	# 				scatter!([cond1_h_m_fp[i][1]],[cond1_h_m_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond1_hm = minimum(Vnullcline(pg_cond1_h_m[2],pg_cond1_h_m))
	# u0_cycle_test_cond1_hm = [-20,max_low_nullcline_cond1_hm-1]
	# u0_stable_test_cond1_hm = [-60,max_low_nullcline_cond1_hm-1]
	#
	# prob_cond1_h_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_hm,tspan_cond,pg_cond1_h_m,callback=cb)
	# prob_cond1_h_m_st = ODEProblem(MQIF!,u0_stable_test_cond1_hm,tspan_cond,pg_cond1_h_m,callback=cb)
	#
	# sol_cond1_h_m_ct = solve(prob_cond1_h_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond1_h_m_st = solve(prob_cond1_h_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond1_h_m_st[1,:],sol_cond1_h_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond1_h_m_ct[1,:],sol_cond1_h_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond1[1]) ; gs = $(cond1[2])")
	md"""Phase plane for (gf,gs)1 regenerative, middle current"""
end

# ╔═╡ 94bdf180-84e7-11eb-3f40-3f78aeefa66e
begin
	# plot_cond2_h_m = plot(cell_potential,V1_cond2_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond2_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond2_h_m[5]!=pg_cond2_h_m[6]
	# 	cond2_h_m_fp1,cond2_h_m_fp2,cond2_h_m_stab=fixedpoints(pg_cond2_h_m)
	# 	cond2_h_m_fp=[cond2_h_m_fp1,cond2_h_m_fp2]
	# 	if length(cond2_h_m_stab)>1
	# 		for i=1:length(cond2_h_m_stab)
	# 			if cond2_h_m_stab[i] > 0
	# 				scatter!([cond2_h_m_fp[i][1]],[cond2_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond2_h_m_stab[i] < 0
	# 				scatter!([cond2_h_m_fp[i][1]],[cond2_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond2_h_m_stab[i] == 0
	# 				scatter!([cond2_h_m_fp[i][1]],[cond2_h_m_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond2_hm = minimum(Vnullcline(pg_cond2_h_m[2],pg_cond2_h_m))
	# u0_cycle_test_cond2_hm = [-20,max_low_nullcline_cond2_hm-1]
	# u0_stable_test_cond2_hm = [-60,max_low_nullcline_cond2_hm-1]
	#
	# prob_cond2_h_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_hm,tspan_cond,pg_cond2_h_m,callback=cb)
	# prob_cond2_h_m_st = ODEProblem(MQIF!,u0_stable_test_cond2_hm,tspan_cond,pg_cond2_h_m,callback=cb)
	#
	# sol_cond2_h_m_ct = solve(prob_cond2_h_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond2_h_m_st = solve(prob_cond2_h_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond2_h_m_st[1,:],sol_cond2_h_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond2_h_m_ct[1,:],sol_cond2_h_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond2[1]) ; gs = $(cond2[2])")
	md"""Phase plane for (gf,gs)2 regenerative, middle current"""
end

# ╔═╡ 40f66220-84e8-11eb-0dd9-452166c81093
begin
	# plot_cond3_h_m = plot(cell_potential,V1_cond3_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond3_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond3_h_m[5]!=pg_cond3_h_m[6]
	# 	cond3_h_m_fp1,cond3_h_m_fp2,cond3_h_m_stab=fixedpoints(pg_cond3_h_m)
	# 	cond3_h_m_fp=[cond3_h_m_fp1,cond3_h_m_fp2]
	# 	if length(cond3_h_m_stab)>1
	# 		for i=1:length(cond3_h_m_stab)
	# 			if cond3_h_m_stab[i] > 0
	# 				scatter!([cond3_h_m_fp[i][1]],[cond3_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond3_h_m_stab[i] < 0
	# 				scatter!([cond3_h_m_fp[i][1]],[cond3_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond3_h_m_stab[i] == 0
	# 				scatter!([cond3_h_m_fp[i][1]],[cond3_h_m_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond3_hm = minimum(Vnullcline(pg_cond3_h_m[2],pg_cond3_h_m))
	# u0_cycle_test_cond3_hm = [-20,max_low_nullcline_cond3_hm-1]
	# u0_stable_test_cond3_hm = [-60,max_low_nullcline_cond3_hm-1]
	#
	# prob_cond3_h_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_hm,tspan_cond,pg_cond3_h_m,callback=cb)
	# prob_cond3_h_m_st = ODEProblem(MQIF!,u0_stable_test_cond3_hm,tspan_cond,pg_cond3_h_m,callback=cb)
	#
	# sol_cond3_h_m_ct = solve(prob_cond3_h_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond3_h_m_st = solve(prob_cond3_h_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond3_h_m_st[1,:],sol_cond3_h_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond3_h_m_ct[1,:],sol_cond3_h_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond3[1]) ; gs = $(cond3[2])")
	md"""Phase plane for (gf,gs)3 regenerative, middle current"""
end

# ╔═╡ 90b80430-84e8-11eb-1c30-a58116c0457f
begin
	# plot_cond4_h_m = plot(cell_potential,V1_cond4_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond4_h_m,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond4_h_m[5]!=pg_cond4_h_m[6]
	# 	cond4_h_m_fp1,cond4_h_m_fp2,cond4_h_m_stab=fixedpoints(pg_cond4_h_m)
	# 	cond4_h_m_fp=[cond4_h_m_fp1,cond4_h_m_fp2]
	# 	if length(cond4_h_m_stab)>1
	# 		for i=1:length(cond4_h_m_stab)
	# 			if cond4_h_m_stab[i] > 0
	# 				scatter!([cond4_h_m_fp[i][1]],[cond4_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond4_h_m_stab[i] < 0
	# 				scatter!([cond4_h_m_fp[i][1]],[cond4_h_m_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond4_h_m_stab[i] == 0
	# 				scatter!([cond4_h_m_fp[i][1]],[cond4_h_m_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond4_hm = minimum(Vnullcline(pg_cond4_h_m[2],pg_cond4_h_m))
	# u0_cycle_test_cond4_hm = [-20,max_low_nullcline_cond4_hm-1]
	# u0_stable_test_cond4_hm = [-60,max_low_nullcline_cond4_hm-1]
	#
	# prob_cond4_h_m_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_hm,tspan_cond,pg_cond4_h_m,callback=cb)
	# prob_cond4_h_m_st = ODEProblem(MQIF!,u0_stable_test_cond4_hm,tspan_cond,pg_cond4_h_m,callback=cb)
	#
	# sol_cond4_h_m_ct = solve(prob_cond4_h_m_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond4_h_m_st = solve(prob_cond4_h_m_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond4_h_m_st[1,:],sol_cond4_h_m_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond4_h_m_ct[1,:],sol_cond4_h_m_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond4[1]) ; gs = $(cond4[2])")
	md"""Phase plane for (gf,gs)4 regenerative, middle current"""
end

# ╔═╡ 40a7fe90-84e9-11eb-2130-a1905ebe5385
md"""Regenerative feedback with high current (vs0 = $(pg_highvs0_highI[3]), I = $(pg_highvs0_highI[1]))"""

# ╔═╡ 3b859a30-84e9-11eb-13db-034f42cd7388
begin
	# plot_cond1_h_h = plot(cell_potential,V1_cond1_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond1_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond1_h_h[5]!=pg_cond1_h_h[6]
	# 	cond1_h_h_fp1,cond1_h_h_fp2,cond1_h_h_stab=fixedpoints(pg_cond1_h_h)
	# 	cond1_h_h_fp=[cond1_h_h_fp1,cond1_h_h_fp2]
	# 	if length(cond1_h_h_stab)>1
	# 		for i=1:length(cond1_h_h_stab)
	# 			if cond1_h_h_stab[i] > 0
	# 				scatter!([cond1_h_h_fp[i][1]],[cond1_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond1_h_h_stab[i] < 0
	# 				scatter!([cond1_h_h_fp[i][1]],[cond1_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond1_h_h_stab[i] == 0
	# 				scatter!([cond1_h_h_fp[i][1]],[cond1_h_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond1_hh = minimum(Vnullcline(pg_cond1_h_h[2],pg_cond1_h_h))
	# u0_cycle_test_cond1_hh = [-20,max_low_nullcline_cond1_hh-1]
	# u0_stable_test_cond1_hh = [-60,max_low_nullcline_cond1_hh-1]
	#
	# prob_cond1_h_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond1_hh,tspan_cond,pg_cond1_h_h,callback=cb)
	# prob_cond1_h_h_st = ODEProblem(MQIF!,u0_stable_test_cond1_hh,tspan_cond,pg_cond1_h_h,callback=cb)
	#
	# sol_cond1_h_h_ct = solve(prob_cond1_h_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond1_h_h_st = solve(prob_cond1_h_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond1_h_h_st[1,:],sol_cond1_h_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond1_h_h_ct[1,:],sol_cond1_h_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond1[1]) ; gs = $(cond1[2])")
	md"""Phase plane for (gf,gs)1 regenerative, high current"""
end

# ╔═╡ 850541a0-84ea-11eb-2ee4-db89807fb982
begin
	# plot_cond2_h_h = plot(cell_potential,V1_cond2_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond2_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond2_h_h[5]!=pg_cond2_h_h[6]
	# 	cond2_h_h_fp1,cond2_h_h_fp2,cond2_h_h_stab=fixedpoints(pg_cond2_h_h)
	# 	cond2_h_h_fp=[cond2_h_h_fp1,cond2_h_h_fp2]
	# 	if length(cond2_h_h_stab)>1
	# 		for i=1:length(cond2_h_h_stab)
	# 			if cond2_h_h_stab[i] > 0
	# 				scatter!([cond2_h_h_fp[i][1]],[cond2_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond2_h_h_stab[i] < 0
	# 				scatter!([cond2_h_h_fp[i][1]],[cond2_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond2_h_h_stab[i] == 0
	# 				scatter!([cond2_h_h_fp[i][1]],[cond2_h_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond2_hh = minimum(Vnullcline(pg_cond2_h_h[2],pg_cond2_h_h))
	# u0_cycle_test_cond2_hh = [-20,max_low_nullcline_cond2_hh-1]
	# u0_stable_test_cond2_hh = [-60,max_low_nullcline_cond2_hh-1]
	#
	# prob_cond2_h_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond2_hh,tspan_cond,pg_cond2_h_h,callback=cb)
	# prob_cond2_h_h_st = ODEProblem(MQIF!,u0_stable_test_cond2_hh,tspan_cond,pg_cond2_h_h,callback=cb)
	#
	# sol_cond2_h_h_ct = solve(prob_cond2_h_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond2_h_h_st = solve(prob_cond2_h_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond2_h_h_st[1,:],sol_cond2_h_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond2_h_h_ct[1,:],sol_cond2_h_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond2[1]) ; gs = $(cond2[2])")
	md"""Phase plane for (gf,gs)2 regenerative, high current"""
end

# ╔═╡ d0bb1f70-84ea-11eb-1597-0dad8f20b183
begin
	# plot_cond3_h_h = plot(cell_potential,V1_cond3_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond3_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond3_h_h[5]!=pg_cond3_h_h[6]
	# 	cond3_h_h_fp1,cond3_h_h_fp2,cond3_h_h_stab=fixedpoints(pg_cond3_h_h)
	# 	cond3_h_h_fp=[cond3_h_h_fp1,cond3_h_h_fp2]
	# 	if length(cond3_h_h_stab)>1
	# 		for i=1:length(cond3_h_h_stab)
	# 			if cond3_h_h_stab[i] > 0
	# 				scatter!([cond3_h_h_fp[i][1]],[cond3_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond3_h_h_stab[i] < 0
	# 				scatter!([cond3_h_h_fp[i][1]],[cond3_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond3_h_h_stab[i] == 0
	# 				scatter!([cond3_h_h_fp[i][1]],[cond3_h_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end
	#
	# #Trajectory
	# max_low_nullcline_cond3_hh = minimum(Vnullcline(pg_cond3_h_h[2],pg_cond3_h_h))
	# u0_cycle_test_cond3_hh = [-20,max_low_nullcline_cond3_hh-1]
	# u0_stable_test_cond3_hh = [-60,max_low_nullcline_cond3_hh-1]
	#
	# prob_cond3_h_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond3_hh,tspan_cond,pg_cond3_h_h,callback=cb)
	# prob_cond3_h_h_st = ODEProblem(MQIF!,u0_stable_test_cond3_hh,tspan_cond,pg_cond3_h_h,callback=cb)
	#
	# sol_cond3_h_h_ct = solve(prob_cond3_h_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond3_h_h_st = solve(prob_cond3_h_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond3_h_h_st[1,:],sol_cond3_h_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond3_h_h_ct[1,:],sol_cond3_h_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond3[1]) ; gs = $(cond3[2])")
	md"""Phase plane for (gf,gs)3 regenerative, high current"""
end

# ╔═╡ 4bf5e800-84eb-11eb-0c40-1980d2947538
begin
	# plot_cond4_h_h = plot(cell_potential,V1_cond4_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,V2_cond4_h_h,linecolor=RGB(0.86,0.06,0.24),linewidth = 1.5,label="V-nullclines")
	# plot!(cell_potential,Vs,linecolor=RGB(0.13,0.55,0.13),linewidth = 1.5,label="Vs-nullcline")
	#
	#
	# if pg_cond4_h_h[5]!=pg_cond4_h_h[6]
	# 	cond4_h_h_fp1,cond4_h_h_fp2,cond4_h_h_stab=fixedpoints(pg_cond4_h_h)
	# 	cond4_h_h_fp=[cond4_h_h_fp1,cond4_h_h_fp2]
	# 	if length(cond4_h_h_stab)>1
	# 		for i=1:length(cond4_h_h_stab)
	# 			if cond4_h_h_stab[i] > 0
	# 				scatter!([cond4_h_h_fp[i][1]],[cond4_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	# 			end
	# 			if cond4_h_h_stab[i] < 0
	# 				scatter!([cond4_h_h_fp[i][1]],[cond4_h_h_fp[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	# 			end
	# 			if cond4_h_h_stab[i] == 0
	# 				scatter!([cond4_h_h_fp[i][1]],[cond4_h_h_fp[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	# 			end
	# 		end
	# 	end
	# end

	#Trajectory
	# max_low_nullcline_cond4_hh = minimum(Vnullcline(pg_cond4_h_h[2],pg_cond4_h_h))
	# u0_cycle_test_cond4_hh = [-20,max_low_nullcline_cond4_hh-1]
	# u0_stable_test_cond4_hh = [-60,max_low_nullcline_cond4_hh-1]
	#
	# prob_cond4_h_h_ct = ODEProblem(MQIF!,u0_cycle_test_cond4_hh,tspan_cond,pg_cond4_h_h,callback=cb)
	# prob_cond4_h_h_st = ODEProblem(MQIF!,u0_stable_test_cond4_hh,tspan_cond,pg_cond4_h_h,callback=cb)
	#
	# sol_cond4_h_h_ct = solve(prob_cond4_h_h_ct,DP5(),reltol=1e-8,abstol=1e-8)
	# sol_cond4_h_h_st = solve(prob_cond4_h_h_st,DP5(),reltol=1e-8,abstol=1e-8)
	# plot!(sol_cond4_h_h_st[1,:],sol_cond4_h_h_st[2,:],linecolor=:orange,linealpha=0.5,label="Trajectory",linewidth=4)
	# plot!(sol_cond4_h_h_ct[1,:],sol_cond4_h_h_ct[2,:],linecolor=:cyan,linealpha=0.5,label="Trajectory",linewidth=4,size=(500,500))
	#
	# xaxis!("V",(-65,-15))
	# yaxis!("Vs",(-65,-15))
	# title!("gf = $(cond4[1]) ; gs = $(cond4[2])")
	md"""Phase plane for (gf,gs)4 regenerative, high current"""
end

# ╔═╡ 9d3d6f4e-8291-11eb-27d2-d9e4d8867b0a
md" ###### Results"

# ╔═╡ fb7856f0-8450-11eb-3539-9d8d1896a63d
md"""Restorative feedback with low current (vs0 = $(pg_lowvs0_lowI[3]), I = $(pg_lowvs0_lowI[1]))"""

# ╔═╡ 5a2f0880-8292-11eb-1d00-37bcb6b6d3c0
begin
	stable_up_l_l,stable_down_l_l,bistable_up_l_l,bistable_down_l_l,cycle_up_l_l,cycle_down_l_l,no_conv_only_up_l_l,no_conv_only_down_l_l,no_conv_cycle_up_l_l,no_conv_cycle_down_l_l,no_conv_stable_up_l_l,no_conv_stable_down_l_l,step_gs_l_l = stability_fill_limits(gf_mesh_lowvs0_lowI,gs_mesh_lowvs0_lowI,stable_only_lowvs0_lowI,cycle_only_lowvs0_lowI,bistable_lowvs0_lowI,no_conv_only_lowvs0_lowI,no_conv_cycle_lowvs0_lowI,no_conv_stable_lowvs0_lowI)

	plot_cond_l_l = plot(legend=:outertopright,size=(1500,1000))
	if sum(stable_up_l_l-stable_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),stable_up_l_l, fill = (stable_down_l_l, 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	end
	if sum(bistable_up_l_l-bistable_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),bistable_up_l_l, fill = (bistable_down_l_l, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	end
	if sum(cycle_up_l_l-cycle_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),cycle_up_l_l, fill = (cycle_down_l_l, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	end
	if sum(no_conv_cycle_up_l_l-no_conv_cycle_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),no_conv_cycle_up_l_l, fill = (no_conv_cycle_down_l_l, 0.5, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
	end
	if sum(no_conv_stable_up_l_l-no_conv_stable_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),no_conv_stable_up_l_l, fill = (no_conv_stable_down_l_l, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	end
	if sum(no_conv_only_up_l_l-no_conv_only_down_l_l)>2
		plot!(gf_mesh_lowvs0_lowI[1:step_gs_l_l]*(100),no_conv_only_up_l_l, fill = (no_conv_only_down_l_l, 0.5, :red),linecolor=:red,linealpha=0,label="Instable")
	end
	xaxis!("gf (%)")
	yaxis!("gs (%)")
	title!("I = $(pg_lowvs0_lowI[1]) ; Vs0 = $(pg_lowvs0_lowI[3])")

	#savefig("BifConductances-45,01-0.1cond-coreected.pdf")
end

# ╔═╡ 8ef756c0-844b-11eb-04e6-cfc34ceb9b85
plot_simu_l_l = plot(plot_cond2_l_l,plot_cond4_l_l,plot_cond1_l_l,plot_cond3_l_l,layout=(2,2),legend = :bottomright,size=(2000,1600))

# ╔═╡ a6437de0-84ff-11eb-0269-c13490ee960c
begin
	plot(plot_cond_l_l,plot_simu_l_l,layout=(1,2),size=(1500,500))
	savefig("BifConductances-45,01-0.1cond.pdf")
end

# ╔═╡ f85b5f80-8450-11eb-256e-f34379563240
md"""Restorative feedback with middle current (vs0 = $(pg_lowvs0_middleI[3]), I = $(pg_lowvs0_middleI[1]))"""

# ╔═╡ c6ab8390-8291-11eb-1ac1-870443c01965
begin
	stable_up_l_m,stable_down_l_m,bistable_up_l_m,bistable_down_l_m,cycle_up_l_m,cycle_down_l_m,no_conv_only_up_l_m,no_conv_only_down_l_m,no_conv_cycle_up_l_m,no_conv_cycle_down_l_m,no_conv_stable_up_l_m,no_conv_stable_down_l_m,step_gs_l_m = stability_fill_limits(gf_mesh_lowvs0_middleI,gs_mesh_lowvs0_middleI,stable_only_lowvs0_middleI,cycle_only_lowvs0_middleI,bistable_lowvs0_middleI,no_conv_only_lowvs0_middleI,no_conv_cycle_lowvs0_middleI,no_conv_stable_lowvs0_middleI)

	plot_cond_l_m = plot(legend=:inside,size=(1500,1000))
	if sum(stable_up_l_m-stable_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),stable_up_l_m, fill = (stable_down_l_m, 0.5, :green),linecolor=:green,linealpha=0.5,label="Stable")
	end
	if sum(bistable_up_l_m-bistable_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),bistable_up_l_m, fill = (bistable_down_l_m, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	end
	if sum(cycle_up_l_m-cycle_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),cycle_up_l_m, fill = (cycle_down_l_m, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	end
	if sum(no_conv_stable_up_l_m-no_conv_stable_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),no_conv_stable_up_l_m, fill = (no_conv_stable_down_l_m, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	end
	if sum(no_conv_cycle_up_l_m-no_conv_cycle_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),no_conv_cycle_up_l_m, fill = (no_conv_cycle_down_l_m, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	end
	if sum(no_conv_only_up_l_m-no_conv_only_down_l_m)>2
		plot!(gf_mesh_lowvs0_middleI[1:step_gs_l_m]*(100),no_conv_only_up_l_m, fill = (no_conv_only_down_l_m, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	end
	xaxis!("gf (%)")
	yaxis!("gs (%)")
	title!("I = $(pg_lowvs0_middleI[1]) ; Vs0 = $(pg_lowvs0_middleI[3])")

	#savefig("BifConductances-45,5-0.1cond-corrected.pdf")
end

# ╔═╡ 303a2fc2-844d-11eb-18cd-07ba350d8651
plot_simu_l_m = plot(plot_cond2_l_m,plot_cond4_l_m,plot_cond1_l_m,plot_cond3_l_m,layout=(2,2),legend = :bottomright,size=(2000,1000))

# ╔═╡ 22e3f0f0-8500-11eb-25eb-1182e8caab4d
begin
	plot(plot_cond_l_m,plot_simu_l_m,layout=(1,2),size=(1500,500))
	savefig("BifConductances-45,5-0.1cond.pdf")
end

# ╔═╡ 08842ef0-8451-11eb-3b92-999a9bb70456
md"""Restorative feedback with high current (vs0 = $(pg_lowvs0_highI[3]), I = $(pg_lowvs0_highI[1]))"""

# ╔═╡ 6eb89a02-8292-11eb-311d-29eca634e7c5
begin
	stable_up_l_h,stable_down_l_h,bistable_up_l_h,bistable_down_l_h,cycle_up_l_h,cycle_down_l_h,no_conv_only_up_l_h,no_conv_only_down_l_h,no_conv_cycle_up_l_h,no_conv_cycle_down_l_h,no_conv_stable_up_l_h,no_conv_stable_down_l_h,step_gs_l_h = stability_fill_limits(gf_mesh_lowvs0_highI,gs_mesh_lowvs0_highI,stable_only_lowvs0_highI,cycle_only_lowvs0_highI,bistable_lowvs0_highI,no_conv_only_lowvs0_highI,no_conv_cycle_lowvs0_highI,no_conv_stable_lowvs0_highI)

	plot_cond_l_h = plot(legend=:outertopright,size=(1500,1000))
	if sum(stable_up_l_h-stable_down_l_h)>2
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),stable_up_l_h, fill = (stable_down_l_h, 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	end
	if sum(bistable_up_l_h-bistable_down_l_h)>2
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),bistable_up_l_h, fill = (bistable_down_l_h, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	end
	if sum(cycle_up_l_h-cycle_down_l_h)>2
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),cycle_up_l_h, fill = (cycle_down_l_h, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	end
	if sum(no_conv_stable_up_l_h-no_conv_stable_down_l_h)>2
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),no_conv_stable_up_l_h, fill = (no_conv_stable_down_l_h, 0.9, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	end
	if sum(no_conv_cycle_up_l_h-no_conv_cycle_down_l_h)>2
		#plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),no_conv_cycle_up_l_h, fill = (no_conv_cycle_down_l_h, 0.9, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),maximum(no_conv_cycle_up_l_h)*ones(size(no_conv_cycle_up_l_h)), fill = (cycle_up_l_h, 0.9, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
	end
	if sum(no_conv_only_up_l_h-no_conv_only_down_l_h)>2
		plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),no_conv_only_up_l_h, fill = (no_conv_only_down_l_h, 0.9, :red),linecolor=:red,linealpha=0,label="Instable",legend=:topleft)
	end
	xaxis!("gf (%)")
	yaxis!("gs (%)")
	title!("I = $(pg_lowvs0_highI[1]) ; Vs0 = $(pg_lowvs0_highI[3])")

	#savefig("BifConductances-45,20-0.1cond-corrected.pdf")
end

# ╔═╡ 39ba1320-8453-11eb-2e3b-99e86470d3d7
plot_simu_l_h = plot(plot_cond2_l_h,plot_cond4_l_h,plot_cond1_l_h,plot_cond3_l_h,layout=(2,2),legend = :bottomright,size=(2000,1000))

# ╔═╡ 7368bbf0-8500-11eb-1685-cb2d1fb50e47
begin
	plot(plot_cond_l_h,plot_simu_l_h,layout=(1,2),size=(1500,500))
	savefig("BifConductances-45,20-0.1cond.pdf")
end

# ╔═╡ 1311c080-8451-11eb-22b0-3b2971f221d3
md"""Regenerative feedback with low current (vs0 = $(pg_highvs0_lowI[3]), I = $(pg_highvs0_lowI[1]))"""

# ╔═╡ d4e6b700-8295-11eb-06bf-535c44212776
begin
	# stable_up_h_l,stable_down_h_l,bistable_up_h_l,bistable_down_h_l,cycle_up_h_l,cycle_down_h_l,no_conv_only_up_h_l,no_conv_only_down_h_l,no_conv_cycle_up_h_l,no_conv_cycle_down_h_l,no_conv_stable_up_h_l,no_conv_stable_down_h_l,step_gs_h_l = stability_fill_limits(gf_mesh_highvs0_lowI,gs_mesh_highvs0_lowI,stable_only_highvs0_lowI,cycle_only_highvs0_lowI,bistable_highvs0_lowI,no_conv_only_highvs0_lowI,no_conv_cycle_highvs0_lowI,no_conv_stable_highvs0_lowI)
	#
	# plot_cond_h_l = plot(legend=:outertopright,size=(1500,1000))
	# if sum(stable_up_h_l-stable_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),stable_up_h_l, fill = (stable_down_h_l, 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	# end
	# if sum(bistable_up_h_l-bistable_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),bistable_up_h_l, fill = (bistable_down_h_l, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	# end
	# if sum(cycle_up_h_l-cycle_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),cycle_up_h_l, fill = (cycle_down_h_l, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	# end
	# if sum(no_conv_cycle_up_h_l-no_conv_cycle_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),no_conv_cycle_up_h_l, fill = (no_conv_cycle_down_h_l, 0.5, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
	# end
	# if sum(no_conv_stable_up_h_l-no_conv_stable_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),no_conv_stable_up_h_l, fill = (no_conv_stable_down_h_l, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	# end
	# if sum(no_conv_only_up_h_l-no_conv_only_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_lowI[1:step_gs_h_l]*(100),no_conv_only_up_h_l, fill = (no_conv_only_down_h_l, 0.5, :red),linecolor=:red,linealpha=0,label="Instable")
	# end
	# xaxis!("gf (%)")
	# yaxis!("gs (%)")
	# title!("I = $(pg_highvs0_lowI[1]) ; Vs0 = $(pg_highvs0_lowI[3])")

	#savefig("BifConductances-35,01-0.001condToVerify.pdf")
end

# ╔═╡ 6e09429e-84c0-11eb-2a52-730ade938fc9
# plot_simu_h_l = plot(plot_cond2_h_l,plot_cond4_h_l,plot_cond1_h_l,plot_cond3_h_l,layout=(2,2),legend=:outertopright,size=(1500,1500))

# ╔═╡ df2d8a50-8500-11eb-15bc-cbf21b2725d0
begin
	# plot(plot_cond_h_l,plot_simu_h_l,layout=(1,2),size=(1500,500))
	#savefig("BifConductances-35,01-0.1cond.pdf")
end

# ╔═╡ 2756a7e0-8451-11eb-105b-91e1f571f79c
md"""Regenerative feedback with middle current (vs0 = $(pg_highvs0_middleI[3]), I = $(pg_highvs0_middleI[1]))"""

# ╔═╡ f7403f70-8438-11eb-3e06-e758dc5bd7c8
begin
	# stable_up_h_m,stable_down_h_m,bistable_up_h_m,bistable_down_h_m,cycle_up_h_m,cycle_down_h_m,no_conv_only_up_h_m,no_conv_only_down_h_m,no_conv_cycle_up_h_m,no_conv_cycle_down_h_m,no_conv_stable_up_h_m,no_conv_stable_down_h_m,step_gs_h_m = stability_fill_limits(gf_mesh_highvs0_middleI,gs_mesh_highvs0_middleI,stable_only_highvs0_middleI,cycle_only_highvs0_middleI,bistable_highvs0_middleI,no_conv_only_highvs0_middleI,no_conv_cycle_highvs0_middleI,no_conv_stable_highvs0_middleI)
	#
	# plot_cond_h_m = plot(legend=:outertopright,size=(1500,1000))
	# if sum(stable_up_h_m-stable_down_h_m)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),stable_up_h_m, fill = (stable_down_h_m, 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	# end
	# if sum(bistable_up_h_m-bistable_down_h_l)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),bistable_up_h_m, fill = (bistable_down_h_m, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	# end
	# if sum(cycle_up_h_m-cycle_down_h_m)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),cycle_up_h_m, fill = (cycle_down_h_m, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	# end
	# if sum(no_conv_cycle_up_h_m-no_conv_cycle_down_h_m)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),no_conv_cycle_up_h_m, fill = (no_conv_cycle_down_h_m, 0.5, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
	# end
	# if sum(no_conv_stable_up_h_m-no_conv_stable_down_h_m)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),no_conv_stable_up_h_m, fill = (no_conv_stable_down_h_m, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	# end
	# if sum(no_conv_only_up_h_m-no_conv_only_down_h_m)>2
	# 	plot!(gf_mesh_highvs0_middleI[1:step_gs_h_m]*(100),no_conv_only_up_h_m, fill = (no_conv_only_down_h_m, 0.5, :red),linecolor=:red,linealpha=0,label="Instable")
	# end
	# xaxis!("gf (%)")
	# yaxis!("gs (%)")
	# title!("I = $(pg_highvs0_middleI[1]) ; Vs0 = $(pg_highvs0_middleI[3])")

	#savefig("BifConductances-35,01-0.001condToVerify.pdf")
end

# ╔═╡ 0516d180-84e9-11eb-2eee-d35a21ae7923
# plot_simu_h_m = plot(plot_cond2_h_m,plot_cond4_h_m,plot_cond1_h_m,plot_cond3_h_m,layout=(2,2),legend=:outertopright,size=(1500,1500))

# ╔═╡ a929c300-8501-11eb-0d94-eb34ede9726e
begin
	# plot(plot_cond_h_m,plot_simu_h_m,layout=(1,2),size=(1500,500))
	#savefig("BifConductances-35,5-0.1cond.pdf")
end

# ╔═╡ 3234f180-8451-11eb-080b-cde94036350f
md"""Regenerative feedback with high current (vs0 = $(pg_highvs0_highI[3]), I = $(pg_highvs0_highI[1]))"""

# ╔═╡ 04fb3140-843b-11eb-1baa-bb5440958fb6
begin
	# stable_up_h_h,stable_down_h_h,bistable_up_h_h,bistable_down_h_h,cycle_up_h_h,cycle_down_h_h,no_conv_only_up_h_h,no_conv_only_down_h_h,no_conv_cycle_up_h_h,no_conv_cycle_down_h_h,no_conv_stable_up_h_h,no_conv_stable_down_h_h,step_gs_h_h = stability_fill_limits(gf_mesh_highvs0_highI,gs_mesh_highvs0_highI,stable_only_highvs0_highI,cycle_only_highvs0_highI,bistable_highvs0_highI,no_conv_only_highvs0_highI,no_conv_cycle_highvs0_highI,no_conv_stable_highvs0_highI)
	#
	# plot_cond_h_h = plot(legend=:outertopright,size=(1000,1000))
	# if sum(stable_up_h_h-stable_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),stable_up_h_h, fill = (stable_down_h_h, 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	# end
	# if sum(bistable_up_h_h-bistable_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),bistable_up_h_h, fill = (bistable_down_h_h, 0.5, :orange),linecolor=:orange,linealpha=0,label="Bistable")
	# end
	# if sum(cycle_up_h_h-cycle_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),cycle_up_h_h, fill = (cycle_down_h_h, 0.9, :pink),linecolor=:pink,linealpha=0,label="Cycle",legend=:outertopright)
	# end
	# if sum(no_conv_cycle_up_h_h-no_conv_cycle_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),no_conv_cycle_up_h_h, fill = (no_conv_cycle_down_h_h, 0.5, :blue),linecolor=:blue,linealpha=0,label="Cycle/Instable")
	# end
	# if sum(no_conv_stable_up_h_h-no_conv_stable_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),no_conv_stable_up_h_h, fill = (no_conv_stable_down_h_h, 0.5, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	# end
	# if sum(no_conv_only_up_h_h-no_conv_only_down_h_h)>2
	# 	plot!(gf_mesh_highvs0_highI[1:step_gs_h_h]*(100),no_conv_only_up_h_h, fill = (no_conv_only_down_h_h, 0.5, :red),linecolor=:red,linealpha=0,label="Instable")
	# end
	# xaxis!("gf (%)")
	# yaxis!("gs (%)")
	# title!("I = $(pg_highvs0_highI[1]) ; Vs0 = $(pg_highvs0_highI[3])")

	#savefig("BifConductances-35,01-0.001condToVerify.pdf")
end

# ╔═╡ b1677320-84eb-11eb-2f97-c5d24260f5fd
# plot_simu_h_h = plot(plot_cond2_h_h,plot_cond4_h_h,plot_cond1_h_h,plot_cond3_h_h,layout=(2,2),legend=:outertopright,size=(1500,1500))

# ╔═╡ 0ad59390-8502-11eb-2f03-0bdc114d3620
begin
	# plot(plot_cond_h_h,plot_simu_h_h,layout=(1,2),size=(1500,500))
	#savefig("BifConductances-35,20-0.1cond.pdf")
end

# ╔═╡ c66ffc6e-9a4a-4ce1-bb52-fefbb76677ab


# ╔═╡ 0f0b344a-f551-4d28-823f-ec81b6ff2f92
function sn_gs(vs0,v0,gf,I)
	gs = gf*I/(I+gf*(v0-vs0)^2)
	return gs
end

# ╔═╡ 59c35498-f56c-4c67-99a5-e2f2a4e2c1ec
function bif_gs_is_gf(gf)
	gs = gf
	return gs
end

# ╔═╡ 7749fd1c-e60b-4c06-9249-7174c110d590
begin
	gf = collect(range(0,stop=1,length=100))
	bif_line = zeros(size(gf))
	sn_gs_gf = zeros(size(gf))
	for i=1:length(gf)
		bif_line[i]=bif_gs_is_gf(gf[i])
		sn_gs_gf[i]=sn_gs(-40,-35,gf[i],20)
	end
	md"""Bifurcations computations for gf;gs plot with I=5, V0 =-40, Vs0=-35"""
end

# ╔═╡ 5fe72b07-4d5b-479b-b0c0-b64740ed2a33
begin
	plot_cond_reg_20 = plot(gf,sn_gs_gf, fill = (0, 0.7, :pink),linecolor=:white,linealpha=0,label="Cycle")
	#plot(gf,bif_line,label="Bifurcation for gf=gs",linewidth=2,linecolor=RGB(0.3,0.4,1))
	#plot!(gf,sn_gs_gf,label="SN",linewidth=3,linecolor=RGB(0.9,0.3,1))
	
	plot!(gf,bif_line, fill = ([sn_gs_gf], 0.7, :orange),linecolor=:white,linealpha=0,label="Bistable")
	plot!(gf,maximum(gf).*ones(size(gf)), fill = ([bif_line], 0.7, :blue),linecolor=:white,linealpha=0,label="Cycle/Instable",legend=:outertopright)
	xaxis!("gf",(0.,1))
	yaxis!("gs",(0.,1))
	title!("I=20, Vs0 = -35")
end

# ╔═╡ 2b725981-83b5-4ca2-b7c1-decf63a95eff
begin
	no_conv_cycle_up_l_h_=zeros(size(no_conv_cycle_up_l_h))
	ind_err = findall(x-> x<=11,no_conv_cycle_up_l_h)
	ind_err_ = findall(x-> x>11,no_conv_cycle_up_l_h)
	no_conv_cycle_up_l_h_[ind_err_].=no_conv_cycle_up_l_h[ind_err_]
	no_conv_cycle_up_l_h_[1]=no_conv_cycle_up_l_h_[2]
	no_conv_cycle_up_l_h_[ind_err].=maximum(no_conv_cycle_up_l_h_)
	
	no_conv_cycle_up_l_h_=no_conv_cycle_up_l_h_./100
end

# ╔═╡ bcb33fa5-682c-4c98-877f-a16bd0afae00
begin
	ind_disp=[]
	for i=1:step_gs_l_h
		if no_conv_cycle_up_l_h_[i]<gf_mesh_lowvs0_highI[i]
			if length(ind_disp)<1
				append!(ind_disp,i)
			end
		end
	end
	ind_disp_gf = findfirst(x->x>=gf_mesh_lowvs0_highI[ind_disp[1]],gf)-1
	
	limcycle = zeros(ind_disp_gf+(step_gs_l_h-ind_disp[1]+1))
	gfcycle = zeros(ind_disp_gf+(step_gs_l_h-ind_disp[1]+1))
	for i=1:ind_disp_gf+(step_gs_l_h-ind_disp[1]+1)
		if i <=ind_disp_gf
			limcycle[i]=gf[i]
			gfcycle[i]=gf[i]
		else
			limcycle[i]=no_conv_cycle_up_l_h_[i-ind_disp_gf+ind_disp[1]-1]
			gfcycle[i]=gf_mesh_lowvs0_highI[i-ind_disp_gf+ind_disp[1]-1]
		end
	end
	
	limcycle_ = zeros(ind_disp[1]+(length(gf)-ind_disp_gf+1))
	gfcycle_ = zeros(ind_disp[1]+(length(gf)-ind_disp_gf+1))
	for i=1:ind_disp[1]+(length(gf)-ind_disp_gf+1)
		if i <=ind_disp[1]
			limcycle_[i]=no_conv_cycle_up_l_h_[i]
			gfcycle_[i]=gf_mesh_lowvs0_highI[i]
		else
			limcycle_[i]=gf[i+ind_disp_gf-ind_disp[1]-1]
			gfcycle_[i]=gf[i+ind_disp_gf-ind_disp[1]-1]
		end
	end
	
	plot_cond_rest_20= plot(gfcycle,limcycle, fill = (0, 0.7, :pink),linecolor=:white,linealpha=0,label="Cycle")
	
	#plot(gf,bif_line,label="Bifurcation for gf=gs",linewidth=2,linecolor=RGB(0.3,0.4,1))
	#plot!(gf_mesh_lowvs0_highI[1:step_gs_l_h],no_conv_cycle_up_l_h_,linecolor=:black,linewidth=3,linealpha=1,label="Instable node")
	
	plot!(gf_mesh_lowvs0_highI[1:ind_disp[1]],no_conv_cycle_up_l_h_[1:ind_disp[1]], fill = ([gf_mesh_lowvs0_highI[1:ind_disp[1]]], 0.7, :blue),linecolor=:white,linealpha=0,label="Cycle/Instable",legend=:outertopright)
	plot!(gf[ind_disp_gf:end],gf[ind_disp_gf:end], fill = ([gf[ind_disp_gf-1].*ones(size(gf[ind_disp_gf:end]))], 0.5, :green),linecolor=:green,linealpha=0,label="Stable")
	plot!(gfcycle_,ones(size(gfcycle_)), fill = ([limcycle_], 0.9, :cyan),linecolor=:cyan,linealpha=0,label="Stable/Instable")
	xaxis!("gf",(0.,1))
	yaxis!("gs",(0.,1))
	title!("I=20, Vs0 = -45")
end

# ╔═╡ cd3c62b7-fdee-4c2b-b32c-5b046c700f0e
begin
	title1 = plot(title = "A. Regenerative feedback", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=16,titlefontcolor=RGB(0,0.4,0.95))
	title2 = plot(title = "B. Restorative feedback", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=16,titlefontcolor=RGB(0,0.4,0.95))
	sub_20 = plot(plot_cond_reg_20,plot_cond_rest_20,size=(850,500))
	plot(title1,title2,plot_cond_reg_20,plot_cond_rest_20,layout = @layout([[A{0.12h} B{0.12h}]; [C D]]),size=(850,500))
	#savefig("MQIF2D-cond_var_.pdf")
end

# ╔═╡ a6f3ad21-c4a0-48fc-aa8a-4a665221cac5
plot(gf_mesh_lowvs0_highI[1:step_gs_l_h]*(100),no_conv_cycle_up_l_h_,linecolor=:blue,linealpha=1,label="Cycle/Instable")

# ╔═╡ Cell order:
# ╟─d2da7b00-828f-11eb-351c-d5e58e150870
# ╠═c0d9fad0-828e-11eb-25b7-55870e15570f
# ╟─c15b4350-828f-11eb-1306-b7c44b018922
# ╟─06214530-828f-11eb-2e9c-096dba628a53
# ╟─080a5082-828f-11eb-11e4-696eafae4279
# ╟─11033850-828f-11eb-2f09-4db9ddbe1d90
# ╟─16c67090-828f-11eb-15c8-d1b6cc6a1865
# ╟─22f18580-828f-11eb-0d06-696d9c03d1b9
# ╟─2ca4bca0-828f-11eb-300e-2749156ae885
# ╟─35c0bcd0-828f-11eb-2f0a-65ddc147eedd
# ╟─39090a00-828f-11eb-08aa-a74e80f4320b
# ╟─1e28c3f0-8ef2-11eb-3567-d9ad51b92cbd
# ╠═49030140-828f-11eb-36ff-5d2834ec2f97
# ╟─4f0899b0-828f-11eb-1585-8b45c1d5ffc0
# ╟─bd92a290-828f-11eb-0bd2-89e017a3729d
# ╟─b866ef60-828f-11eb-2a12-45a3edded36e
# ╟─db8039c0-828f-11eb-1f81-c9adfd61227c
# ╟─e66e13c0-828f-11eb-2710-15b8ca48500c
# ╟─ec024cc0-828f-11eb-0043-df39ff4a0368
# ╟─f361f6f2-828f-11eb-24cf-db17dc893199
# ╠═01d1a4ae-8290-11eb-2acc-bfa854d42c67
# ╠═8ded85a0-828f-11eb-1b86-e9faa599e011
# ╟─ecf3951e-8290-11eb-38d5-07a7d2e1d074
# ╠═43662040-8290-11eb-0049-89099fa90249
# ╟─b327e5d0-8290-11eb-3bcd-491089f35abb
# ╟─34401ad0-8290-11eb-20c9-3b683651b434
# ╠═70186300-8290-11eb-1c6f-75d6cb485173
# ╠═552e6ee2-8290-11eb-2e08-190e44a2318b
# ╟─1acd0f40-8290-11eb-3799-a7b6fceba57d
# ╠═81b948f0-828f-11eb-07fa-278f0310284c
# ╠═55e642a0-828f-11eb-1583-01136f112cf0
# ╠═61d279d0-828f-11eb-001b-eb08bd724932
# ╟─29a59b80-8291-11eb-1fc0-1d4e8fee5cc1
# ╟─03a3fe42-8291-11eb-1e64-873d92c6c4c9
# ╠═a7706122-8eef-11eb-2437-2705f6dec698
# ╟─15f4b6c0-8291-11eb-204a-392f432709f3
# ╟─4a982f20-82ae-11eb-0f18-d142243d24f9
# ╟─ae142db0-8420-11eb-338c-7b9d5123aa59
# ╟─b68fb0d0-82af-11eb-3009-f92bad601035
# ╟─ac79602e-82b1-11eb-3faf-272868f80428
# ╟─6948a750-8291-11eb-1391-17e89878f727
# ╠═1880632e-841c-11eb-1dfd-6306b2bdd5a6
# ╠═db0fb6fe-8437-11eb-0786-f911dacc34a6
# ╠═ce321af0-8437-11eb-249d-2b9dc595f572
# ╠═059fc6ba-6b28-4831-af18-a7f15f6f7ecf
# ╠═a7a811a2-8428-11eb-3b12-31a9eb1bbb4e
# ╠═cbba3d20-8437-11eb-1fb8-dbe2998b1e4c
# ╠═57e45972-8438-11eb-16ec-7f914987e9f2
# ╠═e788c0a2-8291-11eb-0ad7-7f271f14ae4a
# ╠═b3f663f2-8291-11eb-3ec9-27d8e55c112b
# ╠═71be11d2-8292-11eb-15ec-7dd60814a2b1
# ╠═c8503610-8295-11eb-21d4-4980deab4482
# ╠═a4677250-8438-11eb-20e2-11a4f00cd9b5
# ╠═c9a0b680-8438-11eb-0b45-85f6673c3f11
# ╠═fca694c0-84b3-11eb-3872-730076ea9f10
# ╠═10129dd0-83ef-11eb-1c9c-f73b40d2c04d
# ╠═35d6c747-c96c-48f9-925b-67f38ab71df7
# ╠═c3a876b0-8fcc-11eb-0d3e-530aeee079b4
# ╟─6f19b812-8f0f-11eb-12fc-d9ba60894ce9
# ╟─90f164e0-8291-11eb-0f54-151c81ae3bae
# ╠═50c95000-8451-11eb-2f56-9977c8baf8e8
# ╠═ae391d20-914d-11eb-18ca-533973cac20c
# ╠═9b6f2d0e-8293-11eb-1f6f-71301a4774bb
# ╠═0a819880-8444-11eb-06bc-57b6c2bec1e7
# ╠═c8482202-8451-11eb-2192-cfecd9726aa1
# ╟─6b7737fe-84be-11eb-100e-bd880b7eefdd
# ╟─4669dea0-84e6-11eb-1ca9-137c44741630
# ╟─a662a0a0-84e9-11eb-3002-d5db70f50e88
# ╟─1ec2a3c0-844e-11eb-3e3f-eb9d83bd99b7
# ╟─a00e0990-8293-11eb-14f5-958a6f557dab
# ╟─94308d90-8442-11eb-05e8-c7dff6cc828e
# ╟─3e263a20-8443-11eb-378b-a31cd56f3611
# ╟─9522e760-8443-11eb-3c70-97d2ec7e2c79
# ╟─4d034eb0-844e-11eb-058e-9319628a79be
# ╟─4be2f61e-8444-11eb-27a8-318d0f58fec9
# ╟─146538fe-8445-11eb-3810-7fab93bfe878
# ╟─66962930-844c-11eb-0b45-316b1c7671a2
# ╟─bcd1d05e-844c-11eb-228d-1da658340995
# ╟─69bece80-844e-11eb-1039-77245c88f570
# ╟─7f4f1a40-8451-11eb-1a9c-3378bd3385df
# ╟─f6e803f0-8451-11eb-2d23-ff84295dd4a5
# ╟─4530c1a0-8452-11eb-3c24-1905fb6b0fd9
# ╟─e1214210-8452-11eb-2d21-c98ad1ab60be
# ╟─f0ffa2f0-84be-11eb-140e-c3afab500eb8
# ╟─e73ab890-84be-11eb-2e91-71c8509b17f6
# ╟─57d2e0ee-84bf-11eb-3e89-49952681af5c
# ╟─a34b8eb0-84bf-11eb-344c-0dd2bf7f836e
# ╟─1111ab50-84c0-11eb-192d-9f2bf6492011
# ╟─fbad5f7e-84e6-11eb-2195-5f36dad7ead0
# ╟─0d119520-84e7-11eb-1f47-93c154733aa7
# ╟─94bdf180-84e7-11eb-3f40-3f78aeefa66e
# ╟─40f66220-84e8-11eb-0dd9-452166c81093
# ╟─90b80430-84e8-11eb-1c30-a58116c0457f
# ╟─40a7fe90-84e9-11eb-2130-a1905ebe5385
# ╟─3b859a30-84e9-11eb-13db-034f42cd7388
# ╟─850541a0-84ea-11eb-2ee4-db89807fb982
# ╟─d0bb1f70-84ea-11eb-1597-0dad8f20b183
# ╟─4bf5e800-84eb-11eb-0c40-1980d2947538
# ╟─9d3d6f4e-8291-11eb-27d2-d9e4d8867b0a
# ╟─fb7856f0-8450-11eb-3539-9d8d1896a63d
# ╟─5a2f0880-8292-11eb-1d00-37bcb6b6d3c0
# ╟─8ef756c0-844b-11eb-04e6-cfc34ceb9b85
# ╠═a6437de0-84ff-11eb-0269-c13490ee960c
# ╟─f85b5f80-8450-11eb-256e-f34379563240
# ╟─c6ab8390-8291-11eb-1ac1-870443c01965
# ╟─303a2fc2-844d-11eb-18cd-07ba350d8651
# ╠═22e3f0f0-8500-11eb-25eb-1182e8caab4d
# ╟─08842ef0-8451-11eb-3b92-999a9bb70456
# ╠═6eb89a02-8292-11eb-311d-29eca634e7c5
# ╟─39ba1320-8453-11eb-2e3b-99e86470d3d7
# ╠═7368bbf0-8500-11eb-1685-cb2d1fb50e47
# ╟─1311c080-8451-11eb-22b0-3b2971f221d3
# ╟─d4e6b700-8295-11eb-06bf-535c44212776
# ╟─6e09429e-84c0-11eb-2a52-730ade938fc9
# ╟─df2d8a50-8500-11eb-15bc-cbf21b2725d0
# ╟─2756a7e0-8451-11eb-105b-91e1f571f79c
# ╟─f7403f70-8438-11eb-3e06-e758dc5bd7c8
# ╟─0516d180-84e9-11eb-2eee-d35a21ae7923
# ╟─a929c300-8501-11eb-0d94-eb34ede9726e
# ╟─3234f180-8451-11eb-080b-cde94036350f
# ╟─04fb3140-843b-11eb-1baa-bb5440958fb6
# ╟─b1677320-84eb-11eb-2f97-c5d24260f5fd
# ╟─0ad59390-8502-11eb-2f03-0bdc114d3620
# ╠═c66ffc6e-9a4a-4ce1-bb52-fefbb76677ab
# ╠═0f0b344a-f551-4d28-823f-ec81b6ff2f92
# ╠═59c35498-f56c-4c67-99a5-e2f2a4e2c1ec
# ╠═7749fd1c-e60b-4c06-9249-7174c110d590
# ╠═5fe72b07-4d5b-479b-b0c0-b64740ed2a33
# ╠═bcb33fa5-682c-4c98-877f-a16bd0afae00
# ╠═cd3c62b7-fdee-4c2b-b32c-5b046c700f0e
# ╠═2b725981-83b5-4ca2-b7c1-decf63a95eff
# ╠═a6f3ad21-c4a0-48fc-aa8a-4a665221cac5
