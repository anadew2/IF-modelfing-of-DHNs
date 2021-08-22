### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c0d9fad0-828e-11eb-25b7-55870e15570f
begin
	using LinearAlgebra
	using Plots
	using Statistics
	using DifferentialEquations
	using NLsolve
	using Interpolations
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

# ╔═╡ 973e0f10-9150-11eb-3acb-8da1d45e7cbf
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
	Vr = -30
	Vsr = -10
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

# ╔═╡ 10bc3f90-8d98-11eb-1ed5-274d7b9a26bb
md" #### Figure 7 from 'An organizing center in a planar model of neuronal excitability' "

# ╔═╡ 55e46520-8c17-11eb-14f5-3bf63c3d6d06
function n_inf_gd!(v,v0)
	n_inf = 2/(1+exp(-5*(v-v0)))
	return n_inf
end	

# ╔═╡ 5f27272e-8c17-11eb-250c-91229b13e59b
function d_n_inf_gd!(v,v0)
	n_inf = 10 *exp(-5*(v-v0))/(1+exp(-5*(v-v0)))^2
	return n_inf
end

# ╔═╡ d6af2520-8c10-11eb-3100-f13692ab903b
function ogd!(F, x,v0)
    F[1] = x[1]-(x[1]^3)/3-x[2]^2 +2/3
    F[2] = n_inf_gd!(x[1],v0)+x[3]-x[2]
	F[3] = -(1-x[1]^2)+2*x[2]*d_n_inf_gd!(x[1],v0)
end

# ╔═╡ 6fe92d70-8d98-11eb-1e6a-d110c6e12eda
function fp_ogd!(F, x,v0)
    F[1] = x[1]-(x[1]^3)/3-x[2]^2 +2/3
    F[2] = -(1-x[1]^2)+2*x[2]*d_n_inf_gd!(x[1],v0)
end

# ╔═╡ 5b8ebe30-8c17-11eb-0589-670031fde4c8
function d2_n_inf_gd!(v,v0)
	n_inf = -50 *(exp(-5*(v-v0)))*(1-exp(-5*(v-v0)))/(1+exp(-5*(v-v0)))^3
	return n_inf
end	

# ╔═╡ 65f84d00-8c17-11eb-20d7-e95bc1b3fbea
function j_ogd!(J, x,v0)
    J[1, 1] = 1-x[1]^2
    J[1, 2] = -2*x[2]
	J[1, 3] = 0
	J[2, 1] = d_n_inf_gd!(x[1],v0)
    J[2, 2] = -1
	J[2, 3] = 1
    J[3, 1] = 2*x[1]+2*x[2]*d2_n_inf_gd!(x[1],v0)
    J[3, 2] = 2*d_n_inf_gd!(x[1],v0)
	J[3, 3] = 0 
end

# ╔═╡ a6ca3ff0-8d98-11eb-1a05-e7728076f22f
function jfp_ogd!(J, x,v0)
    J[1, 1] = 1-x[1]^2
    J[1, 2] = -2*x[2]
    J[2, 1] = 2*x[1]+2*x[2]*d2_n_inf_gd!(x[1],v0)
    J[2, 2] = 2*d_n_inf_gd!(x[1],v0)
end

# ╔═╡ 2eaabf40-8c26-11eb-30c9-5bc94b04df0c
begin
	v0_mat2=collect(range(-1,stop=1,length=200000))
	tc_gd2 = zeros(size(v0_mat2))
	solu2=[]
	for i=1:length(v0_mat2)
		v0=v0_mat2[i]
		sol=nlsolve((F,x) ->fp_ogd!(F,x,v0),(J,x) ->jfp_ogd!(J,x,v0),[0.0,0.0])
		n0 = sol.zero[2]-n_inf_gd!(sol.zero[1],v0)
		append!(solu2,n0)
		tc_gd2[i]= -n_inf_gd!(-1,v0)
	end
end

# ╔═╡ d4cc5080-8d72-11eb-0802-ffc0d2950c5e
begin
	v0_grid = collect(minimum(v0_mat2):maximum(v0_mat2))
	n0_grid = zeros(size(v0_grid))
	for i=1:length(v0_grid)
		ind = findall(x -> abs(x -v0_grid[i])<0.01,v0_mat2)
		v0 = v0_grid[i]
		min_ind_ = v0_mat2[ind] .-v0
		min_ind_ = abs.(min_ind_)
		min_ind = findfirst(x -> x == minimum(min_ind_),min_ind_)
		n0_grid[i] = solu2[ind[min_ind]]
	end
end

# ╔═╡ 734d1710-8c26-11eb-19a7-df0f80821d16
plot(v0_mat2,solu2,legend=:outertopright)

# ╔═╡ 0fb612a0-8d76-11eb-3fd5-ffcc99e17b45
begin
	itp = interpolate(n0_grid, BSpline(Cubic(Line(OnGrid()))))
end

# ╔═╡ d06a4830-8d6d-11eb-1c25-915f534cf3da
begin 
	#my interpolation of SN curve 
	v0_mat_ = collect(1:0.1:length(v0_grid))
	v0_grid_ = collect(range(minimum(v0_grid),stop=maximum(v0_grid),length=length(v0_mat_)))
	n0_itp = zeros(size(v0_mat_))
	offset = maximum(v0_grid)+1
	for i=1:length(v0_mat_)
		n0_itp[i] = itp(v0_mat_[i])
	end
	plot(v0_mat_.-offset,n0_itp)
	plot!(v0_mat2,tc_gd2)
	xaxis!((-2,2))
	yaxis!((-2,2))
end

# ╔═╡ cc37de90-8d6c-11eb-06df-f7a4d2102de7
begin
	st = 500
	plot(v0_grid_,n0_itp)
	plot!(v0_mat2,solu2)	
	plot!(v0_mat2,tc_gd2,legend=:outertopright)
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
p=(0,-40.0,-40.0,1.0,1.0,0.5,10.0)

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

# ╔═╡ 84e31ab0-8d9d-11eb-3f1a-2d7ad1c0c76d
function sn(vs0,gs,gf,I)
	if I*(gf-gs) >=0 
		if I==0 
			v0 = vs0
		else
			v01 = vs0 + sqrt(I*(gf-gs)/(gf*gs))
			v02 = vs0 - sqrt(I*(gf-gs)/(gf*gs))
			v0 = [v01,v02]
		end
	else
		v0 = [NaN,NaN]
	end
	return v0
end

# ╔═╡ cbcc1702-8e75-11eb-30bb-5957561dbfb4
function sn_1(vs0,gs,gf,I)
	if I*(gf-gs) >=0 
		if I==0 
			v01 = vs0
		else
			v01 = vs0 + sqrt(I*(gf-gs)/(gf*gs))
		end
	else
		v01=NaN
	end
	return v01
end

# ╔═╡ d05102de-8e75-11eb-0ba6-6b3075b7d9e7
function sn_2(vs0,gs,gf,I)
	if I*(gf-gs) >=0 
		if I==0 
			v02 = vs0
		else
			v02 = vs0 - sqrt(I*(gf-gs)/(gf*gs))
		end
	else
		v02 = NaN
	end
	return v02
end

# ╔═╡ 2ca1fe80-912a-11eb-38e6-6f95fd98463a
function sn_gs(vs0,v0,gf,I)
	gs = gf*I/(I+gf*(v0-vs0)^2)
	return gs
end

# ╔═╡ 2c5306c0-8d9e-11eb-394b-37fb44d15668
function tc(vs0,gs,gf)
	return vs0
end

# ╔═╡ e0e6f892-8e60-11eb-346f-fd9642c789a9
function bif_gs_is_gf(gf)
	gs = gf
	return gs
end

# ╔═╡ 8f381aa0-8d9e-11eb-0066-6dd0be3e739c
begin
	vs0 = collect(range(-60,stop=-20,length=100))
	gs = collect(range(0,stop=1,length=100))
	gf = collect(range(0,stop=1,length=100))
	sn_vs0 = zeros(size(vs0))
	tc_vs0 = zeros(size(vs0))
	for i=1:length(vs0)
		sn_vs0[i]=sn(vs0[i],0.5,1,0)
		tc_vs0[i]=tc(vs0[i],0.5,1)
	end
	md"""Saddle-node and transcritical bifurcations computation for I=0 with Vs0"""
end

# ╔═╡ 7cd84390-8fe6-11eb-2a94-11e74ee95578
gr()

# ╔═╡ 043809f0-8d9f-11eb-39af-bd462f5d2a92
begin
	sn_tc_star_plot = plot(vs0,sn_vs0,label="SN",linewidth=3,linecolor=RGB(0,0,0))
	plot!(vs0,tc_vs0,label="TC",linewidth=2,linecolor=RGB(0.9,0.8,0))
	plot!(vs0,tc_vs0, fill = (minimum(vs0).*ones(size(vs0)), 0.7, RGB(0,0.75,0.6)),linecolor=:white,linealpha=0,label="Regenerative")
	plot!(vs0,maximum(vs0).*ones(size(vs0)), fill = (tc_vs0, 0.7, RGB(0.95,0.3,0.6)),linecolor=:white,linealpha=0,label="Restorative",legend=:outertopright)
	xaxis!("Vs0")
	yaxis!("V0")
	title!("SN & TC bifurcations for I* = 0")
end

# ╔═╡ a6359dd0-8e58-11eb-2895-330dc76cebf0
begin
	v_potential = collect(range(p[2]-60,stop=p[2]+60,step=1))
	vs_null = zeros(size(v_potential))
	vnull_star1 = zeros(size(v_potential))
	vnull_star2 = zeros(size(v_potential))
	for i=1:length(v_potential)
		vs_null[i] = Vsnullcline(v_potential[i],p)
		vnull_star1[i],vnull_star2[i]=Vnullcline_(v_potential[i],p)
	end
	
	pp_plot_star = plot(v_potential,vs_null,linecolor=RGB(0,0.7,0.4),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star,fp2_star,stab_star=fixedpoints(p)
	fp_star=[fp1_star,fp2_star]
	if length(stab_star)>1
		i=1 #we have 2 saddles 
		if stab_star[i] == 0
			scatter!([fp_star[i][1]],[fp_star[i][2]],markershape = :star8,markersize = 5,markeralpha = 0.8,markercolor = RGB(0.9,0.7,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Degenerated bifurcation",legend=:outertopright,size=(900,500))#saddle
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Phase Plane, gs=$(p[6]) & gf=$(p[5])")
end

# ╔═╡ 8b0f5570-8e83-11eb-263c-0b79ceb998be
begin
	plot(sn_tc_star_plot,pp_plot_star)
	#savefig("MQIF-bifFig7.pdf")
end

# ╔═╡ 43aa3e10-914a-11eb-1c0b-67c9acf5d861
begin
	p_gs_low = change_gs(0.5,p)
	vnull_star1_gs_low = zeros(size(v_potential))
	vnull_star2_gs_low = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_low[i],vnull_star2_gs_low[i]=Vnullcline_(v_potential[i],p_gs_low)
	end
	
	pp_plot_star_gs_low = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_low,fp2_star_gs_low,stab_star_gs_low=fixedpoints(p_gs_low)
	fp_star_gs_low=[fp1_star_gs_low,fp2_star_gs_low]
	if length(stab_star_gs_low)>1
		 #we have 2 saddles 
		if stab_star_gs_low[i] == 0
			scatter!([fp_star_gs_low[i][1]],[fp_star_gs_low[i][2]],markershape = :star8,markersize = 5,markeralpha = 0.8,markercolor = RGB(0.9,0.7,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Degenerated bifurcation",legend=:outertopright,size=(900,500))#saddle
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0=V0")
	
	md"""Plot Vs0=V0 ; gs<gf"""
end

# ╔═╡ 43e671d2-914c-11eb-281e-fb1f091b64b9
begin
	##HIGH VS0
	p_gs_low_vs0_high = change_vs0(p_gs_low[2]+1,p_gs_low)
	vnull_star1_gs_low_vs0_high = zeros(size(v_potential))
	vnull_star2_gs_low_vs0_high = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_low_vs0_high[i],vnull_star2_gs_low_vs0_high[i]=Vnullcline_(v_potential[i],p_gs_low_vs0_high)
	end
	
	pp_plot_star_gs_low_vs0_high = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_low_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_low_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_low_vs0_high,fp2_star_gs_low_vs0_high,stab_star_gs_low_vs0_high=fixedpoints(p_gs_low_vs0_high)
	fp_star_gs_low_vs0_high=[fp1_star_gs_low_vs0_high,fp2_star_gs_low_vs0_high]
	if length(stab_star_gs_low_vs0_high)>1
		for i=1:length(stab_star_gs_low_vs0_high)
			if stab_star_gs_low_vs0_high[i] == 0
				scatter!([fp_star_gs_low_vs0_high[i][1]],[fp_star_gs_low_vs0_high[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
			end
			if stab_star_gs_low_vs0_high[i] >0 
				scatter!([fp_star_gs_low_vs0_high[i][1]],[fp_star_gs_low_vs0_high[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if stab_star_gs_low_vs0_high[i] <0
				scatter!([fp_star_gs_low_vs0_high[i][1]],[fp_star_gs_low_vs0_high[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0>V0")
	
	##LOW VS0
	p_gs_low_vs0_low = change_vs0(p_gs_low[2]-1,p_gs_low)
	vnull_star1_gs_low_vs0_low = zeros(size(v_potential))
	vnull_star2_gs_low_vs0_low = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_low_vs0_low[i],vnull_star2_gs_low_vs0_low[i]=Vnullcline_(v_potential[i],p_gs_low_vs0_low)
	end
	
	pp_plot_star_gs_low_vs0_low = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_low_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_low_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_low_vs0_low,fp2_star_gs_low_vs0_low,stab_star_gs_low_vs0_low=fixedpoints(p_gs_low_vs0_low)
	fp_star_gs_low_vs0_low=[fp1_star_gs_low_vs0_low,fp2_star_gs_low_vs0_low]
	if length(stab_star_gs_low_vs0_low)>1
		for i=1:length(stab_star_gs_low_vs0_low)
			if stab_star_gs_low_vs0_low[i] == 0
				scatter!([fp_star_gs_low_vs0_low[i][1]],[fp_star_gs_low_vs0_low[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
			end
			if stab_star_gs_low_vs0_low[i] >0 
				scatter!([fp_star_gs_low_vs0_low[i][1]],[fp_star_gs_low_vs0_low[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if stab_star_gs_low_vs0_low[i] <0
				scatter!([fp_star_gs_low_vs0_low[i][1]],[fp_star_gs_low_vs0_low[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0<V0")
	md"""Plot Vs0<V0 & Vs0>V0 ; gs<gf"""
end

# ╔═╡ 54ce40f0-914b-11eb-1ee9-d3e665ab43fc
begin
	p_gs_star = change_gs(p[5],p)
	vnull_star1_gs_star = zeros(size(v_potential))
	vnull_star2_gs_star = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_star[i],vnull_star2_gs_star[i]=Vnullcline_(v_potential[i],p_gs_star)
	end
	
	pp_plot_star_gs_star = plot(vnull_star2_gs_star,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star1_gs_star,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0=V0")
	
	md"""Plot Vs0=V0 ; gs=gf"""
end

# ╔═╡ a75eef40-9150-11eb-00df-630e061d5d6f
begin
	##HIGH VS0
	p_gs_star_vs0_high = change_vs0(p_gs_star[2]+1,p_gs_star)
	vnull_star1_gs_star_vs0_high = zeros(size(v_potential))
	vnull_star2_gs_star_vs0_high = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_star_vs0_high[i],vnull_star2_gs_star_vs0_high[i]=Vnullcline_(v_potential[i],p_gs_star_vs0_high)
	end
	
	pp_plot_star_gs_star_vs0_high = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_star_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_star_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp_star_gs_star_vs0_high,stab_star_gs_star_vs0_high=fixedpoints_gs_is_gf(p_gs_star_vs0_high)
	
	if stab_star_gs_star_vs0_high == 0
		scatter!([fp_star_gs_star_vs0_high[1]],[fp_star_gs_star_vs0_high[2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	end
	if stab_star_gs_star_vs0_high >0 
		scatter!([fp_star_gs_star_vs0_high[1]],[fp_star_gs_star_vs0_high[2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	end
	if stab_star_gs_star_vs0_high <0
		scatter!([fp_star_gs_star_vs0_high[1]],[fp_star_gs_star_vs0_high[2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	end

	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0>V0")
	
	##LOW VS0
	p_gs_star_vs0_low = change_vs0(p_gs_star[2]-1,p_gs_star)
	vnull_star1_gs_star_vs0_low = zeros(size(v_potential))
	vnull_star2_gs_star_vs0_low = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_star_vs0_low[i],vnull_star2_gs_star_vs0_low[i]=Vnullcline_(v_potential[i],p_gs_star_vs0_low)
	end
	
	pp_plot_star_gs_star_vs0_low = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_star_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_star_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp_star_gs_star_vs0_low,stab_star_gs_star_vs0_low=fixedpoints_gs_is_gf(p_gs_star_vs0_low)

	if stab_star_gs_star_vs0_low == 0
		scatter!([fp_star_gs_star_vs0_low[1]],[fp_star_gs_star_vs0_low[2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
	end
	if stab_star_gs_star_vs0_low >0 
		scatter!([fp_star_gs_star_vs0_low[1]],[fp_star_gs_star_vs0_low[2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
	end
	if stab_star_gs_star_vs0_low <0
		scatter!([fp_star_gs_star_vs0_low[1]],[fp_star_gs_star_vs0_low[2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
	end

	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0<V0")
	md"""Plot Vs0<V0 & Vs0>V0 ; gs=gf"""
end

# ╔═╡ 08512612-914c-11eb-3dcc-835c3c2c1d5a
begin
	p_gs_high = change_gs(2,p)
	vnull_star1_gs_high = zeros(size(v_potential))
	vnull_star2_gs_high = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_high[i],vnull_star2_gs_high[i]=Vnullcline_(v_potential[i],p_gs_high)
	end
	
	pp_plot_star_gs_high = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_high,fp2_star_gs_high,stab_star_gs_high=fixedpoints(p_gs_high)
	fp_star_gs_high=[fp1_star_gs_high,fp2_star_gs_high]
	if length(stab_star_gs_high)>1
		 #we have 2 saddles 
		if stab_star_gs_high[i] == 0
			scatter!([fp_star_gs_high[i][1]],[fp_star_gs_high[i][2]],markershape = :star8,markersize = 5,markeralpha = 0.8,markercolor = RGB(0.9,0.7,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Degenerated bifurcation",legend=:outertopright,size=(900,500))#saddle
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0=V0")
	md"""Plot Vs0=V0 ; gs>gf"""
end

# ╔═╡ eb59b1e0-914f-11eb-1be9-8d5937b6c0cc
begin
	##HIGH VS0
	p_gs_high_vs0_high = change_vs0(p_gs_high[2]+1,p_gs_high)
	vnull_star1_gs_high_vs0_high = zeros(size(v_potential))
	vnull_star2_gs_high_vs0_high = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_high_vs0_high[i],vnull_star2_gs_high_vs0_high[i]=Vnullcline_(v_potential[i],p_gs_high_vs0_high)
	end
	
	pp_plot_star_gs_high_vs0_high = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_high_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_high_vs0_high,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_high_vs0_high,fp2_star_gs_high_vs0_high,stab_star_gs_high_vs0_high=fixedpoints(p_gs_high_vs0_high)
	fp_star_gs_high_vs0_high=[fp1_star_gs_high_vs0_high,fp2_star_gs_high_vs0_high]
	if length(stab_star_gs_high_vs0_high)>1
		for i=1:length(stab_star_gs_high_vs0_high)
			if stab_star_gs_high_vs0_high[i] == 0
				scatter!([fp_star_gs_high_vs0_high[i][1]],[fp_star_gs_high_vs0_high[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
			end
			if stab_star_gs_high_vs0_high[i] >0 
				scatter!([fp_star_gs_high_vs0_high[i][1]],[fp_star_gs_high_vs0_high[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if stab_star_gs_high_vs0_high[i] <0
				scatter!([fp_star_gs_high_vs0_high[i][1]],[fp_star_gs_high_vs0_high[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0>V0")
	
	##LOW VS0
	p_gs_high_vs0_low = change_vs0(p_gs_high[2]-1,p_gs_high)
	vnull_star1_gs_high_vs0_low = zeros(size(v_potential))
	vnull_star2_gs_high_vs0_low = zeros(size(v_potential))
	for i=1:length(v_potential)
		vnull_star1_gs_high_vs0_low[i],vnull_star2_gs_high_vs0_low[i]=Vnullcline_(v_potential[i],p_gs_high_vs0_low)
	end
	
	pp_plot_star_gs_high_vs0_low = plot(v_potential,vs_null,linecolor=RGB(0.15,0.7,1),linewidth = 2.5,label="dVs/dt = 0")
	plot!(vnull_star1_gs_high_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	plot!(vnull_star2_gs_high_vs0_low,v_potential,linecolor=RGB(0.8,0.2,0),linewidth = 2.5,label="dV/dt = 0")
	
	fp1_star_gs_high_vs0_low,fp2_star_gs_high_vs0_low,stab_star_gs_high_vs0_low=fixedpoints(p_gs_high_vs0_low)
	fp_star_gs_high_vs0_low=[fp1_star_gs_high_vs0_low,fp2_star_gs_high_vs0_low]
	if length(stab_star_gs_high_vs0_low)>1
		for i=1:length(stab_star_gs_high_vs0_low)
			if stab_star_gs_high_vs0_low[i] == 0
				scatter!([fp_star_gs_high_vs0_low[i][1]],[fp_star_gs_high_vs0_low[i][2]],markershape = :xcross,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0,0.5),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Saddle")#saddle
			end
			if stab_star_gs_high_vs0_low[i] >0 
				scatter!([fp_star_gs_high_vs0_low[i][1]],[fp_star_gs_high_vs0_low[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0.5,0.1,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Unstable")#unstable
			end
			if stab_star_gs_high_vs0_low[i] <0
				scatter!([fp_star_gs_high_vs0_low[i][1]],[fp_star_gs_high_vs0_low[i][2]],markershape = :circle,markersize = 5,markeralpha = 0.6,markercolor = RGB(0,0.9,0),markerstrokewidth = 0.1,markerstrokealpha = 0,markerstrokecolor = :black,markerstrokestyle = :dot,label="Stable")#stable
			end
		end
	end
	
	xaxis!("V",(-50,-30))
	yaxis!("Vs",(-50,-30))
	title!("Vs0<V0")
	md"""Plot Vs0<V0 & Vs0>V0 ; gs>gf"""
end

# ╔═╡ 305071a0-9158-11eb-1e0c-2b59565af016
gr()

# ╔═╡ fc363992-9152-11eb-2343-7f763d2e30aa
begin
	title1 = plot(title = "A.   gs<gs*   with   gs*=gf=1  &  gs=$(p_gs_low[6])", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=25,titlefontcolor=RGB(0,0.4,0.95))
	title2 = plot(title = "B.   gs=gs*   with   gs*=gf=1  &  gs=$(p_gs_star[6])", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=25,titlefontcolor=RGB(0,0.4,0.95))
	title3 = plot(title = "C.   gs>gs*   with   gs*=gf=1  &  gs=$(p_gs_high[6])", grid = false, showaxis = false, axis =false, bottom_margin = -50Plots.px,ticks = nothing,title_location=:left,titlefontsize=25,titlefontcolor=RGB(0,0.4,0.95))
	plot(title1, 
	pp_plot_star_gs_low_vs0_low,pp_plot_star_gs_low,pp_plot_star_gs_low_vs0_high,
	title2,	
	pp_plot_star_gs_star_vs0_low,pp_plot_star_gs_star,pp_plot_star_gs_star_vs0_high,
	title3,
	pp_plot_star_gs_high_vs0_low,pp_plot_star_gs_high,pp_plot_star_gs_high_vs0_high,
	layout = @layout([A{0.04h}; [B C D]; E{0.04h}; [F G H]; I{0.04h}; [J K L]]),size=(1500,1500),legend=:topleft)
	#savefig("MQIF-fig6ofGD.pdf")
end

# ╔═╡ 269a79a0-8e5e-11eb-2655-0589d4f6c3be
begin
	sn_vs0_Im = zeros(length(vs0),2)
	for i=1:length(vs0)
		sn_vs0_Im[i,:]=sn(vs0[i],0.5,1,5)
	end
end

# ╔═╡ b02d67f0-8e6c-11eb-1ed9-e12649f9d465
begin
	gs_05 = collect(range(0.5-0.1,stop=0.5+0.1,step=0.01))
	gs_2 = collect(range(0.01,stop=2,step=0.01))
	sn_gs_Im = zeros(length(gs),2)
	sn_gs_Im_mod_gs_2 = zeros(length(gs_2),2)
	for i=1:length(gs)
		sn_gs_Im[i,:]=sn(-35,gs[i],1,0.1)
	end
	for i=1:length(gs_2)
		sn_gs_Im_mod_gs_2[i,:]=sn(-40,gs_2[i],1,0.1)
	end
end

# ╔═╡ 87223c90-8e5e-11eb-08a8-d5a49e3c8221
begin
	plot(vs0,sn_vs0_Im[:,1],label="SN",linewidth=3,linecolor=RGB(0.45,0.6,0.9))
	plot!(vs0,sn_vs0_Im[:,2],label="SN 2",linewidth=3,linecolor=RGB(0.25,0.6,0.9))
	xaxis!("Vs0")
	yaxis!("V0")
	title!("SN bifurcations for I = 5 & gs = 0.5")
end

# ╔═╡ ce3739b0-8e6c-11eb-3d2c-4f01e9b5e0af
begin
	plot(gs,sn_gs_Im[:,1],label="SN",linewidth=3,linecolor=RGB(0.45,0.6,0.9))
	plot!(gs,sn_gs_Im[:,2],label="SN 2",linewidth=3,linecolor=RGB(0.25,0.6,0.9))
	xaxis!("gs")
	yaxis!("V0")
	title!("SN bifurcations for I = 0.1 & Vs0 = -35 & gf=1")
end

# ╔═╡ 9767a7c8-e99b-4eab-a4a5-b0d78a4c64a9
plotly()

# ╔═╡ 49c37376-595c-4230-bc2f-8e5f31e9a4ce
begin
	plot(gs_2,sn_gs_Im_mod_gs_2[:,1],label="SN, Vs0 = V0",linewidth=3,linecolor=RGB(0.45,0.6,0.9))
	plot!(gs_2,sn_gs_Im_mod_gs_2[:,2],label="SN, Vs0 = V0",linewidth=3,linecolor=RGB(0.25,0.6,0.9))
	#plot!(gs,sn_gs_Im_mod_resto[:,1],label="SN, Vs0 = V0-1",linewidth=3,linecolor=:orange)
	#plot!(gs,sn_gs_Im_mod_resto[:,2],label="SN, Vs0 = V0-1",linewidth=3,linecolor=:orange)
	xaxis!("gs")
	yaxis!("V0",(-45,-35))
	title!("SN bifurcations for I = 0.1 & gf=1")
	#savefig("SN_gs_greater_gf_Ilow.pdf")
end

# ╔═╡ 252ef1dc-db55-46ab-825f-0e57398f2d03
md"""The spike becomes more and more sharp as the current tends to 0"""

# ╔═╡ 54d6332e-912a-11eb-3c39-2d56a730c783
begin
	sn_gs_mI = zeros(length(vs0))
	sn_gs_lI = zeros(length(vs0))
	sn_gs_hI = zeros(length(vs0))
	sn_gs_nlI = zeros(length(vs0))
	for i=1:length(vs0)
		sn_gs_mI[i]=sn_gs(vs0[i],-40,1,5)
		sn_gs_lI[i]=sn_gs(vs0[i],-40,1,0.1)
		sn_gs_hI[i]=sn_gs(vs0[i],-40,1,10)
	end
end

# ╔═╡ 398a91e0-913c-11eb-2d1a-85914a781f74
begin
	plot(vs0,sn_gs_lI,label="I=0.1")
	plot!(vs0,sn_gs_mI,label="I=5")	
	plot!(vs0,sn_gs_hI,label="I=10")
	xaxis!("Vs0")
	yaxis!("gs",(0,1))
	title!("SN bifurcation in the parameter space with gf=1,V0=-40")
end

# ╔═╡ a550fda0-8e6d-11eb-19c2-23ee0b8a4369
begin
	sn_0(a,b) = sn(a,b,1,0)
	sn_1(a,b) = sn_1(a,b,1,5)
	sn_2(a,b) = sn_2(a,b,1,5)
	plot(vs0,gs,sn_0,st=:surface,c=:RdPu_3,colorbar_entry=false,xlabel="Vs0", ylabel="gs",zlabel="V0")
	plot!(vs0,gs,sn_1,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs0", ylabel="gs",zlabel="V0")
	plot!(vs0,gs,sn_2,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs0", ylabel="gs",zlabel="V0")
	title!("SN bifurcation with I = 5")
	#savefig("MQIF_SNbif_vs0_gs.pdf")
end

# ╔═╡ 019f08ae-9146-11eb-373d-f50b19ae6ccc
plotly()

# ╔═╡ bf373790-9145-11eb-14b8-9ffed091aaf8
begin
	sn_gs_n(a,b) = sn_gs(b,a,1,-1)
	sn_gs_p(a,b) = sn_gs(b,a,1,5)
	plot(vs0,vs0,sn_gs_p,st=:surface,c=:PuBu_3,colorbar_entry=false,xlabel="V0", ylabel="Vs0",zlabel="gs")
end

# ╔═╡ 334160f0-9161-11eb-2577-61cb190f8e26
gr()

# ╔═╡ cea0e170-9160-11eb-20d3-d5554e2ca6db
begin
	sn_mod1(a,b) = sn_1(a,b,1,1)
	sn_mod2(a,b) = sn_2(a,b,1,1)
	plot(vs0,gs_2,sn_mod1,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs0", ylabel="gs",zlabel="V0")
	plot!(vs0,gs_2,sn_mod2,st=:surface,c=:YlOrBr_3,colorbar_entry=false,xlabel="Vs0", ylabel="gs",zlabel="V0")
	title!("SN bifurcation with I = 1")
	#savefig("MQIF_SNbif_vs0_gs.pdf")
end

# ╔═╡ ce199600-8e60-11eb-0b3a-15afaf9f6090
begin
	bif_line = zeros(size(gf))
	sn_gs_gf = zeros(size(gf))
	for i=1:length(gf)
		bif_line[i]=bif_gs_is_gf(gf[i])
		sn_gs_gf[i]=sn_gs(-40,-35,gf[i],20)
	end
	md"""Bifurcations computations for gf;gs plot with I=5, V0 =-40, Vs0=-35"""
end

# ╔═╡ edbc4550-8e62-11eb-0e53-156d6b3000a7
begin
	plot(gf,sn_gs_gf,label="SN",linewidth=3,linecolor=RGB(0.9,0.3,1))
	plot!(gf,bif_line,label="Bifurcation for gf=gs",linewidth=2,linecolor=RGB(0.3,0.4,1))
	plot!(gf,sn_gs_gf, fill = (0, 0.7, :pink),linecolor=:white,linealpha=0,label="Cycle")
	plot!(gf,bif_line, fill = ([sn_gs_gf], 0.7, :orange),linecolor=:white,linealpha=0,label="Bistable")
	plot!(gf,maximum(gf).*ones(size(gf)), fill = ([bif_line], 0.7, :blue),linecolor=:white,linealpha=0,label="Cycle/Instable",legend=:outertopright)
	xaxis!("gf",(0.,1))
	yaxis!("gs",(0.,1))
	title!("I=5, Vs0 = -35")
end

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
# ╟─973e0f10-9150-11eb-3acb-8da1d45e7cbf
# ╟─49030140-828f-11eb-36ff-5d2834ec2f97
# ╟─4f0899b0-828f-11eb-1585-8b45c1d5ffc0
# ╟─bd92a290-828f-11eb-0bd2-89e017a3729d
# ╟─b866ef60-828f-11eb-2a12-45a3edded36e
# ╟─db8039c0-828f-11eb-1f81-c9adfd61227c
# ╟─e66e13c0-828f-11eb-2710-15b8ca48500c
# ╟─ec024cc0-828f-11eb-0043-df39ff4a0368
# ╟─f361f6f2-828f-11eb-24cf-db17dc893199
# ╟─01d1a4ae-8290-11eb-2acc-bfa854d42c67
# ╠═8ded85a0-828f-11eb-1b86-e9faa599e011
# ╟─10bc3f90-8d98-11eb-1ed5-274d7b9a26bb
# ╟─d6af2520-8c10-11eb-3100-f13692ab903b
# ╟─65f84d00-8c17-11eb-20d7-e95bc1b3fbea
# ╟─6fe92d70-8d98-11eb-1e6a-d110c6e12eda
# ╟─a6ca3ff0-8d98-11eb-1a05-e7728076f22f
# ╟─55e46520-8c17-11eb-14f5-3bf63c3d6d06
# ╟─5f27272e-8c17-11eb-250c-91229b13e59b
# ╟─5b8ebe30-8c17-11eb-0589-670031fde4c8
# ╠═2eaabf40-8c26-11eb-30c9-5bc94b04df0c
# ╠═d4cc5080-8d72-11eb-0802-ffc0d2950c5e
# ╟─734d1710-8c26-11eb-19a7-df0f80821d16
# ╠═0fb612a0-8d76-11eb-3fd5-ffcc99e17b45
# ╟─d06a4830-8d6d-11eb-1c25-915f534cf3da
# ╟─cc37de90-8d6c-11eb-06df-f7a4d2102de7
# ╟─ecf3951e-8290-11eb-38d5-07a7d2e1d074
# ╠═43662040-8290-11eb-0049-89099fa90249
# ╟─b327e5d0-8290-11eb-3bcd-491089f35abb
# ╟─34401ad0-8290-11eb-20c9-3b683651b434
# ╟─70186300-8290-11eb-1c6f-75d6cb485173
# ╟─552e6ee2-8290-11eb-2e08-190e44a2318b
# ╟─1acd0f40-8290-11eb-3799-a7b6fceba57d
# ╠═81b948f0-828f-11eb-07fa-278f0310284c
# ╠═55e642a0-828f-11eb-1583-01136f112cf0
# ╠═61d279d0-828f-11eb-001b-eb08bd724932
# ╟─29a59b80-8291-11eb-1fc0-1d4e8fee5cc1
# ╟─03a3fe42-8291-11eb-1e64-873d92c6c4c9
# ╟─15f4b6c0-8291-11eb-204a-392f432709f3
# ╟─4a982f20-82ae-11eb-0f18-d142243d24f9
# ╟─ae142db0-8420-11eb-338c-7b9d5123aa59
# ╟─b68fb0d0-82af-11eb-3009-f92bad601035
# ╟─ac79602e-82b1-11eb-3faf-272868f80428
# ╟─6948a750-8291-11eb-1391-17e89878f727
# ╟─84e31ab0-8d9d-11eb-3f1a-2d7ad1c0c76d
# ╟─cbcc1702-8e75-11eb-30bb-5957561dbfb4
# ╟─d05102de-8e75-11eb-0ba6-6b3075b7d9e7
# ╟─2ca1fe80-912a-11eb-38e6-6f95fd98463a
# ╟─2c5306c0-8d9e-11eb-394b-37fb44d15668
# ╟─e0e6f892-8e60-11eb-346f-fd9642c789a9
# ╠═8f381aa0-8d9e-11eb-0066-6dd0be3e739c
# ╠═7cd84390-8fe6-11eb-2a94-11e74ee95578
# ╟─043809f0-8d9f-11eb-39af-bd462f5d2a92
# ╟─a6359dd0-8e58-11eb-2895-330dc76cebf0
# ╟─8b0f5570-8e83-11eb-263c-0b79ceb998be
# ╟─43aa3e10-914a-11eb-1c0b-67c9acf5d861
# ╟─43e671d2-914c-11eb-281e-fb1f091b64b9
# ╟─54ce40f0-914b-11eb-1ee9-d3e665ab43fc
# ╟─a75eef40-9150-11eb-00df-630e061d5d6f
# ╟─08512612-914c-11eb-3dcc-835c3c2c1d5a
# ╟─eb59b1e0-914f-11eb-1be9-8d5937b6c0cc
# ╠═305071a0-9158-11eb-1e0c-2b59565af016
# ╠═fc363992-9152-11eb-2343-7f763d2e30aa
# ╠═269a79a0-8e5e-11eb-2655-0589d4f6c3be
# ╠═b02d67f0-8e6c-11eb-1ed9-e12649f9d465
# ╟─87223c90-8e5e-11eb-08a8-d5a49e3c8221
# ╟─ce3739b0-8e6c-11eb-3d2c-4f01e9b5e0af
# ╠═9767a7c8-e99b-4eab-a4a5-b0d78a4c64a9
# ╠═49c37376-595c-4230-bc2f-8e5f31e9a4ce
# ╟─252ef1dc-db55-46ab-825f-0e57398f2d03
# ╠═54d6332e-912a-11eb-3c39-2d56a730c783
# ╟─398a91e0-913c-11eb-2d1a-85914a781f74
# ╟─a550fda0-8e6d-11eb-19c2-23ee0b8a4369
# ╠═019f08ae-9146-11eb-373d-f50b19ae6ccc
# ╠═bf373790-9145-11eb-14b8-9ffed091aaf8
# ╠═334160f0-9161-11eb-2577-61cb190f8e26
# ╟─cea0e170-9160-11eb-20d3-d5554e2ca6db
# ╠═ce199600-8e60-11eb-0b3a-15afaf9f6090
# ╠═edbc4550-8e62-11eb-0e53-156d6b3000a7
