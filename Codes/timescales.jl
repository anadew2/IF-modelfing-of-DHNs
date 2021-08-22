### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ d5af5020-4479-11eb-0d14-37ed6b209a68
begin
	import Pkg
	Pkg.activate(mktempdir())
	md"""Import packages"""
end

# ╔═╡ 843225f0-447a-11eb-03f6-b99782f3cfd1
Pkg.add("Plots")

# ╔═╡ c73430b0-4479-11eb-0365-575c4023ece4
begin
	Pkg.add("Plots")
	using Plots
end

# ╔═╡ 114d6130-4448-11eb-1fe1-79e9307e281b
md"# Time Scales Exercise
"

# ╔═╡ bb9db670-4480-11eb-0f69-27230b104235
md"## Aguiar model"

# ╔═╡ 06f737b0-4466-11eb-1e44-132d212493f6
celsius = 36

# ╔═╡ f1865370-4465-11eb-36a3-e9fce4a32b1f
tadj = 3.0 ^ ((celsius-36)/ 10 )

# ╔═╡ 8a2e7c20-4465-11eb-07d9-012125275c94
md"### Voltage dependant functions"

# ╔═╡ cf4df32e-4465-11eb-0e63-0b4843e0106e
md"""
HH2 model : iNa,t (proportional to m^3*h) & iK (proportional to n^4)"""

# ╔═╡ a9cc3ee0-4466-11eb-14d1-6bdd94c36b61
vtraub=-55

# ╔═╡ 6ea7b650-4466-11eb-1dbb-8d72f109e1b5
function m_carac_hh2(v)
	v2 = v - vtraub 
	vh = 5
	
	a_m = 0.32 * (13-v2) / ( exp((13-v2)/4) - 1)
	b_m = 0.28 * (v2-40) / ( exp((v2-40)/5) - 1)
	tau_m_hh2 = 1 / (a_m + b_m) / tadj
	m_inf_hh2 = a_m / (a_m + b_m)
	
	return tau_m_hh2,m_inf_hh2
end 

# ╔═╡ 992e8c8e-4471-11eb-3d77-37e98ca55298
function h_carac_hh2(v)
	v2 = v - vtraub 
	vh = 5
	
	a_h = 0.128 * exp((17-v2-vh)/18)
	b_h = 4 / ( 1 + exp((40-v2-vh)/5) )
	tau_h_hh2 = 1 / (a_h + b_h) / tadj
	h_inf_hh2 = a_h / (a_h + b_h)	
	return tau_h_hh2,h_inf_hh2
end 

# ╔═╡ bf8a09a0-4471-11eb-3b66-83b22bb98343
function n_carac_hh2(v)
	v2 = v - vtraub 
	vh = 5
	
	a_n = 0.032 * (15-v2) / ( exp((15-v2)/5) - 1)
	b_n = 0.5 * exp((10-v2)/40)
	tau_n_hh2 = 1 / (a_n + b_n) / tadj
	n_inf_hh2 = a_n / (a_n + b_n)
	return tau_n_hh2,n_inf_hh2
end 

# ╔═╡ 5f4a3a40-4473-11eb-2fc4-67fb906d1425
md"""ICaL model : long-lasting high threshold calcium current (proportional to m^2) """

# ╔═╡ 07211162-4467-11eb-25f1-f51136af6f09
function vtrap(x,c)
	#Traps for 0 in denominator of rate equations
    if (abs(x/c) < 1e-6) 
       vtrap = c + x/2 
    else 
      vtrap = x / (1-exp(-x/c)) 
	end	
	return vtrap	
end

# ╔═╡ 5ca3ab50-4473-11eb-06ae-a7e1ca087792
function m_carac_ical(v)
	a = 1.6 / (1+ exp(-0.072*(v-5)))
	b = 0.02 * vtrap( -(v-1.31), 5.36)

	tau_m = 1/(a+b) / tadj
	m_inf = 1/(1+exp((v+10)/-10))
	return tau_m,m_inf
end

# ╔═╡ aefce830-4473-11eb-3788-9fb999c93b65
md"""INaP model : persistent sodium current (proportinal to m*h)"""

# ╔═╡ 7d5c7b00-4474-11eb-309d-0b82588c3202
vsm = -2

# ╔═╡ 822ceed0-4474-11eb-29b5-b74f94e96e9b
vsh = -5

# ╔═╡ 7d228d00-4479-11eb-332b-0d79441c725f
gamma=0.5

# ╔═╡ 2287eac0-4474-11eb-2199-db797510d042
function m_carac_inap(v)
	v2 = v - vtraub # convert to traub convention

	a = 0.32 * (vsm+13-v2) / ( exp((vsm+13-v2)/4) - 1)
	b = 0.28 * (vsm+v2-40) / ( exp((vsm+v2-40)/5) - 1)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)
	return tau_m,m_inf
end

# ╔═╡ 96bcc1e0-4474-11eb-3295-271e040d0e98
function h_carac_inap(v)
	v2 = v - vtraub # convert to traub convention

	a = 0.128 * exp((vsh+17-v2)/18)
	b = 4 / ( 1 + exp((vsh+40-v2)/5*gamma))
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)
	return tau_h,h_inf
end

# ╔═╡ ba1426b0-4474-11eb-1043-15338a53e007
md"##### Plot time constants"

# ╔═╡ 27174f30-447a-11eb-33e5-2150566f3823
plotly()

# ╔═╡ f6942090-4474-11eb-0880-e9d27bc3aaa2
cell_potential = collect(range(-90,stop=50,length=500))

# ╔═╡ f844cbe0-4476-11eb-276c-c1e5c9aa323b
md"""Initialization"""

# ╔═╡ b011fa50-4476-11eb-2ba8-9f409fe63b5b
begin	
	t_m_hh2 = zeros(1,length(cell_potential))
	m_i_hh2 = zeros(1,length(cell_potential))
	
	t_h_hh2 = zeros(1,length(cell_potential))
	h_i_hh2 = zeros(1,length(cell_potential))
	
	t_n_hh2 = zeros(1,length(cell_potential))
	n_i_hh2 = zeros(1,length(cell_potential))
end 	

# ╔═╡ 12d4c3b0-4478-11eb-04e5-95694ad7bbd0
begin	
	t_m_ical = zeros(1,length(cell_potential))
	m_i_ical = zeros(1,length(cell_potential))
end 

# ╔═╡ 21f2d9e0-4478-11eb-2a77-0980b0ddc553
begin	
	t_m_inap = zeros(1,length(cell_potential))
	m_i_inap = zeros(1,length(cell_potential))
	
	t_h_inap = zeros(1,length(cell_potential))
	h_i_inap = zeros(1,length(cell_potential))
end 

# ╔═╡ 579e58ee-4480-11eb-18de-738026386048
md"""Time constants and gate variables computation"""

# ╔═╡ 1fc33a90-4476-11eb-1566-bd559f113176
for i = 1:length(cell_potential)
	v=cell_potential[i]
	#HH2 characteristics
	t_m_hh2[i],m_i_hh2[i]=m_carac_hh2(v)
	t_h_hh2[i],h_i_hh2[i]=h_carac_hh2(v)
	t_n_hh2[i],n_i_hh2[i]=n_carac_hh2(v)
	#ICaL characteristics
	t_m_ical[i],m_i_ical[i]=m_carac_ical(v)
	#INaP characteristics
	t_m_inap[i],m_i_inap[i]=m_carac_inap(v)
	t_h_inap[i],h_i_inap[i]=h_carac_inap(v)
end

# ╔═╡ 6b74a000-4480-11eb-3ba2-9b12b88b7fd2
md"""Results"""

# ╔═╡ 22fba510-447d-11eb-217e-5f3e38cffd57
begin
	plot(cell_potential,t_m_hh2[:],label="tau_m HH2")
	plot!(cell_potential,t_h_hh2[:],label="tau_h HH2")
	plot!(cell_potential,t_n_hh2[:],label="tau_n HH2")
	plot!(cell_potential,t_m_ical[:],label="tau_m ICaL")
	plot!(cell_potential,t_m_inap[:],label="tau_m INaP")
	plot!(cell_potential,t_h_inap[:],label="tau_h INaP")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Time constant (ms)")
end

# ╔═╡ 9c366270-447e-11eb-28f8-b7b53a13b721
begin
	plot(cell_potential,m_i_hh2[:],label="m_inf HH2")
	plot!(cell_potential,h_i_hh2[:],label="h_inf HH2")
	plot!(cell_potential,n_i_hh2[:],label="n_inf HH2")
	plot!(cell_potential,m_i_ical[:],label="m_inf ICaL")
	plot!(cell_potential,m_i_inap[:],label="m_inf INaP")
	plot!(cell_potential,h_i_inap[:],label="h_inf INaP")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Gate variable at rest")
end

# ╔═╡ e3ba6b80-4480-11eb-3b8a-6d058c216aef
md"""Comments"""

# ╔═╡ f4ba3a00-4480-11eb-2f52-ed18ff095332
md"""ICaL and the n gate variable in HH2 (corresponding to the current iK) seems to be active at a lower time scale that the other voltage dependant current. Indeed, we se that the conductances m for ICaL and n for IK in HH2 vary more slowly than others"""

# ╔═╡ a06dfb00-4465-11eb-2c4e-358aa670830b
md"### Calcium dependant functions"

# ╔═╡ ada102e0-448d-11eb-35e5-bf8a70668c47
taumin	= 0.1	 # (ms)	

# ╔═╡ d1b91510-448c-11eb-2f47-0b687bc0badf
md"""ICaAN : calcium activated nonspecific cationic current (proportional to m^2)"""

# ╔═╡ 4e7e2720-448d-11eb-159b-c3fabc4730f0
cac_ican=0.0005

# ╔═╡ 65e3af70-448d-11eb-06d3-55698316fdfa
beta	= 2.0*10^(-3) #(1/ms)	

# ╔═╡ 88af27f0-448d-11eb-08d9-7f02468a867f
tau_factor = 40  

# ╔═╡ dec8215e-448c-11eb-2b4f-ef8d76a86653
function m_carac_ican(cai)
	alpha2 = beta * (cai/cac_ican)^2
	
	tau_m = tau_factor / (alpha2 + beta) / tadj		
	m_inf = alpha2 / (alpha2 + beta)
	
	if(tau_m < taumin) 
		tau_m = taumin 
	end	
	
	return tau_m,m_inf
end

# ╔═╡ cb7593ce-448d-11eb-2e9d-41d0bb8267a0
md"""IKCa : calcium-dependent potassium current (proportional to m^3)"""

# ╔═╡ eb36d990-448d-11eb-09ea-09e3fc202c2a
cac_ikca = 0.001

# ╔═╡ 6b0e34b0-448e-11eb-2ab2-73fdcc4aa987
function m_carac_ikca(cai)
	car = (cai/cac_ikca)^2
	
    m_inf = car / ( 1 + car )      
    tau_m =  1 / beta / (1 + car) / tadj
	
    if(tau_m < taumin) 
		tau_m = taumin  
	end
	
	return tau_m,m_inf
end

# ╔═╡ b7eba290-448e-11eb-3647-7f00734f6290
md"##### Plot time constants"

# ╔═╡ cbebdfd0-448e-11eb-2479-65c6f1cdae49
cell_Ca = collect(range(10^(-6),stop=0.01,step=10^(-6)))

# ╔═╡ 488bddb0-448f-11eb-1341-f5a1d2b89a61
md"""Initialization"""

# ╔═╡ 58465af0-448f-11eb-1d26-d36bc242e859
begin	
	t_m_ican = zeros(1,length(cell_Ca))
	m_i_ican = zeros(1,length(cell_Ca))
end 

# ╔═╡ 74daa3ae-448f-11eb-0564-29ca38a575c9
begin	
	t_m_ikca = zeros(1,length(cell_Ca))
	m_i_ikca = zeros(1,length(cell_Ca))
end 

# ╔═╡ dc6ac90e-448f-11eb-0516-a7069fa36fd5
for i = 1:length(cell_Ca)
	cai=cell_Ca[i]
	#ICaAN characteristics
	t_m_ican[i],m_i_ican[i]=m_carac_ican(cai)
	#IKCa characteristics
	t_m_ikca[i],m_i_ikca[i]=m_carac_ikca(cai)
end

# ╔═╡ d197a302-448f-11eb-0647-9940ef920759
md"""Time constants and gate variables computation"""

# ╔═╡ 3e5bb8f0-4490-11eb-1e58-37baf741a15d
md"""Results"""

# ╔═╡ 42a4fb60-4490-11eb-1d7a-4369317b36c5
begin
	plot(cell_Ca,t_m_ican[:],label="tau_m ICaAN")
	plot!(cell_Ca,t_m_ikca[:],label="tau_m IKCa")
	plot!(xlabel="Intracellular calcium concentration (mM)")
	plot!(ylabel="Time constant (ms)")
end

# ╔═╡ 767f84a0-4490-11eb-332a-850d0f90eafc
begin
	plot(cell_Ca,m_i_ican[:],label="m_inf ICaAN")
	plot!(cell_Ca,m_i_ikca[:],label="m_inf IKCa")
	plot!(xlabel="Intracellular calcium concentration (mM)")
	plot!(ylabel="Gate variable at rest")
end

# ╔═╡ bc688d00-5d99-11eb-0d12-9febd1d6c3da
md"#### Interpretation"

# ╔═╡ d2bcf370-5d99-11eb-1809-55bcf9ca2308
md""" There are 4 time scales involved in the Aguiar model. """

# ╔═╡ 55f8f600-5d9d-11eb-3ac8-fbe006345392
md"""The first and fastest timescale is composed of the sodium currents INa,t and INa,p. Indeed, each of these two current is porportional to a gate variable which time constant is close to 0. The gating variable change induced in the corresponding ion channels by a small cell voltage perturbation is thus very sensitive. """

# ╔═╡ c79ac630-5d9d-11eb-16fe-1dafe5e85c5d
md""" The second time scale is composed of the potassium current IK and of the calcium current from L type channels ICaL. In fact, each of these two current depends on small but non zero time constants. These two current are thus slower than sodium currents. Also, we can separate these two currents from the other ones that are voltage dependent as the slope of their respecting gating variables at rest are smoother than the others, which makes these channels open slowlier and these current increasing therefore more slowly. """

# ╔═╡ defd4bd0-5d9e-11eb-10ff-9d05a5b09b66
md"""The two last timescales contain calcium dependant currents : IKCa the fastest and ICaAN the slowest. They are separated from each other and from each other timescales because of the value of their time constants, which is significantly higher.
Also, as they are calcium dependant, their response directly depend on the response of the current ICaL, which modifies the intracellular calcium concentration. Consequently, even if the slope of the gating variables at rest seems to be high, one should bare in mind that increasing the calcium concentration is a slow process. When the cell is at rest, the calcium concentration is close to 0. Then, after stimulation, an inward calcium current is created (ICaL). The calcium concentration increases. The time constants of the two calcium dependent current decreases (making them more sensitive to calcium concentration variations). The two considered currents are thus creating progressively an effet on the cell. """

# ╔═╡ ab133a50-60a0-11eb-0a3b-2724907a7773
begin
	pVa = plot(cell_potential,t_m_hh2[:],label="m HH2")
	pVa = plot!(cell_potential,t_h_hh2[:],label="h HH2")
	pVa = plot!(cell_potential,t_n_hh2[:],label="n HH2")
	pVa = plot!(cell_potential,t_m_ical[:],label="m ICaL")
	pVa = plot!(cell_potential,t_m_inap[:],label="m INap")
	pVa = plot!(cell_potential,t_h_inap[:],label="h INap",yaxis=:log,legend=:right)
	xaxis!("Membrane potential (mV)")
	
	pcaa = plot(cell_Ca,t_m_ican[:],label="m ICaAN")
	pcaa = plot!(cell_Ca,t_m_ikca[:],label="m IKCa",legend=:bottomright)
	xaxis!("Calcium concentration (mM)")
	
	pfa = plot(pVa,pcaa,layout=(1,2),linewidth = 1.5,yaxis=:log)
	#size=(1200,600)
	yaxis!("Time constant (ms)",(0.02,14000))
	
	#savefig("AguiarTimeConstants.pdf")
end

# ╔═╡ 8bc739b0-60a2-11eb-19b0-e9535a19e66d
begin
	ppVa = plot(cell_potential,m_i_hh2[:],label="m HH2")
	ppVa = plot!(cell_potential,h_i_hh2[:],label="h HH2")
	ppVa = plot!(cell_potential,n_i_hh2[:],label="n HH2")
	ppVa = plot!(cell_potential,m_i_ical[:],label="m ICaL")
	ppVa = plot!(cell_potential,m_i_inap[:],label="m INap")
	ppVa = plot!(cell_potential,h_i_inap[:],label="h INap",legend=:right)
	xaxis!("Membrane potential (mV)")
	yaxis!(ylabel="Gate variable at rest",(-0.05,1.05))

	ppcaa = plot(cell_Ca,m_i_ican[:],label="m ICaAN")
	ppcaa = plot!(cell_Ca,m_i_ikca[:],label="m IKCa",legend=:right)
	xaxis!("Calcium concentration (mM)")
	yaxis!(ylabel="Gate variable at rest",(-0.05,1.05))
	
	ppfa = plot(ppVa,ppcaa,layout=(1,2),linewidth = 1.5,size=(1200,600))
	#size=(1200,600)
	yaxis!(ylabel="Gate variable at rest",(-0.05,1.05))
	
	savefig("AguiarTimeGates.pdf")
end

# ╔═╡ ee219840-449f-11eb-14be-610f98fdbc83
md"## LeFranc model"

# ╔═╡ eef4edd0-449f-11eb-2411-d9b1b0309238
md"### Voltage dependant functions"

# ╔═╡ efc27700-449f-11eb-1997-93d5df3ffd66
md"""HH2 mdoel : iNa (proportional to m^3*h) & iKDR (proportional to n^4)"""

# ╔═╡ 1d50f380-44a1-11eb-1376-b39d73fe2777
vtraub_lf = -63 #(mV)

# ╔═╡ 5a5deb20-44a1-11eb-07d4-8b9ba18fde66
function h_carac_hh2_lf(v)
	v2 = v - vtraub_lf
	
	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)
	
	return tau_h,h_inf
end 

# ╔═╡ b67ed17e-44a1-11eb-1eee-6dd53ae70ea2
function vtrap_lf(x,y)
	if (abs(x/y) < 1e-6) 
		vtrap = y*(1 - x/y/2)
	else
		vtrap = x/(exp(x/y)-1)
	end
end

# ╔═╡ 3469d7e0-44a0-11eb-09ce-69d2b842c9e3
function m_carac_hh2_lf(v)
	v2 = v - vtraub_lf
	
	a = 0.32 * vtrap_lf(13-v2, 4)
	b = 0.28 * vtrap_lf(v2-40, 5)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)
	
	return tau_m,m_inf
end 

# ╔═╡ 952586a0-44a1-11eb-15de-a31f19d44c61
function n_carac_hh2_lf(v)
	v2 = v - vtraub_lf
	
	a = 0.032 * vtrap_lf(15-v2, 5)
	b = 0.5 * exp((10-v2)/40)
	tau_n = 1 / (a + b) / tadj
	n_inf = a / (a + b)
	
	return tau_n,n_inf
end 

# ╔═╡ fa1bee50-44a1-11eb-3749-732d1fed88c0
md"""IKir model : K+ ion current associated to metabotropic receptors (G protein dependent potassium channel Kir3) (proportional to m)"""

# ╔═╡ 232e4e00-44a2-11eb-0520-19a81bce5059
vhalf=-65 #(mV)

# ╔═╡ 227ead10-44a2-11eb-110c-4907e2b9ae60
zslope=10 #(mV)

# ╔═╡ 49ba29e0-44a2-11eb-1e4d-db081f25bf5b
taukir=10 #(ms)

# ╔═╡ 0ff07660-44a2-11eb-2e63-e5a89b4d2725
function m_carac_ikir(v)
	minf = 1/(1+ exp((v-vhalf)/zslope))
    taum = taukir
	
	return taum,minf 
end

# ╔═╡ 6b25c440-44a2-11eb-3458-75e5693d6904
md"""IfCaL model : fast calcium current from L-type channels characterized by a fast activation and slow inactivation (propotional to oca^2*cca) """ #ILslow2.mod

# ╔═╡ b814eff2-44a3-11eb-30b9-77cf2c48a9ce
zshift=0

# ╔═╡ be93d712-44a3-11eb-24a4-d9f9b6e39265
zpente = 0

# ╔═╡ 7a672470-44a3-11eb-1b0b-011a267dbe89
function cca_carac_ifcal(v)
	cca_ss = 1 / (1+exp((v+(14+zshift))/(zpente+4.03))) 
	cca_tau = 1500
	
	return cca_tau,cca_ss
end

# ╔═╡ ff4e6d60-44a3-11eb-27c7-e7035a6cabef
md"""IsCaL model : slow calcium current from L-type channels (proportional to oca*cca)""" #Ilsuf2.mod

# ╔═╡ 6235a74e-44a8-11eb-071b-27aa3cb81ee2
caslip = 4

# ╔═╡ 96167980-44a5-11eb-3d5d-338b3a909740
function cca_carac_iscal(v)
	cca_ss = 1 / (1+exp((v+(14+zshift))/(zpente+4.03))) 
	cca_tau = 10000
	
	return cca_tau,cca_ss
end

# ╔═╡ 57d96650-44a5-11eb-1368-a1544d547929
function efun(z)
	if (abs(z) < 1e-4) 
		efun = 1 - z/2
	else
		efun = z/(exp(z) - 1)
	end
end

# ╔═╡ 88ff8550-44a2-11eb-0cf5-ad6414e90de5
function oca_carac_ifcal(v)
	oca_ss = -1.12402639e-3+(1.00684407/(1+exp(-(v-(-14.44797322))/3.17052882))^0.54421057)
	
	taufactor=0.5
		
	if (oca_ss>1) 	
		oca_ss=1
	else 
		if (v<=-40) 
			oca_ss=0
		end
	end
	
	v = v+65
	a = 1*efun(.1*(25-v))
	b = 4*exp(-v/18)
	oca_tau = taufactor/(a + b)
	
	return oca_tau,oca_ss
end

# ╔═╡ e83900e0-44a3-11eb-1e9c-ab689c827bc1
function oca_carac_iscal(v)
	oca_ss = -4.77079827e-3+(1.02572113/(1+exp(-(v-(-21.4565198+1))/caslip))^0.47307403)
		
	if (oca_ss<0) 			
		oca_ss=0
	elseif (oca_ss>1)
		oca_ss=1
	elseif (v<=-60) 
		oca_ss=0
	end
	
	taufactor = 160
	
	v = v+65
	a = 1*efun(.1*(25-v))
	b = 4*exp(-v/18)
	oca_tau = (taufactor/(a + b))
	
	
	return oca_tau,oca_ss
end

# ╔═╡ cd1a0820-44a5-11eb-0dac-4dedff66eba6
md"##### Plot time constants"

# ╔═╡ f30955e0-44a5-11eb-0e4f-2523d8c038e7
md"""Initialization"""

# ╔═╡ fa702c00-44a5-11eb-1e87-695643d2f534
begin	
	t_m_hh2_lf = zeros(1,length(cell_potential))
	m_i_hh2_lf = zeros(1,length(cell_potential))
	
	t_h_hh2_lf = zeros(1,length(cell_potential))
	h_i_hh2_lf = zeros(1,length(cell_potential))
	
	t_n_hh2_lf = zeros(1,length(cell_potential))
	n_i_hh2_lf = zeros(1,length(cell_potential))
end 	

# ╔═╡ 2bc5df20-44a6-11eb-18c4-11fb7e77f8ae
begin	
	t_m_ikir_lf = zeros(1,length(cell_potential))
	m_i_ikir_lf = zeros(1,length(cell_potential))
end 

# ╔═╡ 40798de0-44a6-11eb-1b01-636fc5cf76ce
begin	
	t_oca_ifcal_lf = zeros(1,length(cell_potential))
	oca_i_ifcal_lf = zeros(1,length(cell_potential))
	
	t_cca_ifcal_lf = zeros(1,length(cell_potential))
	cca_i_ifcal_lf = zeros(1,length(cell_potential))
end 

# ╔═╡ 6ab97110-44a6-11eb-0731-bb27efbee0c4
begin	
	t_oca_iscal_lf = zeros(1,length(cell_potential))
	oca_i_iscal_lf = zeros(1,length(cell_potential))
	
	t_cca_iscal_lf = zeros(1,length(cell_potential))
	cca_i_iscal_lf = zeros(1,length(cell_potential))
end 

# ╔═╡ 855c4a10-44a6-11eb-316b-91fcbfee76c7
md"""Time constants and gate variables computation"""

# ╔═╡ 8f4d4d80-44a6-11eb-2946-9dc6120506ea
for i = 1:length(cell_potential)
	v=cell_potential[i]
	#HH2 characteristics
	t_m_hh2_lf[i],m_i_hh2_lf[i]=m_carac_hh2_lf(v)
	t_h_hh2_lf[i],h_i_hh2_lf[i]=h_carac_hh2_lf(v)
	t_n_hh2_lf[i],n_i_hh2_lf[i]=n_carac_hh2_lf(v)
	#IKir characteristics
	t_m_ikir_lf[i],m_i_ikir_lf[i]=m_carac_ikir(v)
	#IfCaL characteristics
	t_oca_ifcal_lf[i],oca_i_ifcal_lf[i]=oca_carac_ifcal(v)
	t_cca_ifcal_lf[i],cca_i_ifcal_lf[i]=cca_carac_ifcal(v)
	#IsCaL characteristics
	t_oca_iscal_lf[i],oca_i_iscal_lf[i]=oca_carac_iscal(v)
	t_cca_iscal_lf[i],cca_i_iscal_lf[i]=cca_carac_iscal(v)
end

# ╔═╡ 67744c60-602a-11eb-3ad7-e725cc91a409
gr()

# ╔═╡ eba6ee80-602e-11eb-2dda-e16185ae894e
plotly()

# ╔═╡ 9a572910-44a8-11eb-3019-bb04b52c034c
begin
	plot(cell_potential,t_m_hh2_lf[:],label="tau_m HH2",yaxis=:log)
	plot!(cell_potential,t_h_hh2_lf[:],label="tau_h HH2",yaxis=:log)
	plot!(cell_potential,t_n_hh2_lf[:],label="tau_n HH2",yaxis=:log)
	plot!(cell_potential,t_m_ikir_lf[:],label="tau_m IKir",yaxis=:log)
	plot!(cell_potential,t_oca_ifcal_lf[:],label="tau_oca IfCaL",yaxis=:log)
	plot!(cell_potential,t_cca_ifcal_lf[:],label="tau_cca IfCaL",yaxis=:log)
	plot!(cell_potential,t_oca_iscal_lf[:],label="tau_oca IsCaL",yaxis=:log)
	plot!(cell_potential,t_cca_iscal_lf[:],label="tau_cca IsCaL",yaxis=:log)
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Time constant (ms)")
end

# ╔═╡ 6eb04660-44a9-11eb-0ab0-2db3e39d1802
begin
	plot(cell_potential,t_m_hh2_lf[:],label="tau_m HH2")
	plot!(cell_potential,t_h_hh2[:],label="tau_h HH2")
	plot!(cell_potential,t_h_hh2_lf[:],label="tau_n HH2")
	plot!(cell_potential,t_m_ikir_lf[:],label="tau_m IKir")
	plot!(cell_potential,t_oca_ifcal_lf[:],label="tau_oca IfCaL")
	plot!(cell_potential,t_oca_iscal_lf[:],label="tau_oca IsCaL")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Time constant (ms)")
end

# ╔═╡ 0906bb52-44a9-11eb-258e-e3073a645dd6
begin
	plot(cell_potential,m_i_hh2_lf[:],label="m_inf HH2")
	plot!(cell_potential,h_i_hh2_lf[:],label="h_inf HH2")
	plot!(cell_potential,n_i_hh2_lf[:],label="n_inf HH2")
	plot!(cell_potential,m_i_ikir_lf[:],label="m_inf IKir")
	plot!(cell_potential,oca_i_ifcal_lf[:],label="oca_ss IfCaL")
	plot!(cell_potential,cca_i_ifcal_lf[:],label="cca_ss IfCaL")
	plot!(cell_potential,oca_i_iscal_lf[:],label="oca_ss IsCaL")
	plot!(cell_potential,cca_i_iscal_lf[:],label="cca_ss IsCaL")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Gate variable at rest")
end

# ╔═╡ b0c620de-5c2f-499b-93aa-48ea4e62cecd
gr()

# ╔═╡ 8eae4837-e47e-4e7a-8562-8da0f03ca9f3
begin
	hh_ss = plot(cell_potential,m_i_hh2_lf[:],label="m_inf")
	plot!(cell_potential,h_i_hh2_lf[:],label="h_inf")
	plot!(cell_potential,n_i_hh2_lf[:],label="n_inf")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Gate variable at rest")
	
	hh_tau = plot(cell_potential,t_m_hh2_lf[:],label="tau_m ")
	plot!(cell_potential,t_h_hh2[:],label="tau_h ")
	plot!(cell_potential,t_h_hh2_lf[:],label="tau_n ")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Time constant (ms)")
	
	plot(hh_ss,hh_tau,layout=(1,2),size=(700,350))
	#savefig("HH_ss_tau.pdf")
end

# ╔═╡ d3cf44b2-44a9-11eb-28b7-b9c4724fcbb1
begin
	plot(cell_potential,oca_i_ifcal_lf[:],label="oca_ss IfCaL")
	plot!(cell_potential,cca_i_ifcal_lf[:],label="cca_ss IfCaL")
	plot!(cell_potential,oca_i_iscal_lf[:],label="oca_ss IsCaL")
	plot!(cell_potential,cca_i_iscal_lf[:],label="cca_ss IsCaL")
	plot!(xlabel="Membrane potential (mV)")
	plot!(ylabel="Gate variable at rest")
end

# ╔═╡ 50a11ae0-5cba-11eb-20fd-a73fea2e3647
md"### Calcium dependant functions"

# ╔═╡ 5581bb50-5cba-11eb-37bb-ef8f86d6d52f
md"""ICaAN : calcium activated nonspecific cationic current (proportional to m^2)"""

# ╔═╡ 6a45bdbe-5cba-11eb-1f20-830ac1a98e24
function m_carac_ican_lf(cai)
	cac_ican = 0.9 
	beta = 0.002
	
	alpha2 = beta * (cai/cac_ican)^2
	
    m_inf = alpha2 / ( alpha2 + beta )      
    tau_m =  1000 #value from .ses file
	
	return tau_m,m_inf
end

# ╔═╡ f8636cae-5cba-11eb-049f-156d2214e577
md"""ISK : calcium activated K channels (proportional to o)"""

# ╔═╡ 1953bc90-5cbb-11eb-13ea-f9f3faf51e23
function o_carac_isk(cai)
	k = 0.05
	
    o_inf = cai / ( cai + k )      
    tau_o = 10
	
	return tau_o,o_inf
end

# ╔═╡ 7461b6de-5cbc-11eb-2c80-8d58d6b1e469
md"##### Plot time constants"

# ╔═╡ a88915b0-5cbe-11eb-1055-d18d49d2da93
cell_Ca_lf = collect(range(10^(-6),stop=0.5,step=10^(-4)))

# ╔═╡ 75470dd0-5cbc-11eb-124c-13a7a088aad9
begin	
	t_m_ican_lf = zeros(1,length(cell_Ca_lf))
	m_i_ican_lf = zeros(1,length(cell_Ca_lf))
end 

# ╔═╡ 998688b0-5cbc-11eb-3f57-53ba68d6c5c4
begin	
	t_o_isk = zeros(1,length(cell_Ca_lf))
	o_i_isk = zeros(1,length(cell_Ca_lf))
end 

# ╔═╡ a0b7f942-6025-11eb-1597-55d8b267bd0d
begin
	pV = plot(cell_potential,t_m_hh2_lf[:],label="m HH2")
	pV = plot!(cell_potential,t_h_hh2_lf[:],label="h HH2")
	pV = plot!(cell_potential,t_n_hh2_lf[:],label="n HH2")
	pV = plot!(cell_potential,t_m_ikir_lf[:],label="m IKir")
	pV = plot!(cell_potential,t_oca_ifcal_lf[:],label="oca IfCaL")
	pV = plot!(cell_potential,t_cca_ifcal_lf[:],label="cca IfCaL")
	pV = plot!(cell_potential,t_oca_iscal_lf[:],label="oca IsCaL")
	pV = plot!(cell_potential,t_cca_iscal_lf[:],label="cca IsCaL",yaxis=:log,legend=true)
	xaxis!("Membrane potential (mV)")
	
	pca = plot(cell_Ca_lf,t_m_ican_lf[:],label="m ICaAN")
	pca = plot!(cell_Ca_lf,t_o_isk[:],label="o ISK",legend=:bottomright)
	xaxis!("Calcium concentration (mM)")
	
	pf = plot(pV,pca,layout=(1,2),linewidth = 1.5,yaxis=:log,size=(1200,600))
	yaxis!("Time constant (ms)",(0.02,14000))
	
	savefig("LeFrancTimeConstants.pdf")
end

# ╔═╡ 9165e450-6031-11eb-2766-ad938e7ae86b
begin
	pV1 = plot(cell_potential,m_i_hh2_lf[:],label="m HH2")
	pV1 = plot!(cell_potential,h_i_hh2_lf[:],label="h HH2")
	pV1 = plot!(cell_potential,n_i_hh2_lf[:],label="n HH2")
	pV1 = plot!(cell_potential,m_i_ikir_lf[:],label="m IKir")
	pV1 = plot!(cell_potential,oca_i_ifcal_lf[:],label="oca IfCaL")
	pV1 = plot!(cell_potential,cca_i_ifcal_lf[:],label="cca IfCaL")
	pV1 = plot!(cell_potential,oca_i_iscal_lf[:],label="oca IsCaL")
	pV1 = plot!(cell_potential,cca_i_iscal_lf[:],label="cca IsCaL")
	pV1 = plot!(xlabel="Membrane potential (mV)")
	pV1 = plot!(ylabel="Gate variable at rest",legend=:right,size=(400,600))
	yaxis!((-0.05,1.05))

	pca1 = plot(cell_Ca_lf,m_i_ican_lf[:],label="mf ICaAN")
	pca1 = plot!(cell_Ca_lf,o_i_isk[:],label="o ISK")
	pca1 = plot!(xlabel="Calcium concentration (mM)")
	pca1 = plot!(ylabel="Gate variable at rest",legend=:right)
	yaxis!((-0.05,1.05))
	
	pf1 = plot(pV1,pca1,layout=(1,2),linewidth = 1.5,size=(1200,600))
	
	savefig("LeFrancGate.pdf")
end

# ╔═╡ 9a6e50a0-5cbc-11eb-153f-8977b402aa2e
for i = 1:length(cell_Ca_lf)
	cai=cell_Ca_lf[i]
	#ICaAN characteristics
	t_m_ican_lf[i],m_i_ican_lf[i]=m_carac_ican_lf(cai)
	#ISK characteristics
	t_o_isk[i],o_i_isk[i]=o_carac_isk(cai)
end

# ╔═╡ f1031450-5cbc-11eb-3317-1f40ea067833
begin
	plot(cell_Ca_lf,t_m_ican_lf[:],label="tau_m ICaAN")
	plot!(cell_Ca_lf,t_o_isk[:],label="tau_o ISK")
	plot!(xlabel="Intracellular calcium concentration (mM)")
	plot!(ylabel="Time constant (ms)")
end

# ╔═╡ 1809c580-5cbd-11eb-0c4e-6171dedd4b55
begin
	plot(cell_Ca_lf,m_i_ican_lf[:],label="m_inf ICaAN")
	plot!(cell_Ca_lf,o_i_isk[:],label="o_inf ISK")
	plot!(xlabel="Intracellular calcium concentration (mM)")
	plot!(ylabel="Gate variable at rest")
end

# ╔═╡ e0fdd9e0-5e66-11eb-0ed2-51ad9550e89b
md"#### Calculations"

# ╔═╡ edaf35d0-5e66-11eb-0fff-d5f2e1786452
md"""To compare the time constant of ion currents with multiple gate variables, we will calculate for each of them the maximum of the product of the gate variables time constants. Indeed, if we mathematically define a new gate variable (with less biological meaning than the ones used initially), each current can be compared based on a single value of gate variable time constant. """

# ╔═╡ 31632f50-5fc2-11eb-3795-43129912495f
gr()

# ╔═╡ 205ae9de-5f55-11eb-2dd7-8903f2852577
time=collect(range(0,5,length=1000))

# ╔═╡ 8c36ad90-5f53-11eb-27e3-11ba5de25c5c
function m_time_HH2(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	m=zeros(size(t))
	dt=time[2]-time[1]
	tau_m,m_inf=m_carac_hh2_lf(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			m[i]=0
		end
		if i>1
			m_exp=1-exp(-time[i]/tau_m)
			m[i]=m[i]+m_exp*(m_inf-m[i])	
		end
	end
	
		
	return m
end

# ╔═╡ 36456390-5f58-11eb-2a12-b9a2c4062643
function h_time_HH2(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	h=zeros(size(t))
	dt=time[2]-time[1]
	tau_h,h_inf=h_carac_hh2_lf(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			h[i]=0
		end
		if i>1
			h_exp=1-exp(-time[i]/tau_h)
			h[i]=h[i]+h_exp*(h_inf-h[i])	
		end
	end
	
		
	return h
end

# ╔═╡ 64e84830-5f57-11eb-018f-b9f95f04aa1b
function n_time_HH2(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	n=zeros(size(t))
	dt=time[2]-time[1]
	tau_n,n_inf=n_carac_hh2_lf(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			n[i]=0
		end
		if i>1
			n_exp=1-exp(-time[i]/tau_n)
			n[i]=n[i]+n_exp*(n_inf-n[i])	
		end
	end
	
		
	return n
end

# ╔═╡ 1410dbd0-5fc4-11eb-3131-554175fe60a2
function m_time_IKIR(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	m=zeros(size(t))
	dt=time[2]-time[1]
	tau_m,m_inf=m_carac_ikir(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			m[i]=0
		end
		if i>1
			m_exp=1-exp(-time[i]/tau_m)
			m[i]=m[i]+m_exp*(m_inf-m[i])	
		end
	end
	
		
	return m
end

# ╔═╡ b7e6bef0-5fc9-11eb-20aa-fbad0211854c
function oca_time_IFCAL(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	oca=zeros(size(t))
	dt=time[2]-time[1]
	tau_oca,oca_inf=oca_carac_ifcal(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			oca[i]=0
		end
		if i>1
			oca_exp=1-exp(-time[i]/tau_oca)
			oca[i]=oca[i]+oca_exp*(oca_inf-oca[i])	
		end
	end
	
		
	return oca
end

# ╔═╡ 15e0f930-5fca-11eb-0342-f5d08e6bdfa8
function cca_time_IFCAL(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	cca=zeros(size(t))
	dt=time[2]-time[1]
	tau_cca,cca_inf=cca_carac_ifcal(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			cca[i]=0
		end
		if i>1
			cca_exp=1-exp(-time[i]/tau_cca)
			cca[i]=cca[i]+cca_exp*(cca_inf-cca[i])	
		end
	end
	
		
	return cca
end

# ╔═╡ f58f64e2-5fca-11eb-3ea6-85929a773b20
function oca_time_ISCAL(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	oca=zeros(size(t))
	dt=time[2]-time[1]
	tau_oca,oca_inf=oca_carac_iscal(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			oca[i]=0
		end
		if i>1
			oca_exp=1-exp(-time[i]/tau_oca)
			oca[i]=oca[i]+oca_exp*(oca_inf-oca[i])	
		end
	end
	
		
	return oca
end

# ╔═╡ 00813682-5fcb-11eb-347b-c9bb2db603c7
function cca_time_ISCAL(t,v)
	#time evolution of the gate variable m at time t and constant voltage v
	cca=zeros(size(t))
	dt=time[2]-time[1]
	tau_cca,cca_inf=cca_carac_iscal(v)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			cca[i]=0
		end
		if i>1
			cca_exp=1-exp(-time[i]/tau_cca)
			cca[i]=cca[i]+cca_exp*(cca_inf-cca[i])	
		end
	end
	
		
	return cca
end

# ╔═╡ a8fb46f0-5fc8-11eb-3204-d7e9eda60a41
function m_time_ICAN(t,Ca)
	#time evolution of the gate variable m at time t and constant voltage v
	m=zeros(size(t))
	dt=time[2]-time[1]
	tau_m,m_inf=m_carac_ican_lf(Ca)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			m[i]=0
		end
		if i>1
			m_exp=1-exp(-time[i]/tau_m)
			m[i]=m[i]+m_exp*(m_inf-m[i])	
		end
	end
	
		
	return m
end

# ╔═╡ 814622c0-5fc7-11eb-35f4-03ff614c201c
function o_time_ISK(t,Ca)
	#time evolution of the gate variable m at time t and constant voltage v
	o=zeros(size(t))
	dt=time[2]-time[1]
	tau_o,o_inf=o_carac_isk(Ca)
	
	for i=1:length(t)
		if t[i]==0 | i==1
			o[i]=0
		end
		if i>1
			o_exp=1-exp(-time[i]/tau_o)
			o[i]=o[i]+o_exp*(o_inf-o[i])	
		end
	end
	
		
	return o
end

# ╔═╡ f1e7a250-5fb9-11eb-1d5b-cd8dae74b163
begin 
	v_high = 49
	v_low = -9
	md"""Reference voltage to measure the differential conductance"""
end

# ╔═╡ d2ca5c60-5fc7-11eb-1fa6-7b215a39e099
begin 
	Ca_high = 0.1 #highest cell voltage
	Ca_low = 0.01 #lowest cell voltage
	md"""Reference calcium to measure the differential conductance"""
end

# ╔═╡ 500f8950-5e67-11eb-2419-fb2bc5be6028
md"""INa, proportional to m^3*h :""" 

# ╔═╡ bfe5ab40-5f59-11eb-00c3-a9ed51762f2d
begin 
	m_t_h = m_time_HH2(time,v_high)
	m_t_l = m_time_HH2(time,v_low)	
	h_t_h = h_time_HH2(time,v_high)
	h_t_l = h_time_HH2(time,v_low)
end

# ╔═╡ 784fbee0-5fba-11eb-2a8c-ad2851e762af
begin
	p1 = plot(time, m_t_h, label="V = 49mV")
	p1 = plot!(time, m_t_l, label="V = -9mV",legend=true,title="m(t) HH2")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	p2 = plot(time, h_t_h, label="V = 49mV")
	p2 = plot!(time, h_t_l, label="V = -9mV",legend=true,title="h(t) HH2")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	plot(p1,p2,layout=(1,2),legend = :right)
end

# ╔═╡ b9464bc0-5fc5-11eb-3c92-a1f3282014d7
plotly()

# ╔═╡ 2dcae222-5e69-11eb-2a43-638891bd2bef
md"""IKDR, proportional to n^4:"""

# ╔═╡ 2bc8d570-5fba-11eb-098a-cd7712a2932d
begin 
	n_t_h = n_time_HH2(time,v_high)
	n_t_l = n_time_HH2(time,v_low)
end

# ╔═╡ d0d9e860-5f56-11eb-1131-034e08395038
begin 
	plot(time,n_t_h, label = "n HH2 high")
	plot!(time,n_t_l, label = "n HH2 low")
end

# ╔═╡ 6819a420-5e69-11eb-25d0-2dfdbe7a03b4
md"""IKir, proportional to m"""

# ╔═╡ 32eb4e90-5fc5-11eb-160b-f9f8fa2458fb
begin
	m_t_IKIR_h = m_time_IKIR(time, v_high) 
	m_t_IKIR_l = m_time_IKIR(time, v_low)
end

# ╔═╡ 74263870-5fc5-11eb-0ff4-413089617660
begin 
	plot(time,m_t_IKIR_h, label = "m IKIR high")
	plot!(time,m_t_IKIR_l, label = "m IKIR low")
end

# ╔═╡ b5f77590-5e6a-11eb-1c97-677c0faef0ff
md"""ISK, proportional to o : """

# ╔═╡ 5eeaf840-5fc7-11eb-32d8-8d87361b576d
begin
	o_t_ISK_h = o_time_ISK(time, Ca_high) 
	o_t_ISK_l = o_time_ISK(time, Ca_low)
end

# ╔═╡ 60227160-5fc8-11eb-08cc-5700d1adf57d
begin 
	plot(time,o_t_ISK_h, label = "o ISK high")
	plot!(time,o_t_ISK_l, label = "o ISK low")
end

# ╔═╡ bcd5fe52-5e69-11eb-224e-6ba8a25eaa2f
md"""IfCaL, propotional to oca^2*cca : """ 

# ╔═╡ 44e76520-5fca-11eb-2edf-3d51d7a29d06
begin
	oca_t_IFCAL_h = oca_time_IFCAL(time, v_high) 
	oca_t_IFCAL_l = oca_time_IFCAL(time, v_low)
	cca_t_IFCAL_h = cca_time_IFCAL(time, v_high) 
	cca_t_IFCAL_l = cca_time_IFCAL(time, v_low)
end

# ╔═╡ 74a53a2e-5fca-11eb-1c6f-259d12b71aed
gr()

# ╔═╡ 775df192-5fca-11eb-0d72-5303a0e1e88a
begin
	p1_ifcal = plot(time, oca_t_IFCAL_h, label="V = 49mV")
	p1_ifcal = plot!(time, oca_t_IFCAL_l, label="V = -9mV",legend=true,title="oca(t) IfCaL")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	p2_ifcal = plot(time, cca_t_IFCAL_h, label="V = 49mV")
	p2_ifcal = plot!(time, cca_t_IFCAL_l, label="V = -9mV",legend=true,title="cca(t) IfCaL")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	plot(p1_ifcal,p2_ifcal,layout=(1,2),legend = :right)
end

# ╔═╡ e69d19a0-5fca-11eb-2948-01ea3494996e
plotly()

# ╔═╡ 5fe21b10-5e6a-11eb-17bd-05697da8b58a
md"""ICaAN, proportional to m^2 : """

# ╔═╡ 5ed191f0-5fc9-11eb-26a1-69ac525fbcba
begin
	m_t_ICAN_h = m_time_ICAN(time, Ca_high) 
	m_t_ICAN_l = m_time_ICAN(time, Ca_low)
end

# ╔═╡ 765836d0-5fc9-11eb-10e0-cda40cc6d2dc
begin 
	plot(time,m_t_ICAN_h, label = "m ICAN high")
	plot!(time,m_t_ICAN_l, label = "m ICAN low")
end

# ╔═╡ 2c2af850-5e6a-11eb-3303-8da2e6203569
md"""IsCaL, proportional to oca*cca :""" 

# ╔═╡ 104eeda0-5fcb-11eb-3d83-fdb9b22110ec
begin
	oca_t_ISCAL_h = oca_time_ISCAL(time, v_high) 
	oca_t_ISCAL_l = oca_time_ISCAL(time, v_low)
	cca_t_ISCAL_h = cca_time_ISCAL(time, v_high) 
	cca_t_ISCAL_l = cca_time_ISCAL(time, v_low)
end

# ╔═╡ 20224a10-5fcb-11eb-36a1-19649b383b81
gr()

# ╔═╡ 24721c2e-5fcb-11eb-2714-b72129b0d2d4
begin
	p1_iscal = plot(time, oca_t_ISCAL_h, label="V = 49mV")
	p1_iscal = plot!(time, oca_t_ISCAL_l, label="V = -9mV",legend=true,title="oca(t) IsCaL")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	p2_iscal = plot(time, cca_t_ISCAL_h, label="V = 49mV")
	p2_iscal = plot!(time, cca_t_ISCAL_l, label="V = -9mV",legend=true,title="cca(t) IsCaL")
	xaxis!("Time (ms)")
	yaxis!("Gate")
	
	plot(p1_iscal,p2_iscal,layout=(1,2),legend = :right)
end

# ╔═╡ fcd35440-5fcb-11eb-3aba-2d10224d3bba
plotly()

# ╔═╡ bb57b500-5e64-11eb-0089-4d4cdc95ab12
md"#### Interpretation"

# ╔═╡ bf9b2b10-5e64-11eb-30bb-4720e4206bc6
md"""For the model developped by LeFranc & al., we still have 4 time scales, that differ in multiple order of magnitude. Indeed, the slowest timescales considered here are much slowlier than the one derived from Aguiar's model. """

# ╔═╡ Cell order:
# ╟─114d6130-4448-11eb-1fe1-79e9307e281b
# ╟─bb9db670-4480-11eb-0f69-27230b104235
# ╟─06f737b0-4466-11eb-1e44-132d212493f6
# ╟─f1865370-4465-11eb-36a3-e9fce4a32b1f
# ╟─8a2e7c20-4465-11eb-07d9-012125275c94
# ╟─cf4df32e-4465-11eb-0e63-0b4843e0106e
# ╟─a9cc3ee0-4466-11eb-14d1-6bdd94c36b61
# ╟─6ea7b650-4466-11eb-1dbb-8d72f109e1b5
# ╟─992e8c8e-4471-11eb-3d77-37e98ca55298
# ╟─bf8a09a0-4471-11eb-3b66-83b22bb98343
# ╟─5f4a3a40-4473-11eb-2fc4-67fb906d1425
# ╟─07211162-4467-11eb-25f1-f51136af6f09
# ╟─5ca3ab50-4473-11eb-06ae-a7e1ca087792
# ╟─aefce830-4473-11eb-3788-9fb999c93b65
# ╟─7d5c7b00-4474-11eb-309d-0b82588c3202
# ╟─822ceed0-4474-11eb-29b5-b74f94e96e9b
# ╟─7d228d00-4479-11eb-332b-0d79441c725f
# ╟─2287eac0-4474-11eb-2199-db797510d042
# ╟─96bcc1e0-4474-11eb-3295-271e040d0e98
# ╟─ba1426b0-4474-11eb-1043-15338a53e007
# ╟─d5af5020-4479-11eb-0d14-37ed6b209a68
# ╠═843225f0-447a-11eb-03f6-b99782f3cfd1
# ╠═c73430b0-4479-11eb-0365-575c4023ece4
# ╠═27174f30-447a-11eb-33e5-2150566f3823
# ╠═f6942090-4474-11eb-0880-e9d27bc3aaa2
# ╟─f844cbe0-4476-11eb-276c-c1e5c9aa323b
# ╠═b011fa50-4476-11eb-2ba8-9f409fe63b5b
# ╟─12d4c3b0-4478-11eb-04e5-95694ad7bbd0
# ╟─21f2d9e0-4478-11eb-2a77-0980b0ddc553
# ╟─579e58ee-4480-11eb-18de-738026386048
# ╠═1fc33a90-4476-11eb-1566-bd559f113176
# ╟─6b74a000-4480-11eb-3ba2-9b12b88b7fd2
# ╠═22fba510-447d-11eb-217e-5f3e38cffd57
# ╟─9c366270-447e-11eb-28f8-b7b53a13b721
# ╟─e3ba6b80-4480-11eb-3b8a-6d058c216aef
# ╟─f4ba3a00-4480-11eb-2f52-ed18ff095332
# ╟─a06dfb00-4465-11eb-2c4e-358aa670830b
# ╟─ada102e0-448d-11eb-35e5-bf8a70668c47
# ╟─d1b91510-448c-11eb-2f47-0b687bc0badf
# ╟─4e7e2720-448d-11eb-159b-c3fabc4730f0
# ╟─65e3af70-448d-11eb-06d3-55698316fdfa
# ╟─88af27f0-448d-11eb-08d9-7f02468a867f
# ╟─dec8215e-448c-11eb-2b4f-ef8d76a86653
# ╟─cb7593ce-448d-11eb-2e9d-41d0bb8267a0
# ╟─eb36d990-448d-11eb-09ea-09e3fc202c2a
# ╟─6b0e34b0-448e-11eb-2ab2-73fdcc4aa987
# ╟─b7eba290-448e-11eb-3647-7f00734f6290
# ╟─cbebdfd0-448e-11eb-2479-65c6f1cdae49
# ╟─488bddb0-448f-11eb-1341-f5a1d2b89a61
# ╟─58465af0-448f-11eb-1d26-d36bc242e859
# ╟─74daa3ae-448f-11eb-0564-29ca38a575c9
# ╠═dc6ac90e-448f-11eb-0516-a7069fa36fd5
# ╟─d197a302-448f-11eb-0647-9940ef920759
# ╟─3e5bb8f0-4490-11eb-1e58-37baf741a15d
# ╠═42a4fb60-4490-11eb-1d7a-4369317b36c5
# ╠═767f84a0-4490-11eb-332a-850d0f90eafc
# ╟─bc688d00-5d99-11eb-0d12-9febd1d6c3da
# ╟─d2bcf370-5d99-11eb-1809-55bcf9ca2308
# ╟─55f8f600-5d9d-11eb-3ac8-fbe006345392
# ╟─c79ac630-5d9d-11eb-16fe-1dafe5e85c5d
# ╟─defd4bd0-5d9e-11eb-10ff-9d05a5b09b66
# ╠═ab133a50-60a0-11eb-0a3b-2724907a7773
# ╠═8bc739b0-60a2-11eb-19b0-e9535a19e66d
# ╟─ee219840-449f-11eb-14be-610f98fdbc83
# ╟─eef4edd0-449f-11eb-2411-d9b1b0309238
# ╠═efc27700-449f-11eb-1997-93d5df3ffd66
# ╟─1d50f380-44a1-11eb-1376-b39d73fe2777
# ╟─3469d7e0-44a0-11eb-09ce-69d2b842c9e3
# ╟─5a5deb20-44a1-11eb-07d4-8b9ba18fde66
# ╟─952586a0-44a1-11eb-15de-a31f19d44c61
# ╟─b67ed17e-44a1-11eb-1eee-6dd53ae70ea2
# ╟─8c36ad90-5f53-11eb-27e3-11ba5de25c5c
# ╟─36456390-5f58-11eb-2a12-b9a2c4062643
# ╟─64e84830-5f57-11eb-018f-b9f95f04aa1b
# ╟─fa1bee50-44a1-11eb-3749-732d1fed88c0
# ╟─232e4e00-44a2-11eb-0520-19a81bce5059
# ╟─227ead10-44a2-11eb-110c-4907e2b9ae60
# ╟─49ba29e0-44a2-11eb-1e4d-db081f25bf5b
# ╠═0ff07660-44a2-11eb-2e63-e5a89b4d2725
# ╟─1410dbd0-5fc4-11eb-3131-554175fe60a2
# ╟─6b25c440-44a2-11eb-3458-75e5693d6904
# ╟─b814eff2-44a3-11eb-30b9-77cf2c48a9ce
# ╟─be93d712-44a3-11eb-24a4-d9f9b6e39265
# ╠═88ff8550-44a2-11eb-0cf5-ad6414e90de5
# ╟─7a672470-44a3-11eb-1b0b-011a267dbe89
# ╟─b7e6bef0-5fc9-11eb-20aa-fbad0211854c
# ╟─15e0f930-5fca-11eb-0342-f5d08e6bdfa8
# ╟─ff4e6d60-44a3-11eb-27c7-e7035a6cabef
# ╟─6235a74e-44a8-11eb-071b-27aa3cb81ee2
# ╠═e83900e0-44a3-11eb-1e9c-ab689c827bc1
# ╟─96167980-44a5-11eb-3d5d-338b3a909740
# ╟─57d96650-44a5-11eb-1368-a1544d547929
# ╟─f58f64e2-5fca-11eb-3ea6-85929a773b20
# ╟─00813682-5fcb-11eb-347b-c9bb2db603c7
# ╟─cd1a0820-44a5-11eb-0dac-4dedff66eba6
# ╟─f30955e0-44a5-11eb-0e4f-2523d8c038e7
# ╟─fa702c00-44a5-11eb-1e87-695643d2f534
# ╟─2bc5df20-44a6-11eb-18c4-11fb7e77f8ae
# ╟─40798de0-44a6-11eb-1b01-636fc5cf76ce
# ╟─6ab97110-44a6-11eb-0731-bb27efbee0c4
# ╟─855c4a10-44a6-11eb-316b-91fcbfee76c7
# ╠═8f4d4d80-44a6-11eb-2946-9dc6120506ea
# ╠═67744c60-602a-11eb-3ad7-e725cc91a409
# ╠═a0b7f942-6025-11eb-1597-55d8b267bd0d
# ╠═9165e450-6031-11eb-2766-ad938e7ae86b
# ╠═eba6ee80-602e-11eb-2dda-e16185ae894e
# ╠═9a572910-44a8-11eb-3019-bb04b52c034c
# ╠═6eb04660-44a9-11eb-0ab0-2db3e39d1802
# ╠═0906bb52-44a9-11eb-258e-e3073a645dd6
# ╠═b0c620de-5c2f-499b-93aa-48ea4e62cecd
# ╠═8eae4837-e47e-4e7a-8562-8da0f03ca9f3
# ╟─d3cf44b2-44a9-11eb-28b7-b9c4724fcbb1
# ╟─50a11ae0-5cba-11eb-20fd-a73fea2e3647
# ╟─5581bb50-5cba-11eb-37bb-ef8f86d6d52f
# ╟─6a45bdbe-5cba-11eb-1f20-830ac1a98e24
# ╟─a8fb46f0-5fc8-11eb-3204-d7e9eda60a41
# ╟─f8636cae-5cba-11eb-049f-156d2214e577
# ╟─1953bc90-5cbb-11eb-13ea-f9f3faf51e23
# ╟─814622c0-5fc7-11eb-35f4-03ff614c201c
# ╟─7461b6de-5cbc-11eb-2c80-8d58d6b1e469
# ╟─a88915b0-5cbe-11eb-1055-d18d49d2da93
# ╟─75470dd0-5cbc-11eb-124c-13a7a088aad9
# ╟─998688b0-5cbc-11eb-3f57-53ba68d6c5c4
# ╠═9a6e50a0-5cbc-11eb-153f-8977b402aa2e
# ╠═f1031450-5cbc-11eb-3317-1f40ea067833
# ╠═1809c580-5cbd-11eb-0c4e-6171dedd4b55
# ╟─e0fdd9e0-5e66-11eb-0ed2-51ad9550e89b
# ╟─edaf35d0-5e66-11eb-0fff-d5f2e1786452
# ╠═31632f50-5fc2-11eb-3795-43129912495f
# ╟─205ae9de-5f55-11eb-2dd7-8903f2852577
# ╟─f1e7a250-5fb9-11eb-1d5b-cd8dae74b163
# ╟─d2ca5c60-5fc7-11eb-1fa6-7b215a39e099
# ╟─500f8950-5e67-11eb-2419-fb2bc5be6028
# ╟─bfe5ab40-5f59-11eb-00c3-a9ed51762f2d
# ╠═784fbee0-5fba-11eb-2a8c-ad2851e762af
# ╠═b9464bc0-5fc5-11eb-3c92-a1f3282014d7
# ╟─2dcae222-5e69-11eb-2a43-638891bd2bef
# ╟─2bc8d570-5fba-11eb-098a-cd7712a2932d
# ╟─d0d9e860-5f56-11eb-1131-034e08395038
# ╟─6819a420-5e69-11eb-25d0-2dfdbe7a03b4
# ╠═32eb4e90-5fc5-11eb-160b-f9f8fa2458fb
# ╠═74263870-5fc5-11eb-0ff4-413089617660
# ╟─b5f77590-5e6a-11eb-1c97-677c0faef0ff
# ╟─5eeaf840-5fc7-11eb-32d8-8d87361b576d
# ╟─60227160-5fc8-11eb-08cc-5700d1adf57d
# ╟─bcd5fe52-5e69-11eb-224e-6ba8a25eaa2f
# ╟─44e76520-5fca-11eb-2edf-3d51d7a29d06
# ╠═74a53a2e-5fca-11eb-1c6f-259d12b71aed
# ╟─775df192-5fca-11eb-0d72-5303a0e1e88a
# ╠═e69d19a0-5fca-11eb-2948-01ea3494996e
# ╟─5fe21b10-5e6a-11eb-17bd-05697da8b58a
# ╠═5ed191f0-5fc9-11eb-26a1-69ac525fbcba
# ╟─765836d0-5fc9-11eb-10e0-cda40cc6d2dc
# ╟─2c2af850-5e6a-11eb-3303-8da2e6203569
# ╠═104eeda0-5fcb-11eb-3d83-fdb9b22110ec
# ╠═20224a10-5fcb-11eb-36a1-19649b383b81
# ╟─24721c2e-5fcb-11eb-2714-b72129b0d2d4
# ╠═fcd35440-5fcb-11eb-3aba-2d10224d3bba
# ╟─bb57b500-5e64-11eb-0089-4d4cdc95ab12
# ╟─bf9b2b10-5e64-11eb-30bb-4720e4206bc6
