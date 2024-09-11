# This Julia file plots the figures used in "Selective advantages of redirected helping in a viscous population". Here we compute the critical cost-benefits ratios associated with survival and fecundity benefits when the number of brood failure per patch is either fixed or variable. We first scan the range of parameters of interest and compute the critical cost-benefits ratios in the four cases (fixed and survival benefits, fixed and fecundity benefits, variable and survival benefits and variable and fecundity benefits). We obtain the different figures plotted in the article with the exception of Figure 1 as it only plots an illustration and does not bring any specific results.

#This Julia has been run on Julia 1.10.1, on a Windows 10, 64 bits. Only one library is used: PlotlyJS. If the library needs to be installed, in your Julia terminal, run "import Pkg" and "Pkg.add("PlotlyJS")".


## Libraries

using PlotlyJS

## Parameters

# I set the values of the parameters so functions can be defined. These values will not be used after unless N.

p = 0.0 #Brood failure probability
N = 1000 #Number of offspring per breed
n = 3 #Number of breeding spot per patch
k = 2

#Dispersal parameters
d = 0.5
c = 0.5

#Survival rates
s0 = 0.7
s1 = 0.7

#Dictionary: I store the parameters into a dictionary
parameters = Dict("p" => p, "N" => N, "n" => n, "d" => d, "c" => c, "s0" => s0,  "s1" => s1, "k" => k)


# Functions

#Fixed number of brood failure per patch

### Reproductive value of an offspring
function Psi_fixed(parameters)
    n = parameters["n"]
    d = parameters["d"]
    c = parameters["c"]
    k = parameters["k"]
    s0 = parameters["s0"]
    s1 = parameters["s1"]
    psi = (k*(1-s0)+(n-k)*(1-s1))/(n-k)
    return psi
end

#Variable number of brood failure per patch

### Probability that a patch has k brood failure on n nest
P(k,p,n) = (factorial(n)/(factorial(k)*factorial(n-k)))*p^k*(1-p)^(n-k)

### Probability that a given spot on a patch where k brood failures happen has been won by a local offspring
function h(k,parameters)
    d = parameters["d"]
    c = parameters["c"]
    n = parameters["n"]
    p = parameters["p"]
    h = (1-d)*(n-k)/((1-d)*(n-k)+d*(1-c)*n*(1-p))
    return h
end

### Reproductive value of an offspring where k brood fail
function psi(k)
    s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]
    p = parameters["p"]
    c = parameters["c"]
    return (k*(1-s0)+(n-k)*(1-s1))/((n-k)*(1-d)+n*(1-p)*d*(1-c))
end

### Reproductive value of an offspring
function Psi(parameters)
    s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]
    p = parameters["p"]
    c = parameters["c"]
Psi = (p*(1-s0)+(1-p)*(1-s1))/(1-p)
return Psi
end

# Coefficients of consanguinity

## Fixed number of brood failures per patch

### alpha value for coefficient of consanguinity bar{r} (fixed number of brood failures per patch)
function alpha_fix(m,parameters)
	s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]

    c = parameters["c"]
    if m==n
        alpha = 0
    else
        hm = (1-d)/(1-c*d)
        ss = 2*s1*(1-s1)*hm+(1-s1)^2*hm^2
        sl = s1*(1-s0)*hm+(1-s1)*(1-s0)*hm^2
        ll = (1-s0)^2*hm^2
        alpha = (n-m)/n*(n-m-1)/(n-1)*ss+2*(n-m)/n*m/(n-1)*sl+m/n*(m-1)/(n-1)*ll
        alpha /= n-m
    end
	return alpha
end

### beta value for coefficient of consanguinity bar{r} (fixed number of brood failures per patch)

function beta_fix(m,parameters)
	s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]
    p = parameters["p"]
    c = parameters["c"]
    r_old=1
    if m == n
        beta = s0^2
    else
        hm = (1-d)/(1-c*d)
        ss = s1^2+2*s1*(1-s1)*hm*(n-m-1)/(n-m)+(1-s1)^2*hm^2*(n-m-1)/(n-m)
        sl = s0*s1+s0*(1-s1)*hm+s1*(1-s0)*hm*(n-m-1)/(n-m)+(1-s1)*(1-s0)*hm^2*(n-m-1)/(n-m)
        ll = s0^2+2*s0*(1-s0)*hm+(1-s0)^2*hm^2*(n-m-1)/(n-m)
        beta = (n-m)/n*(n-m-1)/(n-1)*ss+2*(n-m)/n*m/(n-1)*sl+m/n*(m-1)/(n-1)*ll

    end
return beta
end

### coefficient of consanguinity \bar{r} for a fixed number k of brood failures per patch
function r_fix(parameters)
    m = parameters["k"]
    alpha_m = alpha_fix(m,parameters)
    beta_m = beta_fix(m,parameters)
    r = alpha_m/(1-beta_m)
    return r
end

## Variable number of brood failures per patch

### alpha value for the coefficient of consanguinity for a variable number of brood failure

function alpha(m,parameters)
	s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]
    p = parameters["p"]
    c = parameters["c"]
	if m==n
        alpha = 0
    else
        hm = h(m,parameters)
        ss = 2*s1*(1-s1)*hm+(1-s1)^2*hm^2
        sl = s1*(1-s0)*hm+(1-s1)*(1-s0)*hm^2
        ll = (1-s0)^2*hm^2
        alpha = (n-m)/n*(n-m-1)/(n-1)*ss+2*(n-m)/n*m/(n-1)*sl+m/n*(m-1)/(n-1)*ll
        alpha /= n-m
    end
	return alpha
end

### beta value for the coefficient of consanguinity for a variable number of brood failure

function beta(m,parameters)
	s0 = parameters["s0"]
    s1 = parameters["s1"]
    n = parameters["n"]
    d = parameters["d"]
    p = parameters["p"]
    c = parameters["c"]
    hm = h(m,parameters)
    if m == n
        beta = s0^2
    else
        ss = s1^2+2*s1*(1-s1)*hm*(n-m-1)/(n-m)+(1-s1)^2*hm^2*(n-m-1)/(n-m)
        sl = s0*s1+s0*(1-s1)*hm+s1*(1-s0)*hm*(n-m-1)/(n-m)+(1-s1)*(1-s0)*hm^2*(n-m-1)/(n-m)
        ll = s0^2+2*s0*(1-s0)*hm+(1-s0)^2*hm^2*(n-m-1)/(n-m)
        beta = (n-m)/n*(n-m-1)/(n-1)*ss+2*(n-m)/n*m/(n-1)*sl+m/n*(m-1)/(n-1)*ll

    end
	return beta
end

### Coefficient of consanguinity \bar{r} for a variable number of brood failures per patch.
function r_var(parameters)
    n = parameters["n"]
    p = parameters["p"]
    alpha_m = 0
    beta_m = P(n,p,n)*beta(n,parameters)
    for m in range(0,n-1)
        alpha_m += P(m,p,n)*alpha(m,parameters)
        beta_m += P(m,p,n)*beta(m,parameters)
    end

    r = alpha_m/(1-beta_m)
    return r
end

# Relatedness coefficients R_s and R_f

## Fixed number of brood failures per patch

### For survival benefits

function R_s_fixed(parameters)
    Rs_n = r_fix(parameters)
    Rs_d = 1
    n = parameters["n"]
    d = parameters["d"]
    c = parameters["c"]
    hm=(1-d)/(1-c*d)
    Rs_n = r_fix(parameters)*(1-hm)
    Rs_d = 1-r_fix(parameters)*hm
    Rs = Rs_n/Rs_d
    return Rs
end

### For fecundity benefits

function R_f_fixed(parameters)
    d = parameters["d"]
    c = parameters["c"]
    hk = (1-d)/(1-c*d)
    Rf_n = r_fix(parameters)*(1-hk^2)
    Rf_d = 1-r_fix(parameters)*hk
    Rf = Rf_n/Rf_d
    return Rf
end

## Variable number of brood failures per patch

# For survival benefits

function R_s(parameters)
    Rs_n = r_var(parameters)
    Rs_d = 1
    n = parameters["n"]
    p = parameters["p"]
    for k in 1:n-1
        Rs_n -= r_var(parameters)*P(k,p,n)*k*h(k,parameters)/(n*p*(1-p^(n-1)))
        Rs_d -= r_var(parameters)*P(k,p,n)*k*h(k,parameters)/(n*p*(1-p^(n-1)))
    end

    Rs = Rs_n/Rs_d
    return Rs
end

#For fecundity benefits

function R_f(parameters)
    n = parameters["n"]
    d = parameters["d"]
    c = parameters["c"]
    p = parameters["p"]
    Rf_n_disp = r_var(parameters)*d*(1-c)*P(0,p,n)*psi(0)+r_var(parameters)*d*(1-c)*P(n,p,n)*psi(n)
    Rf_n_stay = 0
    Rf_d = 1

    for k in 1:n-1
        Rf_n_disp += r_var(parameters)*d*(1-c)*P(k,p,n)*psi(k)
        Rf_n_stay += r_var(parameters)*(1-d)*P(k,p,n)*k*(1-h(k,parameters))*psi(k)/(n*p*(1-p^(n-1)))
        Rf_d -= r_var(parameters)*P(k,p,n)*k*h(k,parameters)/(n*p*(1-p^(n-1)))
    end

    Rf_n = (Rf_n_disp+Rf_n_stay)/Psi(parameters)
    Rf = Rf_n/Rf_d
    return Rf
end

# Computations
##  We scan four different parameters: adult's survival (s_0=s_1=s), dispersal rate (d), brood failure rate (k=np for fixed and p for variable case), and the patch size (n). For every set of parameters, we compute the critical cost-benefit ratio associated with survival or fecundity benefits when the number of brood failures per patch is either fixed or variable and we store these values in 4 arrays. See manuscript for the formulas of the relatedness coefficients.


M = 100
i=1

parameters["c"] = 0.1

# Arrays storing values for the fixed number of brood failures case

cbs_fixed = zeros(M,M,4,3)
cbf_fixed = zeros(M,M,4,3)


# Arrays storing values for the fixed number of brood failures case

cbs = zeros(M,M,4,3)
cbf = zeros(M,M,4,3)

global i = 0
global j = 0
global k = 0
global l = 0
for s in range(1e-3,1-1e-3,M)
    parameters["s0"] = s
    parameters["s1"] = s
    global i += 1
    for d in range(1e-3,1-1e-3,M)
        parameters["d"] = d
        global j+=1
        for p in [0.05,0.25,0.5,0.75]
            parameters["p"] = p
            global k+=1
            for n in [2,4,8]
                parameters["n"] = n
                global l+=1

		# Fixed case
  		parameters["k"] = parameters["n"]*parameters["p"]
                cbs_fixed[i,j,k,l] = R_s_fixed(parameters)
                cbf_fixed[i,j,k,l] = Psi_fixed(parameters)*R_f_fixed(parameters)/parameters["s0"]

		# Variable case
                cbs[i,j,k,l] = R_s(parameters)
                cbf[i,j,k,l] = Psi(parameters)*R_f(parameters)/parameters["s0"]

            end
            global l = 0
        end
        global k = 0
    end
    global j = 0
end

# Figures

## Fixed number of brood failures
### Figure 2: Survival benefits
#### n = 2

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=2",x=0.5,y=0.8))

ps2 = surface(z=cbs_fixed[:,:,3,1],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)

p = plot(ps2,layout)
display(p)

#### n = 4

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=4",x=0.5,y=0.8))

ps1 = surface(z=cbs_fixed[:,:,2,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),showscale=false,colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
ps2 = surface(z=cbs_fixed[:,:,3,2],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
ps3 =surface(z=cbs_fixed[:,:,4,2],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)

p = plot([ps1,ps2,ps3],layout)
display(p)
#### n = 8

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=8",x=0.5,y=0.8))

ps1 = surface(z=cbs_fixed[:,:,2,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),showscale=false,colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
ps2 = surface(z=cbs_fixed[:,:,3,3],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
ps3 =surface(z=cbs_fixed[:,:,4,3],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)

p = plot([ps1,ps2,ps3],layout)
display(p)
## Figure 3: Fecundity benefits
### n = 2
layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=2",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf_fixed[:,:,3,1],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p = plot(pf2,layout)
display(p) 

### n = 4

layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=4",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf1 = surface(z=cbf_fixed[:,:,2,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf_fixed[:,:,3,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)
pf3 =surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),z=cbf_fixed[:,:,4,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p = plot([pf1,pf2,pf3],layout)
display(p)

### n = 8
layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=8",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf1 = surface(z=cbf_fixed[:,:,2,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf_fixed[:,:,3,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)
pf3 =surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),z=cbf_fixed[:,:,4,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p=plot([pf1,pf2,pf3],layout)
display(p)

## Figure 4: Comparison of the benefits

layout = Layout(yaxis_title="Survival s",xaxis_title="Backward migration rate m",zaxis_title="C:B",font_size=30,zaxis_range=[0,2])

p1 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs_fixed[:,:,2,1]-cbf_fixed[:,:,2,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="red",width=3),showlegend=false)
p2 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs_fixed[:,:,3,1]-cbf_fixed[:,:,3,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="orange",width=3),showlegend=false)
p3 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs_fixed[:,:,4,1]-cbf_fixed[:,:,4,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="green",width=3),showlegend=false)

p = plot([p1,p2,p3],layout)
display(p)


# Variable number of brood failures per patch
## Figure 5: Survival benefits

#### n = 2

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=2",x=0.5,y=0.8))

ps1 = surface(z=cbs[:,:,2,1],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),showscale=false,colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
ps2 = surface(z=cbs[:,:,3,1],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
ps3 =surface(z=cbs[:,:,4,1],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)

p=plot([ps1,ps2,ps3],layout)
display(p)

#### n = 4

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=4",x=0.5,y=0.8))

ps1 = surface(z=cbs[:,:,2,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),showscale=false,colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
ps2 = surface(z=cbs[:,:,3,2],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
ps3 =surface(z=cbs[:,:,4,2],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)

p = plot([ps1,ps2,ps3],layout)
display(p)

#### n = 8

layout = Layout(scene=attr(xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (survival)",fontsize=12),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=8",x=0.5,y=0.8))

ps1 = surface(z=cbs[:,:,2,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),showscale=false,colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
ps2 = surface(z=cbs[:,:,3,3],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
ps3 =surface(z=cbs[:,:,4,3],contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0,showscale=false)
p = plot([ps1,ps2,ps3],layout)
display(p)

## Figure 6: Fecundity benefits
### n = 2
layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=1",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf1 = surface(z=cbf[:,:,2,1],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf[:,:,3,1],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)
pf3 =surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),z=cbf[:,:,4,1],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p = plot([pf1,pf2,pf3],layout)
display(p)

### n = 4

layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=4",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf1 = surface(z=cbf[:,:,2,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf[:,:,3,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)
pf3 =surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),z=cbf[:,:,4,2],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p = plot([pf1,pf2,pf3],layout)
display(p)

### n = 8
layout = Layout(scene=attr(zaxis=attr(range=[0,2]),xaxis_title="← Survival s",yaxis_title="Dispersal d →",zaxis_title="C:B (Fecundity)",fontsize=10),scene_aspectratio=attr(x=1, y=1, z=1),scene_camera_eye=attr(x=1.5, y=1.5, z=1.5),legend=false,title=attr(text="n=8",x=0.5,y=0.8),zaxis=attr(range=[0,2]))

pf1 = surface(z=cbf[:,:,2,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),colorscale = [[0, "rgb(255, 255, 204)"], [1, "rgb(173, 216, 230)"]])
pf2 = surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="black"),x_end=1-1e-3,y=attr(show=true,start=1e-3,color="black",size=0.05),y_end=1-1e-3),z=cbf[:,:,3,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)
pf3 =surface(contours=attr(x=attr(show=true,start=1e-3,size=0.05, color="red"),x_end=1,y=attr(show=true,start=1e-3,color="red",size=0.05),y_end=1-1e-3),z=cbf[:,:,4,3],x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),opacity=0)

p = plot([pf1,pf2,pf3],layout)
display(p)
## Figure 7: Comparison of the benefits
### n = 2

layout = Layout(yaxis_title="Survival s",xaxis_title="Dispersal d",zaxis_title="C:B",font_size=30,zaxis_range=[0,2],title=attr(text="n=2",x=0.5))

p0 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,1,1]-cbf[:,:,1,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="red",width=3),showlegend=false)
p1 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,2,1]-cbf[:,:,2,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="orange",width=3),showlegend=false)
p2 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,3,1]-cbf[:,:,3,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="yellow",width=3),showlegend=false)
p3 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,4,1]-cbf[:,:,4,1],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="green",width=3),showlegend=false)

p = plot([p0,p1,p2,p3],layout)
display(p)
### n = 4

layout = Layout(yaxis_title="Survival s",xaxis_title="Dispersal d",zaxis_title="C:B",font_size=30,zaxis_range=[0,2],title=attr(text="n=4",x=0.5))

p0 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,1,2]-cbf[:,:,1,2],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="red",width=3),showlegend=false)
p1 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,2,2]-cbf[:,:,2,2],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="orange",width=3),showlegend=false)
p2 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,3,2]-cbf[:,:,3,2],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="yellow",width=3),showlegend=false)
p3 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,4,2]-cbf[:,:,4,2],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="green",width=3),showlegend=false)

p = plot([p0,p1,p2,p3],layout)
display(p)

### n = 8

layout = Layout(yaxis_title="Survival s",xaxis_title="Dispersal d",zaxis_title="C:B",font_size=30,zaxis_range=[0,2],title=attr(text="n=8",x=0.5))

p0 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,1,3]-cbf[:,:,1,3],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="red",width=3),showlegend=false)
p1 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,2,3]-cbf[:,:,2,3],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="orange",width=3),showlegend=false)
p2 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,3,3]-cbf[:,:,3,3],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="yellow",width=3),showlegend=false)
p3 = contour(x=range(1e-3,1-1e-3,M),y=range(1e-3,1-1e-3,M),z=cbs[:,:,4,3]-cbf[:,:,4,3],contours=attr(coloring="none",start=0),showscale=false,contours_end=0,line=attr(color="green",width=3),showlegend=false)

p = plot([p0,p1,p2,p3],layout)
display(p)
