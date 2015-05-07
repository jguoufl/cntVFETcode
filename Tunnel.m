%%% Fuunction for Calculating Tunenling Generation Rate%%%
%%%%%%%%%%%%%%%Wenchao Chen, March,2012%%%%%%%%%%%%%%%%%%%

function [G_Ne Tr ind_c J_TE]=Tunnel(Ec_old,Fn_old,R_old,phis)
%%%Ec is band profile along nonlocal mesh line
%%%Fn is quasi Fermi Level along nonlocal mesh line
%%%phis is shottky barrier height

T=300;
q=1.6e-19;
kB=1.38e-23;
hbar=1.055e-34;
m0=9.11e-31;
Effectivemass=0.5;
Mass=Effectivemass*m0;
Richardson=4*pi*q*Effectivemass*m0*kB^2/(2*pi*hbar)^3;

N=length(Ec_old);

if isempty(find(Ec_old==phis))==1
    ind_c=1;
    Ec(ind_c+1:N)=Ec_old(ind_c:N-ind_c);
    Ec(ind_c)=phis;
    Fn(ind_c+1:N)=Fn_old(ind_c:N-ind_c);
    Fn(ind_c)=0;
    R=R_old(ind_c:N);
else
    ind_c=max(find(Ec_old==phis));                        %%Find the first point in the channel
    Ec=Ec_old(ind_c:N);
    Fn=Fn_old(ind_c:N);
    R=R_old(ind_c:N);
end

Veff=-1*Fn_old(N);                                        %%Effective bias at the end of nonlocal mesh line

%%%Calculate electric field%%%
Np=length(R);
Field=zeros(1,Np);
Field(1:Np-1)=(Ec(2:Np)-Ec(1:(Np-1)))./(R(2:Np)-R(1:(Np-1)));
F_Br=Field(Np-1);
Field(Np)=F_Br;

%%%Tr is tunneling probility%%%
for ii=1:Np
    if ii==1
        Tr(ii)=1;
    else
        Tr(ii)=exp(-1*abs(2/hbar*trapz(R(1:ii),abs(sqrt(2*Mass*q*(phis-Ec(1:ii)))))));     %%WKB for Tunneling Probility
    end    
    Gtunnel(ii)=Richardson*T/kB*abs(Field(ii))*Tr(ii)*(log((1+exp((-Ec(ii))*q/kB/T)))-log(1+exp((-Veff-Ec(ii))*q/kB/T)));  %% in unit of /m^3
end
G_Ne=abs(Gtunnel);

N_grid=200;
EnergyGrid=linspace(0,0.5,N_grid);
dEnergy=EnergyGrid(2)-EnergyGrid(1);

TE_bridge=q*Richardson*T/kB*(log((1+exp((-phis-EnergyGrid)*q/kB/T)))-log(1+exp((-1-phis-EnergyGrid)*q/kB/T)));
J_TE=abs(dEnergy*sum(TE_bridge));


