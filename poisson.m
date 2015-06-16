%%%%%%%%% A general purpose 2D non-linear Poisson equation
%% by Jing Guo, Sep 2006
function [Ec_bias, Ne_bias]=poisson(Fn,Ec_old,Laplace_flag, L, Ab,Bb,vol, ind_ch, Tv, ind_C_cnt, phib0)
%%%%%%% Input: 
%% Fn: quasi-Fermi level
%% Ec_old: the old band profile
%% Laplace_flag: 1 for Laplace solution, otherwise for Poisson solution
%% L: the CHANNEL Laplace operator on the L.H.S. of the equation 
%% Ab: the boundary node submatrix on the L.H.S. of the equation
%% ind_ch: the vector that defines channel node index (where charge exists)
%% vol: all element volumns for the finite volumn method 
%%%%%%% Output:
%% Ec_bias: the band profile
%% Ne_bias: the charge density at equil. computed by Ne0*fermi((Fn-Ec)/kBT,1,1/2),

m0=9.11e-31;
kBT=0.0259;
q=1.6e-19;
hbar=1.055e-34;
epso=8.854e-12;
Ne0=1e27;%Ne0=(m0*kBT*q/(2*pi*hbar^2))^(3/2); % the charge density constant for 3D bulk material
nu_tot=length(L);
Ne_old=zeros(nu_tot,1);    % initialize the chage density, only meaningful @ equil.
%d_cnt=10e-9;               % in m, the CNT diameter
%delta=0.42/(d_cnt*1e9)/kBT; % semiconducting CNT, the half CNT Eg normalized by kBT
delta=0;                   % metallic CNT

acc=1.42e-10;
Vpp=-3;
Ne_cnt=8/(3*pi*acc*abs(Vpp))*kBT;

if Laplace_flag==1
    %%%%%% initial guess as Laplace solution
    A=Ab+L; % L.H.S.: boundary condition+Laplace operator
    B=Bb+0; % R.H.S.: boundary condition +zero vector 
    Ec_bias=real(A\B);     % a direct method to solve AX=B
else
    error_inner=1; 
    criterion_inner=1e-3;
    while error_inner>criterion_inner
        zeta=(0-(Ec_old(ind_C_cnt)-phib0))/kBT; % the Fermi level = 0, for CNT
        Ne_old=zeros(nu_tot,1);
        Ne_old(ind_ch)=Ne0*fermi((Fn(ind_ch)-Ec_old(ind_ch))./kBT,1,1/2);        
        Ne_old=Ne_old+Ne_cnt*fermi_cnt(zeta,delta,0)*Tv;    % to be filled 
               
        dummy_prime=zeros(nu_tot,1);        % initialization
        dummy_prime(ind_ch)=-(Ne0/kBT)*fermi((Fn(ind_ch)-Ec_old(ind_ch))./kBT,1,-1/2);
        dummy_prime=dummy_prime-(Ne_cnt/kBT)*fermi_cnt(zeta,delta,1)*Tv;    % to be filled
        
        %%%%% form 1 step of Newton-Ralphson as AX=B
        A=Ab+L-(-q/epso)*spdiags(vol.*dummy_prime,0,nu_tot,nu_tot);     % vol for the Finite volumn method       
        B=Bb+(-q/epso)*sparse(vol.*(Ne_old-dummy_prime.*Ec_old));         
        Ec_bias=real(A\B);
        error_inner=max(abs(full(Ec_bias-Ec_old)))
        Ec_old=Ec_bias;        
    end
end
Ne_bias=Ne_old;
