% 2D, unipolar drift-diffusion solver
% Jing Guo, Purdue CELAB, 2006-09

function [Ne_bias,Fn_bias,Id,Id_x,Id_y,Ge2D,J_TE1,J_TE2, Ge2D_T]=charge(Ec,Fn_bias_old,XI,YI,Vd,ind_s,mu,xg,yg,Ec_bias,Fn_bias1,phis,nu_row,nu_col,ny1,Rad,nx1)
%%%%%%%% Input: 
%% Ec(XI,YI): the 2D band profile in the channel region: the organic film
%% in cluding the source defined by ind_s in the film
%% Fn_ch_old: the quasi Fermi level, which is used only to compute
%% degneracy factor
%% Vd: drain voltage
%% ind_s: the index of the source
%%%%%%%% Output:
%% Ne_bias(XI,YI): the charge density
%% Id: the source drain voltage
%% Fn_bias: the quasi Fermi level computed using Ec and Ne_bias

[GeDensity Ge2D J_TE1 J_TE2 Ge2D_T]=GTunnel(xg,yg,Ec_bias,Fn_bias1,phis,nu_row,nu_col,ny1,Rad,nx1);

% GeDensity=100*GeDensity;
% Ge2D=100*Ge2D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNDAMENTAL physical constants
kBT=0.0259; 
q=1.6e-19; 
m0=9.1e-31;
hbar=1.055e-34;
Ne0=1e27;%2*(m0*kBT*q/(2*pi*hbar^2))^(3/2); % the charge density constant for 3D bulk material
sx=diff(XI);    sy=diff(YI);
%Ef=abs((Ec(1:(Np-1))-Ec(2:Np))./sx);
%mu=1./(1/mu0+(1/vs)*Ef);   % field dependent mobility
Nx=length(XI);
Ny=length(YI);

%%%%% set the boundary conditions
Ne_s=Ne0*fermi((0-Ec(ind_s(1)))/kBT,1,1/2);
Ne_d=Ne0*fermi((-Vd-Ec(Nx*Ny))/kBT,1,1/2);

%%%% Initializeation
Ne_bias=zeros(Nx*Ny,1);
deg_fac_x=zeros(Nx-1,Ny); cx1=zeros(Nx-1,Ny); cx2=zeros(Nx-1,Ny);
deg_fac_y=zeros(Nx,Ny-1); cy1=zeros(Nx,Ny-1); cy2=zeros(Nx,Ny-1);
    
%%%%% compute the electron density using DD model %%%%%% 
zetan=(Fn_bias_old-Ec)./kBT; 
deg_fac_I=fermi(zetan,1,1/2)./fermi(zetan,1,-1/2);    % treat degeneratcy and non-parabolic band structure.
Ec=reshape(Ec,Nx,Ny);       % get Ec(XI,YI)
for ii_y=1:Ny   % degeneracy factor in the x intervals
    ind=((ii_y-1)*Nx+1):(ii_y*Nx);
    deg_fac_x(:,ii_y)=interp1(XI,deg_fac_I(ind),XI(1:Nx-1)+sx*0.5,'cubic'); 
    %% electron flux in x direction is computed by I/q=C1i*n(i+1,j)-C2i*n(i,j)
    cx1(:,ii_y)=mu*(kBT*deg_fac_x(:,ii_y)./sx).*Bern((Ec(1:(Nx-1),ii_y)-Ec(2:Nx,ii_y))./(deg_fac_x(:,ii_y)*kBT));    % coefficient in Schaff-Gummel
    cx2(:,ii_y)=mu*(kBT*deg_fac_x(:,ii_y)./sx).*Bern((Ec(2:Nx,ii_y)-Ec(1:(Nx-1),ii_y))./(deg_fac_x(:,ii_y)*kBT));    % coefficient in Schaff_gunnel
end
for ii_x=1:Nx   % degeneracy factor in the x intervals
    ind=ii_x:Nx:((Ny-1)*Nx+ii_x);
    deg_fac_y(ii_x,:)=interp1(YI,deg_fac_I(ind),YI(1:Ny-1)+sy*0.5,'cubic'); 
    %% electron flux in y direction is computed by I/q=C1i*n(i+1,j)-C2i*n(i,j)
    cy1(ii_x,:)=mu*(kBT*deg_fac_y(ii_x,:)./sy').*Bern((Ec(ii_x,1:(Ny-1))-Ec(ii_x,2:Ny))./(deg_fac_y(ii_x,:)*kBT))';    % coefficient in Schaff-Gummel
    cy2(ii_x,:)=mu*(kBT*deg_fac_y(ii_x,:)./sy').*Bern((Ec(ii_x,2:Ny)-Ec(ii_x,1:(Ny-1)))./(deg_fac_y(ii_x,:)*kBT))';    % coefficient in Schaff_gunnel
end

%% set up the current continuity eqn: AX=B
AA=sparse(Nx*Ny,Nx*Ny);
BB=sparse(Nx*Ny,1);
coef=1;     % the length normalization constant for the derivitive of current 1e-10
for ii_y=2:(Ny-1)
    for ii_x=2:(Nx-1)
        ind=(ii_y-1)*Nx+ii_x;
        AA(ind,ind-Nx)=cy2(ii_x,ii_y-1)*2*coef/(sy(ii_y)+sy(ii_y-1));
        AA(ind,ind-1)=cx2(ii_x-1,ii_y)*2*coef/(sx(ii_x)+sx(ii_x-1));
        AA(ind,ind)=(-cx2(ii_x,ii_y)-cx1(ii_x-1,ii_y))*2*coef/(sx(ii_x)+sx(ii_x-1))+...
            (-cy2(ii_x,ii_y)-cy1(ii_x,ii_y-1))*2*coef/(sy(ii_y)+sy(ii_y-1));
        AA(ind,ind+1)=cx1(ii_x,ii_y)*2*coef/(sx(ii_x)+sx(ii_x-1));
        AA(ind,ind+Nx)=cy1(ii_x,ii_y)*2*coef/(sy(ii_y)+sy(ii_y-1));
    end
end

%%%%%% set the bondary condition
%% the left and right boundary, Neumann b.c. for Jx
ii_l=1; ii_r=Nx;
for ii_y=2:(Ny-1)
    ind_l=(ii_y-1)*Nx+ii_l;     % the left boundary
    AA(ind_l,ind_l)=-cx2(ii_l,ii_y);
    AA(ind_l,ind_l+1)=cx1(ii_l,ii_y);
    GeDensity(ind_l)=0;
    ind_r=(ii_y-1)*Nx+ii_r;     % the right boundary
    AA(ind_r,ind_r)=cx1(ii_r-1,ii_y);
    AA(ind_r,ind_r-1)=-cx2(ii_r-1,ii_y);
    GeDensity(ind_r)=0;
end
%% the top boundary (drain electrode): Direchlet b.c. for charge density
for ind=(Nx*(Ny-1)+1):(Nx*Ny)
    AA(ind,ind)=1;
    BB(ind)=Ne_d;
    GeDensity(ind)=0;
end
%% the bottom boundary (gate insulator): Neumann b.c. for Jy
for ind=1:Nx
    AA(ind,ind)=-cy2(ind,1);
    AA(ind,ind+Nx)=cy1(ind,1);
    GeDensity(ind)=0;
end

%% the CNT boundary (source electrode)
for ii=1:length(ind_s)
    AA(ind_s(ii),:)=sparse(1,Nx*Ny);
    AA(ind_s(ii),ind_s(ii))=1;
    BB(ind_s(ii))=Ne_s;
    GeDensity(ind_s(ii))=0;
end

%%%%%%%% solve for the electron density     
Ne_bias=real(full(AA\(BB-GeDensity))); 
Id=0;
Ne_2D=reshape(Ne_bias,Nx,Ny);       % convert to Ne_2D(x,y)
Id_x=-q*(Ne_2D(2:Nx,1:Ny).*cx1-Ne_2D(1:(Nx-1),1:Ny).*cx2);
Id_y=-q*(Ne_2D(1:Nx,2:Ny).*cy1-Ne_2D(1:Nx,1:(Ny-1)).*cy2);
% Id=trapz(XI,Id_y(1:Nx,Ny-1));  
Id=sum(Id_y(1:Nx,Ny-1));
% the minus sign because the source-drain current is along -x direction
Ec=reshape(Ec,Nx*Ny,1);       % get Ec(XI,YI)
Fn_bias=kBT*anti_dummy(Ne_bias./Ne0,1,1/2)+Ec;
