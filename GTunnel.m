%%%Calculation of Tunneling Induced Carrier Generation%%%%
%%%%%%%%%%%%%%%%%%%%%Wenchao Chen%%%%%%%%%%%%%%%%%%%%%%%%%
function [GeDensity Ge2D J_TE1 J_TE2 Ge2D_T]=GTunnel(xg,yg,Ec_bias1,Fn_bias1,phis,nu_row,nu_col,ny1,Rad,nx1)

[XX, YY]=meshgrid(xg,yg);
Ec_bias=reshape(full(Ec_bias1),nu_col,nu_row)';
Fn_bias=reshape(full(Fn_bias1),nu_col,nu_row)';

q=1.6e-19; 
R_Local=30e-9;                                 %%Radius of Nonlocal Mesh
St_Ag=0.01;                                    %%Starting angel of polar coordinate
N_line=300;                                    %%number of nonlocal lines
Theta=linspace(St_Ag,pi-St_Ag,N_line);
Del_Ag=Theta(2)-Theta(1);                      %%angle step
N_Mesh=300;                                    %%number of points per line
R_Mesh=linspace(0,R_Local,N_Mesh);             %%R in polar coordinate
Del_R=R_Mesh(2)-R_Mesh(1);                     %%R step

for kk=1:length(Theta)
    X_Mesh(:,kk)=cos(Theta(kk))*R_Mesh;                          
    Y_Mesh(:,kk)=sin(Theta(kk))*R_Mesh+Rad;
    Ec_LM(:,kk)=interp2(XX,YY,Ec_bias,X_Mesh(:,kk),Y_Mesh(:,kk)); %%Get Band profile by interpolation
    Fn_LM(:,kk)=interp2(XX,YY,Fn_bias,X_Mesh(:,kk),Y_Mesh(:,kk));
    Ec_temp=Ec_LM(:,kk);
    Fn_temp=Fn_LM(:,kk);
    
    GNe(:,kk)=zeros(N_Mesh,1);
    GNe_T(:,kk)=zeros(N_Mesh,1);
    
    Tunel_p(:,kk)=zeros(N_Mesh,1);                                %%Tunneling probility
    [GRate,Tr,ind_c,J_TE]=Tunnel(Ec_temp',Fn_temp',R_Mesh,phis);
    
    J_TE1(kk)=J_TE;
    G_TE1(kk)=J_TE/q/((1-ind_c/N_Mesh)*R_Local);
    Tunel_p(ind_c:N_Mesh,kk)=Tr;
    
    GNe_T(ind_c:N_Mesh,kk)=GRate;
    GNe(ind_c:N_Mesh,kk)=GRate+G_TE1(kk);                                   %%Carrier Generation in polar coordinate
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Mapping Tunneling Induced Generation From Polar To XY%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ge2D=zeros(nu_row-ny1,nu_col);                                   %%The generation density in XY
Ge2D_T=zeros(nu_row-ny1,nu_col);                                 %%The Tunneling generation density in XY

for jj=(ny1+nx1+1):nu_row                                        %%above the center of the CNT
    
    for ii=1:nu_col
        R_der=sqrt(XX(jj,ii)^2+(YY(jj,ii)-Rad)^2);
        
        if R_der>R_Local                             %%outside of the nonlocal mesh
            Ge2D(jj-ny1,ii)=0;
            Ge2D_T(jj-ny1,ii)=0;
        else
            angle=acos(XX(jj,ii)/R_der);
            Ngrid_angle=round(angle/Del_Ag);
            Ngrid_r=round(R_der/Del_R+1);
            if angle>(pi-St_Ag)||angle<St_Ag
                Ge2D(jj-ny1,ii)=0;
                Ge2D_T(jj-ny1,ii)=0;
            else
                Ge2D(jj-ny1,ii)=GNe(Ngrid_r,Ngrid_angle);
                Ge2D_T(jj-ny1,ii)=GNe_T(Ngrid_r,Ngrid_angle);
            end
        end
        R_der1=sqrt(XX(jj,ii)^2+(YY(jj,ii)-Rad)^2); %%distance to the center of CNT
        if R_der1<Rad                               %%inside the CNT region
            Ge2D(jj-ny1,ii)=0;
            Ge2D_T(jj-ny1,ii)=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Below The Center of the CNT in Channel%%%%%%%%%%%%%%%%%%
xxx=xg(ceil(nu_col/2):length(xg));
xxx=xxx(find(xxx<R_Local&xxx>0));
for zz=2:nx1
    X_Mesh1(:,zz)=R_Mesh;
    Y_Mesh1(:,zz)=(zz-1)*Rad/nx1*ones(length(R_Mesh),1);
    Ec_LM1(:,zz)=interp2(XX,YY,Ec_bias,X_Mesh1(:,zz),Y_Mesh1(:,zz)); %%Get Band profile by interpolation
    Fn_LM1(:,zz)=interp2(XX,YY,Fn_bias,X_Mesh1(:,zz),Y_Mesh1(:,zz));
    Ec_temp=Ec_LM1(:,zz);
    Fn_temp=Fn_LM1(:,zz);
    
    GNe1(:,zz)=zeros(N_Mesh,1);
    GNe1_T(:,zz)=zeros(N_Mesh,1);
    
    Tunel_p1(:,zz)=zeros(N_Mesh,1);                                %%Tunneling probility
    [GRate1,Tr1,ind_c1,J_TE]=Tunnel(Ec_temp',Fn_temp',R_Mesh,phis);
    J_TE2(zz)=J_TE;
    
    G_TE2(zz)=J_TE/q/((1-ind_c/N_Mesh)*R_Local);
    Tunel_p1(ind_c1:N_Mesh,zz)=Tr1;
    
    GNe1_T(ind_c1:N_Mesh,zz)=GRate1;
    GNe1(ind_c1:N_Mesh,zz)=GRate1+G_TE2(zz);
    
    
    Gn_temp=interp1(X_Mesh1(:,zz),GNe1(:,zz),xxx);
    Gn_temp_T=interp1(X_Mesh1(:,zz),GNe1_T(:,zz),xxx);
    
    Ge2D(zz,ceil(nu_col/2):(ceil(nu_col/2)+length(Gn_temp)-1))=Gn_temp';           
    Ge2D(zz,ceil(nu_col/2):-1:(ceil(nu_col/2)-length(Gn_temp)+1))=Gn_temp';
    
    Ge2D_T(zz,ceil(nu_col/2):(ceil(nu_col/2)+length(Gn_temp_T)-1))=Gn_temp_T';        %%%Tunneling Generation
    Ge2D_T(zz,ceil(nu_col/2):-1:(ceil(nu_col/2)-length(Gn_temp_T)+1))=Gn_temp_T';

end


GeDensity=reshape(Ge2D',(nu_row-ny1)*nu_col,1);
                