%%%%%% The 2D DD-Poisson solver for transistors 
%%%%%% Jing Guo Mar.,2007
%%%%%% DD solved in the thin film channel + the source defined by ind_dd
%%%%%% Poisson: solved in the whole channel+oxide region
%%%%%% Poisson: charge exists in the region defined by ind_ch
clc
clear all
close all

inp;    % set up input parameters
init;   % initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ec2D=zeros(nu_row,nu_col,Ng_step+1,Nd_step+1);  % initialize Ec
Ne2D=zeros(nu_row,nu_col,Ng_step+1,Nd_step+1);  % initialize Ne
Id=zeros(Ng_step+1,Nd_step+1);                  % initialize Id

%%%%%% grid spacing for 2D Ec --> 1D Ec
sxp=interp1(xg(1:nu_col-1)+a./2,a,xg,'linear','extrap');    
Ych=yg(find(yg>=0));        % the channel position
ind_dd=bound(3):bound(6);       % the simulation region for 2D drift-diffusion transport
criterion_outer=2e-4;

for ii_vg=1:(Ng_step+1)     % gate voltage sweep
    Vg_bias=Vg0+(ii_vg-1)*Vg_step;
    Ecg=phig-Vg_bias;
    
    phis=phis_mod(ii_vg) %%phis=0.5+0.15*(-Vg_bias+0.6)
    Ecs=phis-Vs;                           %%%%BarrierModulation
    
    for ii_vd=1:(Nd_step+1) % drain voltage sweep
        Vd_bias=Vd0+(ii_vd-1)*Vd_step;
        Ecd=phid-Vd_bias;
        %%% initialization
        Ne_bias=zeros(nu_tot,1);
        Fn_bias=zeros(nu_tot,1); 
        
        %%%%%%%%%%%%%%%%%%%SetTheBoundaryWithRespectToTheSBHModulation%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_cnt=zeros(bound(6),1);  % the position indices of the CNT
for ii=1:nu_col
    for jj=1:nu_row
        i_node=ii+nu_col*(jj-1);
        if (sqrt((xg(ii)-C_x)^2+(yg(jj)-C_y)^2)<=Rad)
            L(i_node,:)=0; 
            Ab(i_node,i_node)=1;
            Bb(i_node,1)=Ecs;
            ind_cnt(i_node)=1;
        end            
    end
end
ind_cnt=find(ind_cnt);
%%% ind_cnt for solving charge.m in the thin film region only, index should
%%% refer to (bound(3)-1)
ind_cnt=ind_cnt-(bound(3)-1);   % the index of the CNT with reference to the 1st node of the channel

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% the initial guess as Laplace solution
        Bb(bound(1):bound(2))=Ecg;      % set the gate boundary condition
        Bb(bound(5):bound(6))=Ecd;      % set the drain boundary condition
        [Ec_bias]=poisson(0,0,1, L, Ab,Bb,vol,ind_ch); % Laplace solution 
        % initial guess for Fn: linear drop over y position
        Fn_bias(ind_dd)=interp1([0 max(yg)], [-Vs -Vd_bias], yg(ceil(ind_dd./bound(2)))); % for 2D Poisson
        error_outer=1;
        %%%%%%%%%%%%% self-consistent iteration
        while error_outer>criterion_outer
            %%%%% solve the non linear Poisson equation
            Ec_bias_old=Ec_bias;
            [Ec_bias]=poisson(Fn_bias,Ec_bias_old,0, L, Ab,Bb,vol,ind_ch);

            %%%%% quasi-Fermi level
            [Ne_bias(ind_dd),Fn_bias(ind_dd),Id_bias,Id_x,Id_y,Ge2D,J_TE1,J_TE2,Ge2D_T]=charge(Ec_bias(ind_dd),...
                Fn_bias(ind_dd),xg,Ych,Vd_bias,ind_cnt,mu,xg,yg,Ec_bias,Fn_bias,phis,nu_row,nu_col,ny1,Rad,nx1);  % solve 2D DD equation
            error_outer=max(abs(full(Ec_bias_old-Ec_bias)))
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Id(ii_vg,ii_vd)=sum(Id_bias)/length(Id_bias)/length(xg)/10;                      %%%in unit of mA/cm^2
        Ec2D(:,:,ii_vg,ii_vd)=reshape(full(Ec_bias),nu_col,nu_row)';
        Ne2D(:,:,ii_vg,ii_vd)=reshape(full(Ne_bias),nu_col,nu_row)';
        Fn2D(:,:,ii_vg,ii_vd)=reshape(full(Fn_bias),nu_col,nu_row)';
        
        %%%calculate current component ratio%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delta_y=Ych(2:length(Ych))-Ych(1:(length(Ych)-1));
        delta_x=xg(2:length(xg))-xg(1:(length(xg)-1));
        SizeofG=size(Ge2D);
        rownumber=SizeofG(1);
        coloumnumber=SizeofG(2);
        for ii=1:rownumber-1
            sumGe_x(ii)=sum(Ge2D(ii,1:coloumnumber-1).*delta_x');          %%Total generation Tunneling+Thermionic
            sumGe_x_T(ii)=sum(Ge2D_T(ii,1:coloumnumber-1).*delta_x');      %%Tunneling generation
        end
        sumGe=sum(sumGe_x.*delta_y');
        sumGe_T=sum(sumGe_x_T.*delta_y');
        
        Id_TE(ii_vg,ii_vd)=(1-sumGe_T/sumGe)*Id(ii_vg,ii_vd);
        Id_T(ii_vg,ii_vd)=(sumGe_T/sumGe)*Id(ii_vg,ii_vd);        
        %%%End of calculate current component ratio%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
    end % end of the drain voltage loop
end     % end of the gate voltage loop

%%%%% visualization
draw