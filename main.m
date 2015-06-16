%%%%%% The 2D DD-Poisson solver for transistors 
%%%%%% Jing Guo Mar.,2007
%%%%%% DD solved in the thin film channel + the source defined by ind_dd
%%%%%% Poisson: solved in the whole channel+oxide region
%%%%%% Poisson: charge exists in the region defined by ind_ch

clear all
%close all

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
    for ii_vd=1:(Nd_step+1) % drain voltage sweep
        Vd_bias=Vd0+(ii_vd-1)*Vd_step;
        Ecd=phid-Vd_bias;
        %%% initialization
        Ne_bias=zeros(nu_tot,1);
        Fn_bias=zeros(nu_tot,1); 
        
        %%% the initial guess as Laplace solution
        Bb(bound(1):bound(2))=Ecg;      % set the gate boundary condition
        Bb(bound(5):bound(6))=Ecd;      % set the drain boundary condition
        if (ii_vg==1) & (ii_vd==1)
            [Ec_bias]=poisson(0,0,1, L, Ab,Bb,vol,ind_ch, Tv, ind_C_cnt, phib0); % Laplace solution 
        end
        % initial guess for Fn: linear drop over y position
        Fn_bias(ind_dd)=interp1([0 max(yg)], [-Vs -Vd_bias], yg(ceil(ind_dd./bound(2)))); % for 2D Poisson
        %error_outer=1e-10;
        %%%%%%%%%%%%% self-consistent iteration
        %while error_outer>criterion_outer
            %%%%% solve the non linear Poisson equation
            Ec_bias_old=Ec_bias;
            [Ec_bias]=poisson(Fn_bias,Ec_bias_old,0, L, Ab,Bb,vol,ind_ch, Tv, ind_C_cnt, phib0);

            %%%%% quasi-Fermi level
            Ne_bias(ind_ch)=Ne0*fermi((Fn_bias(ind_ch)-Ec_bias(ind_ch))./kBT,1,1/2);
        %    [Ne_bias(ind_dd),Fn_bias(ind_dd),Id_bias,Id_x,Id_y]=charge(Ec_bias(ind_dd),...
        %        Fn_bias(ind_dd),xg,Ych,Vd_bias,[],mu);  % solve 2D DD equation
        %    error_outer=max(abs(full(Ec_bias_old-Ec_bias)))
        %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Id(ii_vg,ii_vd)=sum(Id_bias)/length(Id_bias);
        %Ec_bias(ind_cnt)=Ec_bias(ind_cnt)-(phib0-Eg_cnt/2);
        Ev_bias=Ec_bias-Eg_poly;
        %Ev_bias(ind_cnt)=Ec_bias(ind_cnt)-Eg_cnt;
        Ec2D(:,:,ii_vg,ii_vd)=reshape(full(Ec_bias),nu_col,nu_row)';
        Ev2D(:,:,ii_vg,ii_vd)=reshape(full(Ev_bias),nu_col,nu_row)';
        Ne2D(:,:,ii_vg,ii_vd)=reshape(full(Ne_bias),nu_col,nu_row)';
    end % end of the drain voltage loop
end     % end of the gate voltage loop

%%%%% visualization
draw