%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%yg=yg-tins1;            % the gate insulator-channel junction as y=0
%Visulize the simulated results
[XX, YY]=meshgrid(xg,yg);
[XId,YId]=meshgrid(xg,Ych);
%[X_Idx,Y_Idx]=meshgrid(xg(1:nu_col-1)+a*0.5,Ych);
%[X_Idy,Y_Idy]=meshgrid(xg,Ych(1:length(Ych)-1)+diff(Ych)*0.5);
%[Idx_node]=interp2(X_Idx,Y_Idx,Id_x',XId,YId,'linear');
%[Idy_node]=interp2(X_Idy,Y_Idy,Id_y',XId,YId,'linear');
Ec_bias=Ec2D(:,:,Ng_step+1);
Ev_bias=Ev2D(:,:,Ng_step+1);
Ne_bias=Ne2D(:,:,Ng_step+1);
EHOMO=5.4;          % the HOMO level of the polymer channel.

%%%%%% the vacuum level along the vertical direction
figure(1)
for ii_vg=1:(Ng_step+1)
    for ii_vd=1:(Nd_step+1)
        Ec_bias=Ec2D(:,:,ii_vg, ii_vd);   
        plot(yg*1e9,-Ec_bias(:,ceil(nu_col/2)),'-','linewidth',[2]); hold on
    end
end
xlim([min(yg) max(yg)]*1e9)
xlabel('y [nm]','fontsize',[28]);
ylabel('E_c [eV]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

%%%%%% the conduction band edge in the pcolor plot
figure(11)
gmapja(XX*1e9,YY*1e9,-Ec_bias); % the negative sign for p-type
hold on
%quiver(XX,YY,-Exn,-Eyn,1)
%hold off;
%axis image;
xlabel('x [nm]','fontsize',[28]);
ylabel('y [nm]','fontsize',[28]);

[Exn Eyn]=gradient(Ec_bias,xg,yg);
figure(14)
xzoom=10e-9;
yzoom=15e-9;
ind_x_zoom=find(abs(xg)<xzoom);
ind_y_zoom=find((yg>-5e-9)&(yg<yzoom));
gmapja(XX(ind_y_zoom,ind_x_zoom)*1e9,YY(ind_y_zoom,ind_x_zoom)*1e9,...
    -Ec_bias(ind_y_zoom,ind_x_zoom)+EHOMO);
%ind_x_zoom=ind_x_zoom(1:3:length(ind_x_zoom));
%ind_y_zoom=ind_y_zoom(1:3:length(ind_y_zoom));
%quiver(XX(ind_y_zoom,ind_x_zoom)*1e9,YY(ind_y_zoom,ind_x_zoom)*1e9,...
%    Exn(ind_y_zoom,ind_x_zoom),Eyn(ind_y_zoom,ind_x_zoom),'g','linewidth',[2])
axis image;

hold on;
plot((C_x+Dia/2*cos(0:0.01:2*pi))*1e9, (C_y+Dia/2*sin(0:0.01:2*pi))*1e9,...
    'k','linewidth',[2]);
xlabel('x [nm]','fontsize',[28]);
ylabel('y [nm]','fontsize',[28]);

%%%%%%%%% the vacuum level along the x-direction
figure(12)
Ecx=interp2(XX,YY,Ec_bias,xg,1e-9*ones(nu_col,1));
Evx=interp2(XX,YY,Ev_bias,xg,1e-9*ones(nu_col,1));
plot(xg*1e9,-Ecx,'linewidth',[2]); hold on; % negative sign to flip n to p type
%plot(xg*1e9,-Evx,'r','linewidth',[2]); hold on;
plot(xg*1e9,zeros(length(xg),1),'k:','linewidth',[2]); hold on
xlim(40*[-1 1])
xlabel('x [nm]','fontsize',[28]);
ylabel('E_c & E_v [eV]','fontsize',[28]);
%ylim([-1 0.5]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

%%%%%% plot the charge density
figure(2)
for ii_vg=1:(Ng_step+1)
    for ii_vd=1:(Nd_step+1)
        Ne_bias=Ne2D(:,:,ii_vg,ii_vd);
        semilogy(yg*1e9,Ne_bias(:,ceil(nu_col/2)),'-','linewidth',[2]); hold on
    end
end
xlim([min(yg) max(yg)]*1e9)
xlabel('y [nm]','fontsize',[28]);
ylabel('Ne [/cm^3]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);
%%%%%% the conduction band edge in the pcolor plot
% surfU(figure(21),XX*1e9,YY*1e9,Ne_bias);
% xlabel('x [nm]','fontsize',[28]);
% ylabel('y [nm]','fontsize',[28]);

