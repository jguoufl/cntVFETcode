%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%yg=yg-tins1;            % the gate insulator-channel junction as y=0
%Visulize the simulated results
[XX, YY]=meshgrid(xg,yg);
[XId,YId]=meshgrid(xg,Ych);
[X_Idx,Y_Idx]=meshgrid(xg(1:nu_col-1)+a*0.5,Ych);
[X_Idy,Y_Idy]=meshgrid(xg,Ych(1:length(Ych)-1)+diff(Ych)*0.5);
[Idx_node]=interp2(X_Idx,Y_Idx,Id_x',XId,YId,'linear');
[Idy_node]=interp2(X_Idy,Y_Idy,Id_y',XId,YId,'linear');
Ec_bias=Ec2D(:,:,Ng_step+1);
Ne_bias=Ne2D(:,:,Ng_step+1);
Fn_bias=Fn2D(:,:,Ng_step+1);


figure(71)
semilogy(Ych*1e9,1e-6*Ge2D(:,ceil(nu_col/2)),'r','linewidth',[2])        %%%Generation Rate in unit of 1/cm^3/s, vertical direction
xlabel('y (nm)','fontsize',[20])
ylabel('G (1/cm^3/s)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
axis([ 0 50 1e5 1e28])
%%%%%% the vacuum level along the vertical direction
figure(1)
for ii_vg=1:(Ng_step+1)
    for ii_vd=1:(Nd_step+1)
        Ec_bias=Ec2D(:,:,ii_vg, ii_vd);  
        Fn_bias=Fn2D(:,:,ii_vg, ii_vd);
        plot(yg*1e9,Ec_bias(:,ceil(nu_col/2)),'-','linewidth',[2]); hold on
        plot(yg*1e9,Fn_bias(:,ceil(nu_col/2)),'r-','linewidth',[2]); hold on
    end
end
xlim([0 max(yg)]*1e9)
ylim([-1.2 1.5])
xlabel('y [nm]','fontsize',[28]);
ylabel('E_c [eV]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

figure(91)
for ii_vg=1:(Ng_step+1)
    for ii_vd=1:(Nd_step+1)
        Ec_bias=Ec2D(:,:,ii_vg, ii_vd);  
        Fn_bias=Fn2D(:,:,ii_vg, ii_vd);
        plot(xg*1e9,Ec_bias(ny1+5,:),'y','linewidth',[2]); hold on
        plot(xg*1e9,Fn_bias(ny1+5,:),'y.','linewidth',[2]); hold on
        
        plot(yg*1e9-10,Ec_bias(:,ceil(nu_col/2)),'g','linewidth',[2]); hold on
        plot(yg*1e9-10,Fn_bias(:,ceil(nu_col/2)),'g.','linewidth',[2]); hold on
    end
end
xlim([0 20])
ylim([-0.8 0.6])
xlabel('Position [nm]','fontsize',[28]);
ylabel('E_c [eV]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

%%%%%% the conduction band edge in the pcolor plot
figure(11)
gmapja(XX*1e9,YY*1e9,Ec_bias);
hold on
%quiver(XX,YY,-Exn,-Eyn,1)
%hold off;
%axis image;
xlabel('x [nm]','fontsize',[28]);
ylabel('y [nm]','fontsize',[28]);


%%%%%%%%% the vacuum level along the x-direction
figure(12)
Ecx=interp2(XX,YY,Ec_bias,xg,Rad*ones(nu_col,1));
plot(xg,Ecx,'linewidth',[2]);
xlabel('x [nm]','fontsize',[28]);
ylabel('E_c [eV]','fontsize',[28]);
ylim([-1 0.5]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

%%%%%% plot the charge density
figure(2)
for ii_vg=1:(Ng_step+1)
    for ii_vd=1:(Nd_step+1)
        Ne_bias=Ne2D(:,:,ii_vg,ii_vd);
        semilogy(yg*1e9,1e-6*Ne_bias(:,ceil(nu_col/2)),'-','linewidth',[2]); hold on
    end
end
xlim([min(yg) max(yg)]*1e9)
xlabel('y [nm]','fontsize',[28]);
ylabel('Ne [/cm^3]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);
%%%%%% the conduction band edge in the pcolor plot
surfU(figure(21),XX*1e9,YY*1e9,Ne_bias);
xlabel('x [nm]','fontsize',[28]);
ylabel('y [nm]','fontsize',[28]);

%%%%%%% plot position-resolved curretn density
surfU(figure(3),XId*1e9,YId*1e9,Idx_node);
surfU(figure(31),XId*1e9,YId*1e9,Idy_node);

%%%%%%% plot I-V characteristics
Vg=Vg0:Vg_step:(Vg0+Ng_step*Vg_step);
Vd=Vd0:Vd_step:(Vd0+Nd_step*Vd_step);
figure(4)
for ii_vd=1:(Nd_step+1)
    plot(Vg,abs(Id(:,ii_vd)),'linewidth',[2]);
    hold on
end
xlabel('V_G [V]','fontsize',[28]);
ylabel('I_D [mA/cm^2]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);

figure(5)
for ii_vg=1:(Ng_step+1)
    plot(Vd,Id(ii_vg,:),'linewidth',[2]);
    hold on;
end
xlabel('V_D [V]','fontsize',[28]);
ylabel('I_D [A/m]','fontsize',[28]);
set(gca, 'fontsize',[20], 'linewidth',[2]);
set(gca,'position',[0.15 0.20 0.74 0.70]);


