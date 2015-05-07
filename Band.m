%%Calculation of electron drift current%%
%%%%%%%Wenchao Chen, June 2011%%%%%%%%%%%
clc
% clear all;
q=1.6e-19;
global NNx NNy
NNx=201; NNy=301;
N_x=(NNx-1)/2;
N_y=(NNy-1)/2;

Emtr_On=load('Band_minus4V.dat');           %%Load extration file from Sentaurus

%%Band
X_position1=Emtr_On(:,1);            %%X coordinates, um
Y_position1=Emtr_On(:,2);            %%Y coordinates, um
Conduction1=Emtr_On(:,3);            %%Conduction Band Edge
Valance1=Emtr_On(:,4);               %%Valance Band Edge
Htunneling1=Emtr_On(:,5);
QF1=Emtr_On(:,6); 
%%2D Plot
figure(61)
vis(:,1)=[X_position1]';
vis(:,2)=[Y_position1]';
vis(:,3)=Valance1;
[xlin_e ylin_e V2D1]=pr(vis);

figure(61)
vis(:,1)=[X_position1]';
vis(:,2)=[Y_position1]';
vis(:,3)=Conduction1;
[xlin_e ylin_e C2D1]=pr(vis);

figure(61)
vis(:,1)=[X_position1]';
vis(:,2)=[Y_position1]';
vis(:,3)=QF1;
[xlin_e ylin_e QF2D1]=pr(vis);

figure(61)
vis(:,1)=[X_position1]';
vis(:,2)=[Y_position1]';
vis(:,3)=Htunneling1;
[xlin_e ylin_e Ht2D1]=pr(vis);

hold on
figure(1)
% plot(-ylin_e,-C2D(:,N_x),'b','linewidth',2)
% hold on
% plot(-ylin_e,-V2D(:,N_x),'r','linewidth',2)
% hold on
% plot(-ylin_e,QF2D(:,N_x),'r','linewidth',2)
% hold on
plot(-ylin_e*1000,-V2D1(:,N_x),'g.','linewidth',2)
hold on
plot(-ylin_e*1000,QF2D1(:,N_x),'g.','linewidth',2)
xlabel('y (nm)','fontsize',[20])
ylabel('E (eV)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
axis([ 0 300 -1.2 0.9])

hold on
figure(71)
semilogy(-ylin_e*1000,Ht2D1(:,N_x),'g','linewidth',2)
xlabel('y (nm)','fontsize',[20])
ylabel('G (1/cm^3/s)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
axis([ 0 50 1e5 1e28])

hold on
figure(72)
semilogy(xlin_e*1000-100,Ht2D1(NNy,:),'g','linewidth',2)
xlabel('x (nm)','fontsize',[20])
ylabel('G (1/cm^3/s)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
axis([ -30 30 1e5 1e28])

hold on
figure(73)
plot(xlin_e*1000-100,-V2D1(NNy,:),'g','linewidth',2)
hold on
plot(xlin_e*1000-100,QF2D1(NNy,:),'r','linewidth',2)
xlabel('x (nm)','fontsize',[20])
ylabel('E (eV)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
axis([ 0 100 -1.2 0.9])
hold on
plot(-ylin_e*1000-10,-V2D1(:,N_x),'linewidth',2)
hold on
plot(-ylin_e*1000-10,QF2D1(:,N_x),'linewidth',2)

Jcurrent=1e6*1*1*1e-18*sum(sum(Ht2D1))*q/(200e-9)


figure(101)

for rrr=1:NNx
    Jax(rrr)=sum(Ht2D1(:,rrr))*q*1e6*1e-9*0.1;     %%in unit of mA/cm^2
end
hold on;
plot(xlin_e*1000-100,Jax,'b','linewidth',[2])     
xlabel('x (nm)','fontsize',[20])
ylabel('JA (mA/cm^2)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])
xlim(-100,100)



