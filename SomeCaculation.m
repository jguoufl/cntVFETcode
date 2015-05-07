delta_y=Ych(2:length(Ych))-Ych(1:(length(Ych)-1));
delta_x=xg(2:length(xg))-xg(1:(length(xg)-1));

SizeofG=size(Ge2D);
rownumber=SizeofG(1);
coloumnumber=SizeofG(2);

for ii=1:rownumber-1
    sumGe_x(ii)=sum(Ge2D(ii,1:coloumnumber-1).*delta_x');
end
sumGe=sum(sumGe_x.*delta_y')*q/(200e-9)

for jj=1:coloumnumber-1
    sumGe_y(jj)=sum(Ge2D(1:rownumber-1,jj).*delta_y);
end
sumGe_y=q*sumGe_y;

hold on;
figure(101)
plot(xg(1:length(xg)-1)*1e9,sumGe_y*0.1,'r','linewidth',[2])     
xlabel('x (nm)','fontsize',[20])
ylabel('JA (mA/cm^2)','fontsize',[20])
set(gca,'linewidth',[2],'fontsize',[20])