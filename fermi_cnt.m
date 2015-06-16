% This function computes the carrier statistical integral for CNT square root E(k)
% input: zeta=(mu-Em)/kT; delta=Egh/kT; der_flag=0 for integral, =1 for dI/d(zeta)
function [F]=fermi_cnt(zeta,delta,der_flag)
Nband=2;
tail_up=15;
lim_up=max(max(abs(zeta)),delta)+tail_up;

F=0;
if der_flag==0
    for ii=1:Nband
        deltaii=delta*ii;
        %ff=inline('1./(1+exp(sqrt(x^2+delta^2)-zeta))','x','zeta','delta');
        ff=inline('1./(1+exp(sqrt(x^2+delta^2)-zeta))-1./(1+exp(sqrt(x^2+delta^2)+zeta))','x','zeta','delta');
        F=F+myquad(ff,0,lim_up,1e-7,[],zeta,deltaii);
    end
elseif der_flag==1
    for ii=1:Nband
        deltaii=delta*ii;
        %ff=inline('0.25./cosh((sqrt(x^2+delta^2)-zeta)./2).^2','x','zeta','delta');
        ff=inline('0.25./cosh((sqrt(x^2+delta^2)-zeta)./2).^2+0.25./cosh((sqrt(x^2+delta^2)+zeta)./2).^2','x','zeta','delta');
        F=F+myquad(ff,0,lim_up,1e-7,[],zeta,deltaii);
    end
end


%function [F]=dummy(zeta,delta,der_flag)
%tail_up=15;
%lim_up=max(max(abs(zeta)),delta)+tail_up;

%if der_flag==0
%    ff=inline('1./(1+exp(sqrt(x^2+delta^2)-zeta))-1./(1+exp(sqrt(x^2+delta^2)+zeta))','x','zeta','delta');
%    F=myquad(ff,0,lim_up,1e-6,[],zeta,delta);
%elseif der_flag==1
%    ff=inline('0.25./cosh((sqrt(x^2+delta^2)-zeta)./2).^2+0.25./cosh((sqrt(x^2+delta^2)+zeta)./2).^2','x','zeta','delta');
%    F=myquad(ff,0,lim_up,1e-6,[],zeta,delta);
%end
