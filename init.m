%%%%%%% the initialization of 2D Poisson-DD solver
%%%% generate the non-uniform rectangular grid
%%% xg (grid position), a (grid spacing) along X direction
%%% yg (grid position), b (grid spacing) along X direction
y1=tins1;
y2=tins1+Dia+infs;
y3=tins1+Lch;
if alphay1==1
    ny1=round(y1/sy0);
else
    ny1=round(log(1-y1*(1-alphay1)/(sy0*alphay1))/log(alphay1));    
end
if alphay2==1
    ny3=round((y3-y2)/sy0);
else
    ny3=round(log(1-(y3-y2)*(1-alphay2)/(sy0*alphay2))/log(alphay2));
end
ny2=round((y2-y1)/sy0);
nu_row=ny1+ny2+ny3+1;
b=zeros(nu_row-1,1);
b((ny1+1):(ny1+ny2))=sy0*ones(ny2,1);
for ii=(ny1+ny2+1):(nu_row-1)
    b(ii)=sy0*alphay2^(ii-(ny1+ny2));
end
for ii=1:ny1
    b(ii)=sy0*alphay1^(ny1-ii+1);
end
yg=zeros(nu_row,1);
yg(1)=0;
for ii=2:nu_row
    yg(ii)=sum(b(1:(ii-1)));
end
t_bot=yg(ny1+1);            % the actual bottom oxide thickness
t_top=yg(nu_row)-t_bot;     % the actual channel length (top layer thickness)
yg=yg-t_bot;                % the oxide/channel interface as y=0
C_x=0;                      % x positiono f the center of the CNT
C_y=Dia/2+infs;             % y position of the center of the CNT, nm.   

x1=Dia/2;
x2=Lx/2;
nx1=round(x1/sx0);
if alphax==1
    nx2=round((x2-x1)/sx0);
else
    nx2=round(log(1-(x2-x1)*(1-alphax)/(sx0*alphax))/log(alphax));
end
nu_col=2*(nx1+nx2)+1;
a=zeros(nu_col-1,1);
a(1:nx1)=sx0*ones(nx1,1);
for ii=(nx1+1):(nx1+nx2)
    a(ii)=sx0*alphax^(ii-nx1);
end
a((nx1+nx2+1):2*(nx1+nx2))=a(1:(nx1+nx2));
a(1:(nx1+nx2))=a(2*(nx1+nx2):-1:(nx1+nx2+1));
xg=zeros(nu_col,1);
xg(1)=0;
for ii=2:nu_col
    xg(ii)=sum(a(1:(ii-1)));
end
xg=xg-xg(nu_col)/2;
%%%%%%%%%%%%%%% end of grid generation %%%%%%%%%%%%%%%%%

nu_tot=nu_row*nu_col;
bound(1)=1;
bound(2)=nu_col;
bound(3)=nu_col*ny1+1;
bound(4)=nu_col*(ny1+1);
bound(5)=nu_tot-nu_col+1;
bound(6)=nu_tot;
%%%%%%%%%% compute the area of each node for the finite volumn method

%%%%%%%%% set up the parameters for the polymer channel
vol=zeros(nu_tot,1);    % initialization
count=0;                % initialization
spacing=3e-10;          % the spacing between the CNT and the polymer
for ii=2:(nu_col-1)
    for jj=(ceil(bound(4)/nu_col)):(nu_row-1)
        ii_node=(jj-1)*nu_col+ii;
        vol(ii_node)=1/4*(a(ii-1)+a(ii))*(b(jj-1)+b(jj)); % volumn/tranverse length
        if (sqrt((xg(ii)-C_x)^2+(yg(jj)-C_y)^2)>(Rad+spacing))    % exclude the CNT
            count=count+1;
            ind_ch(count)=ii_node;  % record the channel node index            
        end
    end
end

%%%%%%%% Set up the linear part of the poisson equation AV=B
%%% set up the Laplace operator for Poisson Eq. in only the channel region
L=sparse(nu_tot,nu_tot);        % Laplace operator
for i_node=(bound(2)+1):(bound(5)-1)
    if (i_node>bound(2) & i_node<bound(3)) %Gate oxide region
        jy=ceil(i_node/nu_col);
        ix=i_node-(jy-1)*nu_col;
        if (ix>1 & ix<nu_col)
            L(i_node,i_node-nu_col)=(epso1/2)*(a(ix-1)+a(ix))/b(jy-1);
            L(i_node,i_node-1)=(epso1/2)*(b(jy-1)+b(jy))/a(ix-1);
            L(i_node,i_node)=-epso1/2*((a(ix-1)+a(ix))*(1/b(jy-1)+1/b(jy))+(b(jy-1)+b(jy))*(1/a(ix-1)+1/a(ix)));
            L(i_node,i_node+1)=(epso1/2)*(b(jy-1)+b(jy))/a(ix);
            L(i_node,i_node+nu_col)=(epso1/2)*(a(ix-1)+a(ix))/b(jy);
        end
    elseif (i_node>bound(3) & i_node<bound(4)) %Interface
        jy=ceil(i_node/nu_col);
        ix=i_node-(jy-1)*nu_col;
        if (ix>1 & ix<nu_col)
            L(i_node,i_node-nu_col-1)=1/8*(b(jy-1)/a(ix-1))*(epso1);
            L(i_node,i_node-nu_col)=((a(ix-1)+a(ix))/(2*b(jy-1))-1/8*b(jy-1)*(1/a(ix-1)+1/a(ix)))*(epso1);
            L(i_node,i_node-nu_col+1)=1/8*(b(jy-1)/a(ix))*(epso1);
            L(i_node,i_node-1)=3/8*(epso1*(b(jy-1)/a(ix-1))+epso2*(b(jy)/a(ix-1)));
            L(i_node,i_node)=-((a(ix-1)+a(ix))/(2*b(jy-1))+3/8*(b(jy-1)/a(ix-1)+b(jy-1)/a(ix)))*epso1...
                             -((a(ix-1)+a(ix))/(2*b(jy))+3/8*(b(jy)/a(ix-1)+b(jy)/a(ix)))*epso2;
            L(i_node,i_node+1)=3/8*(epso1*(b(jy-1)/a(ix))+epso2*(b(jy)/a(ix)));
            L(i_node,i_node+nu_col-1)=1/8*(b(jy)/a(ix-1))*(epso2);
            L(i_node,i_node+nu_col)=((a(ix-1)+a(ix))/(2*b(jy))-1/8*b(jy)*(1/a(ix-1)+1/a(ix)))*(epso2);
            L(i_node,i_node+nu_col+1)=1/8*(b(jy)/a(ix))*(epso2);
        end
    elseif (i_node>bound(4) & i_node<bound(5))
        jy=ceil(i_node/nu_col);
        ix=i_node-(jy-1)*nu_col;
        if (ix>1 & ix<nu_col)
            L(i_node,i_node-nu_col)=(epso2/2)*(a(ix-1)+a(ix))/b(jy-1);
            L(i_node,i_node-1)=(epso2/2)*(b(jy-1)+b(jy))/a(ix-1);
            L(i_node,i_node)=-epso2/2*((a(ix-1)+a(ix))*(1/b(jy-1)+1/b(jy))+(b(jy-1)+b(jy))*(1/a(ix-1)+1/a(ix)));
            L(i_node,i_node+1)=(epso2/2)*(b(jy-1)+b(jy))/a(ix);
            L(i_node,i_node+nu_col)=(epso2/2)*(a(ix-1)+a(ix))/b(jy);
        end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% apply the boundary condition, set up the boundary part of AX=B
%%%% Ab is the submatrix of A and Bb is the subvector of B for boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ecs=phis-Vs;
Ab=sparse(nu_tot, nu_tot); 
Bb=sparse(nu_tot,1);

%The left and right boudary
i_node_l_0=1;
i_node_r_0=nu_col;
for j=2:nu_row
    i_node_l=(j-1)*nu_col+i_node_l_0;
    i_node_r=(j-1)*nu_col+i_node_r_0;
    Ab(i_node_l,:)=0; Ab(i_node_r,:)=0;
    if i_node_l==bound(5)        
        Ab(i_node_l,i_node_l)=2;
        Ab(i_node_l,i_node_l+1)=-1;
        Ab(i_node_l,i_node_l-nu_col)=-1;
    else
        Ab(i_node_l,i_node_l)=1;
        Ab(i_node_l,i_node_l+1)=-1;
    end
    
    if i_node_r==bound(6)        
        Ab(i_node_r,i_node_r)=2;
        Ab(i_node_r,i_node_r-1)=-1;
        Ab(i_node_r,i_node_r-nu_col)=-1;
    else
        Ab(i_node_r,i_node_r)=1;
        Ab(i_node_r,i_node_r-1)=-1;
    end
end
    
%top boundary
for i_node=bound(5):bound(6)
    Ab(i_node,:)=0;
    Ab(i_node,i_node)=1; 
    %Bb(i_node,1)=Ecd;  % drain b.c. set in the self cons. loop in main.m
end

%%%%% set up the vector that distribute the charge of the CNT to the grids
Tv=zeros(nu_tot,1);
N_pc=100;
delt_sita=2*pi/N_pc;   
for ii_pc=1:N_pc
    sita=ii_pc*delt_sita;
    x_pc=C_x+Dia/2*sin(sita); y_pc=C_y-Dia/2*cos(sita); %The position of the point charge
    ind_x=find(hist(x_pc,xg)); ind_y=find(hist(y_pc,yg)); %The bot-left corner indecies
    ind=nu_col*(ind_y-1)+ind_x;     % the vector index
    Tv(ind)=Tv(ind)+(1/N_pc)/vol(ind);
end  

ind_C_cnt=nu_col*(find(hist(C_y,yg))-1)+find(hist(C_x,xg));

%%%%%%%% the Dirichlet boundary conditions imposed inside the metallic tube 
ind_cnt=zeros(bound(6),1);  % the position indices of the CNT
for ii=1:nu_col
    for jj=1:nu_row
        i_node=ii+nu_col*(jj-1);
        if (sqrt((xg(ii)-C_x)^2+(yg(jj)-C_y)^2)<=Rad)
%            L(i_node,:)=0; 
%            Ab(i_node,i_node)=1;
%            Bb(i_node,1)=Ecs;
            ind_cnt(i_node)=1;
        end            
    end
end
ind_cnt=find(ind_cnt);
%%% ind_cnt for solving charge.m in the thin film region only, index should
%%% refer to (bound(3)-1)
%ind_cnt=ind_cnt-(bound(3)-1);   % the index of the CNT with reference to the 1st node of the channel

%%%%% the bottom gate node
for i_node=bound(1):bound(2) 
    Ab(i_node,i_node)=1; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of boundary condition, end of initialization for Poisson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%