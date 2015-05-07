%%%%%% the input parameters for the 2D Poisson-DD solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of geometric parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Physical constant
epso=8.854e-12;
m0=9.11e-31;
kBT=0.0259;
hbar=1.055e-34;
q=1.6e-19;
Ne0=1e27; %2*(m0*kBT*q/(2*pi*hbar^2))^(3/2); % the charge density constant for 3D bulk material
mu=1.5e-4; %1.5e-4;            % m^2/Vs

Vs=0;               % the source (CNT) voltage
Vg0=-1;              % the bottom gate voltage
Vg_step=0.5;
Ng_step=0;
% phis_mod=[0.6054 0.5486 0.5 0.4595 0.4529 0.4461 0.4431 0.4414]; %%For Diameter 10nm
phis_mod=[0.6137 0.5579 0.5 0.4544 0.4437 0.4425 0.4412 0.4401]; %%For Diameter 5nm
% phis_mod=[0.6060 0.5569 0.5 0.4634 0.4487 0.4420 0.4386 0.4368]; %%For Diameter 2nm

Vd0=1;
Vd_step=0.05;
Nd_step=0;          % the drain (the top electrode) voltage

%%%%%% input the contact work function - the semiconductor affinity
phig=0;             % Gate work function minus semiconductor affinity
phis=0.5;           % The souce-channel SB height
phid=0.5;             % The souce-channel SB height

%Geometric parameters
tins1=20*1e-9;      % bottom insulator thickness, nominal value in m
epso1=4;            % bottom insulator dielectric constant
Lch=300*1e-9;       % top insulator thickness
epso2=6;            % top insulator dielectric contant
top_flag=1;         % top gate flag
Dia=5*1e-9;        % tube diameter in nm
Rad=Dia/2;          % nm
infs=1e-12;

%%%% the simulated region and the grid spacing
Lx=200*1e-9;   % in m
alphax=1.1;     % the grid ratio
sx0=0.1*1e-9    % in m
sy0=0.1*1e-9    % in m
alphay1=1.1;    % the grid ratio in the gate insulator
alphay2=1.1;    % the grid ratio in the thin film channel

% end of input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%