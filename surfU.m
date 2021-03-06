function createaxes(parent1, xdata1, ydata1, zdata1)
%CREATEAXES(PARENT1,XDATA1,YDATA1,ZDATA1)
%  PARENT1:  axes parent
%  XDATA1:  surf xdata
%  YDATA1:  surf ydata
%  ZDATA1:  surf zdata
 
%  Auto-generated by MATLAB on 20-Sep-2006 11:50:13
 
%% Create axes
axes1 = axes('FontSize',20,'Parent',parent1);
ylim(axes1,[0 max(max(ydata1))]);
%zlim(axes1,[-1 0]);
xlabel(axes1,'Gate [nm]');
ylabel(axes1,'Channel [nm]');
zlabel(axes1,'U [eV] or Ne [/cm^3]');
view(axes1,[-37.5 30]);
grid(axes1,'on');
hold(axes1,'all');
 
%% Create surf
surf1 = surf(xdata1,...
  ydata1,zdata1,...
  'Parent',axes1,...
  'ZDataSource','Ec_bias');
 
