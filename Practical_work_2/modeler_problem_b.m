% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
% 1) Export the required variables from pdetool and create a MATLAB script
%    to perform operations on these.
% 2) Define the problem completely using a MATLAB script. See
%    https://www.mathworks.com/help/pde/examples.html for examples
%    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[7 5 1]);
set(ax,'XLimMode','auto');
set(ax,'YLimMode','auto');
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([-10 10 0 -10],'R1');
pdeellip(0,0,10,10,...
0,'C1');
pderect([3.0000000000000036 13.000000000000004 10 0],'R2');
pderect([-4.0999999999999996 -3.9999999999999996 2.1000000000000001 0],'R3');
pderect([-5.7000000000000002 -5.6000000000000005 0 2.1000000000000001],'R4');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','C1-R1-R2-R3-R4')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(12,...
'neu',...
1,...
'-1i*2*pi*566/343',...
'0')
pdesetbd(11,...
'neu',...
1,...
'-1i*2*pi*566/343',...
'0')
pdesetbd(10,...
'neu',...
1,...
'0',...
'0')
pdesetbd(9,...
'neu',...
1,...
'0',...
'0')
pdesetbd(8,...
'neu',...
1,...
'0',...
'0')
pdesetbd(7,...
'neu',...
1,...
'0',...
'0')
pdesetbd(6,...
'neu',...
1,...
'0',...
'0')
pdesetbd(5,...
'neu',...
1,...
'0',...
'0')
pdesetbd(4,...
'neu',...
1,...
'0',...
'0')
pdesetbd(3,...
'neu',...
1,...
'0',...
'0')
pdesetbd(2,...
'neu',...
1,...
'0',...
'0')
pdesetbd(1,...
'neu',...
1,...
'0',...
'0')

% Mesh generation:
setappdata(pde_fig,'trisize',0.075800000000000006);
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')

% PDE coefficients:
pdeseteq(1,...
'-1',...
'(2*pi*566/343)^2',...
'0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['-1              ';...
'(2*pi*566/343)^2';...
'0               ';...
'1.0             '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','96606','10','pdeadworst',...
'0.5','longest','0','1e-4','','fixed','inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','real(u)');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')