%% FINITE ELEMENT METHOD - ACOUSTICS (Frequency domain) with BC
% Viktoriia Boichenko, 2022
% Department of Civil Engineering, University of Coimbra
%-----------------------------------------------------
% This code analyses a 2D problem using FEM with triangular (3 node)
% elements. Allows prescribing Dirichlet, Neumann or Robin boundary
% conditions, which may be frequency-dependent
%% Load mesh and prepare information about model
% don't change this part, except for the name of the file
clear
close all
load test_pw2_fem.mat;
tri=t(1:3,:)';
mat=t(4,:);
nmat=max(mat);
rho=zeros(1,nmat);
c=zeros(1,nmat);
x=p(1,:);
y=p(2,:);
% Detect boundaries
n0=e(1,:); % initial node of segment
n1=e(2,:); % final node of segment
nm1=e(6,:); % element on left
nm2=e(7,:); % element on right
nb=e(5,:);  % boundary code od segment
nn=find(nm1==0 | nm2==0); % boundary segments
nm_bound=max(nm1(nn),nm2(nn)); % material of the element at the boundary

%% Frequencies to be calculated
% will create a figure with colormaps of results for each frequency with value "true".
frequencies=100:10:1000;
plot_results=frequencies*0;

frequencies_third=[125 315 566];
frequencies = [frequencies frequencies_third];

plot_results(1,numel(frequencies)-3)=true;
plot_results(1,numel(frequencies)-2)=true;
plot_results(1,numel(frequencies))=true;
%% Properties of the involved media
c(:)=343.0;
rho(:)=1.22;

%% select groups of nodes to store results
nodes_select={find(x<min(x)+.2);find(x>max(x)-0.2)};
nodes_results={};


%% Boundary conditions definition
%0 - Neumann vn (value in BCValue)
%1 - Dirichlet p (value in BCValue)
%2 - Robin Z (value in BCValue)
%3 - Robin (value in BCValue) + Neumann (unit velocity)
%4 - Robin alfa sound absorption (Real impedance calculated as
%    Z=rho*c*(1+sqrt(1-alfa))/(1-sqrt(1-alfa)) , alfa in BCValue)- can be
%    variable with frequency
BCboundary=unique(nb(nn));
BCType=zeros(1,numel(BCboundary));
BCValue=zeros(1,numel(BCboundary));
% BCType(1:end)=0;
% BCValue(1:end)=0;

% Implementing boundary conditions from the picture we got
% the ones we not initialize are by default Neumann zero and rigid surfaces.
BCType([11,12])=4;
BCValue([11,12])=1; % semicircle

A1=6;
A2=6;
freq=500+A1*10+A2;
alpha1=(0.002*freq+0.04*log10(freq)-(0.0008*freq)^2)/(2*sqrt(A1+10));
alpha2=(0.002*freq+0.04*log10(freq)-(0.0008*freq)^2)/(2*sqrt(A1+5));
alpha3=(0.002*freq+0.04*log10(freq)-(0.0008*freq)^2)/(2*sqrt(A1+1));

BCType(10)=4;
BCValue(10)=alpha1;

BCType(7)=4;
BCValue(7)=alpha2;

BCType(4)=4;
BCValue(4)=alpha3;

%% Internal source points
% Allows considering a source somewhere inside the domain. If no source
% exists, consider Amplitude_source=0.0.
Amplitude_source=1.0;
xs=-7.3;
ys=0.5;


%% Boundary conditions and materials plot
% Plot materials with different colors
figure (1)
clf
pdegplot(g,'SubdomainLabels','on','EdgeLabels','on')
set(1,'Name','Problem geometry and materials')
ax=patch(x(tri)',y(tri)',mat,'linestyle','none');
set(ax,'FaceAlpha',0.6)
% hold on
% for ii=1:nmat
%     nt=find(mat==ii);
%     ti=nt(ceil(numel(nt)/2));
%     xm=mean(x(tri(ti,:)));
%     ym=mean(y(tri(ti,:)));
%     text(xm,ym,strcat('M',num2str(ii)),'fontweight','bold','fontsize',14)
% end
% hold off
colorbar
hold on;
is_on_list=[];
for ii=1:numel(BCboundary)
    if(BCType(ii)==0)
        clr='-g';
        if(BCValue(ii)==0)
            clr=':g';
        end
    elseif(BCType(ii)==1)
        clr='b';
    else
        clr='r';
    end
    try
        nn1=find(nb==BCboundary(ii));
        plot([x(n0(nn1));x(n1(nn1))],[y(n0(nn1));y(n1(nn1))],clr,'linewidth',2)
        if(ismember(BCboundary(ii),is_on_list)==0)
            refel=ceil(numel(nn1)/2);
%             text(x(n0(nn1(refel))),y(n0(nn1(refel))),num2str(ii),'fontweight','bold')
            is_on_list=[is_on_list BCboundary(ii)];
        end
    catch
    end
end
hold off
view([0 90])
axis equal
answer = questdlg("Continue?","Question.","Yes","No","Yes");
if(answer=="No")
    return
end

%% Global Matrices - frequency independent for constant c and rho
nnodes=numel(x);
nelems=size(t,2);

I=zeros(nelems,9);
J=zeros(nelems,9);
X=zeros(nelems,9);
XM=zeros(nelems,9);
% uses sparse matrix assembly for higher performance
for ii=1:nelems     
    gn=tri(ii,:);
    fxe=x(gn);
    fye=y(gn);
    [ke,me]=keme2D(fxe,fye,rho(mat(ii)),c(mat(ii)));
    I0=repmat(gn,3,1);
    I(ii,:)=I0(:);
    J(ii,:)=repmat(gn,1,3);
    X(ii,:)=ke(:);
    XM(ii,:)=me(:);  
end
K=sparse(I(:),J(:),X(:),nnodes,nnodes);
M=sparse(I(:),J(:),XM(:),nnodes,nnodes);
clear I J X XM;

%% Loop through frequencies
for ifreq=1:numel(frequencies)
    frequency=frequencies(ifreq);
    w=2*pi*frequency;
        
    %% Boundary conditions value redifinition if frequency-dependent
    % this section allows overriding the previously defined values, making
    % them frequency dependent
    %BCValue(1:4)=1;
%     BCValue(4)=frequency/1000;
%     BCValue(4)=(0.002*frequency+0.04*log10(frequency)*(0.0008*frequency)^2)/(2*sqrt(8)+5);
% 
    %% Boundary conditions vectors and matrix C - can be frequency dependent
    % Robin conditions
    C=spalloc(nnodes,nnodes,nnodes*20);
    for ii=1:numel(nn)
        nodes=[n0(nn(ii));n1(nn(ii))];
        nbound=find(BCboundary==nb(nn(ii)));
        try
            bt=BCType(nbound);
            if(bt==2 | bt==3)
                bv=BCValue(nbound); % Z
                L=sqrt((x(nodes(2))-x(nodes(1)))^2+(y(nodes(2))-y(nodes(1)))^2);
                C(nodes,nodes)=C(nodes,nodes)+[L/3 L/6;L/6 L/3]*1/bv;
            elseif(bt==4)
                bv=BCValue(nbound); % alfa
                matb=nm_bound(nn(ii));
                Z=rho(matb)*c(matb)*(1+sqrt(1-bv))/(1-sqrt(1-bv));
                L=sqrt((x(nodes(2))-x(nodes(1)))^2+(y(nodes(2))-y(nodes(1)))^2);
                C(nodes,nodes)=C(nodes,nodes)+[L/3 L/6;L/6 L/3]*1/Z;
            end
        catch
        end
    end
    
    % Neumann conditions
    F=zeros(nnodes,1);
    for ii=1:numel(nn)
        nodes=[n0(nn(ii));n1(nn(ii))];
        nbound=find(BCboundary==nb(nn(ii)));
        
        try
            bt=BCType(nbound);
            if(bt==0 | bt==3)
                bv=BCValue(nbound); % Vn
                if (bt==3)
                    bv=1;   %consider unit vn if Robin+Neumann condition
                end                
                L=sqrt((x(nodes(2))-x(nodes(1)))^2+(y(nodes(2))-y(nodes(1)))^2);
                F(nodes(1))=F(nodes(1))-1i*w*bv*L/2;
                F(nodes(2))=F(nodes(2))-1i*w*bv*L/2;
            end
        catch
        end
    end   
    
    %% Internal source points
    % Allows considering a source somewhere inside the domain. If no source
    % exists, consider Amplitude_source=0.0.    
    TP=triangulation(tri,p');
    ts=pointLocation(TP,[xs,ys]);
    if(~isnan(ts))
        area_t=polyarea(x(tri(ts,:)),y(tri(ts,:)));
        F(tri(ts))=F(tri(ts))+Amplitude_source/area_t/3; %remove rho if you have
    end
    
    %% Global matrix for frequency w
    Kw=K+1i*w*C-w^2*M;
    
    %% Dirichlet conditions
    % if a node is subject to Neumann and Dirichlet conditions, Dirichlet
    % overrides Neumann condition.
    for ii=1:numel(nn)
        nodes=[n0(nn(ii));n1(nn(ii))];
        nbound=find(BCboundary==nb(nn(ii)));
        
        try
            bt=BCType(nbound);
            if(bt==1)
                bv=BCValue(nbound); % Z
                Kw(n0(nn(ii)),:)=0;
                Kw(n0(nn(ii)),n0(nn(ii)))=1;
                F(n0(nn(ii)))=bv;
                Kw(n1(nn(ii)),:)=0;
                Kw(n1(nn(ii)),n1(nn(ii)))=1;
                F(n1(nn(ii)))=bv;
            end
        catch
        end
    end      
    
    %% Solve
    P=Kw\F;
    
    %% save solution for frequency
    filename=strcat('results_f_',num2str(frequency));
    save(filename,'x','y','tri','P');
    
    %% Plot results
    if (plot_results(ifreq)==true)
        ax=figure;
        set(ax,'Name',strcat('Freq.(Hz)=',num2str(frequency)))
        subplot(2,2,1)
        trisurf(tri,x,y,real(P))
        view([0 90]);
        axis equal
        shading interp
        title('Real')
        subplot(2,2,2)
        trisurf(tri,x,y,imag(P))
        view([0 90]);
        axis equal
        shading interp
        title('Imag')
        subplot(2,2,3)
        trisurf(tri,x,y,abs(P))
        view([0 90]);
        axis equal
        shading interp
        title('Abs(P)')
        subplot(2,2,4)
        trisurf(tri,x,y,20*log10(abs(P)))
        view([0 90]);
        axis equal
        shading interp
        title('SPL(dB ref 1Pa)')
    end
    for ilines=1:2
        nodes_results{ilines,ifreq}=P(nodes_select{ilines});
    end
end

stop
%%
load results_f_500.mat

ax=figure;
set(ax,'Name',strcat('Time animation'))
f=500;
c=343;
w=2*pi*f;
T=1/f;
dt=T/50;
t=0;
limP=max(abs(P));
while 1
    trisurf(tri,x,y,real(P*exp(1i*w*t)))
    view([0 90]);shading interp;
    caxis([-limP limP]);
    colorbar
    axis equal
    drawnow    
    t=t+dt;
end

