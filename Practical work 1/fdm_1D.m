%% ------------------------------------------------------------------------
% Practical work 1 code to calculate with FDM
function P=fdm_1D(f,c1,c2,L1,L2, p,v,npts,rho)
     % f-frequency
     % c1,c2 – propagation velocities
     % L1,L2 – length of each médium
     % x – position of receivers
     w=2*pi*f;
     dx=(L1+L2)/(npts-1);
     x=0:dx:dx*(npts-1);
     K1=w/c1;
     K2=w/c2;
     pressure_left=p; %Dirichlet condition
     velocity_right=v; %neumann condition

     A=zeros(npts);
     B=zeros(npts,1);
%      boundary condition (Dirichlet, prescribed pressure) on the left-most node:
     A(1,1)=1;
     B(1)=pressure_left;
%      on the right-most node, the Neumann boundary condition  needs to be imposed, approximating the first spatial derivative by backward finite differences
     A(end,end-1:end)=-1/1i/rho/w*[-1/dx 1/dx];
     B(end)=velocity_right;
% the corresponding coefficients are now precalculated as:
     k1=1/dx^2;
     k2=-2/dx^2+K1^2;
     k2_2=-2/dx^2+K2^2;
     k3=1/dx^2;

%      filling matrix A
     for ii=2:npts-1
         if(x(ii)<=L1)
            A(ii,ii-1:ii+1)=[k1 k2 k3];
         else
            A(ii,ii-1:ii+1)=[k1 k2_2 k3];
         end
     end
% the independent term vector is completed:
     B(2:end-1)=0;
%      Solution of the equation system  yields the temperatures throughout the domain:
     P=A\B;
     return 
end