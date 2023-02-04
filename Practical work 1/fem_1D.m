%% ------------------------------------------------------------------------
% Practical work 1 code to calculate with FEM
function P=fem_1D(f,c1,c2,L1,L2,n1,n2,rho)
    % f-frequency
    % c1,c2 – propagation velocities
    % L1,L2 – length of each medium
    % n1,n2 – number of nodes in each medium
    % rho - density of both mediums
    w=2*pi*f;
    dx1=L1/n1;dx2=L2/n2;
    K=zeros(n1+n2+1); % global stiffness matrix
    M=zeros(n1+n2+1);% global mass matrix
    ke=zeros(2,2); % element stiffness matrix
    me=zeros(2,2); % element mass matrix
    node_i=1:n1+n2; % left node of each element -- list of nodes
    node_f=2:n1+n2+1; % right node of each element
    
%     For each element, the corresponding stiffness matrix (ke) must be computed, 
% and then assembled into the global stiffness matrix (K) and mass matrix (M).
    for ii=1:(n1+n2)
        if(ii<=n1)
            Le=dx1;
            c=c1;
        else
            Le=dx2;
            c=c2;
        end
        ke=1/rho*[ 1/Le  -1/Le
                -1/Le   1/Le];
        me=1/(c^2*rho)*[ Le/3  Le/6
                Le/6   Le/3];
        K(node_i(ii):node_f(ii),node_i(ii):node_f(ii))=K(node_i(ii):node_f(ii),node_i(ii):node_f(ii))+ke;
        M(node_i(ii):node_f(ii),node_i(ii):node_f(ii))=M(node_i(ii):node_f(ii),node_i(ii):node_f(ii))+me;
    end
    KM=K-w^2*M;
    KM(1,:)=0;
    KM(1,1)=1;

    B=zeros(n1+n2+1,1);
%     B(1)=1i*w*1; % introducing conditions
    B(1)=1;
%   B(1)=1i*w*0;
    B(n1+n2+1)=0;
    
    P=KM\B;
    return 
end