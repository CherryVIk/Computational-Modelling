%% ------------------------------------------------------------------------
% FEM code to calculate the element stiffness and mass matrices
function [ke,me]=keme2D(fx, fy, rho, c)
    
	x1=fx(1);y1=fy(1);
    x2=fx(2);y2=fy(2);
    x3=fx(3);y3=fy(3);
    
    b1=y2-y3;
    b2=y3-y1;    
    b3=y1-y2;
    c1=x3-x2;   
    c2=x1-x3;   
    c3=x2-x1;   
    
    y32=y3-y2;
    y13=y1-y3;
    y21=y2-y1;

    x23=x2-x3;
    x31=x3-x1;
    x12=x1-x2;

    Ael=abs(x2*y3+x3*y1+x1*y2-x2*y1-x3*y2-x1*y3)/2;
    ke=1/rho*1/4/Ael*[ b1^2+c1^2       b1*b2+c1*c2     b1*b3+c1*c3
                        b1*b2+c1*c2     b2^2+c2^2       b2*b3+c2*c3
                        b1*b3+c1*c3     b2*b3+c2*c3     b3^2+c3^2];
    me=Ael/12/rho/c^2*[2 1 1
                     1 2 1
                     1 1 2];
    %     me=1*Ael/3*eye(3)/rho/c^2;  % uncomment for lumped mass
end