%% ------------------------------------------------------------------------
% Practical work 1 code to calculate ...
function P=solution_1D(f,c1,c2,L1,L2,x)
     % f-frequency
     % c1,c2 – propagation velocities
     % L1,L2 – length of each médium
     % x – position of receivers
     w=2*pi*f;
     k1=w/c1;
     
     k2=w/c2;
     A=complex(zeros(4));
     A(1,:)=[1 exp(-1i*k1*L1) 0 0];
     A(2,:)=[exp(-1i*k1*L1)  exp(0)  -exp(0)  -exp(-1i*k2*(L2))];
     A(3,:)=[(-1i*k1)*exp(-1i*k1*L1)  -(-1i*k1)*exp(0)  -(-1i*k2)*exp(0)  (-1i*k2)*exp(-1i*k2*(L2))];
     A(4,:)=[0 0 (-1i*k2)*exp(-1i*k2*L2) -1*(-1i*k2)];
     B=[1;0;0;0];
     X=A\B;
     P=zeros(numel(x),1);
     for ii=1:numel(x)
         if(x(ii)<=L1 & x(ii)>=0)
             P(ii)=X(1)*exp(-1i*k1*x(ii))+X(2)*exp(-1i*k1*abs(x(ii)-L1));
         elseif(x(ii)<=L1+L2)
             P(ii)=X(3)*exp(-1i*k2*abs(x(ii)-L1))+X(4)*exp(-1i*k2*abs(x(ii)-L1-L2));
         else
             P(ii)=NaN;
         end 
     end
     return 
end