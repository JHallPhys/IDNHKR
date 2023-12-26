% Construct the Floquet operator of size(NxN) for the kicked rotor.

function [U_out,time_out]=UMatrixPartial(U_in,N,N_1,K_class,T,hole_lower,hole_upper,R)

%==========================================================================
% One way
%==========================================================================

K_s = K_class/T; % Quantum Kicking
j_i=linspace(1,N,N);
j_i=j_i-1;
j_i=transpose(j_i);
l1=linspace(1,N,N);
RS=ones(1,N);
RS(1,hole_lower:hole_upper)=sqrt(R);
tic
for k = 1:N


                b=exp(-1i*T*0.5*j_i.^2+2*pi*1i*j_i.*(k-l1)/N-1i*0.5*K_s*cos(2*pi*l1/N));  
                U_in(k,:) =  RS.*(exp(-0.5*1i*K_s*cos(2*pi*k/N))*sum(b))/N;  
       
             
end



U_out=U_in;
time_out=toc;


%==========================================================================
% The Other Way
%==========================================================================

% 
% K_s = K_class/T; % Quantum Kicking
% j_i=linspace(1,N,N);
% j_i=j_i-1;
% j_i=transpose(j_i);
% l1=linspace(1,N,N);
% tic
% for k = 1:N
% 
% 
%            
% 
% 
%                 b=exp(-1i*T*0.5*j_i.^2+2*pi*1i*j_i.*(k-l1)/N-1i*0.5*K_s*cos(2*pi*l1/N));   
%                 U_in(k,:) =  (exp(-0.5*1i*K_s*cos(2*pi*k/N))*sum(b))/N; 
%              
% end
% 
% 
% Rop=eye(N);
% Rop(hole_lower:hole_upper,hole_lower:hole_upper)=(R)
% U_out=U_in*sqrtm(Rop);
% time_out=toc;
% 
% 
% end

%==========================================================================
% Other Map
%==========================================================================


