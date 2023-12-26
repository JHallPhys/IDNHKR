clear all
close all

% parallel.gpu.enableCUDAForwardCompatibility(true) % RTX 3090 ;)


%==========================================================================
%   Quantum System and Schur Parameters
%==========================================================================

%==========================================================================
% g4=0.003
%==========================================================================
% System Parameters
N_1=125; 
% N_1=500; 
% N_1=2000; 
N = 2*N_1; % Hilbert space dimension
K_class =10; % Classical Kicking 
T=2*pi/N; % Effective hbar
kick = K_class/T; % Quantum Kicking
str_ext='.mat'
hole_lower=0.3;
hole_upper=0.6;
hole_lower=N*hole_lower;
hole_upper=N*hole_upper;
R=0.2% This is the reflection strength R=exp(-gamma)*I_{MxM: M<N}
gamma=0.22 % gamma_nat
% gamma=0.35
% gamma=0.48 % gamma_typ
% gamma=0.75
% gamma=0.88 % gamma_inv
dgamma=0.005
% dgamma=0.02
%==========================================================================
%   Matrix Construction and Schur decompesition
%==========================================================================


U=zeros(N,N); % Initialise Flouqet matrix
tic
U=UCheck(N,N_1,K_class,T,R,hole_lower,hole_upper,str_ext);% Check if matrix exists, if it does load it, else make and save it
[psiS,En]=ECheck(U,N,N_1,K_class,T,R,str_ext);% Check if matrix exists, if it does load it, else make and save it
toc


%==========================================================================
%   Project onto Subspace of stability
%==========================================================================

% [psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values

Es=diag(En);

% [psi_2,n_efn]=Psi_lifetime(psiS,Es,eps,set_stability);
%==========================================================================
% Find eigenfuction with decay rate in a given window gamma \pm dgamma
%==========================================================================

% figure(1)
% clf
% hold on
% g1=plot((real(Es)),imag(Es),'k.','Markersize',5);
% 
% title('\epsilon_n=\theta_n-i\gamma_n')


E=1i*log(diag(En));
% return

% figure(2)
% clf
% hold on
% g1=plot((real(E)+pi)./(2*pi),-2*imag(E),'k.','Markersize',10);
% g2=plot(linspace(0,1,100),gamma+dgamma,'b.-'); 
% g2=plot(linspace(0,1,100),gamma-dgamma,'b.-'); 
% xlabel('\theta_n')
% ylabel('-\gamma_n')
% title('\epsilon_n=\theta_n-i\gamma_n')
% axis([0 1 0 1.5])


figure
clf
hold on
g1=plot((wrapTo2Pi(real(E))./(2*pi)),-2*imag(E),'k.','Markersize',10);
g2=plot(linspace(0,1,100),0.88,'g.-'); 
g2=plot(linspace(0,1,100),0.48,'b.-');
g3=plot(linspace(0,1,100),0.22,'r.-'); 
xlabel('\theta_n')
ylabel('-\gamma_n')
title('\epsilon_n=\theta_n-i\gamma_n')
axis([0 1 0 1.1])
set ( gca, 'xdir', 'reverse' )
% return

psi_gamma=get_lifetime_window(psiS,E,gamma,dgamma);

% return

n_efn=size(psi_gamma,2)
% return
% return
% Discrete Phase Space Grid
% n_efn=1
Hus=zeros(N,N,n_efn);
q = linspace(0,1,N); % q interval
p = linspace(0,1,N); % p interval
[qmesh,pmesh]=meshgrid(q,p); 
z = (qmesh+1i*(pmesh)); 
Hus = zeros(N,N); % Average Husimi function array
cs=zeros(N,N);
norm_cs= (2/N)^0.25; % Normalisation constant for the coherent state
tic
% return
for itt = 1:N-1
    itt
   
    cs=Cs_create_component(itt,norm_cs,N,q,z,cs); 
    

    
    phi2(1,1,:)=psi_gamma(itt,:);
    Hus=Hus+conj(cs).*phi2;
    

    cs(:,:)=0;
end
toc


Hus_av=zeros(N,N);

for itt=1:n_efn

Hus_av=Hus_av+abs(Hus(:,:,itt)).^2;

figure
clf
imagesc(q,p,abs(Hus(:,:,itt)).^2)
% colorbar
% caxis([0 1])
colormap(flipud(inferno))
title(strcat('state',num2str(itt)))
set(gca,'YDir','normal')

end


return

viridis=viridis();

figure
clf
imagesc(q,p,Hus_av)
colorbar
colormap(flipud(inferno))
set(gca,'YDir','normal')

Hus_Entropy=-Hus_av.*log(Hus_av);

% figure
% clf
% imagesc(q,p,Hus_Entropy)
% colorbar
% % caxis([0 1])
% colormap(viridis)
% set(gca,'YDir','normal')
return
fname=fname_husimi_single_efn_special(K_class,N,gamma,1,str_ext);
parent_d = cd;  
cd './Husimi_dat' % Directory where matrix is stored
save(fname,'Hus_Entropy'); % save it 
cd(parent_d)


















