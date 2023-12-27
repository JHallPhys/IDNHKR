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
% N_1=125; 
N_1=500; 
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
% gamma=0.22 % gamma_nat
% gamma=0.35
gamma=0.48 % gamma_typ
% gamma=0.75
% gamma=0.88 % gamma_inv
dgamma=0.0118
% dgamma=0.02

%==========================================================================
%  Spectral Stuff 
%==========================================================================

% Diagonalisation
U=zeros(N,N); % Initialise Flouqet matrix
U=UCheck(N,N_1,K_class,T,R,hole_lower,hole_upper,str_ext);% Check if matrix exists, if it does load it, else make and save it
[psiS,En]=ECheck(U,N,N_1,K_class,T,R,str_ext);% Check if matrix exists, if it does load it, else make and save it
toc
Es=diag(En);
E=1i*log(diag(En));
"sort out this bit"
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

%==========================================================================
%  Husimi Stuff 
%==========================================================================

psi_gamma=get_lifetime_window(psiS,E,gamma,dgamma);
n_efn=size(psi_gamma,2)
pause(1)
[q,p,z,dz]=get_husimi_grid(N);
Hus=get_husimi(N,n_efn,q,p,z,psi_gamma);
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


str_title='Average_nefn_:';
str_title_num=num2str(n_efn);
str_gamma='_dgamma_';
str_gamma_num= strrep(num2str(dgamma),'.','p');
str_average_title=strcat(str_title,str_title_num,str_gamma,str_gamma_num)

figure
clf
imagesc(q,p,Hus_av./n_efn)
% colorbar
% caxis([0 1])
colormap(flipud(inferno))
title(str_average_title)
set(gca,'YDir','normal')

%==========================================================================
%  Entropy Stuff 
%==========================================================================

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


















