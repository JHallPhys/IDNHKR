% Function to get the states in a given lifetime window.
function [psi_gamma]=get_lifetime_window(psiS,E,gamma,dgamma)


E=-2*imag(E); % Get the imginary parts of the quasienregies
Ewindow=find(E>gamma-dgamma & E<gamma+dgamma); % Find the energies in the interval
size(Ewindow);
psi_gamma=psiS(:,Ewindow); % Get the closest eigenfunctions insorted
EW=E(Ewindow); 
[delta_gamma,index_close]=sort(abs(EW-gamma));
delta_gamma
psi_gamma=psi_gamma(:,index_close);


end
