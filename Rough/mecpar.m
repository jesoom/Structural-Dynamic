%
% Function mecpar: defines mechanical parameters
%
function [dPar]=mecpar

 % Matrix dPar:
 % collects Young's modulus and Poisson's coefficient for the material;
 % dPar=[E, ni]

  dE=4e9; %N/m2
  dni=.28;
  dPar=[dE, dni];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
