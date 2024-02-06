  
  % function masspar: defines density of the elements

  % Typically, in an engineering context we choose as fundamental units
  % a length and a force, e.g., mm and kN.  In this case the unit of mass
  % is the mass that a kN accelerates by 1 mm/s², that is 10^6 kg,
  % and the mass density is referred to a unit volume,
  % that is 1 mm³ = 10^-9 m³.  Before you proceed, what is the mass density
  % of water when the fundamental units are s, mm and kN?
  % Or, use m, N and kg.
 
 function [el_density]=masspar(nElements)
 
  % If an element is not associated with a material, 
  % then its unit mass remains(reasonably?) equal to zero.
  el_density = zeros(nElements,1);
  
  % (OTOH, if neighbouring elements have zero mass some DoF
  %  may have a zero mass, hence M is no more definite positive,
  %  so please be careful)
  
  % The number of different unit masses/material densities
  density = 21/9.80665*1000; %Kg/m3
  % must be equal to the number of groups of elements
  element_groups = {[1:nElements]};
  %                 [8 9 10 11],
  %                 [12 13 14 15 16 18 19 20]};
  % because they are going to be indexed by the same variable.
  
  for i = 1:length(density)
    el_density(element_groups{i}) = density(i);
  end
  
 end
