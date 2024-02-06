%
% Function stiffm: generates the stiffness matrix for each element
%
function [dKne]=stiffm(dPar,dXY14,thickness,nGtot,dCsiEtaG,dWG)

  % Nodal coordinates
  dXnodes=dXY14(:,1);
  dYnodes=dXY14(:,2);
  
  % Stiffness matrix
  dKne=zeros([8,8]);
  for ng=1:nGtot
    dxg=dCsiEtaG(ng,1);
    dyg=dCsiEtaG(ng,2);
  % dPhi=[(1-dxg)*(1-dyg); (1+dxg)*(1-dyg); (1+dxg)*(1+dyg); (1-dxg)*(1+dyg)]/4;
    dPhidCsi=[-(1-dyg);  (1-dyg); (1+dyg); -(1+dyg)]/4;
    dPhidEta=[-(1-dxg); -(1+dxg); (1+dxg);  (1-dxg)]/4;

    dQmat=dPhidCsi*dPhidEta'-dPhidEta*dPhidCsi';
    ddJ=dXnodes'*dQmat*dYnodes;
     
    dE=dPar(1);
    dni=dPar(2);
    
%     %Plane Stress
%     Ep=dE; nip=dni;
    
    
%Plane Strain
    Ep=dE/(1-dni^2); nip=dni/(1-dni);
    
    dEmat=inv((1/Ep)*[      1, -nip,         0;
                         -nip,    1,         0;
                            0,    0, 2*(1+nip)]);
    
    dBne=zeros([3,8]);
    dBne(1,1:2:end)=-dYnodes'*dQmat;   
    dBne(2,2:2:end)= dXnodes'*dQmat;   
    dBne(3,1:2:end)=dBne(2,2:2:end);   
    dBne(3,2:2:end)=dBne(1,1:2:end);   
    dBne=dBne/ddJ;
    
    dKne=dKne+dWG(ng)*thickness*dBne'*dEmat*dBne*abs(ddJ);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
