%
% Function stress: computes stresses at Gauss points
%
function [dSigma]=stress(du,dPar,nInc,nElements,dXY,nGtot,dCsiEtaG)

  dE=dPar(1);
  dni=dPar(2);
 % %Plane Stress
 %  Ep=dE; nip=dni;
  
    %Plane Strain
    Ep=dE/(1-dni^2); nip=dni/(1-dni);
    
    dEmat=inv((1/Ep)*[      1, -nip,         0;
                         -nip,    1,         0;
                            0,    0, 2*(1+nip)]);
 
  dSigma=zeros([nElements,4*nGtot]);

  for ne=1:nElements
    n14=nInc(ne,1:4);
    dXnodes=dXY(n14,1);
    dYnodes=dXY(n14,2);

    nVne=nInc(ne,5:12);
    dune=du(nVne);
  
    for ng=1:nGtot
      dxg=dCsiEtaG(ng,1);
      dyg=dCsiEtaG(ng,2);
    % dPhi=[(1-dxg)*(1-dyg); (1+dxg)*(1-dyg); (1+dxg)*(1+dyg); (1-dxg)*(1+dyg)]/4;
      dPhidCsi=[-(1-dyg);  (1-dyg); (1+dyg); -(1+dyg)]/4;
      dPhidEta=[-(1-dxg); -(1+dxg); (1+dxg);  (1-dxg)]/4;

      dQmat=dPhidCsi*dPhidEta'-dPhidEta*dPhidCsi';
      ddJ=dXnodes'*dQmat*dYnodes;
     
      dBne=zeros([3,8]);
      dBne(1,1:2:end)=-dYnodes'*dQmat;   
      dBne(2,2:2:end)= dXnodes'*dQmat;   
      dBne(3,1:2:end)=dBne(2,2:2:end);   
      dBne(3,2:2:end)=dBne(1,1:2:end);   
      dBne=dBne/ddJ;

      dSigma(ne,4*ng-3:4*ng-1)=(dEmat*dBne*dune)';    
      
      % Out of Plane Stress (\Sigma zz)
      dSigma(ne,4*ng)=dni*(dSigma(ne,4*ng-3)+dSigma(ne,4*ng-2));
      
      % Vin Mises Stress Definition
%      dSigma(ne,4*ng)=sqrt(dSigma(ne,4*ng-3:4*ng-2)*dSigma(ne,4*ng-3:4*ng-2)'...
%           -dSigma(ne,4*ng-3:4*ng-3)*dSigma(ne,4*ng-2:4*ng-2)+3*dSigma(ne,4*ng-1:4*ng-1)^2);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
