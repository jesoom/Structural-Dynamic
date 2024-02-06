%
% Function stressNodes: stress extrapolations at nodes
%
function [dSigmaNSF,dSigmaNav]=stressNodes(dSigma,nInc,nElements,dXY,nNodes,nGtot,dCsiEtaG,dWG)

   % Recovering Poisson Ratio From dSigma (forth column is \Sigma zz)
   dni=dSigma(1,4)/(dSigma(1,1)+dSigma(1,2));

  % Extrapolation at nodes with shape functions
  if (nGtot > 1)
    dH=zeros([nNodes,nNodes]);
    dP=zeros([nNodes,3]);
  
    for ne=1:nElements
      n14=nInc(ne,1:4);
      dXnodes=dXY(n14,1);
      dYnodes=dXY(n14,2);
      for ng=1:nGtot
        dxg=dCsiEtaG(ng,1);
        dyg=dCsiEtaG(ng,2);
        dPhi=[(1-dxg)*(1-dyg); (1+dxg)*(1-dyg); (1+dxg)*(1+dyg); (1-dxg)*(1+dyg)]/4;
        dPhidCsi=[-(1-dyg);  (1-dyg); (1+dyg); -(1+dyg)]/4;
        dPhidEta=[-(1-dxg); -(1+dxg); (1+dxg);  (1-dxg)]/4;

        dQmat=dPhidCsi*dPhidEta'-dPhidEta*dPhidCsi';
        ddJ=dXnodes'*dQmat*dYnodes;
     
        dH(n14,n14)=dH(n14,n14)+dWG(ng)*dPhi*dPhi'*abs(ddJ);
        dP(n14,:)=dP(n14,:)+dWG(ng)*dPhi*dSigma(ne,4*ng-3:4*ng-1)*abs(ddJ);
      end
    end
    dSigmaNSF=dH\dP;

    % Out of Plane Stress (\Sigma zz)
    dSigmaNSF=[dSigmaNSF, dni*(dSigmaNSF(:,1)+dSigmaNSF(:,2))];
    
    % Von Mises Stress    
%     dSigmaNSF=[dSigmaNSF,sqrt(dSigmaNSF(:,1).^2+dSigmaNSF(:,2).^2-dSigmaNSF(:,1).*dSigmaNSF(:,2)+3*dSigmaNSF(:,3).^2)];
  else
    dSigmaNSF=zeros([nNodes,4*nGtot]);
  end

  
  % Extrapolation at nodes with nodal average
  dSigmaNav=zeros(size(dSigmaNSF));
  if (nGtot > 1)
    for ne=1:nElements
      n14=nInc(ne,1:4);
      dXnodes=dXY(n14,1);
      dYnodes=dXY(n14,2);
      dH=zeros([4,4]);
      dP=zeros([4,3]);
      for ng=1:nGtot
        dxg=dCsiEtaG(ng,1);
        dyg=dCsiEtaG(ng,2);
        dPhi=[(1-dxg)*(1-dyg); (1+dxg)*(1-dyg); (1+dxg)*(1+dyg); (1-dxg)*(1+dyg)]/4;
        dPhidCsi=[-(1-dyg);  (1-dyg); (1+dyg); -(1+dyg)]/4;
        dPhidEta=[-(1-dxg); -(1+dxg); (1+dxg);  (1-dxg)]/4;

        dQmat=dPhidCsi*dPhidEta'-dPhidEta*dPhidCsi';
        ddJ=dXnodes'*dQmat*dYnodes;
     
        dH=dH+dWG(ng)*dPhi*dPhi'*abs(ddJ);
        dP=dP+dWG(ng)*dPhi*dSigma(ne,4*ng-3:4*ng-1)*abs(ddJ);
      end
      dSigmaNav(n14,1:3)=dSigmaNav(n14,1:3)+dH\dP;
    end
  elseif (nGtot==1)
    for ne=1:nElements
      n14=nInc(ne,1:4);
      dSigmaNav(n14,1:3)=dSigmaNav(n14,1:3)+ones([4,1])*dSigma(ne,1:3);
    end
  else
    disp('Fatal error')
    STOP
  end
  
  for nn=1:nNodes
    nElsNods=find(nn==nInc(:,1:4));
    nElsNods=length(nElsNods);
    dSigmaNav(nn,1:3)=dSigmaNav(nn,1:3)/nElsNods;
  end
  
  % Out of Plane Stress (\Sigma zz)
    dSigmaNav(:,4)=dni*(dSigmaNav(:,1)+dSigmaNav(:,2));
    
    % Von Mises Stress    
%   dSigmaNav(:,4)=sqrt(dSigmaNav(:,1).^2+dSigmaNav(:,2).^2-dSigmaNav(:,1).*dSigmaNav(:,2)+3*dSigmaNav(:,3).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
