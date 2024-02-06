%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Fem2D.m                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               Static elastic analysis for plane problems                %
%                      with 4-node finite elements                        %
%                                                                         %
%                       Author: Giuseppe COCCHETTI                        %
%                                                                         %
%                            version 19.05.2011                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       THIS SOFTWARE IS TO BE USED FOR LEARNING PURPOSES ONLY            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is made of a core routine transferring data to the following  %
% functions:                                                              %
%  - geotop: defines geometry and topology data;                          %
%  - mecpar: defines mechanical parameters;    
%  - masspar: defines density of the elements;
%  - locons: defines boundary conditions (loads and constraints);         %
%  - GaussCW: defines the coordinates and weights for Gauss integration;  %
%  - stiffm: generates the stiffness matrix for each element; 
%  - Melement: generates the mass matrix for each element;
%  - assilc: assigns boundary conditions (nodal loads and constraints);   %
%  - syssol: solves the equation system
%            computes eigenvalues and eigenvectors;                       %
%  - stress: computes stresses at Gauss points;                           %
%  - figcre: creates the window of a figure;                              %
%  - nodeconf: displays the nodes of the mesh;                            %
%  - memconf: displays the finite element mesh;                           %
%  - stressNodes: extrapolates stresses at nodes both with shape function %
%                          interpolation and the average of nodal values; %
%  - drawstress: plot stress contours;                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
format short e
disp('Static elastic analysis for plane problems')

%%%%%%%%%%%%%%%%%%%
% Initializations %
%%%%%%%%%%%%%%%%%%%

for only_for_closing=1:1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: acquisition of the structural data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geometry and topology data
[nInc,nElements,dXY,nNodes]=geotop;

% Mechanical parameters
[dPar]=mecpar;
[dMas]=masspar(nElements);

% Thickness
thickness=1;

% Total number of Gauss points [1(=1x1),4(=2x2),9(=3x3),16(=4x4)]
nGtot=4;
[dCsiEtaG,dWG]=GaussCW(nGtot);

% Boundary conditions: loads and constraints
[nCons,dC,nForce,dF,npq,dpq]=locons(dXY,thickness);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the solving system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total number of dofs (before imposing boundary conditions)
nDofTot=max(max(nInc(:,5:12)));

% Global stiffness and mass matrix
dK=zeros([nDofTot,nDofTot]); 
dM=zeros([nDofTot,nDofTot]);

% Set up and assembly of the stiffness and mass matrices
% of each element
for ne=1:nElements
    % Recover the connection numbers
    % n1=nInc(ne,1);  % First node of ne-th element
    % n2=nInc(ne,2);  % Second node of ne-th element
    % n3=nInc(ne,3);  % Third node of ne-th element
    % n4=nInc(ne,4);  % Fourth node of ne-th element
    % n14=[n1,n2,n3,n4];
    % All together:
    n14=nInc(ne,1:4);
    
    % Recover the nodal coordinates for n-th element
    % dXn1=dXY(n1,1);
    % dYn1=dXY(n1,2);
    % dXn2=dXY(n2,1);
    % dYn2=dXY(n2,2);
    % dXn3=dXY(n3,1);
    % dYn3=dXY(n3,2);
    % dXn4=dXY(n4,1);
    % dYn4=dXY(n4,2);
    % dXY14=[dXn1,dYn1; 
    %        dXn2,dYn2; 
    %        dXn3,dYn3; 
    %        dXn4,dYn4];
    % All together:
    dXY14=dXY(n14,:);
    
    
    % Stiffness matrix dKne for the ne-th element
    [dKne]=stiffm(dPar,dXY14,thickness,nGtot,dCsiEtaG,dWG);
    % Mass matrix dMne --- nb dMne is 4x4 because m_xy = 0
    [dMne] = Melement(dMas(ne),dXY14,thickness,nGtot,dCsiEtaG,dWG);
    
    if ne==1
        Kmat_1 =dKne;
    elseif ne==2
        Kmat_2=dKne;
    elseif ne==3
        Kmat_3=dKne;
    elseif ne==4
        Kmat_4=dKne;
    elseif ne==5
        Kmat_5=dKne;
     elseif ne==6
        Kmat_6=dKne;
    elseif ne==7
        Kmat_7=dKne;
    elseif ne==8
        Kmat_8=dKne;
    elseif ne==9
        Kmat_9=dKne;    
    elseif ne==10
        Kmat_10=dKne;
    elseif ne==11
        Kmat_11=dKne;        
    end
    
    % Assembly of the overall stiffness matrix
    nVne=[nInc(ne,5:12)]; % Recovers the Dofs of the ne-th element
    dK(nVne,nVne)=dK(nVne,nVne)+dKne; % Global stiffness matrix dK; Element stiffness matrix dKne
    
    % Assembly of the overall mass matrix
    nVnx = nVne(1:2:end); % x oriented element's DoF
    dM(nVnx,nVnx)=dM(nVnx,nVnx)+dMne;
    nVny = nVne(2:2:end); % y oriented element's DoF
    dM(nVny,nVny)=dM(nVny,nVny)+dMne;
    
end

% Boundary conditions: nodal constraints and loads
[nUs,dUs,nUu,dT]=assilc(nInc,nForce,dF,nCons,dC,npq,dpq,dXY,thickness,nDofTot,nGtot,dCsiEtaG,dWG);

% Solution of linear system
[du,dR, evals, evecs]=syssol(dK,dM,dT,nUu,nUs,dUs,nDofTot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT 1: computational results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Stress computation at Gauss points
% % dSigma(ne,:)=[dsx_1,dsy_1,dtxy_1,dsMises_1,dsx_2,dsy_2,dtxy_2,dsMises_2,...,dsx_nGtot,dsy_nGtot,dtxy_nGtot,dsMises_nGtot]
% [dSigma]=stress(du,dPar,nInc,nElements,dXY,nGtot,dCsiEtaG);
% 
% % Stress extrapolations at nodes
% [dSigmaNSF,dSigmaNav]=stressNodes(dSigma,nInc,nElements,dXY,nNodes,nGtot,dCsiEtaG,dWG);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT 2: graphical representation of the computational results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Characteristic structural dimensions
dXmin=min(dXY(:,1));
dXmax=max(dXY(:,1));
dYmin=min(dXY(:,2));
dYmax=max(dXY(:,2));
dimXmax=dXmax-dXmin;
if (dimXmax==0)
    dimXmax=1;
end
dimYmax=dYmax-dYmin;
if (dimYmax==0)
    dimYmax=1;
end

% % FIGURE 1: UNDEFORMED STRUCTURE
% New window
% nF=Figure number
% vCoordFig =  horizontal & vertical coordinate (in pixel) of the lower left window corner;
%              horizontal & vertical coordinate (in pixel) of the upper right window corner;
% vCoordAxes = horizontal & vertical coordinate (in pixel) of the lower left figure corner; 
%              horizontal & vertical coordinate (in pixel) of the upper right figure corner;
vCoordFig=[10,140,520,320]; 
vCoordAxes=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10];
cTit='Structural Scheme';
figcre(1,vCoordFig,vCoordAxes,cTit);

% Hinge radius (N.B.: graphical convention only!)
dr=max(dimXmax,dimYmax)/100;

% Drawing of nodes
nodeconf(nNodes,dXY,dr);

% Drawing of members
bElemLabels=0;
memconf(bElemLabels,nElements,nInc,dXY,'b-');

% FIGURE 2: DEFORMED AND UNDEFORMED STRUCTURE
% New window
vCoordFig=[500,10,520,320]; 
vCoordAxes=[dXmin-dimXmax/8,dXmax+dimXmax/8,dYmin-dimYmax/8,dYmax+dimYmax/8]; 
cTit='Deformed Configuration';
figcre(2,vCoordFig,vCoordAxes,cTit);

bElemLabels=0;
% Drawing of undeformed structure
memconf(bElemLabels,nElements,nInc,dXY,'k:');

% Amplification factor for the graphical representation of displacements 
dSmax=max(dimXmax,dimYmax)/20;
dUmax=norm(du,inf);
dAmplif=10^(ceil(log10(dSmax/dUmax)))/4;
text(dXmin,dYmax+dimYmax/5,sprintf('Displacement amplification factor: %0.5g',dAmplif))

% Drawing of deformed structure
dXYa=[dXY(:,1)+dAmplif*du(1:2:end,1),dXY(:,2)+dAmplif*du(2:2:end,1)];
memconf(bElemLabels,nElements,nInc,dXYa,'b-');

% Drawing of modal shapes
vCoordFig=[510,20,520,320]; 
for mode = 1:6
    devec = evecs(:,mode);
    figcre(20+mode,vCoordFig,vCoordAxes,sprintf(['\\fontname{Times New Roman}\\fontsize{23} \\it Mode %2d,  T_{\\fontsize{12}%2d} = %.2f ms, \\omega_{\\fontsize{12}%2d}^{\\fontsize{12}2} = %.2f rad^{\\fontsize{12}2}/s^{\\fontsize{12}2} , \\omega_{\\fontsize{12}%2d} = %.2f rad/s, f_{\\fontsize{12}%2d} = %.2f Hz'], mode,mode,(2000*pi)/sqrt(evals(mode)), mode,(evals(mode)) ,mode, sqrt(evals(mode)), mode, sqrt(evals(mode))/(2*pi)));
    vCoordFig(1:2) = vCoordFig(1:2)+10*ones(1,2); 
    dSmax=max(dimXmax,dimYmax)/20;
    dUmax=norm(devec,inf);
    dAmplif=4^ceil(log10(dSmax/dUmax));
    %text(dXmin,dYmax+dimYmax/5,sprintf('Displacement amplification factor: %0.5g',dAmplif))
    
    % Drawing of undeformed structure and modal shape
    memconf(bElemLabels,nElements,nInc,dXY,'k:');
    dXYe1=[dXY(:,1)+dAmplif*devec(1:2:end,1),dXY(:,2)+dAmplif*devec(2:2:end,1)];
    memconf(bElemLabels,nElements,nInc,dXYe1,'m-');
    xlabel(sprintf('\\fontname{Times New Roman} \\fontsize{18} X [m]'))
    ylabel(sprintf('\\fontname{Times New Roman} \\fontsize{18} Y [m]'))
end
