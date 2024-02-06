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

% % FIGURES 3,4,5,6: Stress contour plots
% % New window
% vCoordFig=[500,360,520,320]; 
% vCoordAxes=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10]; 
% bDrawMesh=1;
% bDrawStressSF=0; % bDrawStressSF=1 means stresses extrapolated at nodes with shape functions
%                  % bDrawStressSF=0 means stresses averaged at nodes
% 
% disp(' ..:: For Static Solition: ::..')                
% if (nGtot == 1)
%     disp('Stress extrapolation will be averaged at nodes')
%     bDrawStressSF=0;
% else
%     if (bDrawStressSF == 1)
%         disp('Stress extrapolation at nodes is done with shape functions')
%     else
%         disp('Stress extrapolation at nodes is done by the average of nodal values')
%     end
% end
% 
% for bNStress=1:4
%     if (bNStress == 1)
%       cTit='\sigma_{X}';
%     elseif (bNStress == 2)
%       cTit='\sigma_{Y}';
%     elseif (bNStress == 3)
%       cTit='\tau_{XY}';
%     elseif (bNStress == 4)
%       cTit='\sigma_{Z}';
%     end
%     figcre(2+bNStress,vCoordFig,vCoordAxes,cTit);
%     drawstress(bDrawMesh,bDrawStressSF,bNStress,dSigmaNSF,dSigmaNav,nInc,nElements,dXY);
% end
% 
% % modal stresses
% for mode = [1:6]
%     evc = evecs(:,mode);
%     [eSigma]=stress(evc,dPar,nInc,nElements,dXY,nGtot,dCsiEtaG);
%     [eSigmaNSF,eSigmaNav]=stressNodes(eSigma,nInc, nElements,dXY,nNodes,nGtot,dCsiEtaG,dWG);                             
%     vCoordFig=[550,300,520,320]; 
%     vCoordAxes=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10]; 
%     bDrawMesh=1;
%     bDrawStressSF=1;
%     
%     disp(' ..:: For Dynamic Solition (Modal Stress): ::..') 
%     if (nGtot == 1)
%         disp('Stress extrapolation will be averaged at nodes')
%         bDrawStressSF=0;
%     else
%         if (bDrawStressSF == 1)
%             disp('Stress extrapolation at nodes is done with shape functions')
%         else
%             disp('Stress extrapolation at nodes is done by the average of nodal values')
%         end
%     end
%     for bNStress=1:1
%         if (bNStress == 1)
%            eTit = sprintf(['\\sigma_{X}' , ' for mode %d, T = %.3f ms'],mode, 2000*pi/sqrt(evals(mode)));
%         elseif (bNStress == 2)            
%            eTit = sprintf(['\\sigma_{Y}' , ' for mode %d, T = %.3f ms'],mode, 2000*pi/sqrt(evals(mode)));
%         elseif (bNStress == 3)
%            eTit = sprintf(['\\tau_{XY}' , ' for mode %d, T = %.3f ms'],mode, 2000*pi/sqrt(evals(mode)));
%         elseif (bNStress == 4)               
%            eTit = sprintf(['\\sigma_{Z}' , ' for mode %d, T = %.3f ms'],mode, 2000*pi/sqrt(evals(mode)));
%         end     
%         figcre(200+mode,vCoordFig,vCoordAxes,eTit);
%         drawstress(bDrawMesh,bDrawStressSF,bNStress,eSigmaNSF,eSigmaNav,nInc,nElements,dXY);                 
%     end
% end

%%%%%%%%%%%%%
% The end %%
%%%%%%%%%%%%%

end

%***********************************************************************
%% ************ Lab2 : Modal ContriBution Factors ****************
%***********************************************************************

for only_for_closing=1:1

r_vertic=ones(nDofTot,1);
r_vertic(1:2:end)=0;

r_horiz=ones(nDofTot,1);
r_horiz(2:2:end)=0;

%....:::: Total Static Responses ::::....
% Horizontal Excitation - total static response
dT_tot_horiz=dM*r_horiz;
[du_tot_horiz,dR_tot_horiz, ~, ~]=syssol(dK,dM,dT_tot_horiz,nUu,nUs,dUs,nDofTot);

% Vertical Excitation - total static response
dT_tot_vertic=dM*r_vertic;
[du_tot_vertic,dR_tot_vertic, ~, ~]=syssol(dK,dM,dT_tot_vertic,nUu,nUs,dUs,nDofTot);


%.....:::: Modal Contrib. Responses ::::....
[~,ncol]=size(evecs);
M_star=diag(evecs'*dM*evecs).*eye(ncol);

GammaVec_horiz=zeros(ncol,1);
GammaVec_vertic=zeros(ncol,1);
for i=1:ncol
    GammaVec_horiz(i)=(evecs(:,i)'*dM*r_horiz)/M_star(i,i);
    GammaVec_vertic(i)=(evecs(:,i)'*dM*r_vertic)/M_star(i,i);
end

GammaVec_horiz=GammaVec_horiz.*eye(ncol);
GammaVec_vertic=GammaVec_vertic.*eye(ncol);

R_horiz =dM*evecs*GammaVec_horiz;
R_vertic=dM*evecs*GammaVec_vertic;

ncol=15; %only for 6 modes
horizdisp_72=zeros(ncol,1);
BS_horiz=zeros(ncol,1);
BendM_horiz=zeros(ncol,1);

for i=1:ncol % for mode i
    
    % Horiz Excitation: Horizontal Disp of top Node (Node: 133, dof: 265)
    [du_horiz,dR_horiz, ~, ~]=syssol(dK,dM,R_horiz(:,i),nUu,nUs,dUs,nDofTot);
    horizdisp_72(i,1)=du_horiz(143,1);
    
    % Horiz Excitation: Base Shear (Nodes:1,2,3,4, 273,274,275,276)
    BS_horiz(i)=dR_horiz(1)+dR_horiz(3)+dR_horiz(5)+dR_horiz(7)+...
        dR_horiz(9)+dR_horiz(11);
    
    %Vertical Excitation: Bending Moment in top cross section
    [du_horiz,~, ~, ~]=syssol(dK,dM,R_horiz(:,i),nUu,nUs,dUs,nDofTot);
    duV_1=  du_horiz([1	2	3	4	15	16	13	14],1);
    duV_2= du_horiz([3	4	5	6	17	18	15	16],1);
    duV_3= du_horiz([5	6	7	8	19	20	17	18],1);
    duV_4= du_horiz([7	8	9	10	21	22	19	20],1);
    duV_5= du_horiz([9	10	11	12	23	24	21	22],1);
    
    dR1 =Kmat_1*duV_1 ;
    dR2=Kmat_2 *duV_2;
    dR3=Kmat_3 *duV_3;
    dR4=Kmat_4 *duV_4;
    dR5=Kmat_5 *duV_5;
    

    
    dArm=dXY(2:6,1)';
    
    dRvec=[dR1(4)+dR2(2) , dR2(4)+dR3(2) , dR3(4)+dR4(2) , dR4(4)+dR5(2), dR5(4)];

    BendM_horiz(i)=sum(dArm.*dRvec);

end

% Modal Contrib. Factor : Horiz dips top node
MCF_hor_disp_72=horizdisp_72/du_tot_horiz(143,1);

% Modal Contrib. Factor : Base Shear
MCF_BS=BS_horiz./(dR_tot_horiz(1)+dR_tot_horiz(3)+dR_tot_horiz(5)+dR_tot_horiz(7)+...
    dR_tot_horiz(9) + dR_tot_horiz(11));

% Modal Contrib. Factor : Bending Moment
duV_1=  du_tot_horiz([1	2	3	4	15	16	13	14],1);
duV_2= du_tot_horiz([3	4	5	6	17	18	15	16],1);
duV_3= du_tot_horiz([5	6	7	8	19	20	17	18],1);
duV_4= du_tot_horiz([7	8	9	10	21	22	19	20],1);
duV_5= du_tot_horiz([9	10	11	12	23	24	21	22],1);

dR1 =Kmat_1*duV_1 ;
dR2=Kmat_2 *duV_2;
dR3=Kmat_3 *duV_3;
dR4=Kmat_4 *duV_4;
dR5=Kmat_5 *duV_5;
    
dArm=dXY(2:6,1)';

dRvec=[dR1(4)+dR2(2) , dR2(4)+dR3(2) , dR3(4)+dR4(2) , dR4(4)+dR5(2), dR5(4)];
    
BendM_tot_horiz=sum(dArm.*dRvec);

MCF_BendM=BendM_horiz./BendM_tot_horiz;


format long
disp('..:: Sum MCF Horiz disp top Node ::..')
disp(sum(MCF_hor_disp_72))

disp('..:: Sum MCF Base Shear ::..')
disp(sum(MCF_BS))

disp('..:: Sum MCF Bneding Moment at Top Cross section ::..')
disp(sum(MCF_BendM))

end

cumMCF_HorizDisp=zeros(1,ncol);
cumMCF_Shear=zeros(1,ncol);
cumMCF_Bend=zeros(1,ncol);

for iii=1:ncol
    cumMCF_HorizDisp(iii)=sum(MCF_hor_disp_72(1:iii));
    cumMCF_Shear(iii)=sum(MCF_BS(1:iii));
    cumMCF_Bend(iii)=sum(MCF_BendM(1:iii));

end

myModes=1:ncol;
figure(123)
cumMCF_HorizDisp_p=cumMCF_HorizDisp*100;
plot(myModes(1,1:6),cumMCF_HorizDisp_p(1,1:6),myModes(1,6:15),cumMCF_HorizDisp_p(1,6:15),'--', 'LineWidth', 2.5);
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% Set x-label and y-label
xlabel('Number of Modal', 'FontName', 'Times New Roman', 'FontSize', 18);
ylabel('Comulative Contribution Factor (%)', 'FontName', 'Times New Roman', 'FontSize', 18);
title('Modal Contribution Factor of Horizontal Displacement Node-A', 'FontName', 'Times New Roman', 'FontSize', 18);
hold on
plot(myModes(1,6),cumMCF_HorizDisp_p(1,6), 'ro', 'MarkerSize', 10)
text(myModes(1,6),cumMCF_HorizDisp_p(1,6), num2str(cumMCF_HorizDisp_p(1,6)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
plot(myModes(1,15),cumMCF_HorizDisp_p(1,15), 'ro', 'MarkerSize', 10)
text(myModes(1,15),cumMCF_HorizDisp_p(1,15), num2str(cumMCF_HorizDisp_p(1,15)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
xticks(myModes);
%  %--------------------------------------%
figure(1234)
cumMCF_Shear_p=cumMCF_Shear*100;
plot(myModes(1,1:6),cumMCF_Shear_p(1,1:6),myModes(1,6:15),cumMCF_Shear_p(1,6:15),'--', 'LineWidth', 2.5);
% % add percentage to the y-values
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% Set x-label and y-label
xlabel('Number of Modal', 'FontName', 'Times New Roman', 'FontSize', 18);
ylabel('Comulative Contribution Factor (%)', 'FontName', 'Times New Roman', 'FontSize', 18);
title('Modal Contribution Factor of Shear Force at The Base of The Dam Slice', 'FontName', 'Times New Roman', 'FontSize', 18);
hold on
plot(myModes(1,6),cumMCF_Shear_p(1,6), 'ro', 'MarkerSize', 10)
text(myModes(1,6),cumMCF_Shear_p(1,6), num2str(cumMCF_Shear_p(1,6)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
plot(myModes(1,15),cumMCF_Shear_p(1,15), 'ro', 'MarkerSize', 10)
text(myModes(1,15),cumMCF_Shear_p(1,15), num2str(cumMCF_Shear_p(1,15)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
xticks(myModes);
%--------------------------------------%
figure(12345)
cumMCF_Bend_p=cumMCF_Bend*100;
plot(myModes(1,1:6),cumMCF_Bend_p(1,1:6),myModes(1,6:15),cumMCF_Bend_p(1,6:15),'--', 'LineWidth', 2.5);
% % add percentage to the y-values
yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels 
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
% Set x-label and y-label
xlabel('Number of Modal', 'FontName', 'Times New Roman', 'FontSize', 18);
ylabel('Cumulative Contribution Factor (%)', 'FontName', 'Times New Roman', 'FontSize', 18);
title('Modal Contribution Factor of Overturning Moment at The Base of The Dam Slice', 'FontName', 'Times New Roman', 'FontSize', 18);
hold on
plot(myModes(1,6),cumMCF_Bend_p(1,6), 'ro', 'MarkerSize', 10)
text(myModes(1,6),cumMCF_Bend_p(1,6), num2str(cumMCF_Bend_p(1,6)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
hold on
plot(myModes(1,15),cumMCF_Bend_p(1,15), 'ro', 'MarkerSize', 10)
text(myModes(1,15),cumMCF_Bend_p(1,15), num2str(cumMCF_Bend_p(1,15)), 'VerticalAlignment', 'top','FontName', 'Times New Roman', 'FontSize', 18)
xticks(myModes);
