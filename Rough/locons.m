%
% Function locons: defines boundary conditions (loads and constraints)
%
function [nCons,dC,nForce,dF,npq,dpq]=locons(dXY,thickness)

 % Load matrix dF: 
 % the i-th row of dF collects: the number of the loaded node; 
 %                                      the force direction; 
 %                                      the force intensity.
 % dF(i,1)=node number;
 % dF(i,2)=loaded direction ("1" along "x", "2" along "y");
 % dF(i,3)=action intensity.

 % dF=[1,1,10000];
 % Used for Surface Loads

%  % Soil Pressure Calculation (as Surface Loades)
%  GammaS=18; %kN/m3
%  phi_rad=33*pi/180;
%  
%  Ka=tan( pi/4 - phi_rad/2)^2;
% 
%   maxfx_left=Ka*GammaS*1.95;
%   maxfx_right=Ka*GammaS*5;
%   
%   left_x   =29:-4:1;
%   right_x  =205:4:273;
%   right_y =205:4:237;
%   
%   % Active Thrust Load Distrib x Left
%   LoadDist_left=zeros(1,length(left_x));
%   LoadDist_left(1)=0; LoadDist_left(end)=maxfx_left;
% 
%   for j=2:length(left_x)-1
%       LoadDist_left(j)=maxfx_left*(1 - dXY(left_x(j),2)/1.95);
%   end
%   
%   % Active Thrust Load Distrib x Right
%   LoadDist_right_x=zeros(1,length(right_x));
%   LoadDist_right_x(1)=0; LoadDist_right_x(end)=maxfx_right;
% 
%   for j=2:length(right_x)-1
%       LoadDist_right_x(j)=maxfx_right*(1 - dXY(right_x(j),2)/5);
%   end
%   
%   % Soil Weight Load Distrib y Right
%   LoadDist_right_y=zeros(1,length(right_y));
%   LoadDist_right_y(1)=0;  
%   
%   for j=2:length(right_y)
%       LoadDist_right_y(j)=GammaS*(5-dXY(right_y(j),2));
%   end
%   
%    dF=[0 0 0];
%   
%   % Active Thrust Load Contributaion in x lEFT
%   for i=2:length(left_x)
%       dF=[dF;
%           left_x(i-1) , 1 , (thickness/6)*abs(dXY(left_x(i),2)-dXY(left_x(i-1),2))*(2*LoadDist_left(i-1)+LoadDist_left(i));
%           left_x(i)   , 1 , (thickness/6)*abs(dXY(left_x(i),2)-dXY(left_x(i-1),2))*(2*LoadDist_left(i)+LoadDist_left(i-1));
%         ];
%   end
%   
%     % Active Thrust Load Contributaion in x Right
%   for i=2:length(right_x)
%       dF=[dF;
%           right_x(i-1) , 1 , -(thickness/6)*abs(dXY(right_x(i),2)-dXY(right_x(i-1),2))*(2*LoadDist_right_x(i-1)+LoadDist_right_x(i));
%           right_x(i)   , 1 , -(thickness/6)*abs(dXY(right_x(i),2)-dXY(right_x(i-1),2))*(2*LoadDist_right_x(i)+LoadDist_right_x(i-1));
%         ];
%   end
%   
%   % Soil Weight Load Contributaion in y Right
%   counter=1;
%   for i=2:length(right_y)
%       dF=[dF;
%           right_y(i-1) , 2 , -(thickness/6)*abs(dXY(right_y(i),1)-dXY(right_y(i-1),1))*(2*LoadDist_right_y(i-1)+LoadDist_right_y(i));
%           right_y(i)   , 2 , -(thickness/6)*abs(dXY(right_y(i),1)-dXY(right_y(i-1),1))*(2*LoadDist_right_y(i)+LoadDist_right_y(i-1));
%         ];
%     counter=counter+1;
%   end
%   dF(1,:)=[];
%   
%   [nForce,~]=size(dF);  % nForce=total number of considered loads

nForce=0;
dF=0;

 % Constraint matrix dC: 
 % the i-th row in dC collects: the number of the constrained node; 
 %                                      the direction of the constrained dof; 
 %                                      the magnitude of the imposed displacement.
 % dC(i,1)=node number;
 % dC(i,2)=constrained dof ("1" along "x", "2" along "y");
 % dC(i,3)=magnitude of the imposed displacement.
  dC=[ 1,1,0;
       1,2,0;
       2,1,0;
       2,2,0;
       3,1,0;
       3,2,0;
       4,1,0;
       4,2,0;
       5,1,0;
       5,2,0;
       6,1,0;
       6,2,0;
       
       ];


  [nCons,~]=size(dC);  % nCons=total number of constrained dofs
  
 % Distributed load matrix dpq: 
 % the i-th row of dpq collects: the number of the loaded element; 
 %                               the uniform load in x direction; 
 %                               the uniform load in y direction.
 % Used for Volume Forces (Self-Weight) 
  nElem_sw=55; % number of elements having self-weight
  dpq=[(1:nElem_sw)' , zeros(nElem_sw,1) , -21*ones(nElem_sw,1)];
% dpq=[ 1,0,-.001;
%       2,0,-.001]
  
  [npq,~]=size(dpq);  % npq=total number of considered loads

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
