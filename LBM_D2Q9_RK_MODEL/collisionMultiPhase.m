%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2012, Sébastien Leclaire                                  %
% All rights reserved.                                                    %
%                                                                         %
% Redistribution and use in source and binary forms, with or without      %
% modification, are permitted provided that the following conditions      %
% are met:                                                                %
%                                                                         %
% 1) Redistributions of source code must retain the above copyright       %
%    notice, this list of conditions and the following disclaimer.        %
% 2) Redistributions in binary form must reproduce the above copyright    %
%    notice, this list of conditions and the following disclaimer in the  %
%    documentation and/or other materials provided with the distribution. %
%                                                                         %
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS %
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED   %
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A         %
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      %
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT        %
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   %
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY   %
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT     %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE   %
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R B]=collisionMultiPhase(R,B,siteNeighbor,w,c,upwindOrientation,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COLLISION PROCEDURE FOR MULTIPHASE

% first and second moment
rhoR = R(:,1)+R(:,2)+R(:,3)+R(:,4)+R(:,5)+R(:,6)+R(:,7)+R(:,8)+R(:,9);
rhoB = B(:,1)+B(:,2)+B(:,3)+B(:,4)+B(:,5)+B(:,6)+B(:,7)+B(:,8)+B(:,9);
rho = rhoR+rhoB;
ux  = (R(:,3)+R(:,4)+R(:,5)-R(:,7)-R(:,8)-R(:,9) + B(:,3)+B(:,4)+B(:,5)-B(:,7)-B(:,8)-B(:,9))./rho;
uy  = (R(:,2)+R(:,3)+R(:,9)-R(:,5)-R(:,6)-R(:,7) + B(:,2)+B(:,3)+B(:,9)-B(:,5)-B(:,6)-B(:,7))./rho;

% equilibrium

FEQ_R=get_FeqMultiPhase(rhoR,ux,uy,options.alphaR,w,c);
FEQ_B=get_FeqMultiPhase(rhoB,ux,uy,options.alphaB,w,c);

%% LBGK collision

R = R - options.omegaR*(R-FEQ_R);
B = B - options.omegaB*(B-FEQ_B);

%% evaluate the color gradient

if(strcmp(options.gradient,'isotropic'))
    gradientWeight=[0 1/3 1/12 1/3 1/12 1/3 1/12 1/3 1/12];
else
    gradientWeight=[0 1 1 1 1 1 1 1 1];
end

colorGrad_x=0;
colorGrad_y=0;
rhoDiff=rhoR-rhoB;
for i=2:9
    colorGrad_x = colorGrad_x +gradientWeight(i)*c(1,i)*rhoDiff(siteNeighbor(:,i));
    colorGrad_y = colorGrad_y +gradientWeight(i)*c(2,i)*rhoDiff(siteNeighbor(:,i));
end

%% Two phase perturbation

b=[-4/27 2/27 5/108 2/27 5/108 2/27 5/108 2/27 5/108];

colorGradNorm                                       = sqrt(colorGrad_x.^2+colorGrad_y.^2);
colorGradNorm_squ                                   = colorGradNorm.^2;
 
for i=1:9
    c_dot_colorGrad=c(1,i)*colorGrad_x+c(2,i)*colorGrad_y;
    
    temp = colorGradNorm .* ( w(i).*c_dot_colorGrad.^2./colorGradNorm_squ - b(i) );
    temp(~isfinite(temp))=0; % if there were a division by 0 in the previous line it mean that temp should be 0.
    R(:,i) = R(:,i) + options.aR/2*temp;
    B(:,i) = B(:,i) + options.aB/2*temp;
end

%% Recoloring 

rho = rhoR+rhoB;

cweight=[1 5 20 5 20 5 20 5 20];
cweightR=(1-options.alphaR)./cweight;
cweightB=(1-options.alphaB)./cweight;

FEQ_N=zeros(size(rhoR,1),9);
FEQ_N(:,1) = rhoR.*options.alphaR+rhoB.*options.alphaB;
for k=2:9
   FEQ_N(:,k) = rhoR.*cweightR(k) +  rhoB.*cweightB(k);
end

for i=1:9
    c_dot_colorGrad=c(1,i)*colorGrad_x+c(2,i)*colorGrad_y;
    
    temp=sqrt(c(1,i)^2+c(2,i)^2).*colorGradNorm;
    cosPHI= c_dot_colorGrad./temp;
    cosPHI(~isfinite(cosPHI))=0; % if there were a division by 0 in the previous line it mean that cosPHI should be 0.
        
    temp=options.beta*rhoR.*rhoB.*(cosPHI)./rho.^2;
    RBsum=R(:,i)+B(:,i); % needed to avoid the update of B using R of the previous line
    R(:,i) = RBsum.*rhoR./rho+temp.*FEQ_N(:,i);
    B(:,i) = RBsum.*rhoB./rho-temp.*FEQ_N(:,i);
end

%% swapping R and B in order to do a correct streaming
R=R(:,upwindOrientation);
B=B(:,upwindOrientation);