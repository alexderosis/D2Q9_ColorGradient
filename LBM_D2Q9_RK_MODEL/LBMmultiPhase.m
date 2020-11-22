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

function options = LBMmultiPhase(options)

%% lattice information (periodic boundary condition)

siteNeighbor=get_neighbor(1:options.Nx*options.Ny,options.Nx,options.Ny);

%% WE DEFINE SOME D2Q9 CONSTANT
w=[4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36];
c(1,:)=[0,0,1,1, 1, 0,-1,-1,-1];
c(2,:)=[0,1,1,0,-1,-1,-1, 0, 1];
upwindOrientation=[1 6:9 2:5];

%% INITIALISATION OF THE DENSITY DISTRIBUTION FUNCTION

R=get_FeqMultiPhase(options.rhoR*ones(options.Nx*options.Ny,1),zeros(options.Nx*options.Ny,1),zeros(options.Nx*options.Ny,1),options.alphaR,w,c);
B=get_FeqMultiPhase(options.rhoR/options.gamma*ones(options.Nx*options.Ny,1),zeros(options.Nx*options.Ny,1),zeros(options.Nx*options.Ny,1),options.alphaB,w,c);

[x,y] = ind2sub([options.Nx options.Ny],1:options.Nx*options.Ny);
test=(x-64.5).^2+(y-64.5).^2<=options.R^2;
R = bsxfun(@times,R,test');
B = bsxfun(@times,B,~test');

%plotclr(x,y,(sum(R,2)-sum(B,2))./(sum(R,2)+sum(B,2)),'square');

%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%

disp('start main loop...');tic;

maxDiff=10^11;
maxDiffOld=10^10;
t=0;
tOld=0;
while(t<options.tmax && maxDiff>options.maxDiffTol)  
    
    %% REMEMBER R and B IN ORDER TO CHECK STEADY STATE CONVERGENCE
    Rold=R;
    Bold=B;
    
    %% COLLISION
   
    [R B]=collisionMultiPhase(R,B,siteNeighbor,w,c,upwindOrientation,options);
    
    %% STREAMING
    R=streaming(R,siteNeighbor);
    B=streaming(B,siteNeighbor);

    %% INCREMENT THE TIME STEP
    t=t+1;
        
    %% CHECK STEADY STATE CONVERGENCE AND SIMULATION VERBOSITY
    if mod(t,options.ncheck)==0    
    
        % Steady state convergence check
        fSum        = sum(sum(R+B,1)); % faster to operate on the column first
        maxDiff     = max(max(max(abs(R(:)-Rold(:)),[],1),[],2),max(max(abs(B(:)-Bold(:)),[],1),[],2)); % faster to operate on the column first
        meanDiff    = max(mean(mean(abs(R(:)-Rold(:)),1),2),mean(mean(abs(B(:)-Bold(:)),1),2)); % faster to operate on the column first
                   
        % first and second moment
        rhoR = R(:,1)+R(:,2)+R(:,3)+R(:,4)+R(:,5)+R(:,6)+R(:,7)+R(:,8)+R(:,9);
        rhoB = B(:,1)+B(:,2)+B(:,3)+B(:,4)+B(:,5)+B(:,6)+B(:,7)+B(:,8)+B(:,9);
        rho = rhoR+rhoB;
        ux  = (R(:,3)+R(:,4)+R(:,5)-R(:,7)-R(:,8)-R(:,9) + B(:,3)+B(:,4)+B(:,5)-B(:,7)-B(:,8)-B(:,9))./rho;
        uy  = (R(:,2)+R(:,3)+R(:,9)-R(:,5)-R(:,6)-R(:,7) + B(:,2)+B(:,3)+B(:,9)-B(:,5)-B(:,6)-B(:,7))./rho;
        
        phi=0.9999999999;
        colorFieldNorm=(rhoR-rhoB)./(rhoR+rhoB);
        test=1;
        while (1)       
            inIndex  = colorFieldNorm>=phi;
            outIndex = colorFieldNorm<=-phi;
            if nnz(inIndex)<=0 || nnz(outIndex)<=0
                phi=phi-0.0000000009*test;
                test=10*test;
            else
               break; 
            end
        end

        rhoRin = mean(sum(R(inIndex,:),2));
        Pin=3/5*rhoRin*(1-options.alphaR);

        rhoBin = mean(sum(B(outIndex,:),2));
        Pout=3/5*rhoBin*(1-options.alphaB);  
        
        nuMean=0.5*(options.nuR+options.nuB);
        omegaMean=1/(3*nuMean+0.5);
        if(strcmp(options.gradient,'isotropic'))
            sigmaPredict=2/9*((1+1/options.gamma)*options.rhoR/2)/omegaMean*(options.aR+options.aB);
        else
            sigmaPredict=4/3*((1+1/options.gamma)*options.rhoR/2)/omegaMean*(options.aR+options.aB);
        end
        sigmaLaplace=(Pin-Pout)*options.R;
        
        maxVel=max(sqrt(ux.^2+uy.^2));
        
        % Simulation verbosity
        if mod(t/options.ncheck,30)==0 || t==options.ncheck
            disp(['|   time step |          toc |         fSum |      maxDiff |      meanDiff |        Gamma |            A |         rhoR |            R | sigmaPredict |       maxVel |    Error (%) | '])
        end
        disp(sprintf('| %11.u | %12.5e | %12.5e | %12.5e |  %12.5e | %12.5e | %12.5e | %12.5e | %12.5e | %12.5e | %12.5e | %12.5e |',t,toc,fSum,maxDiff,meanDiff,options.gamma,options.aR,options.rhoR,options.R,sigmaPredict,maxVel,abs(sigmaLaplace-sigmaPredict)/sigmaPredict*100))
        
%         plotclr(x,y,(sum(R,2)-sum(B,2))./(sum(R,2)+sum(B,2)),'square');

    end
    
end

disp('stop main loop...')

fid=fopen(['results','.txt'],'a');
fprintf(fid,'%11.u %11.u %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\r\n',options.ID,t,toc,fSum,maxDiff,meanDiff,options.gamma,options.aR,options.rhoR,options.R,sigmaPredict,maxVel,abs(sigmaLaplace-sigmaPredict)/sigmaPredict*100);
fclose(fid);