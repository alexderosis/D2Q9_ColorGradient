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

%%
% This m-code reproduces some of the results of Table 3 in section 3.2 
% "Steady Bubble" of the following reference:
% Leclaire, S., Reggio, M. and Trépanier, J.-Y. (2012). 
% Numerical evaluation of two recoloring operators for an immiscible 
% two-phase flow lattice Boltzmann model, 
% Applied Mathematical Modelling 36(5): 2237-2252.
%
% The results are stored in the file "results.txt"
%
% Moreover, instead of using the second order anisotropic color gradient,
% it is possible to use the isotropic color 
% gradient of the following reference:
% Leclaire, S., Reggio, M. and Trépanier, J.-Y. (2011).
% Isotropic color gradient for simulating very high-density ratios with a
% two-phase flow lattice Boltzmann model,
% Computers and Fluids 48(1): 98-112.
%
% In average, results are better with the isotropic discretisation.
%
% If you use this code in your work, please cite:
% 1)	the webpage from which you download this m-code and
% 2)	the above two references accordingly. 
%%

clear all
close all

defaultStream = RandStream.getDefaultStream;
reset(defaultStream);
addpath('plotclr');

format long e

temp_R=[20 40];
temp_A=[0.0004 0.01];
temp_rhoR=[1,6];
temp_gamma=[1,6];

ID=0;
for e=1:length(temp_gamma)
for c=1:length(temp_A)
for d=1:length(temp_rhoR)
for h=1:length(temp_R)    
    ID=ID+1;
    A(ID)=temp_A(c);
    rhoR(ID)=temp_rhoR(d);
    gamma(ID)=temp_gamma(e);
    R(ID)=temp_R(h);
end
end
end
end

for i=1:ID
    %% multiphase
       
    options{i}.Nx=127+1;
    options{i}.Ny=127+1;
    options{i}.tmax=200000;
    options{i}.maxDiffTol=10^-10;
    options{i}.ncheck=500;
    options{i}.beta=0.99;
    options{i}.gamma=gamma(i);
    options{i}.aR=A(i);
    options{i}.rhoR=rhoR(i);    
    options{i}.R=R(i);
    
	%% The type of gradient discretisation is changed here.
    options{i}.gradient='anisotropic';%'isotropic';
    
    if(strcmp(options{i}.gradient,'isotropic'))
        options{i}.aR=options{i}.aR*6;
    end
    options{i}.aB=options{i}.aR;
    
    options{i}.alphaB = 0.6;
    options{i}.alphaR = 1-(1-options{i}.alphaB)/options{i}.gamma;

    options{i}.nuR=1/6;
    options{i}.nuB=1/6;
    options{i}.omegaR=1/(3*options{i}.nuR+0.5);
    options{i}.omegaB=1/(3*options{i}.nuB+0.5);
    options{i}.ID=i;
end

% matlabpool 2
for i=1:ID

    LBMmultiPhase(options{i});

end
% matlabpool close