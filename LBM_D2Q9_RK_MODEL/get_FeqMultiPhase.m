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

function FEQ=get_FeqMultiPhase(rho,ux,uy,alpha,w,c)

cweight=[NaN 5 20 5 20 5 20 5 20];

FEQ=zeros(size(rho,1),9);

u_squ=ux.^2+uy.^2;

FEQ(:,1) =  rho.*(alpha-2/3*u_squ);
for k=2:9
    c_dot_u=c(1,k)*ux+c(2,k)*uy;
    FEQ(:,k) = rho.*( (1-alpha)./cweight(k) +w(k)*( 3*c_dot_u + 9/2*c_dot_u.^2 - 3/2*u_squ ));
end