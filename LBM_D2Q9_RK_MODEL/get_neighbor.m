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

function [neighbor,i,j]=get_neighbor(address,Nx,Ny)

% address can be a row vector, in this case, neighbor consist of
% [numel(address),9] size matrix.

[i,j]=ind2sub([Nx Ny],address); % transform the linear index into a two dimensional index

iCond= (i-1==0); % conditional operator to apply periodic boundary condition
jCond= (j-1==0); % conditional operator to apply periodic boundary condition

iE = mod(i,Nx)+1; % periodic boundary condition in the two dimensional plane index
iW = iCond*Nx+(1-iCond).*(i-1); % periodic boundary condition in the two dimensional plane index
jN = mod(j,Ny)+1; % periodic boundary condition in the two dimensional plane index
jS = jCond*Ny+(1-jCond).*(j-1); % periodic boundary condition in the two dimensional plane index

neighbor=zeros(numel(address),9,'uint32');
neighbor(:,1) = uint32(address);
neighbor(:,2) = uint32(sub2ind([Nx Ny],i,jN));      % North neighbor in linear index with periodic boundary condition
neighbor(:,3) = uint32(sub2ind([Nx Ny],iE,jN));     % North East neighbor in linear index with periodic boundary condition
neighbor(:,4) = uint32(sub2ind([Nx Ny],iE,j));      % East neighbor in linear index with periodic boundary condition
neighbor(:,5) = uint32(sub2ind([Nx Ny],iE,jS));     % South East neighbor in linear index with periodic boundary condition
neighbor(:,6) = uint32(sub2ind([Nx Ny],i,jS));      % South neighbor in linear index with periodic boundary condition
neighbor(:,7) = uint32(sub2ind([Nx Ny],iW,jS));     % South West neighbor in linear index with periodic boundary condition
neighbor(:,8) = uint32(sub2ind([Nx Ny],iW,j));      % West neighbor in linear index with periodic boundary condition
neighbor(:,9) = uint32(sub2ind([Nx Ny],iW,jN));     % North West neighbor in linear index with periodic boundary condition

