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

function F=streaming(F,siteNeighbor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STREAMING PROCEDURE
% for D2Q9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp=F(:,2);F(:,2)=F(siteNeighbor(:,6),6);F(siteNeighbor(:,6),6)=temp;
temp=F(:,3);F(:,3)=F(siteNeighbor(:,7),7);F(siteNeighbor(:,7),7)=temp;
temp=F(:,4);F(:,4)=F(siteNeighbor(:,8),8);F(siteNeighbor(:,8),8)=temp;
temp=F(:,5);F(:,5)=F(siteNeighbor(:,9),9);F(siteNeighbor(:,9),9)=temp;
