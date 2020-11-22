%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2010, Stephanie Contardo
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%    * Redistributions of source code must retain the above copyright 
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function h = plotclr(x,y,v,marker,vlim)
% plots the values of v colour coded
% at the positions specified by x and y.
% A colourbar is added on the right side of the figure.
%
% The colourbar strectches from the minimum value of v to its
% maximum.
%
% 'marker' is optional to define the marker being used. The
% default is a point. To use a different marker (such as circles, ...) send
% its symbol to the function (which must be enclosed in '; see example).
%
% 'vlim' is optional, to define the limits of the colourbar.
% v values outside vlim are not plotted
%
% modified by Stephanie Contardo, CSIRO, 2009
% from 'plotc' by Uli Theune, University of Alberta, 2004
%

function plotclr(x,y,v,marker,vlim)

figure(1);
clf
hold on        

map=jet(64);
if nargin >7
    miv = vlim(1) ;
    mav = vlim(2) ;
else
    miv=min(v);
    mav=max(v);
end
clrstep = (mav-miv)/size(map,1) ;

% Plot the points
for nc=1:size(map,1)
    iv = find(v>=miv+(nc-1)*clrstep & v<=miv+nc*clrstep+1000*eps) ;
    if ~isempty(iv)
        plot3(x(iv),y(iv),v(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:))
    end
end

% caxis([miv mav])
caxis([-1 1])
axis([min(x) max(x) min(y) max(y)])
colorbar;
drawnow
pause(0.1)
hold off

% % Re-format the colorbar
% h=colorbar;
% 
% %set(h,'ylim',[1 length(map)]);
% yal=linspace(1,length(map),10);
% set(h,'ytick',yal);
% % Create the yticklabels
% ytl=linspace(miv,mav,10);
% s=char(10,4);
% for i=1:10
%     if min(abs(ytl)) >= 0.001
%         B=sprintf('%-4.3f',ytl(i));
%     else
%         B=sprintf('%-3.1E',ytl(i));
%     end
%     s(i,1:length(B))=B;
% end
% set(h,'yticklabel',s);
% grid on
% view(2)

