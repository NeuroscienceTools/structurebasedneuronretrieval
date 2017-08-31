% Copyright (C) 2014 VRVis.
% All rights reserved.
% Contact: VRVis Forschungs-GmbH (office@vrvis.at)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the VRVis Forschungs-GmbH.
% 4. Neither the name of the VRVis Forschungs-GmbH nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%

function [topIndex] = getTop(top,thresh)

    topBW=top;
       topBW(topBW<thresh)=false;
       topBW(topBW>=thresh)=true;
       topBW=logical(topBW);
       topComp=bwconncomp(topBW);
       
        
        maxpix=0;
        for i=1:size(topComp.PixelIdxList,2)
            if size(topComp.PixelIdxList{i},1)>maxpix
                maxpix=size(topComp.PixelIdxList{i},1);
                topIndex=topComp.PixelIdxList{i};
            end
        end
        
     
end