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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the amout of true positives, precision and recall
%           annotatedResults  n x 2 matrix. The first column are the ordered
%                             indices of the images, the second column is 1
%                             if its positive, 0 if its a negative result
% Parameter
%           tp               amout of true positives
%           precision
%           recall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tp,precision,recall] = getClassificationError(annotatedResults)
   
   tp=sum(annotatedResults(1:sum(annotatedResults(:,2)),2));
   precision=[];
   recall=[];
   
   for border=1:size(annotatedResults,1)
        precision=[precision,(sum(annotatedResults(1:(border),2)))/border];
        recall=[recall,(sum(annotatedResults(1:(border),2)))/(sum(annotatedResults(:,2)))];
   end
end