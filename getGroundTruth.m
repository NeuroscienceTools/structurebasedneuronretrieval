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
% This function returns the ground truth for the 6 testruns
% Return         
%           anno              n x 2 matrix. First column is the index of
%                             the image (fly brain). Index is the same
%                             which is used for pre processing. Second
%                             column is 0 if the image does NOT contains the pattern
%                             image of the testrun, 1 if it contains the
%                             pattern image
% Parameter
%           run               testrun of the paper between 1-6.
%           readDir           Directory which contains the original images.
%                             Only used to get the image index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [anno] = getGroundTruth(run,readDir)
    [anno1,anno2]=textread(sprintf('GroundTruth/annoRun%d_Laszlo_names.txt',run),'%s %d');
    listing = dir(readDir);
    
    anno=zeros(size(listing,1)-2,2);
    for i=3:size(listing,1)
     for j=size(anno1,1):-1:1
         if  strcmp(anno1(j),listing(i).name)
            anno(i-2,1)=i-2;
            anno(i-2,2)=anno2(j);
         end
     end
    end
end
