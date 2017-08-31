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
% Preprocessing of image data to vector volumes
% Return --> nothing, but saves files with a counting number to saveDir in
% this format:
%           fx    m x n x k matrix of the x coordinate of every vector    
%           fy    m x n x k matrix of the y coordinate of every vector  
%           fz    m x n x k matrix of the z coordinate of every vector  
%           top   m x n x k matrix of topology values. Eigenvalue of the third eigenvector of
%                 the structure tensor. Toplogy value is between 0 and 1. A
%                 large topology value means, that it is more likely that
%                 on this position is a structure 
%
% Parameter
%           readDir           The directory which contains the amira files
%           saveDir           Output directory
%           value             m x n x k matrix. Just an ordinary volume matrix.
%           window size       size of the window which is used to compute
%                             the gradients for the tensors
%           tensorresolution  scaling factor for the tensor matrix (usually
%                             1)
%           gamma             Scaling factor for vector normalization.
%                             Performs well between 0.01 and 0.00001
%           redSize           Scaling factor for size reduction of the
%                             input volumes. 
%           smoothRange and   Additional Gaussian smoothing for the input  
%           smoothSigma       volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = generateVectorVolumes(readDir,saveDir,windowsize,tensorresolution,redSize,gamma,smoothRange,smoothSigma)
    matlabpool open 4;
    
    listing = dir(readDir);

    parfor index = 1:size(listing,1)-2
        
      if not(exist(sprintf('%s/%d.mat',saveDir,(index)), 'file'))
            try
                sprintf('Load: %s',sprintf('%s/%s',readDir,listing(index+2).name))
                volume = load_am_data(sprintf('%s/%s',readDir,listing(index+2).name));
            catch err
                sprintf('ERROR: Load: %s',sprintf('%s/%s',readDir,listing(index+2).name))
               volume=zeros(768,768,165);
            end
        
        volume = volume(1:redSize:end,1:redSize:end,1:redSize:end); 
        volume = smooth3(volume,'gaussian',smoothRange,smoothSigma);    
        [fx,fy,fz,top] = getGradients(volume,windowsize,tensorresolution,gamma);
        vectors=[];
        vectors{1}=fx;
        vectors{2}=fy;
        vectors{3}=fz;
        vectors{4}=top;
        
        sprintf('> save to %s/%d.mat',saveDir,(index))
        iSaveX(sprintf('%s/%d.mat',saveDir,(index)),vectors);
      end
    end
    matlabpool close;
end

function iSaveX( fname, vectors )
  save( fname, 'vectors' );
end
