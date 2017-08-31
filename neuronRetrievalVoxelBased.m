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
% This function performs a neuron retrieval on pre-processed drosophila
% brain images. Multiple retrievals (combination of pattern image index and testrun)
% are possible to reduce the times of loading all the data
% Return     
%           precisions        cell list of precision arrays for every retrieval
%           recalls           cell list of recall arrays for every retrieval
%           tps               cell list of true positive amounts for every
%                             retrieval
% Parameter
%           patImgs           N x 2 matrix for N image retrievals. Every
%                             row of the matrix is a combination of the
%                             image index of the pattern (index of the file in inputDir)
%                             and the testrun (1-6).
%           inputDir          Directory which contains the preprocessed
%                             vector fields.
%           readDir           Directory which contains the original images.
%                             Only used to get the image index
%           maxSRange         Is needed for the applyGVF function: 
%                             "Scale factor of the bouning box arround the query
%                             pattern defined by "mask". Is necessary because
%                             the Gradient Vector Flow (GVF) will expand/spread the
%                             query pattern, so the bounding box would be too
%                             small. Usually a value of 1.4 is ok, if gvf2 is <
%                             200"
%           redSize           Needed for the scaling of the masc inScaling factor for size reduction of the
%                             input volumes. 
%           gvf1              Parameter mü for GVF. A value between 0.03 and
%                             0.05 was suitable for testing. In principle the
%                             magnitude of field deformation in one iteration
%           gvf2              Amount of GVF iterations
%           lts               largest topology segmentation.
%                             0: manual segmentation. Just using the mask
%                             1: automatic segmentation within a bounding
%                             box. The bounding box will be created by the
%                             mask.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [precisions,recalls,tps] = neuronRetrievalVoxelBased(patImgs,inputDir,readDir,redSize)
    
    listing = dir(readDir);

    %this part is for loading multiple pattern images. This allows to
    %perform N retrievals while loading every preprocessed image only once
    %instead of N times
    patternBins=[];
   
    for i=1:size(patImgs,1)
        
            %only to get the binary mask
            sprintf('getTestSearchPattern run %d patImgIndex %d',patImgs(i,2),patImgs(i,1))
            [fx,fy,fz,mask] = getTestSearchPattern(patImgs(i,2),patImgs(i,1),2,inputDir,redSize);

            %increase search area (because no GVF)
            mask = imdilate(mask,ones(15,15,15));
            
            pattern = load_am_data(sprintf('%s/%s',readDir,listing(patImgs(i,1)+2).name)); 

            pattern = pattern(1:redSize:end,1:redSize:end,1:redSize:end);

            patternWindow=logical(mask);
            bb=regionprops(patternWindow,'BoundingBox');
            bb=bb(1);
            wx=floor(bb.BoundingBox(2))+1;
            
            wy=floor(bb.BoundingBox(1))+1;
            wz=floor(bb.BoundingBox(3))+1;
            ww=floor(bb.BoundingBox(5));
            wh=floor(bb.BoundingBox(4));
            wl=floor(bb.BoundingBox(6));
            pattern=pattern((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
            mask=mask((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
            patternBin = pattern>(graythresh(pattern)*255);
            patternBin(mask==0)=0;

            patternBin=patternBin(:);

            patternBins{i}.bin=patternBin;
            patternBins{i}.wx=wx;
            patternBins{i}.wy=wy;
            patternBins{i}.wz=wz;
            patternBins{i}.ww=ww;
            patternBins{i}.wh=wh;
            patternBins{i}.wl=wl;
            patternBins{i}.mask=mask;
            patternBins{i}.run=patImgs(i,2);
            patternBins{i}.patImg=patImgs(i,1);
  
    end
    

    results=[];

    %calculate the dice coefficient to the query patterns
    for index=3:size(listing,1)
        close('all')
        fclose('all');
        
        sprintf('load %d',(index-2))
        
        volume = load_am_data(sprintf('%s/%s',readDir,listing(index).name));
        volume = volume(1:redSize:end,1:redSize:end,1:redSize:end);
     
        %for multiple query patterns
        for i=1:size(patternBins,2)
             wx=patternBins{i}.wx;
            wy=patternBins{i}.wy;
            wz=patternBins{i}.wz;
            ww=patternBins{i}.ww;
            wh=patternBins{i}.wh;
            wl=patternBins{i}.wl;
            volume1=volume((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
            
            bin=volume1>(graythresh(volume1)*255);
            
            bin=bin(:);

        
            results{i}.res(index-2)=(sum(bin==patternBins{i}.bin)/size(bin,1)); 
            results{i}.run=patternBins{i}.run;
            results{i}.patImg=patternBins{i}.patImg;
       end
    end

   %calculate precision and recall by comparing the result with the ground
   %truth and true positives
   precisions=[];
   recalls=[];
   tps=[];
   
   
   
   for pb=1:size(patternBins,2) 
       try
                result=results{pb}.res';


                [a,ranking]=sort(result);
             

                annotatedResults =zeros(size(result,1),2);
             
                anno=getGroundTruth(results{pb}.run,readDir);
                for i=0:size(ranking,1)-1
                     annotatedResults(i+1,1)=result(ranking(size(ranking,1)-i));
                     annotatedResults(i+1,2)=anno(ranking(size(ranking,1)-i),2);
                end
    
              
                [tp,precision,recall] = getClassificationError(annotatedResults);

                tps{pb}.tp=tp;
                precisions{pb}.pre=precision;
                recalls{pb}.rec=recall;
                
                tps{pb}.run=results{pb}.run;
                precisions{pb}.run=results{pb}.run;
                recalls{pb}.run=results{pb}.run;
                
                tps{pb}.patImg=results{pb}.patImg;
                precisions{pb}.patImg=results{pb}.patImg;
                recalls{pb}.patImg=results{pb}.patImg;
       catch err
           sprints(err);
           tps{pb}.tp=0;    
           precisions{pb}.pre=0;
           recalls{pb}.rec=0;
       end
  end
        

end
