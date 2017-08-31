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
function [precisions,recalls,tps] = neuronRetrieval3D(patImgs,inputDir,readDir,maxSRange,redSize,gvf1,gvf2,lts)
    warning ('off','all');
    listing = dir(inputDir);

    %this part is for loading multiple pattern images. This allows to
    %perform N retrievals while loading every preprocessed image only once
    %instead of N times
    patternBins=[];
    times=[];
    for i=1:size(patImgs,1)
            tic;
            sprintf('getTestSearchPattern run %d patImgIndex %d and applyGVF',patImgs(i,2),patImgs(i,1))
            [fx,fy,fz,mask] = getTestSearchPattern(patImgs(i,2),patImgs(i,1),lts,inputDir,redSize);

            [pfx,pfy,pfz,wx,wy,wz] = applyGVF(fx,fy,fz,mask,maxSRange,gvf1,gvf2);

            ww=size(pfx,1);
            wh=size(pfx,2);
            wl=size(pfx,3);

            patternBins{i}.pfx=pfx;
            patternBins{i}.pfy=pfy;
            patternBins{i}.pfz=pfz;
            patternBins{i}.wx=wx;
            patternBins{i}.wy=wy;
            patternBins{i}.wz=wz;
            patternBins{i}.ww=ww;
            patternBins{i}.wh=wh;
            patternBins{i}.wl=wl;
            patternBins{i}.mask=mask;
            patternBins{i}.run=patImgs(i,2);
            patternBins{i}.patImg=patImgs(i,1);
            times=[times,toc];
            sprintf('getTestSearchPattern run %d patImgIndex %d and applyGVF time: %f',patImgs(i,2),patImgs(i,1),toc)

    end
    
    sprintf('Average pattern generation time: %f',mean(times))
    timesCompare=[];
    timesSOSD=[];
    results=[];

   %calculate the sum of squared differences to the query patterns
   for index=3:size(listing,1)
        close('all')
        fclose('all');
        sprintf('load %d',(index-2))
        loadCompareTic=tic;
        
        error=0;
        try
            pattern = load(sprintf('%s/%d.mat',inputDir,(index-2)));
            error=0;
        catch err
             error=1;
             for i=1:size(patternBins,2)
                results{i}.res(index-2)=Inf;
                results{i}.run=patternBins{i}.run;
                results{i}.patImg=patternBins{i}.patImg;
             end
        end
        
        
        if error==0
            %for multiple query patterns
            for i=1:size(patternBins,2)
                tic;
                try
                    pfx=patternBins{i}.pfx;
                    pfy=patternBins{i}.pfy;
                    pfz=patternBins{i}.pfz;

                    wx=patternBins{i}.wx;
                    wy=patternBins{i}.wy;
                    wz=patternBins{i}.wz;
                    ww=patternBins{i}.ww;
                    wh=patternBins{i}.wh;
                    wl=patternBins{i}.wl;


                    fx=pattern.vectors{1};
                    fy=pattern.vectors{2};
                    fz=pattern.vectors{3};



                    fx=fx(wx:wx+ww-1,wy:wy+wh-1,wz:wz+wl-1);
                    fy=fy(wx:wx+ww-1,wy:wy+wh-1,wz:wz+wl-1);
                    fz=fz(wx:wx+ww-1,wy:wy+wh-1,wz:wz+wl-1);




                    %set fx,fy,fz to zero where the pattern vector field is 0.
                    %This allows a not rectangular query area
                    fx(abs(pfx)+abs(pfy)+abs(pfz)==0)=0;
                    fy(abs(pfx)+abs(pfy)+abs(pfz)==0)=0;
                    fz(abs(pfx)+abs(pfy)+abs(pfz)==0)=0;



                   if (sum(fx(:))+sum(fy(:))+sum(fz(:)))==0 %check if the are is empty
                       results{i}.res(index-2)=Inf;
                   else
                       results{i}.res(index-2)=summedSquaredDifference(pfx,pfy,pfz,fx,fy,fz);
                   end
                catch err
                    results{i}.res(index-2)=Inf;
                end
                timesSOSD=[timesSOSD,toc];
                results{i}.run=patternBins{i}.run;
                results{i}.patImg=patternBins{i}.patImg;
            end
        end
       timesCompare=[timesCompare,toc(loadCompareTic)];
    end

     sprintf('Average compare time: %f',mean(timesCompare))
     sprintf('Average SSOD time: %f',mean(timesSOSD)/(size(patternBins,2)))

    
   %calculate precision and recall by comparing the result with the ground
   %truth and true positives
   precisions=[];
   recalls=[];
   tps=[];
   
   
   
   for pb=1:size(patternBins,2) 
       try
                result=results{pb}.res';


                [a,ranking]=sort(-result);
             

                annotatedResults =zeros(size(result,1),2);
             
                anno=getGroundTruth(results{pb}.run,readDir);
                for i=0:size(ranking,1)-1
                     annotatedResults(i+1,1)=result(ranking(size(ranking,1)-i));
                     annotatedResults(i+1,2)=anno(ranking(size(ranking,1)-i),2);
                end
    
              
                [tp,precision,recall] = getClassificationError(annotatedResults);

                tps{pb}.tp=tp;
                tps{pb}.results=annotatedResults;
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
           tps{pb}.results=[];
           precisions{pb}.pre=0;
           recalls{pb}.rec=0;
       end
   end
        
     warning ('on','all');
end
