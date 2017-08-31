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
% This function returns the query pattern and mask (binary matrix which
% indicates where the position of the query pattern)
% Return     
%           fx,fy,fz represents the search pattern
%           fx    m x n x k matrix of the x coordinate of every vector    
%           fy    m x n x k matrix of the y coordinate of every vector  
%           fz    m x n x k matrix of the z coordinate of every vector  
%           mask  Binar matrix which indicates the position of the query pattern. 1 if on this
%                 coordinate is the query pattern, 0 if not
% Parameter
%           run               testrun of the paper between 1-6.
%           patimg            image index of the pattern (index of the file in inputDir)
%                             and the testrun (1-6)
%           lts               largest topology segmentation.
%                             0: manual segmentation. Just using the mask
%                             1: automatic segmentation within a bounding
%                             box. The bounding box will be created by the
%                             mask.
%                             2: No segmentation within a bounding box.
%           inputDir          Directory which contains the preprocessed
%                             vector fields.
%           redSize           Scaling factor for size reduction of the
%                             input volumes. Needed for the same size
%                             reduction for the mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fx,fy,fz,mask] = getTestSearchPattern(run,patImg,lts,inputDir,redSize)
    
    %if patImg < 0, the standard image with available segmentation is
    %choosen, if patImg is an image (given by its index), no segmentation
    %or an automatic segmentation will be applied
    if run==1 || run==3 || run==4 || run==5 || run==6
        img=983;
    end
    
    if run==2
        img=668;
    end
    
   
    if patImg > 0
        img=patImg;
    end
 
          if run==1
            data = loadMaskFromMesh('Masks/pIP10proj',true);
            pattern =  load(sprintf('%s/%d.mat',inputDir,img));
            
              mask=zeros(768,768,165);
                data=floor(data);

                for i=1:size(data,1)
                    mask((floor(data(i,1)/420.997*768)-2):(floor(data(i,1)/420.997*768)+2),(floor(data(i,2)/420.997*768)-2):(floor(data(i,2)/420.997*768)+2),(floor(data(i,3)/164*165)-2):(floor(data(i,3)/164*165)+2))=1;   
                end
          end
          
          if run==2
              data = loadMaskFromMesh('Masks/pIP10projR',true);
              pattern =  load(sprintf('%s/%d.mat',inputDir,img));
              
                mask=zeros(768,768,165);
                data=floor(data);

                for i=1:size(data,1)
                    mask((floor(data(i,1)/420.997*768)-2):(floor(data(i,1)/420.997*768)+2),(floor(data(i,2)/420.997*768)-2):(floor(data(i,2)/420.997*768)+2),(floor(data(i,3)/164*165)-2):(floor(data(i,3)/164*165)+2))=1;   
                end
          end
          
          if run==3
              mask = loadMaskFromMesh('Masks/pIP10_arbA_Mask',false);
              pattern =  load(sprintf('%s/%d.mat',inputDir,img));
              
          end
          
          if run==4
              mask = loadMaskFromMesh('Masks/pIP10_arbB_Mask',false);
              pattern =  load(sprintf('%s/%d.mat',inputDir,img));

          end
          
          if run==5
              mask = loadMaskFromMesh('Masks/pIP10_arbC_Mask',false);
              pattern =  load(sprintf('%s/%d.mat',inputDir,img));

          end
          
          if run==6
              data = loadMaskFromMesh('Masks/pIP10proj',true);
              pattern =  load(sprintf('%s/%d.mat',inputDir,img));
            
              mask=zeros(768,768,165);
                data=floor(data);

                for i=1:size(data,1)
                    mask((floor(data(i,1)/420.997*768)-2):(floor(data(i,1)/420.997*768)+2),(floor(data(i,2)/420.997*768)-2):(floor(data(i,2)/420.997*768)+2),(floor(data(i,3)/164*165)-2):(floor(data(i,3)/164*165)+2))=1;   
                end
              
              data1 = loadMaskFromMesh('Masks/pIP10_arbA_Mask',false);
               for i=1:size(data1,1)
                  for j=1:size(data1,2)
                      for k=1:size(data1,3)
                            if data1(i,j,k)==1
                                mask(i,j,k)=1;
                            end
                      end
                  end
               end
              clear data1;
              
              data2 = loadMaskFromMesh('Masks/pIP10_arbB_Mask',false);
               for i=1:size(data2,1)
                  for j=1:size(data2,2)
                      for k=1:size(data2,3)
                            if data2(i,j,k)==1
                                mask(i,j,k)=1;
                            end
                      end
                  end
               end
               clear data2;
              
              data3 = loadMaskFromMesh('Masks/pIP10_arbC_Mask',false);
              for i=1:size(data3,1)
                  for j=1:size(data3,2)
                      for k=1:size(data3,3)
                            if data3(i,j,k)==1
                                mask(i,j,k)=1;
                            end
                      end
                  end
              end
              clear data3;
              
          end
        
      

        mask = mask(1:redSize:end,1:redSize:end,1:redSize:end);
        mask = imdilate(mask,ones(3,3,3));
        
        
        
        fx=pattern.vectors{1};
        fy=pattern.vectors{2};
        fz=pattern.vectors{3};
     
        %if lts == 1, a bounding box with an automatic segmentation will be
        %used instead of the mask, if lts == 2, no segmentation will be
        %applied within the bounding box
        if lts==1
            top=pattern.vectors{4};

            patternWindow=logical(mask);
            bb=regionprops(patternWindow,'BoundingBox');
            
            while size(bb,1)>1
                 mask = imdilate(mask,ones(3,3,3));
                 patternWindow=logical(mask);
                 bb=regionprops(patternWindow,'BoundingBox');
            end
            wx=floor(bb.BoundingBox(2));
            wy=floor(bb.BoundingBox(1));
            wz=floor(bb.BoundingBox(3));
            ww=floor(bb.BoundingBox(5));
            wh=floor(bb.BoundingBox(4));
            wl=floor(bb.BoundingBox(6));
            
            top=top(wx:wx+ww-1,wy:wy+wh-1,wz:wz+wl-1);
          
            thresh = graythresh(top);
           
            topIndex=getTop(top,thresh);
            top=zeros(size(top));
            top(topIndex)=1;
            patternWindow((wx):((wx)+size(top,1)-1),(wy):(ceil(wy)+size(top,2)-1),(wz):((wz)+size(top,3)-1))=top;
            mask=patternWindow;
	
        else if lts==2
         
                patternWindow=logical(mask);
                bb=regionprops(patternWindow,'BoundingBox');
                
               while size(bb,1)>1
                 mask = imdilate(mask,ones(3,3,3));
                 patternWindow=logical(mask);
                 bb=regionprops(patternWindow,'BoundingBox');
                end
                
                wx=floor(bb.BoundingBox(2));
                wy=floor(bb.BoundingBox(1));
                wz=floor(bb.BoundingBox(3));
                ww=floor(bb.BoundingBox(5));
                wh=floor(bb.BoundingBox(4));
                wl=floor(bb.BoundingBox(6));

                mask=zeros(size(patternWindow));
                mask((wx):((wx)+ww-1),(wy):(ceil(wy)+wh-1),(wz):((wz)+wl-1))=1;
	
            end
        end
        
        mask = imdilate(mask,ones(3,3,3));

end
