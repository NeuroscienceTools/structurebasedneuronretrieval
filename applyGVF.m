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
%           fx            m x n x k matrix of the x coordinate of every vector    
%           fy            m x n x k matrix of the y coordinate of every vector  
%           fz            m x n x k matrix of the z coordinate of every vector  
%           mask          Binar matrix which indicates the position of the query pattern. 1 if on this
%                         coordinate is the query pattern, 0 if not
%           maxSRange     Scale factor of the bouning box arround the query
%                         pattern defined by "mask". Is necessary because
%                         the Gradient Vector Flow (GVF) will expand/spread the
%                         query pattern, so the bounding box would be too
%                         small. Usually a value of 1.4 is ok, if gvf2 is <
%                         200
%           gvf1          Parameter mü for GVF. A value between 0.03 and
%                         0.05 was suitable for testing. In principle the
%                         magnitude of field deformation in one iteration
%           gvf2          Amount of GVF iterations
% Parameter
%           pfx,pfy,pfz   query pattern after applying Gradient Vector Flow
%                         (GVF). The query pattern is not a m x n x k
%                         matrix, just the query pattern in the bounding
%                         box. 
%           x,y,z         Represents the offset of the bounding box. x,y,z
%                         is the left upper corner of the bounding box
%                         (where pfx,pfy,pfz starts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pfx,pfy,pfz,x,y,z] = applyGVF(fx,fy,fz,mask,maxSRange,gvf1,gvf2)
      
        fx(not(mask))=0;
        fy(not(mask))=0;
        fz(not(mask))=0;

        %Use only area in a bounding box around the search pattern (given
        %by mask), so GVF is computational cheaper to apply
        patternWindow=logical(mask);
        bb=regionprops(patternWindow,'BoundingBox');
        bb=bb(1);
        wx=floor(bb.BoundingBox(2));
        wy=floor(bb.BoundingBox(1));
        wz=floor(bb.BoundingBox(3));
        ww=floor(bb.BoundingBox(5));
        wh=floor(bb.BoundingBox(4));
        wl=floor(bb.BoundingBox(6));

        %increase the area around the bounding box by maxSRange (resize
        %scale), to give GVF more space (because GVF will increase the 
        %search area.
        wx = floor(max(1,wx+(ww/2)-(ww * maxSRange)/2));
        wy = floor(max(1,wy+(wh/2)-(wh * maxSRange)/2));
        wz = floor(max(1,wz+(wl/2)-(wl * maxSRange)/2));

        ww = min(size(patternWindow,1)-wx-1,(ww * maxSRange));
        wh = min(size(patternWindow,2)-wy-1,(wh * maxSRange));
        wl = min(size(patternWindow,3)-wz-1,(wl * maxSRange));


        fx=fx((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
        fy=fy((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
        fz=fz((wx):(wx+ww-1),(wy):(wy+wh-1),(wz):(wz+wl-1));
        

        [pfx,pfy,pfz] = GVF(fx,fy,fz, gvf1,gvf2);

	
	    %Search area will be defined by the expanded field.
        pfx(abs(pfx)+abs(pfy)+abs(pfz)<0.0001)=0;
        pfy(abs(pfx)+abs(pfy)+abs(pfz)<0.0001)=0; 
        pfz(abs(pfx)+abs(pfy)+abs(pfz)<0.0001)=0; 

        x=wx;
        y=wy;
        z=wz;

end
