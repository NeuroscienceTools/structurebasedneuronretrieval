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
% This function returns a vector field (of unit vectors) from a volume (R ^
% (m x n x k)
% Return         
%           fx    m x n x k matrix of the x coordinate of every vector    
%           fy    m x n x k matrix of the y coordinate of every vector  
%           fz    m x n x k matrix of the z coordinate of every vector  
%           top   m x n x k matrix of topology values. Eigenvalue of the third eigenvector of
%                 the structure tensor. Toplogy value is between 0 and 1. A
%                 large topology value means, that it is more likely that
%                 on this position is a structure 
%
% Parameter
%           volume            m x n x k matrix. Just an ordinary volume matrix.
%           window size       size of the window which is used to compute
%                             the gradients for the tensors
%           tensorresolution  scaling factor for the tensor matrix (usually
%                             1)
%           gamma             Scaling factor for vector normalization.
%                             Performs well between 0.01 and 0.00001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx,fy,fz,top] = getGradients(volume,windowsize,tensorresolution,gamma)
      
    epsilon=0.5; %for vector normalization
   
  
    fx=zeros(size(volume,1)/tensorresolution,size(volume,2)/tensorresolution,size(volume,3)/tensorresolution);
    fy=zeros(size(volume,1)/tensorresolution,size(volume,2)/tensorresolution,size(volume,3)/tensorresolution);
    fz=zeros(size(volume,1)/tensorresolution,size(volume,2)/tensorresolution,size(volume,3)/tensorresolution);
    
    top=zeros(size(volume,1)/tensorresolution,size(volume,2)/tensorresolution,size(volume,3)/tensorresolution);
    
    
    
    for i=floor(windowsize/2)+1:size(volume,1)-floor(windowsize/2)
        for j=floor(windowsize/2)+1:size(volume,2)-floor(windowsize/2)
            for k=floor(windowsize/2)+1:size(volume,3)-floor(windowsize/2)
                if mod(i,tensorresolution)==0 && mod(j,tensorresolution)==0 && mod(k,tensorresolution)==0
                    window = volume((i-(floor(windowsize/2))):(i+(floor(windowsize/2))),(j-(floor(windowsize/2))):(j+(floor(windowsize/2))),(k-(floor(windowsize/2))):(k+(floor(windowsize/2))));
                    Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=part_deriv(window,1);
                   
                    Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=part_deriv(window,2);
                   
                    Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=part_deriv(window,3); 
                   
                    Ixyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)= Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Iyzvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)= Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Ixzvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)= Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    
                    Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)=Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution)*Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                   
                    
                end
            end
        end
    end
    
    Ixvol=smooth3(Ixvol,'gaussian',3,5);
    Iyvol=smooth3(Iyvol,'gaussian',3,5);
    Izvol=smooth3(Izvol,'gaussian',3,5);
    Ixyvol=smooth3(Ixyvol,'gaussian',3,5);
    Iyzvol=smooth3(Iyzvol,'gaussian',3,5);
    Ixzvol=smooth3(Ixzvol,'gaussian',3,5);
    
    for i=floor(windowsize/2)+1:size(volume,1)-floor(windowsize/2)
        for j=floor(windowsize/2)+1:size(volume,2)-floor(windowsize/2)
            for k=floor(windowsize/2)+1:size(volume,3)-floor(windowsize/2)
                if mod(i,tensorresolution)==0 && mod(j,tensorresolution)==0 && mod(k,tensorresolution)==0
                   
                    Ix=Ixvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Iy=Iyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Iz=Izvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Ixy=Ixyvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Iyz=Iyzvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    Ixz=Ixzvol(i/tensorresolution,j/tensorresolution,k/tensorresolution);
                    
                    tensor = [Ix,Ixy,Ixz;Ixy,Iy,Iyz;Ixz,Iyz,Iz];
                    
                   if max(abs(tensor(:)))>0
                        [vec,val]=eig(tensor);
                    
          
                        %normalization of vectors --> after that, the vector
                        %along a structure (usually the smalles gradient)
                        %is the largest.
                        val(1,1)=exp(-gamma*abs(val(1,1)))+epsilon;
                        val(2,2)=exp(-gamma*abs(val(2,2)))+epsilon;
                        val(3,3)=exp(-gamma*abs(val(3,3)))+epsilon;
                        
                        
                       dval=diag(val);
                       
                        [dval,ranking]=sort(dval);
                      
                        
                        val(1,1)=dval(1);
                        val(2,2)=dval(2);
                        val(3,3)=dval(3);
                        
                        vec=vec(:,ranking(3));   
                        top(i/tensorresolution,j/tensorresolution,k/tensorresolution)=1-(abs(dval(1)/dval(3)));
                             
                            %let the vectors show in the same direction if
                            %they are in the same hemishpere 
                            %--> vectors of opposite sides of a tube (like
                            %an axon) show no in the same direction instead
                            %of the opposite
                            if vec(1)<0 && vec(2)<0 && vec(3)<0
                                vec=vec.*-1;
                                
                            end
                            
                            if vec(1)<0 && vec(2)<0 && vec(3)>0
                                vec=vec.*-1;
                                
                            end
                            
                            if vec(1)>0 && vec(2)<0 && vec(3)<0
                                vec=vec.*-1;
                                
                            end
                            
                            if vec(1)>0 && vec(2)<0 && vec(3)>0
                                vec=vec.*-1;
                                
                            end
                           
                        fx(i/tensorresolution,j/tensorresolution,k/tensorresolution)=vec(1);
                        fy(i/tensorresolution,j/tensorresolution,k/tensorresolution)=vec(2);
                        fz(i/tensorresolution,j/tensorresolution,k/tensorresolution)=vec(3);
                        
                        
              
                  end
                end
            end
        end
       
    end
end
   
