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

%Partial derivation of a nxnxn 3D Matrix
% mat = matrix
function deriv = part_deriv(mat,dir)
    

    deriv=0;
    for k=1:size(mat,1) 
        for j=1:size(mat,2)
            for i=1:size(mat,3)
                filter=(i-(floor(size(mat,3)/2)+1));
                if dir==1 && not(filter==0)
                   
                    deriv=deriv+((double(mat(i,j,k)))/filter);
                else if dir==2 && not(filter==0)
                        deriv=deriv+((double(mat(j,i,k)))/filter);
                    else if dir==3 && not(filter==0)
                        deriv=deriv+((double(mat(j,k,i)))/filter);
                        end
                    end
                end
            end
        end
    end
    
end
