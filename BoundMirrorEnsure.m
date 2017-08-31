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

function B = BoundMirrorEnsure(A)
% Ensure mirror boundary condition
%
% The number of rows and columns of A must be greater than 2
%
% for example (X means value that is not of interest)
% 
% A = [
%     X  X  X  X  X   X
%     X  1  2  3  11  X
%     X  4  5  6  12  X 
%     X  7  8  9  13  X 
%     X  X  X  X  X   X
%     ]
%
% B = BoundMirrorEnsure(A) will yield
%
%     5  4  5  6  12  6
%     2  1  2  3  11  3
%     5  4  5  6  12  6 
%     8  7  8  9  13  9 
%     5  4  5  6  12  6
%

% Chenyang Xu and Jerry L. Prince, 9/9/1999
% http://iacl.ece.jhu.edu/projects/gvf

[m,n,l] = size(A);

if (m<3 | n<3 | l<3) 
    error('either the number of rows or columns is smaller than 3');
end

yi = 2:m-1;
xi = 2:n-1;
zi = 2:l-1;
B = A;
B([1 m],[1 n],[1 l]) = B([3 m-2],[3 n-2],[3 l-2]);  % mirror corners
B([1 m],xi,zi) = B([3 m-2],xi,zi);          % mirror left and right boundary
B(yi,[1 n],zi) = B(yi,[3 n-2],zi);          % mirror top and bottom boundary
B(yi,xi,[1 l]) = B(yi,xi,[2 l-2]);          % mirror front and back boundary
