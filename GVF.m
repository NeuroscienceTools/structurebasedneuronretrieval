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

%   GVF of an gradient vector map f.  mu is the GVF regularization coefficient
%   and ITER is the number of iterations that will be computed.  

function [u,v,w] = GVF(fx,fy,fz, mu, ITER)

fx = BoundMirrorExpand(fx);  % Take care of boundary condition
fy = BoundMirrorExpand(fy);  % Take care of boundary condition
fz = BoundMirrorExpand(fz);  % Take care of boundary condition

u = fx; v = fy; w  =   fz;        % Initialize GVF to the gradient
SqrMagf = fx.*fx + fy.*fy + fz.*fz; % Squared magnitude of the gradient field

% Iteratively solve for the GVF u,v
for i=1:ITER,
    
  u = BoundMirrorEnsure(u);
  v = BoundMirrorEnsure(v);
  w = BoundMirrorEnsure(w);
  u = u + mu*6*del2(u) - SqrMagf.*(u-fx);
  v = v + mu*6*del2(v) - SqrMagf.*(v-fy);
  w = w + mu*6*del2(w) - SqrMagf.*(w-fz);

end

u = BoundMirrorShrink(u);
v = BoundMirrorShrink(v);
w = BoundMirrorShrink(w);


end
