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

function [ data ] = loadMaskFromMesh(meshFileName,isProj)
    if isProj

        fid = fopen(meshFileName);
        counter = 0;
        data = {};
        while ~feof(fid)
            currentLine = fgetl(fid);

            if isempty(currentLine)
                continue;
            end

            % Skip the header part until we come across the first piece of data
            if counter == 0 && currentLine(1)~='@'
                continue;
            end

           % Check if new piece of data encountered.
            if currentLine(1)=='@'
                counter = counter+1;
                % Get ready to read the data.
                data{counter} = []; 
            else
                % This is a true data line, read it.
                data{counter} = [data{counter}; sscanf(currentLine,'%f')'];
            end
        end

        fclose(fid);
        data=data{1};

    else
         raw = load_am_data(meshFileName);
         data=zeros(size(raw));
         data(raw>0)=1;
    end
end

