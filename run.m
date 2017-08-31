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

windowsize=5;
tensorresolution=1;
redSize=5;
gamma=0.0001;
smoothRange=5;
smoothSigma=2;
gvf1=0.05;
gvf2=15;
maxSRange=1.4;
searchpatternBorder=15;

readDir='neuron_retrival_amira_files';
inputDir='neuron_retrieval_vector_volumes5';


%data preprocessing
generateVectorVolumes(readDir,inputDir,windowsize,tensorresolution,redSize,gamma,smoothRange,smoothSigma)

%image retrieval
%Manual segmentation
patImgs=[-1,1;-1,2;-1,3;-1,4;-1,5;-1,6];
[precisionsM,recallsM,tpsM] = neuronRetrieval3D(patImgs,inputDir,readDir,maxSRange,redSize,gvf1,gvf2,0);


patImgs=[];
for run=1:6
    anno=getGroundTruth(run,readDir);
    patImgsTemp=1:size(anno,1);
    patImgsTemp=patImgsTemp(anno(:,2)==1);
    runTemp=1:size(patImgsTemp,2);
    runTemp(:)=run;
    patImgs=[patImgs;[patImgsTemp;runTemp]'];
    
end
    
%Automatic segmentation
[precisionsA,recallsA,tpsA] = neuronRetrieval3D(patImgs,inputDir,readDir,maxSRange,redSize,gvf1,gvf2,1);

%No segmentation
[precisionsN,recallsN,tpsN] = neuronRetrieval3D(patImgs,inputDir,readDir,maxSRange,redSize,gvf1,gvf2,2);

%Voxel based approach
[precisionsV,recallsV,tpsV] = neuronRetrievalVoxelBased(patImgs,inputDir,readDir,redSize);


%Converts the precisions and recalls in a form that can be plotted
results_prec_ms_single=[];
results_rec_ms_single=[];
results_prec_as_single=[];
results_rec_as_single=[];
results_prec_ns_single=[];
results_rec_ns_single=[];
results_prec_vox_single=[];
results_rec_vox_single=[];

results_prec_as=[];
results_rec_as=[];
results_prec_ns=[];
results_rec_ns=[];
results_prec_vox=[];
results_rec_vox=[];

for run=1:6
    

    for i=1:6
            results_prec_ms_single=[results_prec_ms_single;precisionsM{i}.pre];
            results_rec_ms_single=[results_rec_ms_single;recallsM{i}.rec];
    end
    
    tempPrec=[];
    tempRec=[];
    for i=1:size(patImgs,1)
        if precisionsA{i}.run==run
            tempRec=[tempRec;recallsA{i}.rec];
            tempPrec=[tempPrec;precisionsA{i}.pre];
        end
        
        if run~=2 && precisionsA{i}.patImg==983
            results_prec_as_single=[results_prec_as_single;precisionsA{i}.pre];
            results_rec_as_single=[results_rec_as_single;recallsA{i}.rec];
        end
        
         if run==2 && precisionsA{i}.patImg==668
            results_prec_as_single=[results_prec_as_single;precisionsA{i}.pre];
            results_rec_as_single=[results_rec_as_single;recallsA{i}.rec];
        end
    end
    results_prec_as=[results_prec_as;trimmean(tempPrec,50)];
    results_rec_as=[results_rec_as;trimmean(tempRec,50)];
    
    tempPrec=[];
    tempRec=[];
    for i=1:size(patImgs,1)
        if precisionsN{i}.run==run
            tempRec=[tempRec;recallsN{i}.rec];
            tempPrec=[tempPrec;precisionsN{i}.pre];
        end
        
        if run~=2 && precisionsN{i}.patImg==983
            results_prec_ns_single=[results_prec_ns_single;precisionsN{i}.pre];
            results_rec_ns_single=[results_rec_ns_single;recallsN{i}.rec];
        end
        
         if run==2 && precisionsN{i}.patImg==668
            results_prec_ns_single=[results_prec_ns_single;precisionsN{i}.pre];
            results_rec_ns_single=[results_rec_ns_single;recallsN{i}.rec];
        end
    end
    results_prec_ns=[results_prec_ns;trimmean(tempPrec,50)];
    results_rec_ns=[results_rec_ns;trimmean(tempRec,50)];

    tempPrec=[];
    tempRec=[];
    for i=1:size(patImgs,1)
        if precisionsV{i}.run==run
            tempRec=[tempRec;recallsV{i}.rec];
            tempPrec=[tempPrec;precisionsV{i}.pre];
        end
        
        if run~=2 && precisionsV{i}.patImg==983
            results_prec_vox_single=[results_prec_vox_single;precisionsV{i}.pre];
            results_rec_vox_single=[results_rec_vox_single;recallsV{i}.rec];
        end
        
         if run==2 && precisionsV{i}.patImg==668
            results_prec_vox_single=[results_prec_vox_single;precisionsV{i}.pre];
            results_rec_vox_single=[results_rec_vox_single;recallsV{i}.rec];
        end
    end
    results_prec_vox=[results_prec_vox;trimmean(tempPrec,50)];
    results_rec_vox=[results_rec_vox;trimmean(tempRec,50)];

end

%Precision-Recall plot
p=plot(mean(results_rec_ms_single),mean(results_prec_ms_single));
set(p,'Color','blue','LineWidth',1)
hold on
p=plot(mean(results_rec_as_single),mean(results_prec_as_single));
set(p,'Color','green','LineWidth',1)
hold on
p=plot(mean(results_rec_as),mean(results_prec_as),'--');
set(p,'Color','green','LineWidth',1)
hold on
p=plot(mean(results_rec_ns_single),mean(results_prec_ns_single));
set(p,'Color','magenta','LineWidth',1)
hold on
p=plot(mean(results_rec_ns),mean(results_prec_ns),'--');
set(p,'Color','magenta','LineWidth',1)
hold on
p=plot(mean(results_rec_vox_single),mean(results_prec_vox_single));
set(p,'Color','red','LineWidth',1)
hold on
p=plot(mean(results_rec_vox),mean(results_prec_vox),'--');
set(p,'Color','red','LineWidth',1)
hold off
legend({'ms','as','truncated mean as','ns','truncated mean ns','vox','truncated mean vox'})
ylabel('Precision');
xlabel('Recall');