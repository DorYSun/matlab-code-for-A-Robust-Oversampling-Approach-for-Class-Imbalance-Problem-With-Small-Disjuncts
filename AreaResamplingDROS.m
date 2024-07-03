function [trainMinorityResamplingNNG,noiseNoSampleGenIndex]=AreaResamplingDROS(trainMinority,trainMajority,NNG,project1,K_Maj,cAngle,g)
%% % % % % % % % % % % % % % % % % % % % 
K_Maj=min([K_Maj size(trainMajority,1)]);
%% % % % step 1 distance computation and sorting
[N_min,dim]=size(trainMinority);

[I_sortdistMinToMaj,~]= knnsearch(trainMajority,trainMinority, 'k',K_Maj);



barrierKMaj=trainMajority; 
N_barKMaj=size(barrierKMaj,1);
chosedMajRM=trainMajority; 

%% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % revised in 211209, save memory
RMatrix=single(zeros(N_min));% % % % compute the relationship between minority samples: 1 means related
ValueMaxAngleMajMatrix=single(zeros(N_min));% % % % compute the related project value
IndexMaxAngleMajMatrix=single(zeros(N_min));% % % % record the related index
PL1='ABCDEFGHIJKLMN';
PL2='ABCDEFGHIJKLMN';
%% % % % step 2 relationship computation

for i=1:N_min
% 
%     i=i;
%     tarminI=N_min-i+1;
%     tarminI=412;
    iMin=trainMinority(i,:);
    RMinIndex=i+1:N_min;
    remainedMin=trainMinority(RMinIndex,:);
    N_RMin=size(remainedMin,1);
    drectioniMinTochosedMajRM=-1*bsxfun(@minus,chosedMajRM,iMin);
    drectioniMinTochosedMajRMNorm = bsxfun(@rdivide,drectioniMinTochosedMajRM,sqrt(sum(drectioniMinTochosedMajRM.^2,2)));
    drectioniMinTochosedMajRMNorm(isnan(drectioniMinTochosedMajRMNorm))=0;
    drectioniMinTochosedMajRMNorm(isinf(drectioniMinTochosedMajRMNorm))=0;
    minValue=zeros(1,N_RMin);
    minIndex=zeros(1,N_RMin);
    for j=1:N_RMin
        jRMin=remainedMin(j,:);
        drectioniRMinTochosedMajRM=-1*bsxfun(@minus,chosedMajRM,jRMin);
        drectioniRMinTochosedMajRMNorm = bsxfun(@rdivide,drectioniRMinTochosedMajRM,sqrt(sum(drectioniRMinTochosedMajRM.^2,2)));
        drectioniRMinTochosedMajRMNorm(isnan(drectioniRMinTochosedMajRMNorm))=0;
        drectioniRMinTochosedMajRMNorm(isinf(drectioniRMinTochosedMajRMNorm))=0;
        project=sum(drectioniMinTochosedMajRMNorm.*drectioniRMinTochosedMajRMNorm,2);
        [minValue(1,j), minIndex(1,j)]= min(project,[],1);
    end
    pindex=find(minValue>=project1);

    RMatrix(i,RMinIndex(pindex))=1;
    RMatrix(RMinIndex(pindex),i)=1;    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    ValueMaxAngleMajMatrix(i,RMinIndex(pindex))=minValue(pindex);
    ValueMaxAngleMajMatrix(RMinIndex(pindex),i)=minValue(pindex);
    IndexMaxAngleMajMatrix(i,RMinIndex(pindex))=minIndex(pindex);
    IndexMaxAngleMajMatrix(RMinIndex(pindex),i)=minIndex(pindex);
end

%% % % % step 3-5 searchlight sacnning structure of each minority sample
cellInf=cell(N_min,10);% save the center, central angle, vertex and .ect for each searchlight structure
N_nonoise=0;% % % % number of normal minortiy samples
noiseNoIrelIndex=[];
noiseNoScannedMajIndex=[];
noiseNoPosRIndex=[];
for i=1:N_min
    iMin=trainMinority(i,:);
    irelIndex=find(RMatrix(i,:)==1);
    if(any(size(irelIndex)==0))
% % % % % %\*210316
        noiseNoIrelIndex(end+1)=i;
% % % % % %*/210316
        continue;
    end
    irelMin=trainMinority(irelIndex,:);
%% % % % step 3 related samples and final selected samples
%% % % % % % % % %compute the center of selected sampels and the target one
    nnOneMaj=trainMajority(I_sortdistMinToMaj(i,1),:);
    nnKMaj=trainMajority(I_sortdistMinToMaj(i,1:K_Maj),:);
    meannnKMaj=mean(nnKMaj,1);
    iMinMinusmeannnKMaj=iMin-meannnKMaj;
    iMinMinusmeannnKMajNorm=iMinMinusmeannnKMaj./repmat(sqrt(sum(iMinMinusmeannnKMaj.^2,2)),1,dim);
    iMinMinusmeannnKMajNorm(isnan(iMinMinusmeannnKMajNorm))=0;
    % % % % % % % % % % % % % % average projects as the vertex
    irelMinMinusiMin=irelMin-repmat(iMin,size(irelMin,1),1);
    projects=irelMinMinusiMin*iMinMinusmeannnKMajNorm';
    projectsPos01=projects>0;
    
    meanlength=sum(projects(projectsPos01))/sum(projectsPos01);
    meanlength(isnan(meanlength))=0;
    CircleCenter=meanlength*iMinMinusmeannnKMajNorm+iMin;
    aaa=1;
    finalCIrelMin=irelMin;

    
%% % % % step 5 compute the radius     
    referedDirection=-1*iMinMinusmeannnKMajNorm;
    BarKMajMinuscirCenter=barrierKMaj-repmat(CircleCenter,N_barKMaj,1);
    BarKMajMinuscirCenterNorm=BarKMajMinuscirCenter./repmat(sqrt(sum(BarKMajMinuscirCenter.^2,2)),1,dim);
    BarKMajMinuscirCenterNorm(isnan(BarKMajMinuscirCenterNorm))=0;
    projectSD=referedDirection*BarKMajMinuscirCenterNorm';
    sdbKMajIndex=projectSD>=cAngle;
    fallInScannedAreaMaj=barrierKMaj(sdbKMajIndex,:);
    distimin=pdist2(iMin,CircleCenter,'euclidean');
    distTofallInScannedAreaMaj=pdist2(CircleCenter,fallInScannedAreaMaj,'euclidean');
    [mindist, minIndex]=min(distTofallInScannedAreaMaj,[],2);
    radius=min([distimin mindist])+0.5*(mindist-distimin);
    radius2=0;
    if( size(fallInScannedAreaMaj,1)==0)
        continue;
    end
    if(radius<=0)
        continue;
    end
    if(sum(abs(referedDirection))==0)
        fprintf('zero vector refD; \n');
        continue;
    end
    if(isempty(radius))
        fprintf('empty radius; \n');
        continue;
    end
    N_nonoise=N_nonoise+1;
    cellInf{N_nonoise,1}=radius;
    cellInf{N_nonoise,2}=radius2;
    cellInf{N_nonoise,3}=cAngle;
    cellInf{N_nonoise,4}=CircleCenter;
    cellInf{N_nonoise,5}=referedDirection;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    cellInf{N_nonoise,6}=iMin;
    cellInf{N_nonoise,7}=nnOneMaj;
    cellInf{N_nonoise,8}=nnKMaj;
    cellInf{N_nonoise,9}=finalCIrelMin;
    cellInf{N_nonoise,10}=fallInScannedAreaMaj;
end
noiseNoSampleGenIndex.noiseNoIrelIndex=noiseNoIrelIndex;
noiseNoSampleGenIndex.noiseNoPosRIndex=noiseNoPosRIndex;
noiseNoSampleGenIndex.noiseNoScannedMajIndex=noiseNoScannedMajIndex;
if(N_nonoise==0)
    trainMinorityResamplingNNG=[];
    return;
end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
overGSDN=0;
NForOnce=ceil(NNG/N_nonoise);
for i=1:N_nonoise
    radius=cellInf{i,1};
    radius2=cellInf{i,2};
    cAngle=cellInf{i,3};
    CircleCenter=cellInf{i,4};
    referedDirection=cellInf{i,5};
%% % random Length and normalized rand direction vector 
    randomLength=( rand([1 NForOnce])*(1-g) + g ) *radius;   
    randNumEachTime=NForOnce;
    aRNorm=RandmDirectionForLAO201216(cAngle,referedDirection,randNumEachTime,dim); 
    project_i_aR=aRNorm*referedDirection';
    maxProject=max(project_i_aR);   

    conditionIndex=find(project_i_aR>=cAngle-0.000000001);
    conditionIndexSize=size(conditionIndex,1);
    satisfieCISize=conditionIndexSize;
    satisfiedaRNorm=aRNorm(conditionIndex,:);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    if(satisfieCISize==0)
      fprintf('no one existed rand direction and continue; \n');
      continue;     
    else
        randpermSatisfiedD=randi([1 satisfieCISize],1,NForOnce);
        chosedaRNorm=satisfiedaRNorm(randpermSatisfiedD(1:NForOnce),:);
    end
    chosedNewPoints=diag(randomLength)*chosedaRNorm+repmat(CircleCenter,NForOnce,1);
    trainMinorityResampling(overGSDN+1:overGSDN+NForOnce,:)=chosedNewPoints;
    overGSDN=overGSDN+NForOnce;         
end
randpermRes=randperm(size(trainMinorityResampling,1));
trainMinorityResamplingNNG=trainMinorityResampling(randpermRes(1:NNG),:);
end

function [finalRandDirection]=RandmDirectionForLAO201216(CentralAngle,NormDirectionFromCCToimin,randNumEachTime,dim)
%% % % % % % % % % % % % % % % % % % % % % % % 
randForOneNew=100;
limitMaxCentralAngle=0.99;

if(CentralAngle>=0.99)
    finalRandDirection=NormDirectionFromCCToimin;
    return;
else
    randDirection= randi([1 100],[randForOneNew*randNumEachTime dim])-50;
    aRNorm=randDirection./repmat(sqrt(sum(randDirection.^2,2)),1,dim);
end
project_i_aR=aRNorm*NormDirectionFromCCToimin';
maxProject=max(project_i_aR);
selectedRandDirection=[];
conditionIndex=find(project_i_aR>=CentralAngle);
selectedRandDirection=[selectedRandDirection; aRNorm(conditionIndex,:)];

count2=0;
while(maxProject<limitMaxCentralAngle)
    count2=count2+1;
    randDirection=aRNorm+1*repmat(NormDirectionFromCCToimin,randForOneNew*randNumEachTime,1);
    aRNorm=randDirection./repmat(sqrt(sum(randDirection.^2,2)),1,dim);
    project_i_aR=aRNorm*NormDirectionFromCCToimin';
    maxProject=max(project_i_aR);
    conditionIndex=find(project_i_aR>=CentralAngle);
    selectedRandDirection=[selectedRandDirection; aRNorm(conditionIndex,:)];
end

project_i_aR=selectedRandDirection*NormDirectionFromCCToimin';
index0_1=find(project_i_aR>0 & project_i_aR<=0.1);
index0_2=find(project_i_aR>0.1 & project_i_aR<=0.2);
index0_3=find(project_i_aR>0.2 & project_i_aR<=0.3);
index0_4=find(project_i_aR>0.3 & project_i_aR<=0.4);
index0_5=find(project_i_aR>0.4 & project_i_aR<=0.5);
index0_6=find(project_i_aR>0.5 & project_i_aR<=0.6);
index0_7=find(project_i_aR>0.6 & project_i_aR<=0.7);
index0_8=find(project_i_aR>0.7 & project_i_aR<=0.8);
index0_9=find(project_i_aR>0.8 & project_i_aR<=0.9);
index0_10=find(project_i_aR>0.9 & project_i_aR<=1);
indexV=[{index0_1} {index0_2} {index0_3} {index0_4} {index0_5} {index0_6} {index0_7} {index0_8} {index0_9} {index0_10}];

size0_1=size(index0_1,1);
size0_2=size(index0_2,1);
size0_3=size(index0_3,1);
size0_4=size(index0_4,1);
size0_5=size(index0_5,1);
size0_6=size(index0_6,1);
size0_7=size(index0_7,1);
size0_8=size(index0_8,1);
size0_9=size(index0_9,1);
size0_10=size(index0_10,1);

sizeV=[size0_1 size0_2 size0_3 size0_4 size0_5 size0_6 size0_7 size0_8 size0_9 size0_10];
Non0Ind=find(sizeV>0);
sizeNon0Ind=size(Non0Ind,2);
minSizeNon0=min(sizeV(Non0Ind));

finalIndex=[];
for ii=1:sizeNon0Ind
    size0_ii=sizeV(Non0Ind(ii));
    randpermD0_ii=randperm(size0_ii);
    
    
    index0_ii=indexV{Non0Ind(ii)};
    selectedIndex0_ii=index0_ii(randpermD0_ii(1:minSizeNon0));
    
    finalIndex=[finalIndex selectedIndex0_ii'];
end

secondSelectedRandDirection=selectedRandDirection(finalIndex,:);


selectedSize=size(secondSelectedRandDirection,1);
randpermD=randperm(selectedSize);
finalSize=min([selectedSize randForOneNew*randNumEachTime]);
finalRandDirection=secondSelectedRandDirection(randpermD(1:finalSize),:);
project_i_aR=finalRandDirection*NormDirectionFromCCToimin';

if(any(size(find(project_i_aR>CentralAngle))==0))
    fprintf('no one existed rand direction in RandDForLAO; \n');
    aaa=1;
end
end
