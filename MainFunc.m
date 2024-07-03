clc;
clear;
close all;
%% % % % % %                                                        12     13        14       15       16 
col=[2];
% % % % % % % 1       2    
Allmethod=[{'Ori'} {'DROS'}];
Cmethod=Allmethod(col);
N_Method=size(col,2);
% % % % % % % % % % % % % % % % % % % % % % % % % % 
axisX1X2Y1Y2=[1 7 1 7; 0 1 0 1;-1.5 1.5 -1.5 1.5];
datasetname=[{'dataset'};];
abcdef=['(a) ';'(b) ';'(c) ';'(d) ';'(e) ';'(f) ';'(g) ';'(h) ';'(i) ';'(j) ';'(k) ';'(l) ';'(m) ';'(n) ';];
FontSizeSet=8;
MarkerS=1;
MarkerS_Min=2;
MarkerS2=3;
DS=1;
OneCircleOneRing=importdata( '.\OneCircleOneRing.mat');
dataset=[OneCircleOneRing.Maj;OneCircleOneRing.Min];
[totalSize,size2D]=size(dataset);
minorityIndex=find(dataset(:,size2D)==1);
minorityDataset=dataset(minorityIndex,:);
minoritySize=size(minorityIndex,1);
majorityIndex=find(dataset(:,size2D)==0);
majorityDataset=dataset(majorityIndex,:);
majoritySize=size(majorityIndex,1);
K_fold=5;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
randpermMinIndex=randperm(minoritySize);
randpermMajIndex=randperm(majoritySize);
fold_min_size=floor(minoritySize/K_fold);
fold_maj_size=floor(majoritySize/K_fold); 
% % % % % % % % % %randomly dividing the current block into train set and test set
i=1;
testMinIndex=randpermMinIndex((i-1)*fold_min_size+1:i*fold_min_size);
testMajIndex=randpermMajIndex((i-1)*fold_maj_size+1:i*fold_maj_size);
testIndex=[minorityIndex(testMinIndex);majorityIndex(testMajIndex)];
testSet1=dataset(testIndex,1:size2D-1);
testLabel1=dataset(testIndex,size2D);
% % % %     
trainMinIndex1=randpermMinIndex(1:(i-1)*fold_min_size);
trainMinIndex2=randpermMinIndex(i*fold_min_size+1:minoritySize);
trainMajIndex1=randpermMajIndex(1:(i-1)*fold_maj_size);
trainMajIndex2=randpermMajIndex(i*fold_maj_size+1:majoritySize);
trainIndex=[minorityIndex(trainMinIndex1);minorityIndex(trainMinIndex2);majorityIndex(trainMajIndex1);majorityIndex(trainMajIndex2)];
trainSet1=dataset(trainIndex,:);
trainSet1=dataset;  
trainMinorityIndex=find(trainSet1(:,size2D)==1);
trainMajorityIndex=find(trainSet1(:,size2D)==0);
trainMinority=trainSet1(trainMinorityIndex,1:size2D-1);
trainMajority=trainSet1(trainMajorityIndex,1:size2D-1);    
trainMinorityNumber=size(trainMinority,1);
trainMajorityNumber=size(trainMajority,1);
N=ceil(trainMajorityNumber/trainMinorityNumber);
NNG=trainMajorityNumber-trainMinorityNumber;        
[trainMinorityResampling{2},~]=AreaResamplingDROS(trainMinority,trainMajority,NNG,-0.7660,7,0.5,1);

axisX1=min(dataset(:,1));
axisX2=max(dataset(:,1));
axisY1=min(dataset(:,2));
axisY2=max(dataset(:,2));
figure;
subplot(1,2,1);
h1=plot(trainMajority(:,1),trainMajority(:,2),'k.','MarkerSize',MarkerS);
hold on
h2=plot(trainMinority(:,1),trainMinority(:,2),'ro','MarkerSize',MarkerS_Min,'MarkerFaceColor','r');
xlabel('(a) Original data','FontSize',FontSizeSet);
axis([axisX1 axisX2 axisY1 axisY2]);

nm=1;
trainMinorityRes=trainMinorityResampling{col(nm)};
nmethod=abcdef(nm+1,:);
methodName=Allmethod{col(nm)};
len=size(methodName,2);
nmethod(end+1:end+len)=methodName;
subplot(1,2,nm+1);
plot(trainMajority(:,1),trainMajority(:,2),'k.','MarkerSize',MarkerS);hold on
plot(trainMinority(:,1),trainMinority(:,2),'ro','MarkerSize',MarkerS_Min,'MarkerFaceColor','r');hold on
h3=plot(trainMinorityRes(:,1),trainMinorityRes(:,2),'g+','MarkerSize',MarkerS);
xlabel(nmethod,'FontSize',FontSizeSet);
axis([axisX1 axisX2 axisY1 axisY2]);
legend([h1 h2 h3],'majority','minority','new minority');
aa=1;
