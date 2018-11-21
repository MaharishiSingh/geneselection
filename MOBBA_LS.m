path1=genpath('C:\Users\MAHARISHI SINGH\Desktop\Dissertation Research\multi_objective bat algorithm\MOBBA_LS\data3');
addpath(path1);

disp('path added');

clear
close all
clc

localSearchMethod=0;
set=1;

if set==1
    S='srbct3.mat';
elseif set==2
    S='Leukemia3.mat';
elseif set==3
    S='Prostate3.mat';
    end
load(S);

nPop=20;
numVar=500;
data=data';
[indx scor]=FisherSc(data(itrain,:),label(itrain));
data=data(:,indx(1:numVar));

xTrain=data(itrain,:);
yTrain=label(itrain);

%% Cost Function

CostFunction=@(x,y) fitFunc1(x,y);      % Cost Function


%% MOBBA7 Parameters

VarSize=[1 numVar];%feaure space size

Qmin=0;         % Frequency minimum
Qmax=1;         % Frequency maximum
delta=Qmax-Qmin;

sigma=.7;% for velocity, Equatio3

Vmin=1;% Minimum Velocity
Vmax=[];% No bounadary as we use a sigmoid function 1/(1+e^-x)

Amin=0;%Loadness
Amax=1;
Loudness=Amax;
alpha=.9;%loudness decreasing factor
V=zeros(1,nPop);
pulseTap=.25;% there are four operations: estimate/increase/decrease/injection
MaxIt=30;

injRate=.01;
extRate=.01;

f1NumInt=zeros(1,MaxIt);
f1NumExt=zeros(1,MaxIt);

%% Initialization

empty_bat.Position=[];
empty_bat.Cost=[];
empty_bat.Rank=[];
empty_bat.DominationSet=[];
empty_bat.DominatedCount=[];
empty_bat.CrowdingDistance=[];
empty_bat.Velocity=[];
empty_bat.Frequency=[];
empty_bat.Loudness=[];
empty_bat.PulsNum=[];

pop=repmat(empty_bat,nPop,1);

for i=1:nPop
    
    pop(i).Position=logical(randi(2,VarSize)-1);
    pop(i).Velocity=ones(VarSize)*Vmin;
    pop(i).Frequency=(Qmax-Qmin)*rand;
    
    %pop(i).Loudness=Amax;
    pop(i).PulsNum=1;
    
    pop(i).Cost=CostFunction(xTrain(:,pop(i).Position),yTrain);
    
end

%% compute global best

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Sort Population
F1=pop(F{1});
nF1=numel(F1);
Stopping_Criterion=8;%we can perform post analysis to obtain lower number of genes

%% main loop
it=0;
while true
    if (pop(1).Cost(2)<Stopping_Criterion)
        r=reshape([pop.Cost],[2 numel(pop)]);
        r1=r(1,:);
        f=find(r==min(r1),1);
        if(r(2,f)<Stopping_Criterion)
            break;
        end
    end
    if (pop(1).Cost(2)<Stopping_Criterion)
        break;
    end
    it=it+1;
    %% move bats
    for i=1:nPop
        %% generate new solution using Eq 3,4,7
        
        % use gBest randomly from Fron1
        pop(i).Frequency=delta*rand;
        
        % use gBest randomly from Fron1
        selectedfInd=randi(nF1);
        gBest=F1(selectedfInd);
        
        %obtain velocity
        Velocity = pop(i).Velocity+((pop(i).Position-gBest.Position)*pop(i).Frequency);
        V(i)=sum(Velocity);
        
        newPosition=((1./(1+exp(-Velocity)))>sigma);
        newCost = CostFunction(xTrain(:,newPosition),yTrain);
        
        %replace global best
        if(Dominates(newCost,gBest.Cost))
            F1(selectedfInd).Position=newPosition;
            F1(selectedfInd).Cost=newCost;
            F1(selectedfInd).PulsNum=F1(selectedfInd).PulsNum+1;%gBest=no pulse tap
        end
        
        %% replace new solution
        if(rand<Loudness && Dominates(newCost,pop(i).Cost)) %if obtained a better solution
            % Evaluate each solution using a classifir
            pop(i).Position=newPosition;
            pop(i).Cost=newCost;
            %update loudness
            pop(i).PulsNum=pop(i).PulsNum+pulseTap;
        end
    end
    Loudness=alpha*Loudness;
    
    %% non-dominated sorting
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);
    
    % Store F1
    F1=pop(F{1});
    nF1=numel(F1);
    f1NumInt(it)=nF1;
    
    %% local search
    for i=1:nPop
        if(localSearchMethod==2 || localSearchMethod==0)
            p1=pop(i);
            if(Loudness<rand) %whatever the loudness decreases,
                nn=randi(numel(F1));
            else
                nn=randi(nPop);
            end
            p2=pop(nn);%concentrate on population
            
            %injection
            LInd=randperm(numVar,floor(injRate*sum(p2.Position))+1) ;
            p1.Position(LInd)=p2.Position(LInd);%injection
            
            %investigate
            p1.Cost = CostFunction(xTrain(:,p1.Position),yTrain);
            if(Dominates(p1.Cost,pop(i).Cost))
                pop(i)=p1;
                pop(i).PulsNum=pop(i).PulsNum+pulseTap;
            end
        end
        
        if(localSearchMethod==1 || localSearchMethod==0)
            %% method1 only mutating
            pRate=1-(1/(1+log(1+pop(i).PulsNum)));
            if(pRate>rand)%whatever the pulse rate increases, more concentration puts over leader set
                nn=randi(nPop);
            else
                nn=randi(numel(F1));
            end
            p1=pop(nn);
            
            % first: decrease
            posInd=find(p1.Position);
            npos=numel(posInd);
            pInd=randperm(npos,floor(extRate*npos+1));
            p1.Position(posInd(pInd))=false;
            
            p1.Cost = CostFunction(xTrain(:,p1.Position),yTrain);
            if(Dominates(p1.Cost,pop(nn).Cost))
                pop(nn)=p1;
                pop(nn).PulsNum=pop(nn).PulsNum+pulseTap;
            end
            
            % then: increase
            posInd=find(~p1.Position);
            npos=numel(posInd);
            pInd=randperm(npos,floor(extRate*npos+1));
            p1.Position(posInd(pInd))=true;
            
            p1.Cost = CostFunction(xTrain(:,p1.Position),yTrain);
            if(~Dominates(pop(nn).Cost,p1.Cost))
                pop(nn)=p1;
                pop(nn).PulsNum=pop(nn).PulsNum+pulseTap;
            end
        end
        
    end
    
    %% Obtain Global Best
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);
    
    F1=pop(F{1});
    nF1=numel(F1);
    f1NumExt(it)=numel(F1);
    
    %% visualization
    %FF(it).F1=[F1.Cost];
    if(mod(it,5)==0)
        disp(['Iteration ' num2str(it) '-F1Num:' num2str(numel(F1))...
            '--cost:' num2str(F1(1).Cost(1))...
            '--num:' num2str(F1(1).Cost(2))...
            ]);
    end
end
%% Results

%% evaluate results
C=[];
P=[];
dtrain=data(itrain,:);
dtest=data(itest,:);
ltrain=label(itrain);
ltest=label(itest);

M=[1 2 3];%methods 1 KNN, 2 NBY, 3 SVM , if you want to evaluate only using some classifiers just you need to change the M to for instance M=[1,3]
for rr=1:numel(pop)
    for method=M
        selGenes=pop(rr).Position;
        data0=dtrain(:,selGenes);
        if method==1
            Model= fitcknn(data0,ltrain,'NumNeighbors',3);
        elseif method==2
            Model=fitcnb(data0,ltrain,'Distribution','normal');
        elseif method==3
            Model= fitcecoc(data0,ltrain);%'polynomial','PolynomialOrder',3,'Solver','SMO');
        end
        y=predict(Model,dtest(:,selGenes));% predict the current sample (one out)
        missNum=sum(y~=ltest);
        error=(missNum/numel(ltest));% in percent
        accuracy(rr,method)= 1-error;
    end
end
accuracy= accuracy
pop.Cost

