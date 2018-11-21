function [index score ] = FisherSc(X,Y)
%   returns the most significant features first.
%   Fisher Score Methods written by Ashkan Dashtban dashtban@ieee.org
%   this code is copyright. plz contact me if needed.
%   X, the feature set, each raw is an instance
%   Y, the label
%   
%    Usage:
%    [index score]=FisherSc(X,Y);
%
%   note that label can be nominal and the dataset can be sparse!
%   please cite us if you use this code!
%   copyright @ Dashtban@ieee.org, ashkan0021@gmail.com
%
%
% note the higher the value of score the more valuable the feature.
%
%
% end

if (issparse(X))
    X=full(X);
end

if (size(X,1)~=size(Y,1))
    error('number of samples and Labels is different in X and Y');
end

clabel=unique(Y);%class label
numC=length(clabel);

%[~, numF] = size(X);
[ind, numF] = size(X);
score = zeros(1,numF);

% some statistic for classes
cInd = cell(numC,1);% 
Ni = zeros(numC,1); % cardinality of each class

for j = 1:numC
    
    %convert to a ordinal format from nomina or other cases
    cInd{j} = find(Y(:)==clabel(j));
    
    %store the cardinality of each class
    Ni(j) = length(cInd{j});
    
end

% calculate score for each features
for i = 1:numF
    enum = 0;
    denom = 0;
    fi = X(:,i);
    ui = mean(fi);
    
    for j = 1:numC
        
%mean of the number of samples that fall within the jth class
        uij = mean(fi(cInd{j}));
        
% coresponding variance                
        Vij = var(fi(cInd{j}),1);
        
%num of samples in jth class * (uij-ui)^2 
        enum = enum + (uij-ui)^2; 
        denom = denom +  Vij;
    end
    
    if enum == 0
        score(i) = 0;
    else
        if denom == 0
            score(i) = 1000;
        else
            score(i) = enum/denom;
        end
    end
end

%[~, out.fList] = sort(out.W, 'descend');
[ score1 index1] = sort(score, 'descend');
score=score1;
index=index1;
