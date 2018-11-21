
function cost=fitFunc1(xTrain,yTrain)
% index: an integer vector containing the index of the feature 
% xTrain is the train data completely
% yTrain is the label of train data
%
% 
%
%
%
mcr = crossval('mcr',xTrain,yTrain,'Predfun',@predFunc1);
cost=[mcr size(xTrain,2)];
end
