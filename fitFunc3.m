function cost=fitFunc3(xTrain,yTrain)
% index: an integer vector containing the index of the feature 
% xTrain is the train data completely
% yTrain is the label of train data
%
% 
%
%
%
mcr = crossval('mcr',xTrain,yTrain,'Predfun',@predFunc3);
cost=[mcr size(xTrain,2)];

end
