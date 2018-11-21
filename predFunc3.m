
function y=predFunc3(X,Y,Z)
% X is the train data completely
% Z is the test data completely
% Y is the label of train data
%
%
%
%
Model=fitcnb(X,Y,'Distribution','kernel');% could be normal or mvmn or even mn
y=predict(Model,Z);

end
