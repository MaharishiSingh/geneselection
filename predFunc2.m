
function y=predFunc2(X,Y,Z)
% X is the train data completely
% Z is the test data completely
% Y is the label of train data
%
%
%
%
Model= fitcknn(X,Y,'NumNeighbors',3);
y=predict(Model,Z);

end
