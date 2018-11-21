
% this version is knncls1 but uses 5-fold cross validation to estimate
% fitness of individuals
function y=predFunc1(X,Y,Z)
% X is the train data 
% Y is the label of train data
% Z is the test data
%
%
%% uncmment only one to test other method
% Model= fitcknn(X,Y,'NumNeighbors',1);
% Model= fitcknn(X,Y,'NumNeighbors',3);
% Model=fitcnb(X,Y,'Distribution','kernel');
% Model=fitcnb(X,Y,'Distribution','normal');
Model= fitcecoc(X,Y);

y=predict(Model,Z);

end

%% Distribution
% 'kernel' Kernel smoothing density estimate.
% 'mn'
% Multinomial distribution. If you specify mn, then all features are components of a multinomial distribution. Therefore, you cannot include 'mn' as an element of a cell array of strings. For details, see Algorithms.
% 'mvmn'
% Multivariate multinomial distribution. For details, seeAlgorithms.
% 'normal'
% Normal (Gaussian) distribution.

