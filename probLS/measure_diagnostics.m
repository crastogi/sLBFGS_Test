% Need to make totEpochs and nDataPoints global
clear all
global nDataPoints path nEpochs x_min data batchsize stochIters measureMoves;
global nIncorrectDirs nCorrectDirs truePositive falsePositive trueNegative falseNegative;
% Set parameters
nSamples = 25;
measureMoves = true;

totEpochs = zeros(1, nSamples);
totDist = totEpochs;
totICDirs = totEpochs;
totCDirs = totEpochs;
totTP = totEpochs;
totFP = totEpochs;
totTN = totEpochs;
totFN = totEpochs;
nDataPoints = 0;

for q = 1:nSamples
    run__minimal_example;
    totDist(q) = norm(path(:,end)-x_min{1});
    totEpochs(q) = nEpochs;
    totICDirs(q) = nIncorrectDirs;
    totCDirs(q) = nCorrectDirs;
    totTP(q) = truePositive;
    totFP(q) = falsePositive;
    totTN(q) = trueNegative;
    totFN(q) = falseNegative;
end
disp([num2str(mean(totEpochs)) ' ± ' num2str(sqrt(var(totEpochs))) ' epochs']);
disp([num2str(mean(totDist)) ' ± ' num2str(sqrt(var(totDist))) ' distance to optimum']);

nTotEpochs = sum(totEpochs);
excessEvals = (nDataPoints - ((nTotEpochs+1)*length(data)+2*batchsize*stochIters*nTotEpochs+batchsize*stochIters*nTotEpochs))/batchsize;
nTotEpochs = nTotEpochs/nSamples;
nDataPoints = nDataPoints/(length(data)*nSamples);

disp([num2str(excessEvals/(nSamples*nTotEpochs*stochIters)) ' mean line search evaluations'])
% Correct directions
disp([num2str(100*mean(totCDirs./(totICDirs+totCDirs))) '% ± ' num2str(100*sqrt(var(totCDirs./(totICDirs+totCDirs)))) '% correct directions']);
% FN/TP
disp([num2str(100*mean(totFN./totCDirs)) '% FN/' num2str(100*mean(totTP./totCDirs)) '% TP (± ' num2str(100*sqrt(var(totTP./totCDirs))) '%)']);
% TN/FP
disp([num2str(100*mean(totTN./totICDirs)) '% TN/' num2str(100*mean(totFP./totICDirs)) '% FP (± ' num2str(100*sqrt(var(totFP./totICDirs))) '%)']);
% Precision
disp(['Precision: ' num2str(mean(totTP./(totTP+totFP))) ' ± ' num2str(sqrt(var(totTP./(totTP+totFP))))]);
% Recall
disp(['Recall: ' num2str(mean(totTP./(totTP+totFN))) ' ± ' num2str(sqrt(var(totTP./(totTP+totFN))))]);

% Correct Dirs: False Negative, True Positive
% Incorrect Dirs: True Negative, False Positive