% Need to make totEpochs and nDataPoints global
clear all
global nDataPoints path nEpochs x_min data batchsize stochIters measureMoves;
global CDLCFD CDLCFI CDLIFD CDLIFI IDLCFD IDLCFI IDLIFD IDLIFI;
% Set parameters
nSamples = 10;
measureMoves = true;


totEpochs = zeros(1, nSamples);
metrics = zeros(nSamples, 8);
totDist = totEpochs;
nDataPoints = 0;
for q = 1:nSamples
    run__minimal_example;
    totDist(q) = norm(path(:,end)-x_min{1});
    totEpochs(q) = nEpochs;
    metrics(q,:) = [CDLCFD CDLCFI CDLIFD CDLIFI IDLCFD IDLCFI IDLIFD IDLIFI];
end
disp([num2str(mean(totEpochs)) ' ± ' num2str(sqrt(var(totEpochs))) ' epochs']);
disp([num2str(mean(totDist)) ' ± ' num2str(sqrt(var(totDist))) ' distance to optimum']);

% Replace this once sLBFGS loop is replaced
nTotEpochs = sum(totEpochs);
excessEvals = (nDataPoints - ((nTotEpochs+1)*length(data)+2*batchsize*stochIters*nTotEpochs+batchsize*stochIters*nTotEpochs))/batchsize;
nTotEpochs = nTotEpochs/nSamples;
nDataPoints = nDataPoints/(length(data)*nSamples);
disp([num2str(excessEvals/(nSamples*nTotEpochs*stochIters)) ' mean line search evaluations'])

% Correct/Incorrect directions
totCD = sum(metrics(:,1:4),2);
totICD= sum(metrics(:,5:8),2);
totDirs = totCD + totICD;
% TP, FP, TN, FN
totTP = sum(metrics(:,1:2),2);
totFN = sum(metrics(:,3:4),2);
totTN = sum(metrics(:,5:6),2);
totFP = sum(metrics(:,7:8),2);
% correct direction/correct move ratio
CDFD = (metrics(:,1) + metrics(:,3))./totDirs;
CDFI = (metrics(:,2) + metrics(:,4))./totDirs;
IDFD = (metrics(:,5) + metrics(:,7))./totDirs;
IDFI = (metrics(:,6) + metrics(:,8))./totDirs;

disp([num2str(100*mean(totCD./(totICD+totCD))) '% ± ' num2str(100*sqrt(var(totCD./(totICD+totCD)))) '% correct directions']);
% FN/TP
disp([num2str(100*mean(totFN./totCD)) '% FN/' num2str(100*mean(totTP./totCD)) '% TP (± ' num2str(100*sqrt(var(totTP./totCD))) '%)']);
% TN/FP
disp([num2str(100*mean(totTN./totICD)) '% TN/' num2str(100*mean(totFP./totICD)) '% FP (± ' num2str(100*sqrt(var(totFP./totICD))) '%)']);
% Precision
disp(['Precision: ' num2str(mean(totTP./(totTP+totFP))) ' ± ' num2str(sqrt(var(totTP./(totTP+totFP))))]);
% Recall
disp(['Recall: ' num2str(mean(totTP./(totTP+totFN))) ' ± ' num2str(sqrt(var(totTP./(totTP+totFN))))]);
% Correct/incorrect direction move ratio
disp(['Correct Move/Correct Direction: ' num2str(100*mean(CDFD)) '% ± ' num2str(sqrt(var(CDFD))) ...
    '; Incorrect Move/Correct Direction: ' num2str(100*mean(CDFI)) '% ± ' num2str(sqrt(var(CDFI)))]);
disp(['Correct Move/Incorrect Direction: ' num2str(100*mean(IDFD)) '% ± ' num2str(sqrt(var(IDFD))) ...
    '; Incorrect Move/Incorrect Direction: ' num2str(100*mean(IDFI)) '% ± ' num2str(sqrt(var(IDFI)))]);
    
% Correct Dirs: False Negative, True Positive
% Incorrect Dirs: True Negative, False Positive