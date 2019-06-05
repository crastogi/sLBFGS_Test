% Need to make totEpochs and nDataPoints global
clear all
global totEpochs totDist nDataPoints nSamples path nEpochs x_min data batchsize stochIters;
totEpochs = 0;
totDist = 0;
nDataPoints = 0;
nSamples = 50;

for q = 1:nSamples
    run__minimal_example;
    totDist = totDist + norm(path(:,end)-x_min{1});
    totEpochs = totEpochs + nEpochs;
    %pause(2);
end 
excessEvals = (nDataPoints - ((totEpochs+1)*length(data)+2*batchsize*stochIters*totEpochs+batchsize*stochIters*totEpochs))/batchsize;
totEpochs = totEpochs/nSamples;
nDataPoints = nDataPoints/(length(data)*nSamples);
totDist = totDist/nSamples;
excessEvals/(nSamples*totEpochs*stochIters)
