function [x] = pattern_search(x0, inputfunction)
    x = x0;
    d = 1;
    bestF = inputfunction(x);
    nDims = size(x0,2);
    
    while (true)
        isUpdate = 0;
        for currDim = 1:nDims
            currPos = x;
            currPos(currDim) = currPos(currDim)+d;
            fUp = inputfunction(currPos);
            currPos(currDim) = currPos(currDim)-2*d;
            fDn = inputfunction(currPos);
            if (fUp < bestF)
                x(currDim) = x(currDim)+d;
                bestF = fUp;
                isUpdate = 1;
            elseif (fDn < bestF)
                x(currDim) = x(currDim)-d;
                bestF = fDn;
                isUpdate = 1;
            end
        end
        disp(bestF);
        if (isUpdate == 0) 
            d = d/2;
            disp(['Reduced step size; d = ', num2str(d)]);
            if d < 1E-7
                break;
            end
        end
    end
end