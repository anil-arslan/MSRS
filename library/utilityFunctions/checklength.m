function outputData = checklength(inputData, sizeData, message)
    %checklength Summary of this function goes here
    %   Detailed explanation goes here
    arguments (Input)
        inputData
        sizeData (1, :) {mustBeInteger, mustBePositive} = 1
        message (1, :) char = 'lenght must be'
    end
    arguments (Output)
        outputData
    end
    if isstruct(inputData)
        if sizeData == 1
            sizeData = max(structfun(@numel, inputData));
        end
        for field = fieldnames(inputData).'
            outputData.(char(field)) = checklength(inputData.(char(field)), sizeData);
        end
    else
        if isscalar(inputData)
            outputData = repmat(inputData, 1, sizeData);
        elseif length(inputData) ~= sizeData
            error([message, ' ', num2str(sizeData)]);
        else
            outputData = inputData;
        end
    end
end