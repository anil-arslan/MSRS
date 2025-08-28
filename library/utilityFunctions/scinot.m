function str = scinot(value, significantFigure)
    if nargin < 2
        significantFigure = 1;
    end
    numberOfValues = numel(value);
    significand = 10.^mod(log10(value), 1);
    exponent = floor(log10(value));
    str = strings(1, numberOfValues);
    for valueID = 1 : numberOfValues
        if sign(exponent(valueID)) == -1
            str(valueID) = sprintf(['%0.' num2str(significantFigure) 'f•10^{-%i}'], significand(valueID), abs(exponent(valueID)));
        else
            str(valueID) = sprintf(['%0.' num2str(significantFigure) 'f•10^%i'], significand(valueID), exponent(valueID));
        end
    end
end