function output = T1_fit_function(data, TI)

% Finds estimates of T1, |a|, and |b| using a nonlinear least
% squares approach together with polarity restoration. 
% The model +-|ra + rb*exp(-t/T1)| is used. 

% Output is: [T1, a, b, resid]

% Modifcation of qMRlab code.

[x,y] = size(data);

if y ~= length(TI)
    error('Check input, number of columns must == number of TIs')
end

% To set start options, lets assume T1 is in seconds
scaled = 0;
if max(TI) > 100 % in milliseconds
    TI = TI/1000;
    scaled = 1;
end


ft = fittype( 'a+b*exp(-x/T)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0, 0, -Inf];
opts.StartPoint = [1, 1, 1e-4];
opts.Upper = [10, 10, Inf];
xData = TI(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = zeros(x,4);

for i = 1:x
    yData= data(i,:);
    [fitresult, gof] = fit( xData, yData(:), ft, opts );
    output(i, 1:3) = coeffvalues(fitresult);
    output(i,4) = gof.sse;
end


if scaled % restore units
    output(:,1) = output(:,1)*1000;
end


