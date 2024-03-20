function output = biexp_T1_fit_function(data, TI)

% General model:
%      val(x) = a+b*exp(-x/S) + c*exp(-x/L)


% Output is: [T1L, T1S, a, b, c, resid]

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

ft = fittype( 'a+b*exp(-x/S) + c*exp(-x/L)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.5 0 -Inf -Inf -Inf];
opts.StartPoint = [0.91 0.79 0.51 0.47 0.61];
opts.Upper = [Inf 0.5 Inf Inf Inf];
xData = TI(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = zeros(x,6);

for i = 1:x
    yData= data(i,:);
    [fitresult, gof] = fit( xData, yData(:), ft, opts );
    output(i, 1:5) = coeffvalues(fitresult);
    output(i,6) = gof.sse;
end


if scaled % restore units
    output(:,1:2) = output(:,1:2)*1000;
end


