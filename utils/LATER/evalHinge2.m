function [out,qEarly] = evalHinge2(data,params,qIgnoreEarly,mode)

assert(mode~=2||qIgnoreEarly,'qIgnoreEarly must be true for fit mode 2')

% for early component, only slope varies. intercept is always 0 (in probit
% space, so 50% point)
earlyInt = 0;
% others from input
laterInt    = params(1);
laterSlope  = params(2);
if ~qIgnoreEarly
    earlySlope  = params(3);
end

% for points before intersection of two lines (judged along x) calc
% distance to early line, for rest to later line.
if qIgnoreEarly
    qEarly = false(size(data,1),1);
else
    intersection = (laterInt-earlyInt)/(earlySlope-laterSlope);
    qEarly = data(:,1)>intersection;
end




switch mode
    case 1
        % least squares
        dist( qEarly) = data( qEarly,2)-(earlyInt+earlySlope*data( qEarly,1));
        dist(~qEarly) = data(~qEarly,2)-(laterInt+laterSlope*data(~qEarly,1));
        out = sum(dist.^2);
    case 2
        % komolgorov-smirnov
        mu = -laterInt./laterSlope;
        sd = -1./laterSlope;
        datTest = (data(:,1)-mu)./sd;
        [~,p] = kstest(datTest);
        out = 1-p;
end


if 0
    clf
    r = [min(data(:,1)) max(data(:,1))];
    plot(data(:,1),data(:,2),'r')
    hold on
    if ~qIgnoreEarly
        if any(qEarly)
            idx = find(qEarly); idx(end+1) = idx(end)+1;
            plot(data(idx,1),data(idx,2),'c')
        end
        plot(r,earlySlope*r,'b')
    end
    plot(r,laterSlope*r+laterInt,'k')
    
    pause
end