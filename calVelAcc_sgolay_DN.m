function data = calVelAcc_sgolay_DN(data,ETparams)

% prepare parameters
%--------------------------------------------------------------------------
% span of filter, use minimum length of saccade
span = ceil(ETparams.minSaccadeDurm/1000*ETparams.samplingFreq);
% number of tabs of filter
F = 2*ceil(span)-1;


% store unfiltered data unfiltered data
%--------------------------------------------------------------------------
data.Xpix = data.X;
data.Ypix = data.Y;


% compute eye positions, velocities, and accelerations in pixels, convert
% to degree
%--------------------------------------------------------------------------

if 0 % some debugging, can delete it
    % this have better onset transient response and correct phase
    X2 = sgFilt(X,0,F);
    V2 = sgFilt(X,1,F) * ETparams.samplingFreq;
    A2 = sgFilt(X,2,F) * ETparams.samplingFreq^2;
    
    % getvelacc from ignace
    [V3,A3,X3] = getvelacc(X,span,ETparams.samplingFreq);
    
    ETparams.data(i,j).X = filter(g(:,1),1,X);
    ETparams.data(i,j).Y = filter(g(:,1),1,Y);
    ETparams.data(i,j).velX = filter(g(:,2),1,X) * ETparams.samplingFreq;
    ETparams.data(i,j).accX = filter(g(:,3),1,X) * ETparams.samplingFreq^2;
    plot(X,'k');                          hold on; plot(X2,'r');plot(X3,'g');plot(ETparams.data(i,j).X,'b');
    
    figure
    hold on; plot(X,'k--'); plot(V2,'r'); plot(-V3,'g');plot(ETparams.data(i,j).velX,'b');
    Xs=[X(20:40); X2(20:40); X3(20:40); ETparams.data(i,j).X(20:40)].'
    Vs=[X(20:40); V2(20:40); V3(20:40); ETparams.data(i,j).velX(20:40)].'
    As=[X(20:40); A2(20:40); A3(20:40); ETparams.data(i,j).accX(20:40)].'
    pause
end


% Calculate the filtered position, velocity and acceleration
[tempP,tempV,tempA] = sgFilt([data.X; data.Y],[0 1 2],F);

% convert to degree and calculate derivative magnitudes
data.X      = tempP(1,:) / ETparams.angleInPixelsH;
data.Y      = tempP(2,:) / ETparams.angleInPixelsV;

velX        = tempV(1,:) / ETparams.angleInPixelsH;
velY        = tempV(2,:) / ETparams.angleInPixelsV;
data.vel    = hypot(velX, velY) * ETparams.samplingFreq;

accX        = tempA(1,:) / ETparams.angleInPixelsH;
accY        = tempA(2,:) / ETparams.angleInPixelsV;
data.acc    = hypot(accX, accY) * ETparams.samplingFreq^2;




%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps, [],2);     % along second dimension (column) as that is how data is currently organized
end