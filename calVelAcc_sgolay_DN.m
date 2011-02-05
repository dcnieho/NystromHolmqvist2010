function calVelAcc_sgolay_DN(i,j)
global ETparams

% Lowpass filter window length
smoothInt = ETparams.minSaccadeDur; % in seconds

% Span of filter
span = ceil(smoothInt*ETparams.samplingFreq);

% Calculate how many degrees one pixel spans.
[angleInPixelsH, angleInPixelsV] = degrees2pixels(1,ETparams.viewingDist...
    ,ETparams.screenSz,ETparams.screenDim);

% Calculate unfiltered data
%--------------------------------------------------------------------------
ETparams.data(i,j).Xorg = ETparams.data(i,j).X;
ETparams.data(i,j).Yorg = ETparams.data(i,j).Y;

ETparams.data(i,j).velXorg = [0 diff(ETparams.data(i,j).X)]/angleInPixelsH*ETparams.samplingFreq;
ETparams.data(i,j).velYorg = [0 diff(ETparams.data(i,j).Y)]/angleInPixelsV*ETparams.samplingFreq;
ETparams.data(i,j).velOrg = sqrt(ETparams.data(i,j).velXorg.^2 + ETparams.data(i,j).velYorg.^2);


% Pixel values, velocities, and accelerations
%--------------------------------------------------------------------------
N = 2;                 % Order of polynomial fit
F = 2*ceil(span)-1;    % Window length
% [b,g] = sgolay(N,F);   % Calculate S-G coefficients

% % Extract relevant gaze coordinates for the current trial.
% X = ETparams.data(i,j).X;
% Y = ETparams.data(i,j).Y;

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

    plot(ETparams.data(i,j).velXorg,'k'); hold on; plot(X,'k--'); plot(V2,'r'); plot(-V3,'g');plot(ETparams.data(i,j).velX,'b');
    Xs=[X(20:40); X2(20:40); X3(20:40); ETparams.data(i,j).X(20:40)].'
    Vs=[X(20:40); V2(20:40); V3(20:40); ETparams.data(i,j).velX(20:40)].'
    As=[X(20:40); A2(20:40); A3(20:40); ETparams.data(i,j).accX(20:40)].'
    pause
end


% Calculate the filtered velocity and acceleration
[tempP,tempV,tempA] = sgFilt([ETparams.data(i,j).X; ETparams.data(i,j).Y],[0 1 2],F);

ETparams.data(i,j).X = tempP(1,:);
ETparams.data(i,j).Y = tempP(2,:);

ETparams.data(i,j).velX = tempV(1,:);
ETparams.data(i,j).velY = tempV(2,:);
ETparams.data(i,j).vel  = hypot(ETparams.data(i,j).velX, ETparams.data(i,j).velY)/angleInPixelsH*ETparams.samplingFreq;

ETparams.data(i,j).accX = tempA(1,:);
ETparams.data(i,j).accY = tempA(2,:);
ETparams.data(i,j).acc  = hypot(ETparams.data(i,j).accX, ETparams.data(i,j).accY)/angleInPixelsH*ETparams.samplingFreq^2;


%%% helpers
function varargout = sgFilt(x,difforder,ntaps)
% wrapper for convenient syntax, fitting always using 2nd order polynomial
for p=1:length(difforder)
    varargout{p} = savitzkyGolayFilt(x,2,difforder(p),ntaps, [],2);
end