function calVelAcc_sgolay(i,j)
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
[b,g] = sgolay(N,F);   % Calculate S-G coefficients

% Extract relevant gaze coordinates for the current trial.
X = ETparams.data(i,j).X;
Y = ETparams.data(i,j).Y;

% Calculate the velocity and acceleration
ETparams.data(i,j).X = filter(g(:,1),1,X);
ETparams.data(i,j).Y = filter(g(:,1),1,Y);

ETparams.data(i,j).velX = filter(g(:,2),1,X);
ETparams.data(i,j).velY = filter(g(:,2),1,Y);
ETparams.data(i,j).vel = sqrt(ETparams.data(i,j).velX.^2 + ETparams.data(i,j).velY.^2)/angleInPixelsH*ETparams.samplingFreq;

ETparams.data(i,j).accX = filter(g(:,3),1,X);
ETparams.data(i,j).accY = filter(g(:,3),1,Y);
ETparams.data(i,j).acc = sqrt(ETparams.data(i,j).accX.^2 + ETparams.data(i,j).accY.^2)/angleInPixelsH*ETparams.samplingFreq^2;

