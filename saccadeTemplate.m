function [eye_pos_vel] = saccadeTemplate(amplitude,duration,tstep,order)

% function saccadeTemplate(amplitude,duration,tstep,order)
%
% Outputs eye position (order==0) in degrees or eye velocity (order==1) in
% degrees/s from the model given in Appendix A of Garcia-Perez and Peli
% (2001). Intrasaccadic Perception. Journal of Neuroscience 21(18):
% 7313-7322.
% The duration in milliseconds given as an input argument, and the output
% is sampled with a timstep of tstep milliseconds.
%
% This corresponds to the minimum snap model for horizontal saccades from
% Harwood, Mezey & Harris (1999). The Spectral Main Sequence of Human
% Saccades. Journal of Neuroscience 19(20): 9098-9106. It should be noted
% that they rejected it in this paper as an accurate description of human
% horizontal saccades.
%
% I, DCN, have found that saccadic velocity profiles, including the skew
% induced by the long tail are fit well by a skewed normal, as might be
% expected from Collewijn, Erkelens & Steinman (1988). Binocular
% co-ordination of human vertical saccadic eye movements. Journal of
% Physiology 404, pp 183-197. Note that the 
% code:
% -------------------------------------------------------------------------
%     % note: (nlinfit is from statistic toolbox, but you could use
%     % fminsearch if you set up your own objective function)
%
%     % sp contains the velocity profile of a single saccade, its what were fitting
%     plot(sp)
%     hold on
%
%     % fit normal first for init
%     normfun = @(p,x) p(3).*normpdf(x,p(1),p(2));
%     startingVals = [20 10 max(sp)*20];
%     coefEsts = nlinfit(1:length(sp), sp, normfun, startingVals);  
%     line(1:length(sp), normfun(coefEsts, 1:length(sp)), 'Color','g');
% 
%     % http://azzalini.stat.unipd.it/SN/Intro/intro.html
%     % http://en.wikipedia.org/wiki/Skew_normal_distribution
%     % params:
%     % 1: mu
%     % 2: sigma
%     % 3: shape (alpha)
%     % 4: scale (data isn't normalized)
%     modelFun =  @(p,x) p(4).* (2./p(2) .* normpdf((x-p(1))./p(2)) .* normcdf(p(3).*((x-p(1))./p(2))));
% 
%     startingVals = [coefEsts(1:2) 5 coefEsts(3)];
%     coefEsts = nlinfit(1:length(sp), sp, modelFun, startingVals);
%     line(1:length(sp), modelFun(coefEsts, 1:length(sp)), 'Color','r');
% -------------------------------------------------------------------------
%
% Lastly, for a very accurate model of horizontal saccades, look at Enderle
% & Zhou (2010). Models of Horizontal Eye Movements II; A 3rd-Order Linear
% Saccade Model.

% sample time around middle of saccade
time = (tstep:tstep:duration) - (duration+tstep)/2;

% convert to seconds for garcia-perez and peli equation A5
time     = time ./ 1000;
duration = duration ./ 1000;

switch order
    case 0
        % eye position
        % Garcia-Perez and Peli, 2001, Eq. A5
        eye_pos_vel = amplitude/2 + (amplitude)   .*     ...
                    (                                    ...
                      (35.*time)   ./(16.*duration)    - ...
                      (35.*time.^3)./(4 .*duration.^3) + ...
                      (21.*time.^5)./(    duration.^5) - ...
                      (20.*time.^7)./(    duration.^7)   ...
                    );
    case 1
        % eye velocity - note that this is
        % Garcia-Perez and Peli, 2001, Eq. A6
        eye_pos_vel = 35.*amplitude ./ (16.*duration)  .* ...
                    (                                     ...
                      1-(4.*time.^2)./(duration.^2)       ...
                    ).^3;
    otherwise
        error('No saccade template available for order %d',order)
end