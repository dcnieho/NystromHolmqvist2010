function [data] = interpSacIgnace(data,order,start1,fin1,start2,fin2,smark,qDebug)

if nargin<8
    qDebug = false;
end

if isempty(smark),
    disp('NO SACCADES');
else
    isaccsta        = smark(1:2:end);           % isaccsta is index laatste sample voor de saccade
    isaccend        = smark(2:2:end);           % isaccend is index eerste  sample na   de saccade

    for p=1:length(isaccsta),

        coef1 = mypolyfit(start1,order,fin1,data,qDebug);       % het fitten van traject 1
        coef3 = mypolyfit(start2,order,fin2,data,qDebug);       % het fitten van traject 3

        if order == 3,                         % derde orde interpolatie
            a0 = coef1(4);  
            a1 = coef1(3); 
            a2 = coef1(2);
            a3 = coef1(1);
            c1 = coef3(3);
            c2 = coef3(2);
            c3 = coef3(1);
        end

        if order == 2,                         % tweede orde interpolatie
            a0 = coef1(3);  
            a1 = coef1(2); 
            a2 = coef1(1);
            a3 = 0.0;
            c1 = coef3(2);
            c2 = coef3(1);
            c3 = 0.0;
        end

        if order == 1,                         % eerste orde interpolatie
            a0 = coef1(2);  
            a1 = coef1(1); 
            a2 = 0.0;
            a3 = 0.0;
            c1 = coef3(1);
            c2 = 0.0;
            c3 = 0.0;
        end

        t1 = isaccsta(p);                           % start en eind punt van de interpolatie
        t2 = isaccend(p) - 1;

        % c0 = (3*a0 + a1*t1 - c1*t1 + 2*a1*t2 - 2*c1*t2 + 2*a2*t1*t2 - 2*c2*t1*t2 + a2*t2*t2 - c2*t2*t2 + 3*a3*t1*t2*t2 - 3*c3*t1*t2*t2)/3; % uitgecomment want niet gebruikt in Ignace's code
        b0 = (3*a0*t1*t1 + a1*t1*t1*t1 - c1*t1*t1*t1 - 6*a0*t1*t2 + 2*a2*t1*t1*t1*t2 - 2*c2*t1*t1*t1*t2 + 3*a0*t2*t2 + 3*a3*t1*t1*t1*t2*t2 - 3*c3*t1*t1*t1*t2*t2)/(3*( t1 + t2)*( t1 + t2));
        b1 = (c1*t1*t1 - 2*a1*t1*t2 - 2*a2*t1*t1*t2 + 2*c2*t1*t1*t2 + a1*t2*t2 - 3*a3*t1*t1*t2*t2 + 3*c3*t1*t1*t2*t2)/((t1 - t2)*(t1 - t2));
        b2 = (a1*t1 - c1*t1 + a2*t1*t1 - 2*c2*t1*t2 + a2*t2*t2 + 3*a3*t1*t2*t2 - 3*c3*t1*t2*t2)/((-t1 + t2)*(-t1 + t2));
        b3 = (-a1 + c1 + 3*a3*t1*t1 - 2*a2*t2 + 2*c2*t2 - 6*a3*t1*t2 + 3*c3*t2*t2)/(3*(-t1 + t2)*(-t1 + t2));

        % debug plot, data voor correctie
        if qDebug
            figure(3345);
            idxs = [start1:fin2];
            subplot(2,1,1)
            plot(idxs,data(idxs));
        end
        
        if 1
            % bereken de reconstructed pursuit, let op! nare offsets in
            % deze data.
            % We berekenen een extra sample voor en na de saccade. Dan
            % kunnen we de nieuwe pursuit gladjes in de data hangen door de
            % zorgen dat het verschil tussen de berekende data en de echte
            % data nul is.
            new_t = [isaccsta(p)-1:isaccend(p)];
            newdata = b0 + b1*new_t + b2*new_t.^2 + b3*new_t.^3;
            
            % corrigeer offsets aan begin en einde van saccade
            offstart =  newdata(1  )-data(isaccsta(p)-1);
            offeind  = -newdata(end)+data(isaccend(p))+offstart;    % dit is gelijk aan de saccade amplitude!
            % hang geinterpoleerde data erin op de plaats van de saccade
            data(isaccsta(p):isaccend(p)-1) = newdata(2:end-1)      - offstart;
            data(isaccend(p):end)           = data(isaccend(p):end) - offeind;
        else
            % OUDE CODE: gebruik je dit dan heb je een plateau (twee
            % dezelfde oogposities) net voor en na je saccade
            
            % hang geinterpoleerde data erin op de plaats van de saccade
            new_t = [isaccsta(p):isaccend(p)-1];
            data(isaccsta(p):isaccend(p)-1) = b0 + b1*new_t + b2*new_t.^2 + b3*new_t.^3;
            
            % corrigeer offsets aan begin en einde van saccade
            data(isaccsta(p):end) = data(isaccsta(p):end) - (data(isaccsta(p)) - data(isaccsta(p)-1));
            data(isaccend(p):end) = data(isaccend(p):end) - (data(isaccend(p)) - data(isaccend(p)-1));
        end
        
        % debug plot, data voor correctie
        if qDebug
            subplot(2,1,2)
            plot(idxs,data(idxs));
        end
    end
end




% helpers
function [c,t1,y2] = mypolyfit(start,order,eind,data,qDebug)

if nargin<5
    qDebug = false;
end

if qDebug
    disp(sprintf('order = %d',order));
    disp(sprintf('start = %d',start));
    disp(sprintf('einde = %d',eind));
end

t1  = [start:eind];
y1  = data(start:eind).';

% polyfit, scaled and centered. You get that when requesting the third
% output, which we'll need to undo the transformation of x
[c_sc,~,mu]   = polyfit(t1,y1,order);

% undo scaling and centering
% solve A_new*(X_new)^2 +B_new *(X_new) + C = A*(X)^2 + B*(X) + C, where
% X_new is the scaled and centered X defined as X_new = (X 
% mean(X))/std(X), A_new, B_new etc are the polynomial
% coefficients fit to the scaled data, and A, B etc are the polynomial
% coefficients for unscaled data. It can similarly be derived for other
% orders.
if order == 2
    c(3) = c_sc(1)*mu(1)^2/mu(2)^2 - c_sc(2)*mu(1)/mu(2) + c_sc(2);
    c(2) = (c_sc(2)*mu(2) - 2*c_sc(1)*mu(1))/mu(2)^2;
    c(1) = c_sc(1)/mu(2)^2;
elseif order == 1
else
    error('not implemented for this order')
end