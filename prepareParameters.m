function ETparams = prepareParameters(ETparams)

% calculate extends of screen around the origin (we'll change the data's
% origin in the same way in prepareData)
ETparams.screen.rect.pix = [...
     [0 0]                      - ETparams.screen.subjectStraightAhead, ...
     ETparams.screen.resolution - ETparams.screen.subjectStraightAhead
    ];

% calculate screen extends in degree
pixPerMeter     = ETparams.screen.resolution ./ ETparams.screen.size;
screenExtends   = ETparams.screen.rect.pix ./ [pixPerMeter pixPerMeter];
ETparams.screen.rect.deg = atan2(screenExtends,ETparams.screen.viewingDist)*180/pi;