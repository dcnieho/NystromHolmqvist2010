README-file of the event detection algorithm accompanying the article
Nyström, M. & Holmqvist, K. (in press), "An adaptive algorithm for fixation, saccade, and glissade detection in eye-tracking data". Behavior Research Methods
 
If you have any questions or comments on the paper or thecode, send me an email at
marcus.nystrom@humlab.lu.se. 

Observe that the algorithm is suitable ONLY for data collected from viewers
keeping their heads relatively still while watching static stimuli. It is not 
designed to handle data containing smooth pursuit movements.

The algorithm is implemented in Matlab and tested on version R2008a.
Usage: Run the file: beginEventDetection.m 

The command 'load('1250Hz_3_Participants.mat');' in the file beginEventDetection.m 
loads raw gaze positions from three participants into the variable ETdata, which holds data in a structure
array, sorted after participant and trial.
Notice that the example data contains portions or poor data (where something when wrong during recording). This is to show how the 
algorithms handles noisy data.

Example: ETdata(2,6).X contains all x-coordinates for participant 2 recorded during 
trial 6. Similarily, ETdata(2,6).Y contains corresponding y-coordinates.

Results are stored in the folder 'DetectionResults' in the file 'DetectionResults.mat'. The results can be accessed from the
variable ETparams.
