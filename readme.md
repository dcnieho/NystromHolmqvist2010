The code in this repository is a reimplementation of Nyström, M. & Holmqvist, K. (2010), "An adaptive algorithm for fixation, saccade, and glissade detection in eye-tracking data". Behavior Research Methods 42(1): 188-204. It processes the recorded eye movement data to extract saccades, fixations, and glissades (the latter are now often called post-saccadic oscillations). When using this code, in addition to Nystr, please cite Niehorster, Siu & Li (2015). [See below](#citation).

# Differences from original implementation
First, the internals of the algorithm have been rewritten extensively with an eye on increasing performance. Furthermore, quite a few additions have been made. This is a non-exhaustive list:

  - optionally (and by default) less strict about NaN/missing data during events.
  - calculates a lot of information about classified events.
  - option to detrend the gaze velocity data by median filtering it and then subtracting the output of this median filter from the gaze velocity trace.
  - option to cross-correlate (detrended or original) gaze velocity data with a saccade template for additional noise robustness in saccade classification.
  - can detect blinks by means of thresholding change-of-pupil size signal. This enables also detecting partial blinks.
  - option to merge nearby saccades/glissades
  - option to determine saccade onset by fitting straight line to accelerating flank of velocity profile of saccade, and seeing where this line intersects with vel==0. This was used to get precise saccade onset times in a low-sampling frequency data (see Oliva et al. (2017) below).
  - includes post-processing utilities to merge glissades into their associated saccades, remove saccades from the position and velocity traces, and keep only the saccades in the position and velocity traces.
  - detailed plotting functionality

# Citation
When using this code, please cite Niehorster, Siu & Li (2015). If using ETparams.saccade.onsetRefineMethod=2, please additionally cite Oliva, Niehorster, Jarodzka & Holmqvist (2017). Example citation:
> Saccades were classified using the Niehorster, Siu & Li (2015) implementation of the Nyström & Holmqvist (2010) algorithm, with default settings. In addition, saccade onsets were determined using the method of Oliva, Niehorster, Jarodzka & Holmqvist (2017).

NB: it is probably good to discuss these methods in a few lines each. It is furthermore *important* that if you change settings from their default, you note these changes in your article.

References:
> Nyström, M. & Holmqvist, K. (2010), "An adaptive algorithm for fixation, saccade, and glissade detection in eye-tracking data". Behavior Research Methods 42(1): 188-204. doi: 10.3758/BRM.42.1.188
>
> Niehorster, D.C., Siu, W.W.F., & Li, L. (2015). Manual tracking enhances smooth pursuit eye movements. Journal of Vision 15(15), 11. doi: 10.1167/15.15.11
>
> Oliva, M., Niehorster, D.C., Jarodzka, H., & Holmqvist, K. (2017). Social Presence Influences Saccadic and Manual Responses. I-Perception 8(1). doi: 10.1177/2041669517692814

# Usage
Example data from two experiments is provided with this implementation.

  1. In the folder `NiehorsterSiuLi2015`, data from two participants from one condition from Niehorster Siu & Li (2015) (citation above) is provided. This dataset contains pursuit and saccades. Run `eventClassificationNiehorsterSiuLi2015.m` to run event classification, including creating of data traces containing only pursuit and only saccades.
  2. In the folder `pictureViewing`, data from two participants, 10 trials each, is provided, recorded during an unpublished experiment. This dataset contains fixations and saccades. Run `eventClassificationPictureViewing.m` to run event classification, including implicit classification of fixations in the data.
