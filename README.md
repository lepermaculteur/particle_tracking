Adaptation and extension of original particle tracking code written by Nicholas T. Ouellette, September 2010  

Noticeable add-on (February 2022):  
-Added a selection on the eccentricity of the detected regions to decide whether or not they are real particles  
-Implementation of a "in-focus" threshold : For each particle, a small square sample (centered on the particle center) of 12x12pixels is taken from the background
removed image and sent to a "focus measurement algorithm" (see 'Helmli and Scherer's mean method').
If the value returned is below a certain threshold, the particle is considered as 'out of focus' and 
is rejected. The "focus measurement algorithm" has been selected among about twenty other algorithms presented in
"Pertuz et al. / Pattern Recognition (2013)" after being tested on calibration images.  
-Correction of  minor indexation errors in the initialization and finalization step of PredictiveTracker  
