# fabisearch 0.0.3.6

### Edits ###
1. Reproducibility using doRNG loops. For the refitting and statistical inference procedures of the main detect.cps function, we now utilize doRNG loops which ensures that the results are reproducible.
2. Changed some of the functionality for the est.net function. Now, it includes a changepoints variable which can be used to input the calculated changepoints from the data. This allows multiple stationary networks to be estimated with only a single call of the function.
3. The net.3dplot function has been adjusted with some more easy to use ways of specifying communities. 

### Minor Changes ###
1. Minor changes to some of the text/explanation of how the functions work as well as the included datasets.
2. Included more fMRI datasets and changed fmridata into gordfmri. There is now gordfmri and AALfmri which are the same dataset excep parcellated according to different atlases (Gordon atlas and AAL atlas respectively).
3. Added atlases (gordatlas and AALatlas) for use with the net.3dplot function.
