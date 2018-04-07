# CatecholamineEEG

Code for my behavioural and EEG analysis on the effects of catecholamines, conducted at Monash University, analysed in Trinity College Dublin. 
Includes data processing/analysis files (matlab/eeglab) and Inferential Statistics/data analysis files (R)

**Loughnane**, G.M., Barnes, J.M., Nandam, L.S., O’Connell, R.J. & Bellgrove, M.A. (2018). Catecholamine modulation of evidence accumulation during perceptual decision formation, in review.

## Full Analysis Pipeline:

### In Matlab, run scripts to process raw eeg files:
(note: can skip steps 1 and 2 because they are just for finding noisy eeg channeles, which have already been found listed in 'runafew_BL.m' for step 3)

1. 'get_chanvars_drugs.m' - to extract a measure of electrode noise
2. 'check4badchans_drugs.m' - checks for noisy channeles, you manually list any noisy channeles in 'runafew_BL.m'
3. 'runafew_drugs.m' - calls 'p3b_epochs.m' to extract ERP, Alpha, and Pupil Diameter epochs 
4. 'p3b_epochs.m' - extracts and saves individualised ROI electrodes for each participant's pre-target alpha measurements 
5. 'p3b_plot.m' - plots relevant ERPs and extracts subject-by-subject data (behavioural, EEG/ERP) and saves to .csv for import into R for inferential analysis

### Run Inferential Statistics/analysis in R (with Rstudio):
1. anova_covar_drugs.R
2. lme_drugs.R
3. mediation_drugs.R

## Figures From Study

### Figure 1

![alt text](https://github.com/gerontium/CatecholamineEEG/blob/master/Figure_1_low_res.png)

### Figure 2

![alt text](https://github.com/gerontium/CatecholamineEEG/blob/master/Figure_2_low_res.png)

### Figure 3

![alt text](https://github.com/gerontium/CatecholamineEEG/blob/master/Figure_3_low_res.png)

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/3.0/"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-nc-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/3.0/">Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License</a>.
