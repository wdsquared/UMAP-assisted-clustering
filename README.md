# UMAP-assisted clustering of COVID-19 clinical data

Code for the paper: "Machine learning of COVID-19 clinical data identifies population structures with therapeutic potential". David Greenwood (1,2), Tom Taverner (3), Nicola Adderley (3), Malcolm Price (3,4), Krishna Gokhale (3), Chris Sainsbury (3), Suzy Gallier (5,6), Carly Welch (5,6), Elizabeth Sapey (5,6,7), Duncan Murray (1,5), Hilary Fanning (5), Simon Ball (1,5,7), Krishnarajah Nirantharakumar (3,7), Wayne Croft (1,2), Paul Moss (1,5). Submitted (2021).

1 - Institute of Immunology and Immunotherapy, University of Birmingham, Birmingham, UK. 
2 - The Centre for Computational Biology, University of Birmingham, Birmingham, UK.
3 - Institute of Applied Health Research, University of Birmingham, Birmingham, UK.
4 - NIHR Birmingham Biomedical Research Centre, University Hospitals Birmingham NHS Foundation Trust and University of Birmingham, Birmingham, UK
5 - University Hospitals Birmingham NHS Foundation Trust, Birmingham, UK. 
6 - Institute of Inflammation and Ageing, University of Birmingham, Birmingham, UK. 
7 - Health Data Research UK

# Contact 

Lead - Paul Moss - p 'dot' moss 'at' bham 'dot' ac 'dot' uk

GitHub - David Greenwood - drg707 'at' student 'dot' bham 'dot' ac 'dot' uk, wdsquared 'dot' projects 'at' gmail 'dot' com

# Summary

Clinical outcomes for patients admitted to hospital with acute COVID-19 are heterogeneous and there is interest in defining patient subgroups for prognostic modelling and development of targeted treatment algorithms.

We obtained 28 demographic and laboratory variables obtained on the day of admission in three independent cohorts of patients admitted to hospital with COVID-19. These comprised a training cohort (n= 6099) and two validation cohorts during the first and second waves of the pandemic (n=996; n=1011). Uniform manifold approximation and projection (UMAP) dimension reduction and Gaussian mixture model (GMM) analysis was used to define patient clusters in the training cohort and applied to confirmation cohorts.

GMM analysis defined 29 separate patient clusters within the training cohort and this population structure was directly transferable to patients within the validation cohorts.  Individual clusters were associated with markedly different rates of mortality which were highly predictive within the confirmation datasets. 

Deconvolution of clinical features within each cluster identified unexpected relationships between clinical variables that might potentially be tar-getable for therapeutic management. Integration of large and diverse clinical datasets using UMAP-assisted clustering can therefore identify novel patient subgroups. This can be used to provide detailed prognostic information and has the potential to uncover unexpected interactions between clinical variables. This application of machine learning represents a powerful approach for delineating mechanisms of disease pathogenesis and guiding potential novel therapeutic interventions for disorders of unmet need."

# Content 

The modelling framework, validation procedure and cluster characterisation are implemented in UMAP-clustering-pipeline-clean-March1st-2022.R.

The imputation procedure is implemented in "Cohort-imputation-clean-March-1st-2022.R".

# Data

The clinical data reported in the study cannot be deposited in a public repository because of ethical constraints.
To request access to the development data, contact the Geriatric Medicine Research Collaborative (https://www.gemresearchuk.com/, gemresearch 'dot' uk 'at' gmail.com). 
To request access to the validation data, contact the PIONEER Health Data Research Hub (https://www.pioneerdatahub.co.uk/, pio-neer 'at' uhb 'dot' nhs 'dot' uk).
