\name{pbi} 
\alias{pbi} 
\alias{pbi-class} 
\title{Environment Object with Functions for the Course Practical Bioinformatics} 
\description{ 
The functions within the pbi environment are utility functions used in the  
course Practical Bioinformatics at the University of Potsdam. 
} 
\section{Methods}{ 
\itemize{ 
\item \code{\link[pbi:pbi_clusterSilhouette]{pbi$clusterSilhouette}} - determine clusterr quality using silhouette index 
\item \code{\link[pbi:pbi_clusterSimIndex]{pbi$clusterSimIndex}} - determine clusterr quality using silhouette index 
\item \code{\link[pbi:pbi_cohensD]{pbi$cohensD}} - effect size between two means 
\item \code{\link[pbi:pbi_cohensH]{pbi$cohensH}} - effect size for a 2x2 contingency table 
\item \code{\link[pbi:pbi_cohensW]{pbi$cohensW}} - effect size for a 2x2 o larger contingency tables 
\item \code{\link[pbi:pbi_cv]{pbi$cv}} - coefficient of variation 
\item \code{\link[pbi:pbi_dassoc]{pbi$dassoc}} - assocplot with residual shading 
\item \code{\link[pbi:pbi_dcorr]{pbi$dcorr}} - pairwise correlations and p-values for data frames and matrices 
\item \code{\link[pbi:pbi_dcorrplot]{pbi$dcorrplot}} - visualize a pairwise correlation matrix 
\item \code{\link[pbi:pbi_df2md]{pbi$df2md}} - convert a data frame or matrix into Markdown code 
\item \code{\link[pbi:pbi_dist]{pbi$dist}} - extension of the stats::dist function with more measures for binary data and correlation 
\item \code{\link[pbi:pbi_domainplot]{pbi$domainplot}} - plot protein domains similar to those on the prosite website 
\item \code{\link[pbi:pbi_dpairs]{pbi$dpairs}} - improved pairs plot considering the data types 
\item \code{\link[pbi:pbi_dpairs.legend]{pbi$dpairs.legend}} - adding legends to pairs plots 
\item \code{\link[pbi:pbi_epsilonSquared]{pbi$epsilonSquared}} - effect size for comparing three or more means of non-normal data 
\item \code{\link[pbi:pbi_etaSquared]{pbi$etaSquared}} - effect size for comparing three or more means of normal data 
\item \code{\link[pbi:pbi_file.cat]{pbi$file.cat}} - show a text file within the R-console 
\item \code{\link[pbi:pbi_file.head]{pbi$file.head}} - show the first lines of a text file 
\item \code{\link[pbi:pbi_grect]{pbi$grect}} - add colored background and grid lines to an existing plot 
\item \code{\link[pbi:pbi_impute]{pbi$impute}} - missing value imputation 
\item \code{\link[pbi:pbi_lmplot]{pbi$lmplot}} - xyplot with regression line and correlation coefficient 
\item \code{\link[pbi:pbi_mi]{pbi$mi}} - mutual information 
\item \code{\link[pbi:pbi_modelQuality]{pbi$modelQuality}} - model quality for classification and regression models 
\item \code{\link[pbi:pbi_modus]{pbi$modus}} - most often qualitative level in a variable 
\item \code{\link[pbi:pbi_msa2pwm]{pbi$msa2pwm}} - position weight matrix for an alignment 
\item \code{\link[pbi:pbi_mw]{pbi$mw}} - molecular weight of an protein sequence 
\item \code{\link[pbi:pbi_package.deps]{pbi$package.deps}} - display package dependencies 
\item \code{\link[pbi:pbi_pastel]{pbi$pastel}} - create pastel color scales 
\item \code{\link[pbi:pbi_pca.biplot]{pbi$pca.biplot}} - improved biplot for prcomp objects 
\item \code{\link[pbi:pbi_pca.corplot]{pbi$pca.corplot}} - correlation plot for original variables and PC's 
\item \code{\link[pbi:pbi_pca.pairs]{pbi$pca.pairs}} - improved pairs plot for first PC's 
\item \code{\link[pbi:pbi_pca.plot]{pbi$pca.plot}} - screeplot for prcomp objects 
\item \code{\link[pbi:pbi_pca.variances]{pbi$pca.variances}} - returns absolute variances for the PC's 
\item \code{\link[pbi:pbi_pca.varplot]{pbi$pca.varplot}} - plot variances for the PC's 
\item \code{\link[pbi:pbi_pca.toData]{pbi$toData}} - returns the original data for a prcomp object 
\item \code{\link[pbi:pbi_pcor]{pbi$pcor}} - determine partial correlation 
\item \code{\link[pbi:pbi_pcor.test]{pbi$pcor.test}} - test significance of partial correlations 
\item \code{\link[pbi:pbi_prosite2regex]{pbi$prosite2regex}} - convert a prosite regular expression into a normal regular expression 
\item \code{\link[pbi:pbi_protscale]{pbi$protscale}} - protscale plot for sliding window properties of proteins 
\item \code{\link[pbi:pbi_readFasta]{pbi$readFasta}} - read fasta file into a list of sequences 
\item \code{\link[pbi:pbi_readPepinfo]{pbi$readPepinfo}} - read data from the EMBOSS pepinfo tool 
\item \code{\link[pbi:pbi_report.chisq.test]{pbi$report.chisq.test}} - formatted text for a Chisq test 
\item \code{\link[pbi:pbi_report.conf.int]{pbi$report.conf.int}} - formatted text to report a confidence interval 
\item \code{\link[pbi:pbi_report.pval]{pbi$report.pval}} - formatting a p-value for reports 
\item \code{\link[pbi:pbi_sem]{pbi$sem}} - standard error of the mean 
\item \code{\link[pbi:pbi_searchFasta]{pbi$searchFasta}} - search a FASTA file with a regular expression 
\item \code{\link[pbi:pbi_text2fasta]{pbi$text2fasta}} - convert sequence text into FASTA formatted files 
\item \code{\link[pbi:pbi_tkregex]{pbi$tkregex}} - graphical user interface to test regular expressions 
\item \code{\link[pbi:pbi_wininstall]{pbi$wininstall}} - install a R script as a executable Batch file on windows 
\item \code{\link[pbi:pbi_wordFreq]{pbi$wordFreq}} - determine word frequencies in sequences 
\item \code{\link[pbi:pbi_xyplot]{pbi$xyplot}} - improved xy-plot with grid and correlation coefficient 
} 
} 
\examples{ 
set.seed(124) 
ls(pbi) 
}  
 
