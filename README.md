Analysis code for AIL B6xBFMI

Introduction:

  To identify additional causal gene variants contributing to variance in obesity–associated traits, an additional advanced intercross       line was generate using a different substrain of the BFMI (BFMI861-S1) which is also carrier of the obesity locus on chromosome 3 and in   addition possesses the highest liver triglycerides compared to the other BFMI lines.   
  The identification or the assessment of the importance of candidate genes and their variants that contribute to the differences in body   weight, liver weight and fatty liver will allow us in the future to identify individuals at risk of specific pathological phenotypes       such as obesity, insulin resistance, type 2 diabetes and fatty liver.

Pipeline:

-QTL mapping:
  General linear models were fitted to the trait values including subfamily, litter size and sex as fixed effects. Different statistical     models were used to do the mapping of all phenotypes (y) onto the different genotype models: Body weight: y = sex + mother + marker       genotype; Liver weight: y = sex + mother + marker genotype; Liver triglycerides: y = marker genotype. P-values were corrected for         multiple testing using a Bonferroni correction. Correction for multiple testing was performed using the number of informative SNPs as     the total number of tests performed. A logarithm (base 10) of odds (LOD) score after Bonferroni correction above 4.5 was deemed to be     ‘genome-wide highly significant’ and above 4 was deemed ‘genome-wide significant’ when supported by multiple adjacent SNP markers. QTL     regions were defined by a conservative 1.5 LOD drop from the top marker; region start and end positions are defined by the first markers   upstream and downstream of the top position that have a LOD score 1.5 LOD lower than the top marker. 
  Looking for additional effects outside the major effect QTL (jObes1) located on chromosome 3, we used a multiple QTL mapping approach.11   We adjust our single QTL model to compensate for the effect of the jObes1 locus by including the top marker from the chromosome 3 region   (SNP UNC5048297) as an additional cofactor into the model: Y = sex + mother + UNC5048297 + marker genotype + error.

-Candidate genes prioritization:
  For prioritization, first all coding genes in the identified regions were downloaded using bioMART (R package). As a next step, all SNPs   within these genes were listed using high coverage (BFMI860-12) to medium coverage sequencing (BFMI861-S1 and BFMI861-S2) data of three substrains of the BFMI compared to the reference genome BCFtools. The impact of these SNPs was checked using Ensembl Variant Effect Predictor and candidate genes were ranked according to their impact on protein structure (from highly deleterious to tolerated). 
