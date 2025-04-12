This is the code used to generate the variant peptide database of QuaVaProt

libraries_and_functions.R contains all functions and libraries used in the pipeline.

NCI_GDC.R and COSMIC.R to were used to process their respective mutation data into peptide variants, with evaluations.

database_builder_gdc_cosmic.R was used to process the output tables of the above R files into a .sqlite file, for use in the QuaVaProt app
