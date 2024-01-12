# General workflow

### Inputs

- RNAseq: The read_data.R script in the parent directory must have produced the main RNAseq RDS object.
    - Input in rule `prepare_pheno`
- Genotypes: The raw filtered (but not imputed) vcf form SNP archer
    - Input in rule `beagle`

### Outputs

Main results are the genes with eqtls and gxe eqtl his in the detections folder.

# TODO:

mapping scripts `eQTLmapping/scripts/runGridLMM.R` and `eQTLmapping/scripts/runGxEGridLMM.R` require some manual input of the fixed effect design, and specially the eQTL script needs to be edited to caputure the correct column for the main eQTL effect.