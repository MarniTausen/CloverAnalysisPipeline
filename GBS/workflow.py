from gwf import Workflow

gwf = Workflow()

## Genotyping

## SNPS
gwf.target("Genotyper", outputs = ['Clover_SP_dbSNP_V1.1.final.vcf'],
           cores=2, memory='32g', walltime="240:00:00", account="NChain") << """
./Unifiedgenotyper.sh
"""

## ALL_SITES
gwf.target("Genotyper_all", outputs = ['Clover_SP_dbSNP_V1.1.final.all.vcf'],
           cores=2, memory='32g', walltime="240:00:00", account="NChain") << """
./Unifiedgenotyperall.sh
"""

## SNPS
gwf.target("SelectVariants",
           inputs=['Clover_SP_dbSNP_V1.1.final.vcf'],
           outputs = ['Clover_SP_dbSNP_V1.1.clean.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./Selectvariants.sh
"""

gwf.target("SelectVariantsMAF5",
           inputs=['Clover_SP_dbSNP_V1.1.clean.vcf'],
           outputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./SelectvariantsMAF5.sh
"""

## ALL_SITES
gwf.target("SelectVariants_all",
           inputs=['Clover_SP_dbSNP_V1.1.final.all.vcf'],
           outputs = ['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./Selectvariantsall.sh
"""

gwf.target("SelectVariantsMAF5_all",
           inputs=['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           outputs=['Clover_SP_dbSNP_V1.1.MAF.all.vcf'],
           cores=1, memory='32g', walltime="08:00:00", account="NChain") << """
./SelectvariantsMAF5all.sh
"""

## SNP density

gwf.target("SNPdensity", inputs=['Clover_SP_dbSNP_V1.1.clean.all.vcf'],
           outputs = ['snp_density.csv'],
           cores=1, memory='90g', walltime="12:00:00", account="NChain") << """
./SNPdensity.sh 100
"""

## Site Frequency Spectrum

gwf.target("SFS", inputs=['Clover_SP_dbSNP_V1.1.clean.vcf'],
           outputs = ['sfs.csv'],
           cores=1, memory='8g', walltime="02:00:00", account="NChain") << """
./SFS.sh
"""

## Linkage disequilibrium

gwf.target("LD_vcftools", inputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           outputs = ['LD_10k_window.geno.ld'],
           cores=1, memory='4g', walltime="04:00:00", account="NChain") << """
./LD.sh
"""

gwf.target("LD_figures", inputs=['LD_10k_window.geno.ld'],
           outputs=['chr1.png','chr2.png','chr3.png','chr4.png','chr5.png',
                    'chr6.png','chr7.png','chr8.png','chr9.png','chr10.png',
                    'chr11.png','chr12.png','chr13.png','chr14.png','chr15.png',
                    'chr16.png'],
           cores=1, memory='4g', walltime="04:00:00", account="NChain") << """
./LDR.sh
"""

## Genetic Relationship Matrix

gwf.target("GRM", inputs=['Clover_SP_dbSNP_V1.1.MAF.vcf'],
           outputs=['GRM.csv'], memory='2g', account="NChain") << """
./GRM.sh
"""
