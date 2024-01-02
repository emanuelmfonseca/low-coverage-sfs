# Import necessary libraries
import os
import subprocess
import pandas as pd
import random
from io import StringIO as SIO
import shutil

# Function to extract genotypes from VCF table
def extract_genotypes(vcf_table, position):
    genotypes = vcf_table.split('|')
    genotypes = [int(x) for x in genotypes]
    return genotypes[position]

# Set working directory and parameters
working_directory = '/Users/emanuelmfonseca/projects/LowCoverageGenomes/Analyses/Demographic_models'
model = '1D_bottlegrowth_nu1_0.4_nu2_1_T_0.4'
coverages = [1, 3, 5, 10, 30]
replicates = 10
ns = 20
n_snsps = 20000
inbreed_coefs = [0.1, 0.5, 0.9] 

# Loop through inbreeding coefficients
for inbreed_coef in inbreed_coefs:
    model_ = f'{model}_F_{inbreed_coef}'
    set_wc = os.path.join(working_directory, model_, 'heterogeneous_coverage', f'coverage_{coverages[0]}')
    
    # Create directories if they do not exist
    os.makedirs(set_wc, exist_ok=True)
    
    for replicate in range(replicates):
        save_vcf_file = os.path.join(set_wc, 'VCF_SLiM')
        os.makedirs(save_vcf_file, exist_ok=True)
        
        save_snp_matrix = os.path.join(set_wc, 'Genomic_files', f'Replicate_{replicate+1}', 'Reference_genome')
        os.makedirs(save_snp_matrix, exist_ok=True)
        
        os.chdir(save_vcf_file)
        
        s = 0
        i = 0
        
        # Generate VCF files using SLiM
        while i < n_snsps:
            silent = subprocess.run(f'slim -d Fis={inbreed_coef} -d ns={ns} -d i={replicate+1} /Users/emanuelmfonseca/projects/LowCoverageGenomes/Scripts/Analyses_scripts/Bottlegrowth.slim', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            
            if i == 0:
                vcf_path = os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}_prov.vcf")
                fid = os.path.join(os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}.vcf"))
                fid = open(fid, 'w')
                
                # Process the VCF header
                with open(vcf_path, 'r') as file:
                    for line in file:
                        if line.startswith('#'):
                            if not line.startswith('##'):
                                samples_name = ['tsk_' + str(i + 1) for i in range(len(line.split('\t')[9:]))]
                                line = '\t'.join((line.split('\t')[:9] + samples_name + ['\n']))
                            silent = fid.write(line)
                fid.close()
                
            vcf_path = os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}_prov.vcf")
            with open(vcf_path, 'r') as file:
                lines = [line for line in file if not line.startswith('#')]
            
            # Skip empty files
            if len(lines) == 0:
                continue
            else:
                random.seed(i)
                random_line = random.randint(1, len(lines)) - 1
                line = lines[random_line]
                
                s += 500
                line = '\t'.join([line.split('\t')[0], str(s), *line.split('\t')[2:]])
                
                fid = os.path.join(os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}.vcf"))
                fid = open(fid, 'a')
                silent = fid.write(line)
                fid.close()
                
                i += 1
            
        os.remove(os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}_prov.vcf"))
            
        # Read VCF file and create SNP matrix
        vcf_path = os.path.join(set_wc, 'VCF_SLiM', f"simulation_SLiM_bottlegrowth_rep_F_{inbreed_coef}_replicate{replicate + 1}.vcf")
        with open(vcf_path, 'r') as f:
            lines = [l for l in f if not l.startswith('#')]
            
            vcf_data = pd.read_csv(
                SIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                    'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t', header=None,
            ).rename(columns={'#CHROM': 'CHROM'})
        
        positions = vcf_data.iloc[:,1]
        genotypes0, genotypes1 = vcf_data.iloc[:,9:].applymap(lambda x: extract_genotypes(x, 0)), vcf_data.iloc[:,9:].applymap(lambda x: extract_genotypes(x, 1))
        
        # Concatenate genotypes to create SNP matrix
        for ncol in range(len(genotypes0.columns)):
            if ncol == 0:
                genotypes = pd.concat([genotypes0.iloc[:,ncol], genotypes1.iloc[:,ncol]], axis=1)
            else:
                genotypes = pd.concat([genotypes, pd.concat([genotypes0.iloc[:,ncol], genotypes1.iloc[:,ncol]], axis=1)], axis=1)
        
        snp_matrix = pd.concat([positions, genotypes], axis=1)
        
        # Save SNP matrix to a text file
        snp_matrix.to_csv(os.path.join(save_snp_matrix, f'SNPs_reference_genome_replicate{replicate + 1}.txt'), sep='\t', header=False, index=False)
    
    # Write population file
    nind = 0
    with open(os.path.join(set_wc, 'VCF_SLiM', 'popfile.txt'), "w") as popfile:
        for pop, i in enumerate([ns]):
            for ii in range(i):
                silent = popfile.write('tsk_' + str(nind + 1) + ' ' + 'pop' + str(pop + 1) + '\n')
                nind += 1
    
    # Copy results for different coverages
    coverages_ = coverages[1:]
    src_dir = os.path.join(working_directory, model_, 'heterogeneous_coverage', f'coverage_{coverages[0]}')
    for coverage_ in coverages_:
        dst_dir = os.path.join(working_directory, model_, 'heterogeneous_coverage', f'coverage_{coverage_}')
        if os.path.exists(dst_dir):
            shutil.rmtree(dst_dir)
        shutil.copytree(src_dir,dst_dir)
