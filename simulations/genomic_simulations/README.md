# What does this folder contain?

This directory includes four scripts, with one serving as the primary script ('genomic_script_generator.R'). Additionally, there are three secondary scripts that are invoked by the main script: 'generate_genomes.R', 'simulate_reads_and_processing_genomes.R', and 'variant_calling_multisample.R'. The primary script generates R and Slurm scripts that will be used for execution in a High-Performance Computing (HPC) environment. To customize the script creation process, users are required to update information such as software paths, Slurm credentials, and relevant paths.

For any questions, please contact Emanuel M. Fonseca (emanuelmfonseca@arizona.email.edu or emanuelmfonseca@gmail.com).