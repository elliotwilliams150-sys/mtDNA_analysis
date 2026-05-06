# Print start time 
echo "Start time:" >> test_output.txt
date >> test_output.txt

# actual thing you want to run
bwa-mem2 index /mnt/claw-raid/elliot/genomes/GRCh38/GRCh38.primary_assembly.genome.fa 

# Print end time 
echo "End time:" >> test_output.txt
date >> test_output.txt
