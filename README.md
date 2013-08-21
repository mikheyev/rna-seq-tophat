##Using tophat2 to map transcripts


1.  move.sh
    - Move data from SQC into the working directory, and merges files
2. tophat1.sh
   - Map using tophat2, then extract unmapped reads, and convert them back to fastq.
   - This step is parallelized on the cluser.
   - Output in rnaseq   
3. merge_junctions.sh
   - Merge junctions from first tophat run
4. tophat2.sh
   - Try re-mapping the unmapped reads
   - Output in rnaseq2
5. tophat2.sh
   - Try re-mapping the unmapped reads 
   - Output into rnaseq2
   - extract reads still unmapped into *unmapped*
   - This step is likewise parallelized.
6. trinity.sh
   - assemble unmpped reads using trinity
7. merge_tophat.sh
   - merge all of the tophat assemblies
8. cufflinks.sh
   - find transcripts
9. merge_assembly.sh
   - combine trinity and cufflinks assemblies, add ERCC controls
10. rsem.sh
   - map libraries to assemblies to get counts
