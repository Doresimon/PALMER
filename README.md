Pre-mAsking Long reads for Mobile Element inseRtion 
<<<<<<< HEAD
PALMER is used to detect non-reference MEI events within the masked sequence data. It uses the reference-aligned BAM files from long-read technology as inputs. 
1. It delineates large genome data into small bins (100kb bins), which will be piped into multi-threads and processed separately, to increase the efficiency of program. 
2. The known repeats (L1NE-1, Alus or SVAs in reference) are used to mask the portions of reads that aligned to these repeats. 
3. After obtaining the pre-masked reads, PALMER searches against the insertion sequence library (for L1Hs sequence, GenBank Accession: L19088)  by using Blastn. Then it uses the cigar information in the reads and joint the truncated segments from Blastn into one in a single read. These reads with insertion sequence are considered as supporting reads. 
4. PALMER identifies the candidate TSD motif in 50bp 5� upstream and 3kb 3� downstream of insertion sequence for each read. 
5. It then runs a module for filtering the candidate TSD motif and identifying transduction/polyA sequence. PALMER will define the supporting read with/without valid TSD motif and with/without valid transduction and polyA sequence. The ideal structure of an event having transduction sequence should be 5�-TSD-L1Hs-polyA-TransD-polyA-TSD-3�. 
6. Afterwards, PALMER will cluster all supporting reads in one loci (or supporting one event) of the genome, as well as cluster the TSD motif among all supporting reads for one event and choose the most confident TSD motif with/without transduction/polyA sequence. PALMER will obtain two numbers for one event, the number of supporting reads and the number of supporting reads with predicted TSD motif. 
7. Finally, PALMER will combine all events in each bin and output all candidate non-reference MEIs.
=======
PALMER is used to detect non-reference LINE-1 events within the masked sequence data. It uses the reference-aligned BAM files from long-read technology as inputs. 
1. It delineates large genome data into small bins (100kb bins), which will be piped into multi-threads and processed separately, to increase the efficiency of program. 
2. The known repeats (L1NE-1, Alus or SVAs in reference) are used to mask the portions of reads that aligned to these repeats. 
3. After obtaining the pre-masked reads, PALMER searches against the insertion sequence library (for L1Hs sequence, GenBank Accession: L19088)  by using Blastn. Then it uses the cigar information in the reads and joint the truncated segments from Blastn into one in a single read. These reads with insertion sequence are considered as supporting reads. 
4. PALMER identifies the candidate TSD motif in 50bp 5’ upstream and 3kb 3’ downstream of insertion sequence for each read. 
5. It then runs a module for filtering the candidate TSD motif and identifying transduction/polyA sequence. PALMER will define the supporting read with/without valid TSD motif and with/without valid transduction and polyA sequence. The ideal structure of an event having transduction sequence should be 5’-TSD-L1Hs-polyA-TransD-polyA-TSD-3’. 
6 Afterwards, PALMER will cluster all supporting reads in one loci (or supporting one event) of the genome, as well as cluster the TSD motif among all supporting reads for one event and choose the most confident TSD motif with/without transduction/polyA sequence. PALMER will obtain two numbers for one event, the number of supporting reads and the number of supporting reads with predicted TSD motif. 
7. Finally, PALMER will combine all events in each bin and output all candidate non-reference MEIs.
>>>>>>> 6ac7718ccb1106499c54c11687f591800cba1f9f
