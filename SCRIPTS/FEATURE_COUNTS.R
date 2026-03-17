compact_featureCounts <- function(sample) {
                              featureCounts(sample,
                              
                              # annotation
                              annot.ext="/work/sdz852/WORKING/RNA-seq/H5N1/JAKE/ref/genome.gtf",
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType="exon",
                              GTF.attrType="gene_id",
                              chrAliases=NULL,
                              
                              # level of summarization
                              useMetaFeatures=TRUE,
                              
                              # overlap between reads and features
                              allowMultiOverlap=FALSE,
                              minOverlap=1,
                              largestOverlap=FALSE,
                              readExtension5=0,
                              readExtension3=0,
                              read2pos=NULL,
                              
                              # multi-mapping reads
                              countMultiMappingReads=FALSE,
                              fraction=FALSE,
                              
                              # read filtering
                              minMQS=0,
                              splitOnly=FALSE,
                              nonSplitOnly=FALSE,
                              primaryOnly=TRUE,
                              ignoreDup=FALSE,
                              
                              # strandness
                              strandSpecific=0,
                              
                              # exon-exon junctions
                              juncCounts=FALSE,
                              genome=NULL,
                              
                              # parameters specific to paired end reads
                              isPairedEnd=TRUE,
                              requireBothEndsMapped=TRUE,
                              checkFragLength=FALSE,
                              minFragLength=10,
                              maxFragLength=1000,
                              countChimericFragments=TRUE,
                              autosort=TRUE,

                              # Threads
                              nthreads=80
)}
