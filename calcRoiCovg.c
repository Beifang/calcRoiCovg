///  Date: 05.14.2013 
///  Beifang Niu
///  Modified the code to allow user to customize categories of mutation rates. 
///  E.g. TpC or TpT or any XpX 
//
/// Author: Cyriac Kandoth
/// Date: 01.19.2011
/// Description: Counts bases with sufficient read-depth in regions of interest within two BAMs
//
/// Notes:
/// - If ROIs of the same gene overlap, they will not be merged. Use BEDtools' mergeBed if needed
/// - The totals written at the end count each base only once, even if it is in multiple ROIs
//

#define _GNU_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

#include "calcRoiCovg.h"

pileup_data_t data;

// usage infor
void usage(void)
{

    // usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Version 0.1\n");
    fprintf(stderr, "Usage: calcRoiCovg <bam1> <bam2> <roi_file> <ref_seq_fasta> <output_file>\n\n");

    fprintf(stderr, "        -q INT    filtering reads with mapping quality less than INT [%d]\n", data.min_mapq);
    fprintf(stderr, "        -n INT    minimum reads depth for bam1 [%d]\n", data.min_depth_bam1);
    fprintf(stderr, "        -t INT    minimum reads depth for bam2 [%d]\n", data.min_depth_bam2);
    fprintf(stderr, "        -c STRING bp class types, delimited by comma, default: \"AT,CG,CpG\"\n");
    
    fprintf( stderr, "\n\nROI file should be a tab-delimited list of [chrom, start, stop, annotation]" );
    fprintf( stderr, "\nwhere start and stop are both 1-based chromosomal loci. For example:" );
    fprintf( stderr, "\n\n20\t44429404\t44429608\tELMO2\nMT\t5903\t7445\tMT-CO1\n" );
    fprintf( stderr, "\nROI file *must* be sorted by chromosome/contig names\n\n" );

    exit(1);

}

// get options
int mGetOptions(int rgc, char *rgv[])
{
    int c;

    while ((c = getopt(rgc, rgv, "q:n:t:c:")) >= 0)
    {
        switch (c) {

            case 'q': data.min_mapq = atoi(optarg); break;
            case 'n': data.min_depth_bam1 = atoi(optarg); break;
            case 't': data.min_depth_bam2 = atoi(optarg); break;
            case 'c': data.bp_class_types = optarg; break;

            default: fprintf(stderr, "Unrecognized option '-%c'.\n", c); return 1;
        }

    }

    return(c);
}

int main(int argc, char *argv[])
{

    // shared across functions
    //
    data.ref_id   = -1;
    data.ref_seq  = NULL;
    data.bp_class = NULL;

    // set default min_mapq and 
    // min_depths
    //
    data.min_mapq = 20;
    data.min_depth_bam1 = 6;
    data.min_depth_bam2 = 8;

    // bp class types
    //
    data.bp_class_types = (char*)malloc(MAX_BP_CLASS_TYPES_STRING_LEN);

    // default bp class types
    // change default class types
    //
    //data.bp_class_types
    //= 
    //"CpG,AT,CG,IUB,UNKNOWN";
    //
    data.bp_class_types = "AT,CG,CpG";

    //fprintf(stderr, "bp class types: %s\n", data.bp_class_types);

    mGetOptions(argc, argv);

    // usage
    if ((argc-optind) == 0)
    {
        usage();
    }

    // Open both BAM files and load their index files
    data.sam1 = samopen(argv[optind], "rb", 0);
    if (!data.sam1) fprintf(stderr, "Failed to open BAM file %s\n", argv[optind]);

    bam_index_t *idx1 = bam_index_load(argv[optind]);
    if (!idx1) fprintf(stderr, "BAM index file is not available for %s\n", argv[optind]);

    data.sam2 = samopen(argv[optind+1], "rb", 0);
    if (!data.sam2) fprintf(stderr, "Failed to open BAM file %s\n", argv[optind+1]);

    bam_index_t *idx2 = bam_index_load(argv[optind+1]);
    if (!idx2) fprintf(stderr, "BAM index file is not available for %s\n", argv[optind+1]);

    // Open the file with the annotated regions of interest
    FILE *roiFp = fopen(argv[optind+2], "r");
    if (!roiFp) fprintf(stderr, "Failed to open ROI file %s\n", argv[optind+2]);

    // Load an index to the reference sequence fasta file
    faidx_t *ref_fai = fai_load( argv[optind+3] );
    if (!ref_fai) fprintf(stderr, "Failed to open reference fasta file %s\n", argv[optind+3]);

    // Open the output file to write to
    FILE* outFp = fopen( argv[optind+4], "w" );
    if (!outFp) fprintf(stderr, "Failed to open output file %s\n", argv[optind+4]);

    // Show the user any and all errors they need to 
    // fix above before quitting the program
    if (   !data.sam1 
        || !idx1 
        || !data.sam2 
        || !idx2 
        || !roiFp 
        || !ref_fai 
        || !outFp )
        
        return 1;

    // seperate class type string 
    // instead of using fixed enum()
  
    // Memmory allocation
    data.bp_class_container = (char **)malloc((MAX_BP_CLASS_TYPES+1)*sizeof(char *));

    uint8_t i,j;

    for (i=0; i<=MAX_BP_CLASS_TYPES; i++)
    {
        data.bp_class_container[i] = (char *)malloc((MAX_BP_CLASS)*sizeof(char));
    }

    data.bp_class_number = separateString(data.bp_class_types, ',', data.bp_class_container);

    // IUB && UNKNOWN
    // 
    data.iub     = data.bp_class_number;
    data.unknown = data.iub + 1;

    // init class type to index hash 
    // using Khash 
    // actually don't have to 
    // but maybe useful later
    // when function extension
    //
    khiter_t kIter;
    khash_t(s) *classTypeMap = kh_init(s);

    int ret;
    for (i=0; i< data.bp_class_number; i++)
    {
        kIter = kh_put(s, classTypeMap, data.bp_class_container[i], &ret);
        kh_value(classTypeMap, kIter) = i;
    }
    // clean up
    kh_destroy(s, classTypeMap);

    // Write a header with column titles 
    // for the output file
    fprintf( outFp, "#NOTE: Last line in file shows non-overlapping totals across all ROIs\n" );
    fprintf( outFp, "#Gene\tROI\tLength\tCovered\t" );

    // write output file header
    //
    for (i=0; i< (data.bp_class_number - 1); i++)
    {
        fprintf( outFp, "%ss_Covered\t", data.bp_class_container[i] );
    }

    fprintf( outFp, "%ss_Covered\n", data.bp_class_container[data.bp_class_number - 1] );

    // bp class lengths
    // &&
    // change to upper 
    for (i=0; i<data.bp_class_number; i++)
    {
        data.bp_class_lengths[i] = strlen(data.bp_class_container[i]);

        for (j=0; j<data.bp_class_lengths[i]; j++)
        {
            data.bp_class_container[i][j] = toupper(data.bp_class_container[i][j]);
        }
    }

    // Initialize a header hash to check for valid ref_names in bam1
    bam_init_header_hash(data.sam1->header);
    khash_t(s) *hdr_hash = (khash_t(s)*)data.sam1->header->hash;
    
    // Initialize the counters for the total number of 
    // non-overlapping bases in all ROIs

    data.tot_covd_bases = 0;

    for (i=0; i<=data.bp_class_number; i++)
    {
        data.tot_base_cnt[i] = 0;
    }

    //data.tot_covd_bases = data.tot_base_cnt[AT] 
    //                    = data.tot_base_cnt[CG] 
    //                    = data.tot_base_cnt[CpG] = 0;

    size_t length;

    char ref_name[50];
    char gene_name[100];

    char *line = NULL;
    
    line = (char*)malloc(200);
    
    while (getline(&line, &length, roiFp) != -1)
    {

        if (
                sscanf( line, "%s %lu %lu %s", ref_name, 
                    (unsigned long *)&data.beg, 
                    (unsigned long *)&data.end, gene_name ) == 4
           )
        {
            int ref_id;

            // If this region is valid in bam1,
            // we'll assume it's also valid in bam2
            
            khiter_t iter = kh_get(s, hdr_hash, ref_name);

            if ( 
                 iter == kh_end(hdr_hash) 
                 || 
                 data.beg > data.end 
               )
            {
                fprintf(stderr, "Skipping invalid ROI: %s", line);
            }
            else
            {
                // Make the start locus a 
                // 0-based coordinate
                //
                --data.beg; 

                ref_id = kh_value( hdr_hash, iter );
                
                uint32_t bases = data.end - data.beg;

                // data.covd_bases = data.base_cnt[AT] 
                //                 = data.base_cnt[CG] 
                //                 = data.base_cnt[CpG] = 0;

                // change inition enum part
                // to the following loop
                //
                data.covd_bases = 0;

                for (i=0; i<data.bp_class_number; i++)
                {
                    data.base_cnt[i] = 0 ;
                }
                
                // calloc also sets them to zero
                data.bam1_cvg = (bool*)calloc( bases, sizeof( bool )); 

                // Load this whole chromosome's refseq unless 
                // already loaded for the previous ROI
                //
                if (data.ref_seq == NULL || ref_id != data.ref_id)
                {
                    if (data.ref_seq)  free(data.ref_seq);
                    if (data.bp_class) free(data.bp_class);

                    data.ref_seq = fai_fetch(ref_fai, data.sam1->header->target_name[ref_id], &data.ref_len);
                    data.bp_class = (char*)malloc( data.ref_len * sizeof( char ));

                    //memset(data.bp_class, UNKNOWN, data.ref_len);
                    //set all UNKNOWN
                    //
                    memset(data.bp_class, data.unknown, data.ref_len);

                    data.ref_id = ref_id;
                }

                // If the ROI is at a chromosome tip, edit it so 
                // we can look for CpGs without segfaulting
                //
                if (data.beg == 0) ++data.beg;
                if (data.end == data.ref_len) --data.end;

                // Pileup bam1 and tag all the bases which 
                // have sufficient read depth
                 
                // Initialize pileup 
                bam_plbuf_t *buf1 = bam_plbuf_init(pileup_func_1, &data); 

                bam_fetch(data.sam1->x.bam, idx1, ref_id, data.beg, data.end, buf1, fetch_func);

                bam_plbuf_push(0, buf1);
                bam_plbuf_destroy(buf1);

                // Pileup bam2 and count bases with sufficient 
                // read depth, and tagged earlier in bam1

                // Initialize pileup
                bam_plbuf_t *buf2 = bam_plbuf_init(pileup_func_2, &data); 
                bam_fetch(data.sam2->x.bam, idx2, ref_id, data.beg, data.end, buf2, fetch_func);

                bam_plbuf_push(0, buf2);
                bam_plbuf_destroy(buf2);
                
                //fprintf( outFp, "%s\t%s:%lu-%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
                //        gene_name,
                //        ref_name,
                //        (unsigned long)data.beg+1, 
                //        (unsigned long)data.end, 
                //        (unsigned long)bases,
                //        (unsigned long)data.covd_bases, 
                //        (unsigned long)data.base_cnt[AT],
                //        (unsigned long)data.base_cnt[CG], 
                //        (unsigned long)data.base_cnt[CpG] );
                //

                fprintf(outFp, "%s\t%s:%lu-%lu\t%lu\t%lu\t", gene_name, ref_name,
                        (unsigned long)data.beg+1, 
                        (unsigned long)data.end, 
                        (unsigned long)bases,
                        (unsigned long)data.covd_bases);

                for (j=0; j<(data.bp_class_number - 1); j++)         
                {
                    fprintf(outFp, "%lu\t", (unsigned long)data.base_cnt[j]);
                }

                fprintf(outFp, "%lu\n", (unsigned long)data.base_cnt[data.bp_class_number - 1]);

                free(data.bam1_cvg);

            }
        }
        else
        {

            fprintf(stderr, "Badly formatted ROI: %s", line);
            fprintf(stderr, "\nROI file should be a tab-delimited list of [chrom, start, stop, annotation]");
            fprintf(stderr, "\nwhere start and stop are both 1-based chromosomal loci");
            fprintf(stderr, "\nFor example:\n20\t44429404\t44429608\tELMO2\nMT\t5903\t7445\tMT-CO1\n");
            fprintf(stderr, "\nNOTE: ROI file *must* be sorted by chromosome/contig names\n\n");
            
            return 1;
        }

    }
    // The final line in the file contains the 
    // non-overlapping base counts across all ROIs

    //fprintf( outFp, "#NonOverlappingTotals\t\t\t%lu\t%lu\t%lu\t%lu\n",
    //        (unsigned long)data.tot_covd_bases,
    //        (unsigned long)data.tot_base_cnt[AT],
    //        (unsigned long)data.tot_base_cnt[CG],
    //        (unsigned long)data.tot_base_cnt[CpG] );

    fprintf(outFp, "#NonOverlappingTotals\t\t\t%lu\t", (unsigned long)data.tot_covd_bases);
    
    for (j=0; j<(data.bp_class_number - 1); j++)
    {
        fprintf(outFp, "%lu\t", (unsigned long)data.tot_base_cnt[j]);
    }

    fprintf(outFp, "%lu\n", (unsigned long)data.tot_base_cnt[data.bp_class_number - 1]);

    // Cleanup
    //

    // bp_class container 
    for (i=0; i<=MAX_BP_CLASS_TYPES; i++)
    {
        if (data.bp_class_container[i])
        {
           free(data.bp_class_container[i]); 
        }
    }

    if (data.bp_class_container)
    {
        free(data.bp_class_container);
    }

    if (line) free(line);

    if (data.ref_seq) free(data.ref_seq);
    if (data.bp_class) free(data.bp_class);

    bam_index_destroy( idx1 );
    bam_index_destroy( idx2 );
    
    samclose( data.sam1 );
    samclose( data.sam2 );
    
    fai_destroy( ref_fai );
    
    fclose( roiFp );
    fclose( outFp );

    return 0;

}

