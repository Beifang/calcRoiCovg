/// Date: 05.14.2013 
/// Author: Beifang Niu
/// Modified the code to allow user to customize categories of mutation rates. 
/// E.g. TpC or TpT or any XpX 
//
//
/// Author: Cyriac Kandoth
/// Date: 01.19.2011
/// Description: Counts bases with sufficient read-depth in regions of interest within two BAMs
/// Notes:
/// - If ROIs of the same gene overlap, they will not be merged. Use BEDtools' mergeBed if needed
/// - The totals written at the end count each base only once, even if it is in multiple ROIs
//

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>

#include "sam.h"
#include "faidx.h"
#include "khash.h"

// Set bp class container
// 
// XpY
// XRpYZ
// XRNpYZM
//
// until 5 spectrum
// 0.1 version only consider XpY
// but still set large memory  
// in advance
//
// 11 = 2*5 + 1 
// avoid memory problems
//
#define MAX_BP_CLASS_TYPES_STRING_LEN 500
#define MAX_BP_CLASS_TYPES 100
#define MAX_BP_CLASS 11

KHASH_MAP_INIT_STR(s, int)

// Initializes the header hash in 
// a given bam_header_t.
//
// Defined in samtools/bam_aux.c

void bam_init_header_hash(bam_header_t *header);

// Create some values that represent the different 
// refseq bp-classes

enum bp_class_t { AT, CG, CpG, IUB, UNKNOWN };

typedef struct
{
    // The start and stop of a region of interest
    uint32_t beg;
    uint32_t end;
    
    // Minimum mapping quality of the reads to pileup
    int min_mapq; 
   
    // Minimum read depth required per bam 
    int min_depth_bam1;
    int min_depth_bam2; 

    // Tags bases in a region of bam1 with the 
    // minimum required read-depth
    bool *bam1_cvg; 
    
    // Counts bases in a region that has the 
    // minimum read depth in both bams
    uint32_t covd_bases; 

    // Counts covered bases in an ROI of 4 bp-classes 
    // AT, CG, CpG, IUB
    //uint32_t base_cnt[4]; 
    //
    // extend to 100 class tyes
    uint32_t base_cnt[MAX_BP_CLASS_TYPES];

    // Counts bases in all ROIs that have the 
    // minimum read depth in both bams
    uint32_t tot_covd_bases;

    // Counts covered bases in all ROIs of 4 bp-classes 
    // AT, CG, CpG, IUB
    //uint32_t tot_base_cnt[4];
    //
    // extend to 100 class tyes 
    uint32_t tot_base_cnt[MAX_BP_CLASS_TYPES];

    // Contains the reference sequence for the entire 
    // chromosome a region lies in
    char *ref_seq;

    // This prevents counting the same base twice when ROIs 
    // overlap in a chromosome
    char *bp_class; 

    // user customize bp class types
    char *bp_class_types;
    char **bp_class_container;

    uint8_t bp_class_number;
    uint8_t bp_class_lengths[MAX_BP_CLASS_TYPES];
    uint8_t unknown;
    uint8_t iub;

    // A chromosome's ID in the the BAM header hash, 
    // and it's length
    int ref_id;
    int ref_len; 
    
    // The two bam files that need to be piled-up
    samfile_t *sam1;
    samfile_t *sam2; 

} pileup_data_t;

// get class type 
static bool getClass( char pre, 
                      char mid, 
                      char pos, 
                      char *container, uint8_t length)
{
    // AT or CG
    if (length == 2)
    {
        if ( (mid == container[0]) || (mid == container[1]) )
        {
            return true;
        }
    
    }
    else
    {
        // CpG or XpY
        // will extend to XXXXXpYYYYY
        if ( ((mid == container[0]) && (pos == container[2])) 
             ||
             ((mid == container[2]) && (pre == container[0]))
           )
        {
           return true;
        }

    }

    return false;

}


// Callback for bam_fetch() pushed only alignments 
// that pass the minimum mapping quality
//
static int fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf = (bam_plbuf_t*)data;
    
    // Invokes the callback function specified 
    // in bam_plbuf_init()

    bam_plbuf_push(b, buf); 
    
    return 0;
}

// Callback for bam_plbuf_init() when running a pileup on bam1
static int pileup_func_1(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    pileup_data_t *tmp = (pileup_data_t*)data;
    
    if (pos >= tmp->beg && pos < tmp->end)
    {
        // Count the number of reads that pass the mapping 
        // quality threshold across this base
        int i;
        int mapq_n = 0;
        
        for (i = 0; i < n; ++i)
        {
            const bam_pileup1_t *base = pl + i;
            if (!base->is_del && base->b->core.qual >= tmp->min_mapq)
            {
                mapq_n++;
            }
        }
        
        tmp->bam1_cvg[pos - tmp->beg] = (mapq_n >= tmp->min_depth_bam1);
    }
    
    return 0;

}

// Callback for bam_plbuf_init() when running a pileup on bam2
static int pileup_func_2(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    pileup_data_t *tmp = (pileup_data_t*)data;

    if (pos >= tmp->beg && pos < tmp->end && tmp->bam1_cvg[pos - tmp->beg])
    {
        // Count the number of reads that pass the mapping 
        // quality threshold across this base
        //
        int i;
        int mapq_n = 0;
        
        uint8_t j;

        for (i = 0; i <n; ++i)
        {
            const bam_pileup1_t *base = pl + i;
            if (!base->is_del && base->b->core.qual >= tmp->min_mapq)
            {
                mapq_n++;
            }

        }

        if (mapq_n >= tmp->min_depth_bam2)
        {

            int class = (int)(tmp->bp_class[pos]);

            // default IUB 
            bool notiub = false;
            bool isHit  = false;

            if (class == tmp->unknown)
            {

                char base = toupper(tmp->ref_seq[pos]);
                char prev_base = toupper(tmp->ref_seq[pos-1]);
                char next_base = toupper(tmp->ref_seq[pos+1]);

                for (j=0; j<tmp->iub; ++j)
                {
                    isHit = getClass( prev_base, 
                                           base, 
                                      next_base, 
                                      tmp->bp_class_container[j], 
                                      tmp->bp_class_lengths[j] );

                    if (isHit)
                    {
                        notiub = true;
                        class = j;
                        break;
                    }

                }

                if (!notiub)
                {
                    class = tmp->iub;
                }

                ++tmp->covd_bases;
                ++tmp->base_cnt[class];
                ++tmp->tot_covd_bases;
                ++tmp->tot_base_cnt[class];

                // Tag this as seen and save its class 
                // for an overlapping ROI
                //
                tmp->bp_class[pos] = (char)class;

            }
            else
            {
                ++tmp->covd_bases;
                ++tmp->base_cnt[class];
            }

        }
    }

    return 0;

}

// Separate string based on given 
// char delimiter
//
int separateString(char *string, char c, char **container)
{
    int z = 0;
    int x, y = 0;
    int tokenNum = 0;
    int length = strlen(string);

    char buff1[MAX_BP_CLASS_TYPES_STRING_LEN];
    char array[MAX_BP_CLASS_TYPES_STRING_LEN];

    for (x = 0; x < length; x++)
    {
        if (string[x + 1] == '\0')
        {
            buff1[z] =  string[x];
            buff1[z + 1] =  '\0';

            if (strlen(buff1))
            {             
                strcpy(array, buff1);
                strcpy(container[y], buff1);
                //fprintf(stderr, "%s\n", array);
                tokenNum++;
            }
        }
        else if (string[x] != c)
        {
            buff1[z] =  string[x];
            z++;

            if (string[x + 1] == c)
            {
                buff1[z] = '\0';
                z = 0;

                if (strlen(buff1))
                {
                    strcpy(array, buff1);
                    strcpy(container[y], buff1);
                    y++;
                    //fprintf(stderr, "%s\n", array);
                    tokenNum++;
                }

            }
        }
    }

    return(tokenNum);

}

