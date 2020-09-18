#include "samfaidx.h"

#define DEFAULT_FASTA_LINE_LEN 60

static unsigned char comp_base[256] = {
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
 32, '!', '"', '#', '$', '%', '&', '\'','(', ')', '*', '+', ',', '-', '.', '/',
'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
'@', 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', '[', '\\',']', '^', '_',
'`', 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', '{', '|', '}', '~', 127,
128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

static void reverse_complement(char *str, const hts_pos_t len) {
    char c;
    hts_pos_t i = 0, j = len - 1;

    while (i <= j) {
        c = str[i];
        str[i] = comp_base[(unsigned char)str[j]];
        str[j] = comp_base[(unsigned char)c];
        i++;
        j--;
    }
}

static void reverse(char *str, const hts_pos_t len) {
    char c;
    hts_pos_t i = 0, j = len - 1;

    while (i < j) {
        c = str[i];
        str[i] = str[j];
        str[j] = c;
        i++;
        j--;
    }
}


static int write_line(faidx_t *faid, FILE *file, const char *line, const char *name,
                      const int ignore, const int length, const hts_pos_t seq_len) {
    int id;
    hts_pos_t beg, end;

    if (seq_len < 0) {
        fprintf(stderr, "[faidx] Failed to fetch sequence in %s\n", name);

        if (ignore && seq_len == -2) {
            return EXIT_SUCCESS;
        } else {
            return EXIT_FAILURE;
        }
    } else if (seq_len == 0) {
        fprintf(stderr, "[faidx] Zero length sequence: %s\n", name);
    } else if (fai_parse_region(faid, name, &id, &beg, &end, 0)
               && (end < INT_MAX) && (seq_len != end - beg)) {
        fprintf(stderr, "[faidx] Truncated sequence: %s\n", name);
    }

    hts_pos_t i, seq_sz = seq_len;

    for (i = 0; i < seq_sz; i += length)
    {
        hts_pos_t len = i + length < seq_sz ? length : seq_sz - i;
        if (fwrite(line + i, 1, len, file) < len ||
            fputc('\n', file) == EOF) {
            // print_error_errno("faidx", "failed to write output");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static int write_output(faidx_t *faid, FILE *file, const char *name, const int ignore,
                        const int length, const int rev,
                        const char *pos_strand_name, const char *neg_strand_name,
                        enum fai_format_options format) {
    hts_pos_t seq_len;
    char *seq = fai_fetch64(faid, name, &seq_len);

    if (format == FAI_FASTA) {
        fprintf(file, ">%s%s\n", name, rev ? neg_strand_name : pos_strand_name);
    } else {
        fprintf(file, "@%s%s\n", name, rev ? neg_strand_name : pos_strand_name);
    }

    if (rev && seq_len > 0) {
        reverse_complement(seq, seq_len);
    }

    if (write_line(faid, file, seq, name, ignore, length, seq_len)
        == EXIT_FAILURE) {
        free(seq);
        return EXIT_FAILURE;
    }

    free(seq);

    if (format == FAI_FASTQ) {
        fprintf(file, "+\n");

        char *qual = fai_fetchqual64(faid, name, &seq_len);

        if (rev && seq_len > 0) {
            reverse(qual, seq_len);
        }

        if (write_line(faid, file, qual, name, ignore, length, seq_len)
            == EXIT_FAILURE) {
            free(qual);
            return EXIT_FAILURE;
        }

        free(qual);
    }

    return EXIT_SUCCESS;
}


static int read_regions_from_file(faidx_t *faid, hFILE *in_file, FILE *file, const int ignore,
                                  const int length, const int rev,
                                  const char *pos_strand_name,
                                  const char *neg_strand_name,
                                  enum fai_format_options format) {
    kstring_t line = {0, 0, NULL};
    int ret = EXIT_FAILURE;

    while (line.l = 0, kgetline(&line, (kgets_func *)hgets, in_file) >= 0) {
        if ((ret = write_output(faid, file, line.s, ignore, length, rev, pos_strand_name, neg_strand_name, format)) == EXIT_FAILURE) {
            break;
        }
    }

    free(line.s);

    return ret;
}

SamFaidx::SamFaidx(/* args */)
{
}

SamFaidx::~SamFaidx()
{
}
int SamFaidx::SamFaidxCommand(const char *faFileName, const char *regions, const char* output_file)
{
    int c, ignore_error = 0, rev = 0;
    int line_len = DEFAULT_FASTA_LINE_LEN ;/* fasta line len */
    char *region_file = NULL; // list of regions from file, one per line
    char *pos_strand_name = ""; // Extension to add to name for +ve strand
    char *neg_strand_name = "/rc"; // Extension to add to name for -ve strand
    char *strand_names = NULL; // Used for custom strand annotation
    FILE* file_out = stdout;/* output stream */
    enum fai_format_options format = FAI_FASTA;

    faidx_t *fai = fai_load_format(faFileName, format);

    if ( !fai ) {
        fprintf(stderr, "[faidx] Could not load fai index of %s\n", faFileName);
        return EXIT_FAILURE;
    }

    /** output file provided by user */
    if( output_file != NULL ) {
        file_out = fopen( output_file, "w" );

        if( file_out == NULL) {
            fprintf(stderr,"[faidx] Cannot open \"%s\" for writing :%s.\n", output_file, strerror(errno) );
            return EXIT_FAILURE;
        }
    }

    int exit_status = EXIT_SUCCESS;

    exit_status = write_output(fai, file_out, regions, ignore_error, line_len, rev, pos_strand_name, neg_strand_name, format);

    fai_destroy(fai);

    if (fflush(file_out) == EOF) {
        fprintf(stderr, "[faidx] failed to flush output\n");
        exit_status = EXIT_FAILURE;
    }

    if( output_file != NULL) fclose(file_out);
    free(strand_names);

    return exit_status;
}
