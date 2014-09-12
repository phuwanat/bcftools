/*  reheader.c -- reheader subcommand.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include <htslib/tbx.h>
#include "bcftools.h"
#include "khash_str2str.h"

typedef struct _args_t
{
    htsFile *fp;
    char **argv, *fname, *samples_fname, *header_fname, *contigs_fname;
    int argc, file_type;
}
args_t;

static void read_header_file(char *fname, kstring_t *hdr)
{
    kstring_t tmp = {0,0,0};
    hdr->l = 0;

    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
    {
        kputsn(tmp.s,tmp.l,hdr);
        kputc('\n',hdr);
    }
    if ( hts_close(fp) ) error("Close failed: %s\n", fname);
    free(tmp.s);

    while ( hdr->l>0 && isspace(hdr->s[hdr->l-1]) ) hdr->l--;  // remove trailing newlines
    kputc('\n',hdr);
}

static int set_sample_pairs(char **samples, int nsamples, kstring_t *hdr, int idx)
{
    int i, j, n;

    // Are these samples "old-name new-name" pairs?
    void *hash = khash_str2str_init();
    for (i=0; i<nsamples; i++)
    {
        char *key, *value;
        key = value = samples[i];
        while ( *value && !isspace(*value) ) value++;
        if ( !*value ) break;
        *value = 0; value++;
        while ( isspace(*value) ) value++;
        khash_str2str_set(hash,key,value);
    }
    if ( i!=nsamples )  // not "old-name new-name" pairs
    {
        khash_str2str_destroy(hash);
        return 0;
    }

    while ( hdr->l>0 && isspace(hdr->s[hdr->l-1]) ) hdr->l--;  // remove trailing newlines
    hdr->s[hdr->l] = 0;

    kstring_t tmp = {0,0,0};
    i = j = n = 0;
    while ( hdr->s[idx+i] && hdr->s[idx+i])
    {
        if ( hdr->s[idx+i]=='\t' )
        {
            hdr->s[idx+i] = 0;

            if ( ++n>9 )
            {
                char *ori = khash_str2str_get(hash,hdr->s+idx+j);
                kputs(ori ? ori : hdr->s+idx+j, &tmp);
            }
            else
                kputs(hdr->s+idx+j, &tmp);

            kputc('\t',&tmp);

            j = ++i;
            continue;
        }
        i++;
    }
    char *ori = khash_str2str_get(hash,hdr->s+idx+j);
    kputs(ori ? ori : hdr->s+idx+j, &tmp);

    if ( hash ) khash_str2str_destroy(hash);

    hdr->l = idx;
    kputs(tmp.s, hdr);
    kputc('\n', hdr);
    free(tmp.s);

    return 1;
}

static void set_samples(char **samples, int nsamples, kstring_t *hdr)
{
    // Find the beginning of the #CHROM line
    int i = hdr->l - 2, ncols = 0;
    while ( i>=0 && hdr->s[i]!='\n' )
    {
        if ( hdr->s[i]=='\t' ) ncols++;
        i--;
    }
    if ( i<0 || strncmp(hdr->s+i+1,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",45) ) error("Could not parse the header: %s\n", hdr->s);

    // Are the samples "old-sample new-sample" pairs?
    if ( set_sample_pairs(samples,nsamples,hdr, i+1) ) return;

    // Replace all samples
    if ( ncols!=nsamples+8 )
        fprintf(stderr, "Warning: different number of samples: %d vs %d\n", nsamples,ncols-8);

    ncols = 0;
    while ( ncols!=9 )
    {
        i++;
        if ( hdr->s[i]=='\t' ) ncols++;
    }
    hdr->l = i;

    for (i=0; i<nsamples; i++)
    {
        kputc('\t', hdr);
        kputs(samples[i], hdr);
    }
    kputc('\n', hdr);
}

static void parse_fai_record(char *line, void *chr_dict, kstring_t *out)
{
    char *ss = line, *se = line;
    while ( *se && !isspace(*se) ) se++;
    char tmp = *se; *se = 0;
    if ( chr_dict && !khash_str2int_has_key(chr_dict, ss) ) return;

    kputs("##contig=<ID=", out);
    kputs(ss, out);
    kputs(",length=", out);

    *se = tmp;
    ss = se;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) error("Could not parse .fai: %s\n", line);
    se = ss;
    while ( *se && isdigit(*se) ) se++;
    if ( ss==se ) error("Could not parse .fai: %s\n", line);
    *se = 0;

    kputs(ss, out);
    kputs(">\n", out);
}

static void parse_dict_record(char *line, void *chr_dict, kstring_t *out)
{
    if ( strncmp("@SQ\t",line,4) ) return;

    char *ss, *se, tmp;
    if ( !(ss=strstr(line,"SN:")) ) error("Could not parse .dict file, the SN field not present: %s\n", line);
    ss += 3;
    se = ss;
    while ( *se && !isspace(*se) ) se++;
    tmp = *se; *se = 0;
    if ( chr_dict && !khash_str2int_has_key(chr_dict, ss) ) return;
    kputs("##contig=<ID=", out);
    kputs(ss, out);
    *se = tmp;

    if ( !(ss=strstr(line,"LN:")) ) error("Could not parse .dict file, the LN field not present: %s\n", line);
    ss += 3;
    se = ss;
    while ( *se && !isspace(*se) ) se++;
    tmp = *se; *se = 0;
    kputs(",length=", out);
    kputs(ss, out);
    *se = tmp;

    if ( (ss=strstr(line,"M5:")) )
    {
        ss += 3;
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        tmp = *se; *se = 0;
        kputs(",md5=", out);
        kputs(ss, out);
        *se = tmp;
    }
    if ( (ss=strstr(line,"AS:")) )
    {
        ss += 3;
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        tmp = *se; *se = 0;
        kputs(",assembly=", out);
        kputs(ss, out);
        *se = tmp;
    }
    if ( (ss=strstr(line,"UR:")) )
    {
        ss += 3;
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        tmp = *se; *se = 0;
        kputs(",URL=", out);
        kputs(ss, out);
        *se = tmp;
    }
    kputs(">\n", out);
}

static void set_contigs(args_t *args, kstring_t *hdr)
{
    // Remove any existing contig lines
    char *ss = hdr->s;
    while ( (ss=strstr(ss,"\n##contig=<")) )
    {
        char *se = ss+1;
        while ( *se && *se!='\n' ) se++;
        if ( !*se ) error("Could not parse the header\n");
        memmove(ss+1,se+1,hdr->l - (se - hdr->s));
        hdr->l -= se - ss;
    }

    // If indexed, we add only contigs which are actually used
    void *chr_dict = NULL;
    if ( args->file_type & HTS_GZ )
    {
        tbx_t *idx = tbx_index_load(args->fname);
        if ( idx )
        {
            int i, nseq;
            const char **seq = tbx_seqnames(idx, &nseq);
            chr_dict = khash_str2int_init();
            for (i=0; i<nseq; i++) khash_str2int_inc(chr_dict, strdup(seq[i]));
            tbx_destroy(idx);
            free(seq);
        }
    }

    kstring_t tmp = {0,0,0}, ctgs = {0,0,0};
    int is_fai = -1;
    htsFile *fp = hts_open(args->contigs_fname, "r");
    if ( !fp ) error("Could not read: %s\n", args->contigs_fname);
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
    {
        // first time here, is this a .dict or fai file?
        if ( is_fai==-1 ) is_fai = tmp.s[0]=='@' ? 0 : 1;

        if ( is_fai )
            parse_fai_record(tmp.s, chr_dict, &ctgs);
        else 
            parse_dict_record(tmp.s, chr_dict, &ctgs);
    }
    if ( hts_close(fp) ) error("Close failed: %s\n", args->contigs_fname);
    free(tmp.s);

    if ( chr_dict) khash_str2int_destroy_free(chr_dict);

    // Add the new contig lines
    int i = hdr->l - 2;
    while ( i>=0 && hdr->s[i]!='\n' ) i--;
    if ( i<0 || strncmp(hdr->s+i+1,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",38) ) error("Could not parse the header: %s\n", hdr->s);
    ks_resize(hdr, hdr->l+ctgs.l+1);
    i++;    // after the newline
    memmove(hdr->s+i+ctgs.l, hdr->s+i, hdr->l - i);
    memcpy(hdr->s+i,ctgs.s,ctgs.l);
    hdr->l += ctgs.l;
    hdr->s[hdr->l] = 0;
    free(ctgs.s);
}

BGZF *hts_get_bgzfp(htsFile *fp);
static void reheader_vcf_gz(args_t *args)
{
    BGZF *fp = hts_get_bgzfp(args->fp);
    if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length )
        error("Failed to read %s: %s\n", args->fname, strerror(errno));

    kstring_t hdr = {0,0,0};
    char *buffer = (char*) fp->uncompressed_block;

    // Read the header and find the position of the data block
    if ( buffer[0]!='#' ) error("Could not parse the header, expected '#', found '%c'\n", buffer[0]);

    int skip_until = 1;     // end of the header in the current uncompressed block
    while (1)
    {
        if ( buffer[skip_until]=='\n' )
        {
            skip_until++;
            if ( skip_until>=fp->block_length )
            {
                kputsn(buffer,skip_until,&hdr);
                if ( bgzf_read_block(fp) != 0 || !fp->block_length ) error("FIXME: No body in the file: %s\n", args->fname);
                skip_until = 0;
            }
            // The header has finished
            if ( buffer[skip_until]!='#' )
            {
                kputsn(buffer,skip_until,&hdr);
                break;
            }
        }
        skip_until++;
        if ( skip_until>=fp->block_length )
        {
            kputsn(buffer,fp->block_length,&hdr);
            if (bgzf_read_block(fp) != 0 || !fp->block_length) error("FIXME: No body in the file: %s\n", args->fname);
            skip_until = 0;
        }
    }

    int nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
        samples = hts_readlines(args->samples_fname, &nsamples);
    if ( args->header_fname )
    {
        free(hdr.s); hdr.s = NULL; hdr.l = hdr.m = 0;
        read_header_file(args->header_fname, &hdr);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &hdr);
        int i;
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }
    if ( args->contigs_fname ) set_contigs(args, &hdr);

    // Output the modified header
    BGZF *bgzf_out = bgzf_dopen(fileno(stdout), "w");
    bgzf_write(bgzf_out, hdr.s, hdr.l);
    free(hdr.s);

    // Output all remainig data read with the header block
    if ( fp->block_length - skip_until > 0 )
    {
        if ( bgzf_write(bgzf_out, buffer+skip_until, fp->block_length-skip_until)<0 ) error("Error: %d\n",fp->errcode);
    }
    if ( bgzf_flush(bgzf_out)<0 ) error("Error: %d\n",bgzf_out->errcode);

    // Stream the rest of the file without as it is, without decompressing
    ssize_t nread;
    int page_size = getpagesize();
    char *buf = (char*) valloc(page_size);
    while (1)
    {
        nread = bgzf_raw_read(fp, buf, page_size);
        if ( nread<=0 ) break;

        int count = bgzf_raw_write(bgzf_out, buf, nread);
        if (count != nread) error("Write failed, wrote %d instead of %d bytes.\n", count,(int)nread);
    }
    if (bgzf_close(bgzf_out) < 0) error("Error: %d\n",bgzf_out->errcode);
    free(buf);
}
static void reheader_vcf(args_t *args)
{
    kstring_t hdr = {0,0,0};
    htsFile *fp = args->fp;
    while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )
    {
        kputc('\n',&fp->line);  // hts_getline eats the newline character
        if ( fp->line.s[0]!='#' ) break;
        kputsn(fp->line.s,fp->line.l,&hdr);
    }

    int nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
        samples = hts_readlines(args->samples_fname, &nsamples);
    if ( args->header_fname )
    {
        free(hdr.s); hdr.s = NULL; hdr.l = hdr.m = 0;
        read_header_file(args->header_fname, &hdr);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &hdr);
        int i;
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }
    if ( args->contigs_fname ) set_contigs(args, &hdr);

    int out = STDOUT_FILENO;
    if ( write(out, hdr.s, hdr.l)!=hdr.l ) error("Failed to write %d bytes\n", hdr.l);
    free(hdr.s);
    if ( fp->line.l )
    {
        if ( write(out, fp->line.s, fp->line.l)!=fp->line.l ) error("Failed to write %d bytes\n", fp->line.l);
    }
    while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )   // uncompressed file implies small size, we don't worry about speed
    {
        kputc('\n',&fp->line);
        if ( write(out, fp->line.s, fp->line.l)!=fp->line.l ) error("Failed to write %d bytes\n", fp->line.l);
    }
}

static bcf_hdr_t *strip_header(bcf_hdr_t *src, bcf_hdr_t *dst)
{
    bcf_hrec_t *src_hrec, *dst_hrec, *tmp;
    bcf_hdr_t *out = bcf_hdr_init("r");
    int i;
    for (i=0; i<dst->nhrec; i++)
    {
        // first insert lines which do not code BCF ids, their order does not matter
        dst_hrec = dst->hrec[i];
        if ( dst_hrec->type==BCF_HL_FLT || dst_hrec->type==BCF_HL_INFO || dst_hrec->type==BCF_HL_FMT || dst_hrec->type== BCF_HL_CTG ) continue;
        bcf_hdr_add_hrec(out, bcf_hrec_dup(dst_hrec));
    }
    for (i=0; i<src->nhrec; i++)
    {
        // now transfer header lines which define BCF ids
        src_hrec = src->hrec[i];

        if ( src_hrec->type==BCF_HL_FLT || src_hrec->type==BCF_HL_INFO || src_hrec->type==BCF_HL_FMT || src_hrec->type== BCF_HL_CTG )
        {
            int j = bcf_hrec_find_key(src_hrec, "ID");
            dst_hrec = bcf_hdr_get_hrec(dst, src_hrec->type, "ID", src_hrec->vals[j], NULL);
            if ( !dst_hrec ) continue;

            tmp = bcf_hrec_dup(dst_hrec);

            j = bcf_hrec_find_key(src_hrec, "IDX");
            if ( j>=0 )
            {
                j = atoi(src_hrec->vals[j]);
                hrec_add_idx(tmp, j);
            }
            bcf_hdr_add_hrec(out, tmp);
        }
    }
    bcf_hdr_sync(out);
    for (i=0; i<dst->nhrec; i++)
    {
        // finally add new structured fields
        dst_hrec = dst->hrec[i];
        if ( dst_hrec->type==BCF_HL_FLT || dst_hrec->type==BCF_HL_INFO || dst_hrec->type==BCF_HL_FMT || dst_hrec->type== BCF_HL_CTG )
        {
            int j = bcf_hrec_find_key(dst_hrec, "ID");
            tmp = bcf_hdr_get_hrec(out, dst_hrec->type, "ID", dst_hrec->vals[j], NULL);
            if ( !tmp )
                bcf_hdr_add_hrec(out, bcf_hrec_dup(dst_hrec));
        }
    }
    for (i=0; i<dst->n[BCF_DT_SAMPLE]; i++) bcf_hdr_add_sample(out, dst->samples[i]);
    bcf_hdr_add_sample(out, NULL);
    bcf_hdr_destroy(dst);
    return out;
}

static void reheader_bcf(args_t *args, int is_compressed)
{
    htsFile *fp = args->fp;
    bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to read the header: %s\n", args->fname);
    kstring_t htxt = {0,0,0};
    int hlen;
    htxt.s = bcf_hdr_fmt_text(hdr, 1, &hlen);
    htxt.l = hlen;

    int i, nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
        samples = hts_readlines(args->samples_fname, &nsamples);
    if ( args->header_fname )
    {
        free(htxt.s); htxt.s = NULL; htxt.l = htxt.m = 0;
        read_header_file(args->header_fname, &htxt);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &htxt);
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }

    bcf_hdr_t *hdr_out = bcf_hdr_init("r");
    bcf_hdr_parse(hdr_out, htxt.s);
    if ( args->header_fname ) hdr_out = strip_header(hdr, hdr_out);

    // write the header and the body
    htsFile *fp_out = hts_open("-",is_compressed ? "wb" : "wbu");
    bcf_hdr_write(fp_out, hdr_out);

    bcf1_t *rec = bcf_init();
    while ( bcf_read(fp, hdr, rec)==0 )
    {
        // sanity checking, this slows things down. Make it optional?
        bcf_unpack(rec, BCF_UN_ALL);
        if ( rec->rid >= hdr_out->n[BCF_DT_CTG] || strcmp(bcf_hdr_int2id(hdr,BCF_DT_CTG,rec->rid),bcf_hdr_int2id(hdr_out,BCF_DT_CTG,rec->rid)) )
            error("The CHROM is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_CTG,rec->rid));

        for (i=0; i<rec->d.n_flt; i++)
        {
            int id = rec->d.flt[i];
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_FLT,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken FILTER ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->d.n_flt )
            error("The FILTER is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.flt[i]));

        for (i=0; i<rec->n_info; i++)
        {
            int id = rec->d.info[i].key;
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !hdr_out->id[BCF_DT_ID][id].key ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_INFO,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken INFO ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->n_info )
            error("The INFO tag is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.info[i].key));

        for (i=0; i<rec->n_fmt; i++)
        {
            int id = rec->d.fmt[i].id;
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !hdr_out->id[BCF_DT_ID][id].key ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_FMT,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken FORMAT ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->n_fmt )
            error("The FORMAT tag is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.fmt[i].id));

        bcf_write(fp_out,hdr_out,rec);
    }
    bcf_destroy(rec);

    free(htxt.s);
    hts_close(fp_out);
    bcf_hdr_destroy(hdr_out);
    bcf_hdr_destroy(hdr);
}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Modify header of VCF/BCF files, change sample names.\n");
    fprintf(stderr, "Usage:   bcftools reheader [OPTIONS] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --contigs <file>    set contig lines from .fai or .dict file\n");
    fprintf(stderr, "    -h, --header <file>     new header\n");
    fprintf(stderr, "    -s, --samples <file>    new sample names\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_reheader(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;

    static struct option loptions[] =
    {
        {"contigs",1,0,'c'},
        {"header",1,0,'h'},
        {"samples",1,0,'s'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "s:h:c:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'c': args->contigs_fname = optarg; break;
            case 's': args->samples_fname = optarg; break;
            case 'h': args->header_fname = optarg; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else usage(args);
    }
    else args->fname = argv[optind];

    if ( !args->samples_fname && !args->header_fname && !args->contigs_fname ) usage(args);
    if ( !args->fname ) usage(args);

    args->fp = hts_open(args->fname,"r");
    if ( !args->fp ) error("Failed to open: %s\n", args->fname);
    args->file_type = hts_fd_type(args->fp);

    if ( HTS_FT(args->file_type)==HTS_FT_VCF || HTS_FT(args->file_type)==HTS_FT_TXT )
    {
        if ( args->file_type & HTS_GZ )
            reheader_vcf_gz(args);
        else
            reheader_vcf(args);
    }
    else if ( HTS_FT(args->file_type) == HTS_FT_BCF )
    {
        if ( args->contigs_fname ) error("The -c option works only for VCF, sorry!\n");
        reheader_bcf(args, args->file_type & HTS_GZ);
    }
    else
        error("The file type is not recognised: %s\n", args->fname);

    if ( hts_close(args->fp) ) error("Close failed: %s\n", args->fname);
    free(args);
    return 0;
}
