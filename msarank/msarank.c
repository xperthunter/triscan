/****************************************************************************\
 File: msarank.c
 Author: Ian Korf
 License: Public Domain
\****************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>

#define MAX_LINE 1048576

static int VERBOSE = 0;

typedef struct MultipleSequenceAlignment {
	char *id;
	char *ac;
	char *de;
	int   sq;
	char **sids;
	char **seqs;
} MultipleSequenceAlignment, *msa_t;

typedef struct StockholmReader {
	char *filename;
	gzFile fp;
} StockholmReader, *STOCK;

STOCK stock_reader(char *filename) {
	STOCK sfr = malloc(sizeof(StockholmReader));
	sfr->filename = filename;
	sfr->fp = gzopen(filename, "r");
	if (!sfr->fp) {
		fprintf(stderr, "error opening %s\n", filename);
		exit(1);
	}
	char *line = malloc(32);
	gzgets(sfr->fp, line, 32);
	if (strncmp(line, "# STOCKHOLM", 11) != 0) {
		fprintf(stderr, "not STOCKHOLM format\n");
		exit(1);
	}
	free(line);
	return sfr;
}

msa_t read_stockholm(STOCK sfr) {
	msa_t msa = malloc(sizeof(MultipleSequenceAlignment));
	msa->id = NULL;
	msa->ac = NULL;
	msa->de = NULL;
	msa->sq = 0;
	msa->sids = NULL;
	msa->seqs = NULL;

	char *line = malloc(MAX_LINE);

	int record = 0;
	while (gzgets(sfr->fp, line, MAX_LINE) != NULL) {
		line[strcspn(line, "\r\n")] = '\0';
		if (strncmp(line, "#=GF ID", 7) == 0) {
			msa->id = strdup(line+9);
		} else if (strncmp(line, "#=GF AC", 7) == 0) {
			msa->ac = strdup(line+9);
		} else if (strncmp(line, "#=GF DE", 7) == 0) {
			msa->de = strdup(line+9);
		} else if (strncmp(line, "#=GF SQ", 7) == 0) {
			msa->sq = atoi(line+9);
			msa->sids = malloc(sizeof(char*) * msa->sq);
			msa->seqs = malloc(sizeof(char*) * msa->sq);
		} else if (strcmp(line, "//") == 0) {
			return msa;
		} else if (line[0] == '#') {
			continue;
		} else {
			int seq_start = 0;
			for (int i = 0; i < MAX_LINE; i++) {
				if (line[i] < 33) {
					line[i] = '\0';
					seq_start = i;
					break;
				}
			}
			for (int i = seq_start; i < MAX_LINE; i++) {
				if (line[i] >= 'A') {
					seq_start = i;
					break;
				}
			}

			msa->sids[record] = strdup(line);
			msa->seqs[record] = strdup(line+seq_start);
			record++;
		}
	}

	return NULL;
}

void free_msa(msa_t msa) {
	free(msa->id);
	free(msa->ac);
	free(msa->de);
	for (int i = 0; i < msa->sq; i++) free(msa->sids[i]);
	for (int i = 0; i < msa->sq; i++) free(msa->seqs[i]);
}

double seq_similarity(const char *s1, const char *s2) {
	int slen = strlen(s1);
	int match = 0;
	int total = 0;
	for (int i = 0; i < slen; i++) {
		if (s1[i] == '.') continue;
		if (s2[i] == '.') continue;
		if (s1[i] == s2[i]) match++;
		total++;
	}
	return (double)match/(double)total;
}

static char *help = "\
usage: msarank [options] <stockholm>\n\
  -t <float> threshold (not yet implemented)\n\
  -c <int>   cpus (not yet implemented)\n\
  -v         verbose\n\
";

int main(int argc, char **argv) {
	double  threshold = 0.5;
	int     cpus = 2;

	// CLI
	if (argc == 1) {printf("%s", help); exit(1);}
	int opt;
	while ((opt = getopt(argc, argv, "i:e:o:t:a:v")) != -1) {
		switch (opt) {
			case 't': threshold = atof(optarg); break;
			case 'c': cpus      = atoi(optarg); break;
			case 'v': VERBOSE   = 1;
		}
	}
	char *filename = argv[argc-1];

	if (VERBOSE) printf("%s %d %f\n", filename, cpus, threshold);

	// main loop
	STOCK sfp = stock_reader(filename);
	msa_t msa;
	while ((msa = read_stockholm(sfp)) != NULL) {
		printf("%s %s %s %d\n", msa->ac, msa->id, msa->de, msa->sq);
		// allocate matrix for all pairwise comparisions
		double ** mat = malloc(sizeof(double*) * msa->sq);
		for (int i = 0; i < msa->sq; i++)
			mat[i] = malloc(sizeof(double) * msa->sq);

		// make comparisions (mirror diagonal)
		for (int i = 0; i < msa->sq; i++) {
			for (int j = i+1; j < msa->sq; j++) {
				mat[i][j] = seq_similarity(msa->seqs[i], msa->seqs[j]);
				mat[j][i] = mat[i][j];
			}
		}

		// calculate similarity of each row (skipping self)
		for (int i = 0; i < msa->sq; i++) {
			double sum = 0;
			for (int j = 0; j < msa->sq; j++) {
				if (i == j) continue;
				sum += mat[i][j];
			}
			double average = sum / (msa->sq - 1);
			if (VERBOSE) printf("%s %f\n", msa->sids[i], average);
		}

		// cleanup
		for (int i = 0; i < msa->sq; i++) free(mat[i]);
		free(mat);
		free_msa(msa);
	}

	exit(0);
}
