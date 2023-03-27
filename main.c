#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <htslib/khash.h>

#define BUFSIZE 100000
#define NKMERS_INIT 100
#define KMER_LENGTH 31

// Creating a hash table type linking a k-mer to its index in the input LD matrix
KHASH_MAP_INIT_STR(kmer_map, int);

// Creating a hash table type that represents a set of k-mers
KHASH_SET_INIT_STR(kmer_set);

// A struct that contains the k-mers sequence and presence/absence data of the input file
typedef struct {
	char **kmers;
	int **pav;
	int n_kmers; // Number of rows
	int n_samples; // Number of columns
} kmerData;

// A struct that contains the LD values for a set of k-mers; essentially a square matrix with k-mer identifiers
typedef struct {
	char **kmers;
	double **ld;
	int n_kmers;
} ldMatrix;

// A function that reads the first row from the file and returns the the number of fields
int read_first_row(FILE* input);

// A function that reads the k-mer sequences and presence/absence data from the file and returns a kmerData struct
// input: the file to read the input from (must have been processed with read_first_row before)
// kmer_size: the number of rows to initialize the k-mer sequence and pav struct members with
// n_samples: the number of samples to read per row of the input file (should be known from read_first_row)
kmerData read_pav(FILE* input, int kmer_size, int n_samples);

// A function that prints the contents of the pav member of a kmerData matrix
void print_pav(kmerData pav_data);

// A function that de-allocates memory allocated to a kmerData struct
void free_kmers(kmerData pav_data);

// A function that de-allocates memory allocated to a ldMatrix struct
void free_ld(ldMatrix ld_matrix);

// A function that takes a kmerData struct as input and returns a ldMatrix struct
// of the pairwise LD values between all k-mers
ldMatrix compute_ld(kmerData pav_data);

// A function that computes the sum of an int array
int sum(int *array, int length);

// A function that takes two integer arrays as input and counts the number of
// positions where both values are 1
int joint_sum(int *a, int *b, int length);

// A function that prints the contents of a LD matrix (only the numeric part)
void print_ld(ldMatrix ld_matrix);

// A function that clusters the most similar k-mers in a LD matrix using a greedy algorithm
ldMatrix cluster_ld(ldMatrix ld_matrix);

int main(int argc, char* argv[]) {
	// Opening the k-mers presence/absence matrix for reading
	FILE *pav_file = fopen(argv[1], "r");

	// Creating variables
	int n_columns;
	kmerData kmer_data;
	ldMatrix ld_matrix, clustered_ld;

	// Determining the number of columns and samples from the first row of the file
	n_columns = read_first_row(pav_file);
	fprintf(stderr, "The file %s contains %d columns.\n", argv[1], n_columns);

	// Reading the k-mer sequences and presence/absence data from the file
	kmer_data = read_pav(pav_file, NKMERS_INIT, n_columns - 1);
	fprintf(stderr, "%d k-mers read from %d samples.\n", kmer_data.n_kmers, kmer_data.n_samples);

	// We no longer need the file to be open
	fclose(pav_file);

	// Debug printing
	//print_pav(kmer_data);
	
	// Computing the pairwise LD between all k-mers
	ld_matrix = compute_ld(kmer_data);

	// Debug printing
	// print_ld(ld_matrix);

	// Clustering the k-mers in that matrix using the LD values
	clustered_ld = cluster_ld(ld_matrix);

	// Debug printing
	print_ld(clustered_ld);

	// De-allocating memory
	free_kmers(kmer_data);
	free_ld(ld_matrix);

	return 0;
}

// A function that reads the first row from the file and returns the the number of fields
int read_first_row(FILE* input) {

	// The variable that will be output
	int n_columns = 0;
	char buffer[BUFSIZE], *token;

	// Checking that the file was properly opened
	if(fgets(buffer, BUFSIZE, input) == NULL) {
		fprintf(stderr, "Unable to read input file. Aborting\n");
		exit(1);
	}

	// Removing the newline character from the input string
	buffer[strcspn(buffer, "\n")] = '\0';

	// Splitting the first line into fields in order to count them
	token = strtok(buffer, "\t");

	while(token != NULL) {
		n_columns++;
		token = strtok(NULL, "\t");
	}

	return n_columns;
}

// A function that reads the k-mer sequences and presence/absence data from the file and returns the number of records
kmerData read_pav(FILE* input, int kmer_size, int n_samples) {

	// Initializing some variables
	int n_kmers = 0;
	char buffer[BUFSIZE], *token;
	kmerData pav_data;

	// Initializing the struct members for k-mer sequence and pav data
	pav_data.kmers = (char**) malloc(kmer_size * sizeof(char*));
	pav_data.pav = (int**) malloc(kmer_size * sizeof(int*));

	// Next we need to read the rows one after the other
	while(fgets(buffer, BUFSIZE, input) != NULL)	{
		buffer[strcspn(buffer, "\n")] = '\0';

		// We ensure dynamic memory allocation of the kmers and pav_data arrays
		if(n_kmers >= kmer_size) {
			kmer_size *= 2;
			pav_data.pav = (int**) realloc(pav_data.pav, kmer_size * sizeof(int*));
			pav_data.kmers = (char**) realloc(pav_data.kmers, kmer_size * sizeof(char*));
		}

		// The first value read from the tab-separated column is the k-mer sequence
		token = strtok(buffer, "\t");
		pav_data.kmers[n_kmers] = strdup(token);

		// Allocating memory for the current pav_data array and filling it
		pav_data.pav[n_kmers] = (int*) malloc(n_samples * sizeof(int));

		for(int i = 0; i < n_samples; i++) {
			token = strtok(NULL, "\t");
			pav_data.pav[n_kmers][i] = atoi(token);
		}

		// The next token read should be NULL after that
		token = strtok(NULL, "\t");
		if(token != NULL) {
			fprintf(stderr, "Data fields present beyond expected length. Aborting.\n");
			exit(1);
		}

		n_kmers++;
	}

	pav_data.n_samples = n_samples;
	pav_data.n_kmers = n_kmers;

	return pav_data;
}

// A function that prints the contents of a pav_data matrix
void print_pav(kmerData pav_data) {

	for(int i = 0; i < pav_data.n_kmers; i++) {
		for(int j = 0; j < pav_data.n_samples; j++) {
			printf("%d\t", pav_data.pav[i][j]);
		}
		printf("\n");
	}
}

// A function that de-allocates memory allocated to a kmerData struct
void free_kmers(kmerData pav_data) {

	for(int i = 0; i < pav_data.n_kmers; i++) {
		free(pav_data.kmers[i]);
		free(pav_data.pav[i]);
	}

	free(pav_data.kmers);
	free(pav_data.pav);
}

// A function that takes a kmerData struct as input and returns a ldMatrix struct
// of the pairwise LD values between all k-mers
ldMatrix compute_ld(kmerData pav_data) {
	// Creating the output object
	ldMatrix ld_matrix;

	// Creating double vectors that store the proportion of 1s and 0s for a given k-mer
	double *p1, *p0, d, pij; // d and pij are used in LD calculation
	p1 = (double*) malloc(pav_data.n_kmers * sizeof(double));
	p0 = (double*) malloc(pav_data.n_kmers * sizeof(double));

	// The number of k-mers is the same as in the input pav_data
	ld_matrix.n_kmers = pav_data.n_kmers;

	// The k-mers can also be directly copied from the kmerData input
	ld_matrix.kmers = (char**) malloc(ld_matrix.n_kmers * sizeof(char*));

	for(int i = 0; i < pav_data.n_kmers; i++) {
		ld_matrix.kmers[i] = strdup(pav_data.kmers[i]);
	}

	// Allocating memory for the ld double matrix
	ld_matrix.ld = (double**) malloc(ld_matrix.n_kmers * sizeof(double*));

	for(int i = 0; i < ld_matrix.n_kmers; i++) {
		ld_matrix.ld[i] = (double*) malloc(ld_matrix.n_kmers * sizeof(double));
	}

	// Pre-computing the frequencies of each k-mer across samples to speed up computation
	for(int i = 0; i < pav_data.n_kmers; i++) {
		p1[i] = sum(pav_data.pav[i], pav_data.n_samples) / (double) pav_data.n_samples;
		p0[i] = 1.0 - p1[i];
	}

	// Now actually computing the LD values and filling the matrix
	for(int i = 0; i < ld_matrix.n_kmers - 1; i++) {
		for(int j = i + 1; j < ld_matrix.n_kmers; j++) {
			pij = joint_sum(pav_data.pav[i], pav_data.pav[j], pav_data.n_samples) / (double) pav_data.n_samples;
			// Computing the LD and setting the value directly in the ld matrix
			d = pij - p1[i] * p1[j];
			ld_matrix.ld[i][j] = (d * d) / (p1[i] * p1[j] * p0[i] * p0[j]);
			// Also setting the values at [j][i] because the matrix is symmetric
			ld_matrix.ld[j][i] = ld_matrix.ld[i][j];
		}
	}

	// Also setting all the values along the diagonal to 1
	for(int i = 0; i < ld_matrix.n_kmers; i++) {
		ld_matrix.ld[i][i] = 1.0;
	}

	return ld_matrix;
}

// A function that computes the sum of an int array
int sum(int *array, int length) {
	int output = 0;

	for(int i = 0; i < length; i++) {
		output += array[i];
	}

	return output;
}

// A function that takes two integer arrays as input and counts the number of
// positions where both values are 1
int joint_sum(int *a, int *b, int length) {
	int output = 0;

	for(int i = 0; i < length; i++) {
		output += a[i] && b[i];
	}

	return output;
}

// A function that prints the contents of a LD matrix
// k-mer sequences are printed at the beginning of the line
void print_ld(ldMatrix ld_matrix) {
	for(int i = 0; i < ld_matrix.n_kmers; i++) {
		printf("%s\t", ld_matrix.kmers[i]);
		for(int j = 0; j < ld_matrix.n_kmers; j++) {
			printf("%f", ld_matrix.ld[i][j]);
			if(j < ld_matrix.n_kmers - 1) printf("\t");
		}

		printf("\n");
	}
}

// A function that de-allocates memory allocated to a ldMatrix struct
void free_ld(ldMatrix ld_matrix) {
	for(int i = 0; i < ld_matrix.n_kmers; i++) {
		free(ld_matrix.kmers[i]);
		free(ld_matrix.ld[i]);
	}

	free(ld_matrix.kmers);
	free(ld_matrix.ld);
}

// A function that clusters the most similar k-mers in a LD matrix using a greedy algorithm
ldMatrix cluster_ld(ldMatrix ld_matrix) {

	// Creating the output ldMatrix struct and initializing/allocating memory for the members
	ldMatrix clust_matrix;
	clust_matrix.n_kmers = ld_matrix.n_kmers;
	clust_matrix.kmers = (char**) malloc(clust_matrix.n_kmers * sizeof(char*));
	clust_matrix.ld = (double**) malloc(clust_matrix.n_kmers * sizeof(double*));

	for(int i = 0; i < clust_matrix.n_kmers; i++) {
		clust_matrix.ld[i] = (double*) malloc(clust_matrix.n_kmers * sizeof(double));
	}

	// Creating an array of k-mers that have been clustered
	// The last one in the array is the next one to be used for clustering
	char** clustered_kmers = (char**) malloc(ld_matrix.n_kmers * sizeof(char*));

	// Creating other character and integer variables
	char *i_kmer, *j_kmer, *best_kmer;
	int i_index, j_index, n_clustered = 0;
	double max_ld;

	// Initializing basic khash-related variables
	int kh_return;
	khiter_t k;

	// Initializing a hash table that maps a k-mer to its index in the input ld_matrix
	khash_t(kmer_map) *kmer_index_map = kh_init(kmer_map);

	// Initializing two hash tables that are sets representing k-mers that
	// have been clustered and k-mers that have yet to be clustered
	khash_t(kmer_set) *to_cluster, *clustered;
	to_cluster = kh_init(kmer_set);
	clustered = kh_init(kmer_set);

	// Filling the hash table with its values; also filling the to_cluster set
	for(int i = 0; i < ld_matrix.n_kmers; i++) {
		// Setting the values in the kmer_index_map hash table
		k = kh_put(kmer_map, kmer_index_map, strdup(ld_matrix.kmers[i]), &kh_return);
		kh_value(kmer_index_map, k) = i;

		// Filling the set with k-mers that have yet to be clustered (all except the first)
		if(i > 0) {
			k = kh_put(kmer_set, to_cluster, strdup(ld_matrix.kmers[i]), &kh_return);
		}
	}

	// Adding the first k-mer to the set of clustered k-mers
	k = kh_put(kmer_set, clustered, strdup(ld_matrix.kmers[0]), &kh_return);
	clustered_kmers[0] = strdup(ld_matrix.kmers[0]);
	n_clustered++;
	
	// Debug printing
	//puts("Begin clustering");

	// Adding the k-mer with the highest LD relative to the last one added,
	// until all k-mers have been clustered
	while(kh_size(to_cluster) > 0) {
		//printf("The size of the to_cluster hash set is %d\n", kh_size(to_cluster));

		// Getting the sequence of the last clustered k-mer and its index in the original LD matrix
		i_kmer = clustered_kmers[n_clustered - 1];
		i_index = kh_value(kmer_index_map, kh_get(kmer_map, kmer_index_map, i_kmer));
		//printf("i_kmer: %s ; i_index: %d\n", i_kmer, i_index);

		// Resetting the maximum LD
		max_ld = -1.0;

		// We need to find the sample with the highest_ld relative to the last clustered k-mer
		for(k = kh_begin(to_cluster); k != kh_end(to_cluster); k++) {
			if(kh_exist(to_cluster, k) == 0) continue;

			j_kmer = strdup(kh_key(to_cluster, k));
			j_index = kh_value(kmer_index_map, kh_get(kmer_map, kmer_index_map, j_kmer));
			//printf("j_kmer: %s ; j_index: %d\n", j_kmer, j_index);

			if(ld_matrix.ld[i_index][j_index] > max_ld) {
				max_ld = ld_matrix.ld[i_index][j_index];
				best_kmer = strdup(j_kmer);
			}
		}

		// Now that the best k-mer has been found, we need to:
		// First, add it to the clustered set
		k = kh_put(kmer_set, clustered, strdup(best_kmer), &kh_return);
		// Second, remove it from the to_cluster set
		k = kh_get(kmer_set, to_cluster, best_kmer);
		kh_del(kmer_set, to_cluster, k);
		// Third, add it to the clustered_kmers array and increment the number of clustered k-mers
		clustered_kmers[n_clustered++] = strdup(best_kmer);
	}

	// Now we need to arrange the values in the output matrix according to the order of the k-mers in the cluster
	if(n_clustered != clust_matrix.n_kmers) {
		fprintf(stderr, "The number of expected k-mers (%d) does not match the number of clustered k-mers (%d)",
				clust_matrix.n_kmers, n_clustered );
		exit(1);
	}

	// Setting the k-mer sequence in the output struct
	for(int i = 0; i < n_clustered; i++) {
		clust_matrix.kmers[i] = strdup(clustered_kmers[i]);
	}

	// Looping over all k-mer sequences to set the values in the matrix
	for(int i = 0; i < clust_matrix.n_kmers; i++) {
		for(int j = 0; j < clust_matrix.n_kmers; j++) {
			i_index = kh_value(kmer_index_map, kh_get(kmer_map, kmer_index_map, clust_matrix.kmers[i]));
			j_index = kh_value(kmer_index_map, kh_get(kmer_map, kmer_index_map, clust_matrix.kmers[j]));
			clust_matrix.ld[i][j] = ld_matrix.ld[i_index][j_index];
		}
	}

	return clust_matrix;
}

