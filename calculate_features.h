#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>

int getNumPeptides(const char *filename) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	int numPeptides = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
		printf("File %s not found\n", filename);
        exit(1);
	}

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
		numPeptides ++;
    }

    fclose(fp);
    if (line)
        free(line);
	return numPeptides;
}

void readPeptides(const char *filename, char ***peptides) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	int numPeptides = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
		printf("File %s not found\n", filename);
        exit(1);
	}

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
		(*peptides)[numPeptides] = (char *) malloc (sizeof(char) * (read+1));
		strcpy((*peptides)[numPeptides], line);
		numPeptides ++;
    }

    fclose(fp);
    if (line)
        free(line);
}

void readSupportVectors(const char *filename, float ***support_vectors) {
	printf("Reading support vector file %s ...\n", filename);
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	int row = 0;
	int col = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
		printf("File %s not found\n", filename);
        exit(1);
	}

    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
    	char *pt;
    	pt = strtok (line,"\t");
	
		col = 0;
    	while (pt != NULL) {
			(*support_vectors)[col][row] = atof(pt);
        	pt = strtok (NULL, ",");
			col ++;
    	}
		row ++;
    }

    fclose(fp);
    if (line)
        free(line);
	printf("... Completed reading support vector file %s \n", filename);
}

int getVectorLength(const char *filename) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	int row = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
		printf("File %s not found\n", filename);
        exit(1);
	}

    while ((read = getline(&line, &len, fp)) != -1) {
		row ++;
    }

    fclose(fp);
    if (line)
        free(line);
	return row;
}

void readVector(const char *filename, float **vector) {
	printf("Reading file %s ...\n", filename);
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	int row = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
		printf("File %s not found\n", filename);
        exit(1);
	}

    while ((read = getline(&line, &len, fp)) != -1) {
    	char *pt;
    	pt = strtok (line,"\t");
	
    	while (pt != NULL) {
			(*vector)[row] = atof(pt);
        	pt = strtok (NULL, ",");
    	}
		row ++;
    }

    fclose(fp);
    if (line)
        free(line);
	printf("... Completed reading file %s \n", filename);
}

//Substitute for matlab strfind
//https://www.mathworks.com/help/matlab/ref/strfind.html
int numOfOccurenceOfCharInString(char **peptideAddress, char c) {
	char *peptide = *peptideAddress;
	int len = strlen(peptide);
	int count = 0;
	for (int i=0; i<len; i++) {
		if (c == peptide[i]) {
			count ++;
		}
	}
	return count;
}

int sum_dot(float *v1, float *v2, int len) {
	float sum = 0.0;
	for (int i=0; i<len; i++) {
		sum += v1[i]*v2[i];
	}
	return sum;
}

/////bond_3D.m//////
void get_amino_acid_wgt_mat(char ***peptidesAddress, int numPeptides, double **wgt_mat_addr) {
	char **peptides = (*peptidesAddress);
	double *wgt_mat = (*wgt_mat_addr);
	for (int i=0; i<numPeptides; i++) {
		char *peptide = peptides[i];
		int peptide_len = strlen(peptide) - 1;
		for (int j=0; j<peptide_len-2; j++) {
			int index_i = ((int)peptide[j]) - (96-32);
			int index_j = ((int)peptide[j+1]) - (96-32);
			int index_k = ((int)peptide[j+2]) - (96-32);
			wgt_mat[26*26*index_i + 26*index_j + index_k] += 1.0;
		}
	}
	double total_wgt = 0;
	for (int i=0; i<26*26*26; i++) {
		wgt_mat[i] += 2.2204e-16;
		total_wgt += wgt_mat[i];
	}
	for (int i=0; i<26*26*26; i++) {
		//printf("abcde %lf, %lf, %lf\n", wgt_mat[i], wgt_mat[i]/total_wgt, log(wgt_mat[i]/total_wgt));
		wgt_mat[i] = log(wgt_mat[i]/total_wgt);
	}
}
////////////////////
/////amino_bond_3D.m//////
void get_amino_acid_bond_vector(char ***peptidesIdentifiedAddress, 
		char ***peptidesUnIdentifiedAddress, 
		int numPeptidesIdentified, 
		int numPeptidesUnIdentified, 
		double **amino_acid_feature_vector_address) {
	printf("Creating amino acid vector ...\n");
	char **peptidesIdentified = (*peptidesIdentifiedAddress);
	char **peptidesUnIdentified = (*peptidesUnIdentifiedAddress);
	double *amino_acid_feature_vector = (*amino_acid_feature_vector_address);
	
	double *amino_acid_wgt_mat_identified = (double *) malloc(sizeof(double) * 26 * 26 * 26);
	for (int i=0; i<26 * 26 * 26; i++) {
		amino_acid_wgt_mat_identified[i] = 0.0;
	}
	double *amino_acid_wgt_mat_unidentified = (double *) malloc(sizeof(double) * 26 * 26 * 26);
	for (int i=0; i<26 * 26 * 26; i++) {
		amino_acid_wgt_mat_unidentified[i] = 0.0;
	}
	get_amino_acid_wgt_mat(&peptidesIdentified, numPeptidesIdentified, &amino_acid_wgt_mat_identified);
	get_amino_acid_wgt_mat(&peptidesUnIdentified, numPeptidesUnIdentified, &amino_acid_wgt_mat_unidentified);

	for (int i=0; i<numPeptidesIdentified; i++) {
		double total_wgt_identified = 0.0;
		double total_wgt_unidentified = 0.0;
		char *peptide = peptidesIdentified[i];
		int pep_length = strlen(peptide) - 1;
		for (int j=0; j<pep_length-2; j++) {
			int index_i = ((int)peptide[j]) - (96-32);
			int index_j = ((int)peptide[j+1]) - (96-32);
			int index_k = ((int)peptide[j+2]) - (96-32);
			double weight_for_identified = amino_acid_wgt_mat_identified[26*26*index_i + 26*index_j + index_k];
			double weight_for_unidentified = amino_acid_wgt_mat_unidentified[26*26*index_i + 26*index_j + index_k];
			total_wgt_identified += weight_for_identified;
			total_wgt_unidentified += weight_for_unidentified;
		}
		//printf("i wgti=%lf, wgtj=%lf\n", total_wgt_identified, total_wgt_unidentified);
		amino_acid_feature_vector[i]= total_wgt_identified/total_wgt_unidentified;
	}
	for (int i=0; i<numPeptidesUnIdentified; i++) {
		double total_wgt_identified = 0.0;
		double total_wgt_unidentified = 0.0;
		char *peptide = peptidesUnIdentified[i];
		int pep_length = strlen(peptide) - 1;
		for (int j=0; j<pep_length-2; j++) {
			int index_i = ((int)peptide[j]) - (96-32);
			int index_j = ((int)peptide[j+1]) - (96-32);
			int index_k = ((int)peptide[j+2]) - (96-32);
			double weight_for_identified = amino_acid_wgt_mat_identified[26*26*index_i + 26*index_j + index_k];
			double weight_for_unidentified = amino_acid_wgt_mat_unidentified[26*26*index_i + 26*index_j + index_k];
			total_wgt_identified += weight_for_identified;
			total_wgt_unidentified += weight_for_unidentified;
		}
		//printf("wgti=%lf, wgtj=%lf\n", total_wgt_identified, total_wgt_unidentified);
		amino_acid_feature_vector[i+numPeptidesIdentified] = total_wgt_identified/total_wgt_unidentified;
	}

	free(amino_acid_wgt_mat_identified);
	free(amino_acid_wgt_mat_unidentified);
	printf("... Created amino acid vector\n");
}
////////////////////

void calculate_features(float ***featureMatrix, float **y, int *numOnes, int *numZeros, int readingTest) {

	char *identified_filename;
	if (readingTest == 0) {
		identified_filename = strdup("dataFiles/identifiedYP.txt");
	}
	else {
		identified_filename = strdup("testFiles/identifiedEK12.txt");
	}
	int numPeptidesIdentified = getNumPeptides(identified_filename);
	printf("Num of identified Peptides = %d\n", numPeptidesIdentified);
	*numOnes = numPeptidesIdentified;

	char **peptides_identified = (char **) malloc(numPeptidesIdentified * sizeof(char *));
	readPeptides(identified_filename, &peptides_identified);

	char *unidentified_filename;
	if (readingTest == 0) {
		unidentified_filename = strdup("dataFiles/un_identifiedYP.txt");
	}
	else {
		unidentified_filename = strdup("testFiles/un_identifiedEK12.txt");
	}
	int numPeptidesUnIdentified = getNumPeptides(unidentified_filename);
	printf("Num of unidentified Peptides = %d\n", numPeptidesUnIdentified);
	*numZeros = numPeptidesUnIdentified;

	char **peptides_unidentified = (char **) malloc(numPeptidesUnIdentified * sizeof(char *));
	readPeptides(unidentified_filename, &peptides_unidentified);

	const char *sv_filename = "STEPP_SupportVectors.txt";
	float **support_vectors = (float **) malloc(sizeof(float *) * 35);
	for (int i=0; i<35; i++) {
		support_vectors[i] = (float *) malloc(sizeof(float) * 9653);
	}
	readSupportVectors(sv_filename, &support_vectors);

	const char *wgt_filename = "STEPP_Weights.txt";
	int wgt_vec_length = getVectorLength(wgt_filename);
	float *support_weights = (float *) malloc(sizeof(float) * wgt_vec_length);
	readVector(wgt_filename, &support_weights);

	const char *fv_mw_filename = "Feature_Values_MW.txt";
	int fv_mw_vec_length = getVectorLength(fv_mw_filename);
	float *fv_mw = (float *) malloc(sizeof(float) * fv_mw_vec_length);
	readVector(fv_mw_filename, &fv_mw);

	const char *fv_eh_filename = "Feature_Values_EH.txt";
	int fv_eh_vec_length = getVectorLength(fv_eh_filename);
	float *fv_eh = (float *) malloc(sizeof(float) * fv_eh_vec_length);
	readVector(fv_eh_filename, &fv_eh);

	const char *fv_hh_filename = "Feature_Values_HH.txt";
	int fv_hh_vec_length = getVectorLength(fv_hh_filename);
	float *fv_hh = (float *) malloc(sizeof(float) * fv_hh_vec_length);
	readVector(fv_hh_filename, &fv_hh);

	const char *fv_kh_filename = "Feature_Values_KH.txt";
	int fv_kh_vec_length = getVectorLength(fv_kh_filename);
	float *fv_kh = (float *) malloc(sizeof(float) * fv_kh_vec_length);
	readVector(fv_kh_filename, &fv_kh);

	const char *fv_rm_filename = "Feature_Values_RM.txt";
	int fv_rm_vec_length = getVectorLength(fv_rm_filename);
	float *fv_rm = (float *) malloc(sizeof(float) * fv_rm_vec_length);
	readVector(fv_rm_filename, &fv_rm);

	const char *fv_gp_filename = "Feature_Values_GP.txt";
	int fv_gp_vec_length = getVectorLength(fv_gp_filename);
	float *fv_gp = (float *) malloc(sizeof(float) * fv_gp_vec_length);
	readVector(fv_gp_filename, &fv_gp);

	const char *fv_zp_filename = "Feature_Values_ZP.txt";
	int fv_zp_vec_length = getVectorLength(fv_zp_filename);
	float *fv_zp = (float *) malloc(sizeof(float) * fv_zp_vec_length);
	readVector(fv_zp_filename, &fv_zp);

	const char *fv_bk_filename = "Feature_Values_BK.txt";
	int fv_bk_vec_length = getVectorLength(fv_bk_filename);
	float *fv_bk = (float *) malloc(sizeof(float) * fv_bk_vec_length);
	readVector(fv_bk_filename, &fv_bk);

	const char *nf_mean_filename = "STEPP_NormFactor_mean.txt";
	int nf_mean_vec_length = getVectorLength(nf_mean_filename);
	float *nf_mean = (float *) malloc(sizeof(float) * nf_mean_vec_length);
	readVector(nf_mean_filename, &nf_mean);

	const char *nf_std_filename = "STEPP_NormFactor_std.txt";
	int nf_std_vec_length = getVectorLength(nf_std_filename);
	float *nf_std = (float *) malloc(sizeof(float) * nf_std_vec_length);
	readVector(nf_std_filename, &nf_std);

	float vec_count[20];
	char *pattern = "ACDEFGHIKLMNPQRSTVWY";	

	int numPeptides = numPeptidesIdentified + numPeptidesUnIdentified;
	(*featureMatrix) = (float **) malloc (sizeof(float *) * numPeptides);
	*y = (float *) malloc (sizeof(float) * numPeptides);

	for (int i=0; i<numPeptides; i++) {
		(*featureMatrix)[i] = (float *) malloc (sizeof(float) * (35+1)); //+1 for amino acid feature
		float *feature_vector = (*featureMatrix)[i];
		int peptide_len = 0;
		if (i < numPeptidesIdentified) {
			//Step 1: transform the peptide string into
			for (int j=0; j<20; j++) {
				vec_count[j] = numOfOccurenceOfCharInString(&(peptides_identified[i]), pattern[j]);
			}
			peptide_len = strlen(peptides_identified[i]) - 1;
			(*y)[i] = 1;
		} 
		else {
			//Step 1: transform the peptide string into
			for (int j=0; j<20; j++) {
				vec_count[j] = numOfOccurenceOfCharInString(&(peptides_unidentified[i-numPeptidesIdentified]), pattern[j]);
			}
			peptide_len = strlen(peptides_unidentified[i-numPeptidesIdentified]) - 1;
			(*y)[i] = -1;
		}
		//Step 2: compute the 35 feature values 
		// Feature 1 is length
		feature_vector[0] = peptide_len;
		// Feature 2 is molecular weight
		feature_vector[1] = sum_dot(fv_mw,vec_count,20) - (feature_vector[0] - 1.0)*18.015;
		// Feature 3 is the number of non-polar hydrophobic residues
		feature_vector[2] = vec_count[0] + vec_count[4] + vec_count[5] + vec_count[7] + vec_count[9] + vec_count[10] + vec_count[12] + vec_count[17] + vec_count[18] + vec_count[19];
		// Feature 4 is the number of polar hydrophilic residues
		feature_vector[3] = feature_vector[0] - feature_vector[2];
		// Feature 5 is the number of uncharged polar hydrophilic residues
		feature_vector[4] = vec_count[1] + vec_count[11] + vec_count[13] + vec_count[15] + vec_count[16];
		// Feature 6 is the number of charged polar hydrophilic residues
		feature_vector[5] = vec_count[2] + vec_count[3] + vec_count[6] + vec_count[8] + vec_count[14];
		// Feature 7 is the number of positively charged polar hydrophilic residues
		feature_vector[6] = vec_count[6] + vec_count[8] + vec_count[14];
	        // Feature 8 is the number of negatively charged polar hydrophilic residues
		feature_vector[7] = vec_count[2] + vec_count[3];
	        // Feature 9 is Hydrophobicity_Eisenberg
		feature_vector[8] = sum_dot(fv_eh, vec_count, 20)/feature_vector[0];
	        // Feature 10 is Hydrophilicity Hopp-Woods
		feature_vector[9] = sum_dot(fv_hh, vec_count, 20)/feature_vector[0];
	        // Feature 11 is Hydrophobicity_Kyte-Doolittle
		feature_vector[10] = sum_dot(fv_kh, vec_count, 20)/feature_vector[0];
	        // Feature 12 is Hydropathicity Roseman
		feature_vector[11] = sum_dot(fv_rm, vec_count, 20)/feature_vector[0];
	        // Feature 13 is polarity Grantham
		feature_vector[12] = sum_dot(fv_gp, vec_count, 20)/feature_vector[0];
	        // Feature 14 is polarity zimmerman
		feature_vector[13] = sum_dot(fv_zp, vec_count, 20)/feature_vector[0];
	        // Feature 15 is bulkiness
		feature_vector[14] = sum_dot(fv_bk, vec_count, 20)/feature_vector[0];
	        // Features 16-35 are amino acid counts
		for (int k=15; k<35; k++) {
			feature_vector[k] = vec_count[k-15];
		}
	
	}
	double *amini_acid_feature = (double *) malloc(sizeof(double) * numPeptides);
	get_amino_acid_bond_vector(&peptides_identified, &peptides_unidentified, numPeptidesIdentified, numPeptidesUnIdentified, &amini_acid_feature);
	for (int i=0; i<numPeptides; i++) {
		float *feature_vector = (*featureMatrix)[i];
		feature_vector[35] = amini_acid_feature[i];
	}

	//Free all
	free(amini_acid_feature);
	free(support_weights);
	free(fv_mw);
	free(fv_eh);
	free(fv_hh);
	free(fv_kh);
	free(fv_rm);
	free(fv_gp);
	free(fv_zp);
	free(fv_bk);
	free(nf_mean);
	free(nf_std);
	for (int i=0; i<35; i++) { //+1 for amino acid feature
		free(support_vectors[i]);
	}
	free(support_vectors);

	for (int i=0; i<numPeptidesIdentified; i++) {
		//printf("%d: %s\n", i+1, peptides_identified[i]);
		free(peptides_identified[i]);
	}
	free(peptides_identified);

	for (int i=0; i<numPeptidesUnIdentified; i++) {
		//printf("%d: %s\n", i+1, peptides_unidentified[i]);
		free(peptides_unidentified[i]);
	}
	free(peptides_unidentified);
}

void scale_features(int numPeptides, float ***featureMatrixAddress) {
	float **featureMatrix = (*featureMatrixAddress);
	float *min = (float *) malloc (sizeof(float) * 36);
	float *max = (float *) malloc (sizeof(float) * 36);
	for (int i=0; i<36; i++) {
		min[i] = DBL_MAX;
		max[i] = -DBL_MAX;
	}

	for (int i=0; i<numPeptides; i++) {
		float *featureVector = featureMatrix[i];
		for (int j=0; j<36; j++) {
			if (featureVector[j] < min[j]) {
				min[j] = featureVector[j];
			}
			if (featureVector[j] > max[j]) {
				max[j] = featureVector[j];
			}
		}
	}

	for (int i=0; i<numPeptides; i++) {
		float *featureVector = featureMatrix[i];
		for (int j=0; j<36; j++) {
			if (fabs(max[j] - min[j]) > 1e-8) { //divide by zero 
				featureVector[j] = (featureVector[j] - min[j])/(max[j] - min[j]);
			}
		}
	}

	free(min);
	free(max);
	
}

void print_features(int numPeptides, float ***featureMatrixAddress, float *y) {
	float **featureMatrix = (*featureMatrixAddress);

	for (int i=0; i<numPeptides; i++) {
		float *featureVector = featureMatrix[i];
		for (int j=0; j<35; j++) {
			printf("%0.4f,", featureVector[j]);
		}
		printf("%d\n",(int) y[i]);
	}

}
