#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctime>
#include <algorithm>
#include"svm.h"
#include"calculate_features.h"

#define UNEQUAL_WEIGHT_UNEQUAL_SAMPLES 0
#define EQUAL_WEIGHT_EQUAL_SAMPLES 1
#define NUM_FOLDS 10
#define SCALE_FEATURES 1
#define PRINT_FEATURES 0
#define GRID_SEARCH 0
#define PARAM_NU 0.57
#define PARAM_C 7e4

//Setting algorithm type
//#define VERSION UNEQUAL_WEIGHT_UNEQUAL_SAMPLES
#define VERSION EQUAL_WEIGHT_EQUAL_SAMPLES

//See documentation from: https://github.com/cjlin1/libsvm
void convertToLibsvmFormat(float ***featureMatrixAddress, 
		float **yAddress, 
		int numPeptides,
		int keep1s, int keep0s,
		struct svm_problem *prob) {
	
	float **featureMatrix = *featureMatrixAddress;
	float *y = *yAddress;

	prob->l = keep1s + keep0s;
	prob->y = (double *) malloc (sizeof(double) * prob->l);

	srand(1);
	int *perm = (int *) malloc (sizeof(int) * numPeptides);
	for (int p=0; p<numPeptides; p++) {
		perm[p] = p;
	}
	std::random_shuffle(perm, &perm[numPeptides-1]);

	prob->x = (struct svm_node **) malloc (prob->l * sizeof(struct svm_node *));

	int keep1Counter = 0;
	int keep0Counter = 0;
	int ii=0;
	for (int i=0; i<numPeptides; i++) {
		int p = perm[i];
		//p = i;
		if (fabs(y[p]-1.0) < 1e-8) { //y[p] is equal to 1
			if (keep1Counter >= keep1s) {
				continue;
			}
			keep1Counter ++;
		}
		if (fabs(y[p]+1.0) < 1e-8) { //y[p] is equal to -1
			if (keep0Counter >= keep0s) {
				continue;
			}
			keep0Counter ++;
		}
		prob->y[ii] = y[p];

		int numNonZerosInRow = 0;
		float *featureVector = featureMatrix[p];
		for (int j=0; j<36; j++) {
			if ((j == 6 || j == 29 || j == 0
					|| j == 2 || j == 3 || j == 5
					|| j == 35
					)
					&& fabs(featureVector[j]) > 1e-8) { 
				numNonZerosInRow ++;
			}
		}
		prob->x[ii] = (struct svm_node *) malloc ((numNonZerosInRow + 1) * sizeof(struct svm_node));
		int k = 0;

		if (featureVector[6] > 1e-8) {
			struct svm_node elem;
			elem.index = 1;
			elem.value = featureVector[6];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[29] > 1e-8) {
			struct svm_node elem;
			elem.index = 2;
			elem.value = featureVector[29];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[0] > 1e-8) {
			struct svm_node elem;
			elem.index = 3;
			elem.value = featureVector[0];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[2] > 1e-8) {
			struct svm_node elem;
			elem.index = 4;
			elem.value = featureVector[2];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[3] > 1e-8) {
			struct svm_node elem;
			elem.index = 5;
			elem.value = featureVector[3];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[5] > 1e-8) {
			struct svm_node elem;
			elem.index = 6;
			elem.value = featureVector[5];
			prob->x[ii][k] = elem;
			k++;
		}
		if (featureVector[35] > 1e-8) {
			struct svm_node elem;
			elem.index = 7;
			elem.value = featureVector[35];
			prob->x[ii][k] = elem;
			k++;
		}
		//Setting last element as -1
		struct svm_node elem;
		elem.index = -1;
		elem.value = 0.0;
		prob->x[ii][k] = elem;
		ii ++;
	}
	free(perm);
}

void print_svm_problem_data(svm_problem *prob) {
	printf("Number of rows is %d\n", prob->l);
	printf("Printing first few rows only\n");
	for (int i=0; i<prob->l && i<10; i++) {
		printf("%lf: ", prob->y[i]);

		for (int j=0; ; j++) {
			struct svm_node elem = prob->x[i][j];
			printf("(%d, %0.5lf)  ", elem.index, elem.value);
			if (elem.index == -1) {
				printf("\n");
				break;
			}
		}
	}
}

void print_kernel(int kernel_type) {
	switch (kernel_type) {
		case LINEAR:
		printf("Using LINEAR kernel");
		break;
		case POLY:
		printf("Using POLYNOMIAL kernel");
		break;
		case RBF:
		printf("Using RBF kernel");
		break;
		case SIGMOID:
		printf("Using SIGMOID kernel");
		break;
	}
}

void free_svm_prob(svm_problem *probAddr) {
	svm_problem prob = *probAddr;
	free(prob.y);
	for (int i=0; i<prob.l; i++) {
		free(prob.x[i]);
	}
	free(prob.x);
}

void evaluate_model(int numPeptides, double **yAddr, double **predYAddr,
		double *sensitivity, double *specificity, 
		double *precision, double *recall, 
		double *f1, double *accuracy) {

	double *y = *yAddr;
	double *predY = *predYAddr;
	int tp = 0;
	int tn = 0;
	int fp = 0;
	int fn = 0;
	for (int i=0; i<numPeptides; i++) {
		if (fabs(y[i]-1.0) < 1e-8) { //Actual is 1
			if (fabs(predY[i]-1.0) < 1e-8) { //Pred is 1
				tp ++;
			}
			else { //Pred is 0
				fn ++;
			}
		}
		else { //Actual is 0
			if (fabs(predY[i]-1.0) < 1e-8) { //Pred is 1
				fp ++;
			}
			else { //Pred is 0
				tn ++;
			}
		}
	}
	printf("\t\t------Confusion matrix------\n");
	printf("\t\t\t\t    Predicted\t\t\n");
	printf("\t\t\t\t\tN\t\t\tY\t\t\n");
	printf("\t\tAct(N)\t\t%d\t\t\t%d\t\t\n", tn, fp);
	printf("\t\tAct(Y)\t\t%d\t\t\t%d\t\t\n", fn, tp);
	
	*sensitivity = (tp*1.0)/((tp+fn)*1.0);
	*specificity = (tn*1.0)/((tn+fp)*1.0);
	
	printf("\n\n");
	printf("sensitivity = %f\n", *sensitivity);
	printf("specificity = %f\n", *specificity);

	*precision = (tp*1.0)/((tp+fp)*1.0);
	*recall = (tp*1.0)/((tp+fn)*1.0);

	printf("\n\n");
	printf("precision = %f\n", *precision);
	printf("recall = %f\n", *recall);

	*f1 = 2* (*precision) * (*recall) /((*precision)+(*recall));
	printf("f1-score = %f\n", *f1);

	*accuracy = ((tp+tn)*1.0)/((tn+fp+fn+tp)*1.0);
	printf("accuracy = %f\n", *accuracy);
}

int main(int argc, char **argv) {
	float **featureMatrix;
	float **testFeatureMatrix;
	float *y;
	float *testY;
	int numPeptidesIdentified;
	int numPeptidesUnIdentified;
	int testNumPeptidesIdentified;
	int testNumPeptidesUnIdentified;
	calculate_features(&featureMatrix, &y, &numPeptidesIdentified, &numPeptidesUnIdentified, 0);
	calculate_features(&testFeatureMatrix, &testY, &testNumPeptidesIdentified, &testNumPeptidesUnIdentified, 1);

	int numPeptides = numPeptidesIdentified + numPeptidesUnIdentified;
	int testNumPeptides = testNumPeptidesIdentified + testNumPeptidesUnIdentified;

	if (SCALE_FEATURES) {
		printf("Scaling features....\n");
		scale_features(numPeptides, &featureMatrix);
		scale_features(testNumPeptides, &testFeatureMatrix);
		printf("feature scaling done....\n");
	}
	else {
		printf("Features are NOT scaled\n");
		printf("Set SCALE_FEATURES to 1 in top of main.cpp and recompile\n");
	}
	if (PRINT_FEATURES) {
		printf("Printing features....\n");
		print_features(numPeptides, &featureMatrix, y);
		printf("Printed features.\n");

		printf("Printing features....\n");
		print_features(testNumPeptides, &testFeatureMatrix, testY);
		printf("Printed features.\n");
		exit(1);
	}

	struct svm_problem prob;
	struct svm_problem all_prob;
	struct svm_problem test_prob;
	struct svm_parameter param;

	int keep1s = numPeptidesIdentified;
	int keep0s = numPeptidesUnIdentified;
	convertToLibsvmFormat(&featureMatrix, &y, numPeptides, keep1s, keep0s, &all_prob);

	keep1s = testNumPeptidesIdentified;
	keep0s = testNumPeptidesUnIdentified;
	convertToLibsvmFormat(&testFeatureMatrix, &testY, testNumPeptides, keep1s, keep0s, &test_prob);

	//For making the model, we sample some of the data
	if (VERSION) {
		//Sampling both identified and unidentified equally
		keep1s = numPeptidesIdentified/4;
		keep0s = keep1s;
	}
	else {
		//Sampling unidentified in the same ratio as in the data
		keep1s = numPeptidesIdentified/4;
		keep0s = keep1s * ceil((numPeptidesUnIdentified*1.0)/numPeptidesIdentified);
	}
	convertToLibsvmFormat(&featureMatrix, &y, numPeptides, keep1s, keep0s, &prob);


	//Free featureMatrix
	for (int i=0; i<(numPeptides); i++) {
		free(featureMatrix[i]);
	}
	free(featureMatrix);
	free(y);

	printf("free done\n");

	//Free testFeatureMatrix
	for (int i=0; i<(testNumPeptides); i++) {
		free(testFeatureMatrix[i]);
	}
	free(testFeatureMatrix);
	free(testY);

	printf("test free done\n");

	double *target = NULL;
	double *test_target = NULL;
	if (VERSION) {
	param.svm_type = NU_SVC; 
	}
	else {
	param.svm_type = C_SVC; 
	}
	param.kernel_type = LINEAR;
	print_kernel(param.kernel_type);
	printf("\n");

	param.eps = 1e-1;
	param.shrinking = 0;
	param.cache_size = 1000; //in MBs
	if (VERSION) {
		param.nr_weight = 0; //Sets equal weight for both labels
	}
	else {
		int wgt_lbl[2];
		wgt_lbl[0] = 1; //LABEL for identified
		wgt_lbl[1] = -1; //LABEL for unidentified
		param.weight_label = wgt_lbl;
		double wgt[2];
		//Since there are more Unidentified than identified, giving higher weight to identified
		wgt[0] = (numPeptidesUnIdentified*1.0)/numPeptidesIdentified; 
		wgt[1] = 1.0; 
		param.weight = wgt;
		param.nr_weight = 2; //Two labels
	}

	double sensitivity, specificity, precision, recall, f1, accuracy;

	if (GRID_SEARCH) {
		target = (double *) malloc (prob.l * sizeof(double));
		double possibleCs [] = {6e5, 4e5, 2e5};
		int length_C_tries = 4;

		double possibleNUs [] = {0.46,	0.47,	0.48,	0.49,	0.50,	0.51,	0.52,	0.53,	0.54,	0.55,	0.56,	0.57,	0.58,	0.59,	0.60,	0.61,	0.62,	0.63};
		int length_nu_tries = 20;
		int total_tries = length_C_tries;

		if (VERSION) {
			total_tries = length_nu_tries;
		}

		for (int c=0; c<total_tries; c++) {
			for (int i=0; i<(prob.l); i++) {
				target[i]= -1.0;
			}
			printf("nu....nu....nu....nu....\n");
			if (VERSION) {
			param.nu = possibleNUs[c];
			}
			else {
			param.C = possibleCs[c];
			}
			printf("nu Equals = %f\n\n", param.nu);
	
			int nr_fold = NUM_FOLDS;
			printf("\nChecking svm parameters....\n");
			svm_check_parameter(&prob, &param);
			printf(".... SVM parameters found okay\n");
			printf("\nPrinting features....\n");
			print_svm_problem_data(&prob);
			printf(".... Printed features.\n");
			printf("\nStarting svm now....\n");
			svm_cross_validation(&prob, &param, nr_fold, target);
			printf(".... Ending svm now.\n");
			
			evaluate_model(prob.l,
					&(prob.y), &target,
					&sensitivity, &specificity, 
					&precision, &recall, 
					&f1, &accuracy);
			printf("\nRESULT: ");
			print_kernel(param.kernel_type);
			printf(",");
			if (VERSION) {
			printf("nu Equals = %f,", param.nu);
			}
			else {
			printf("C Equals = %f,", param.C);
			}
			printf("sensitivity=%f,specificity=%f,precision=%f,recall=%f,f1=%f\n", sensitivity, specificity, precision, recall, f1);
			if (VERSION) {
			printf("Table,%f,%f,%f,%f,%f,%f\n",param.nu,sensitivity,specificity,precision,recall,f1);
			}
			else {
			printf("Table,%f,%f,%f,%f,%f,%f\n",param.C,sensitivity,specificity,precision,recall,f1);
			}
	
			printf("DONE_nu....DONE_nu....DONE_nu....DONE_nu....\n\n");
		}
	
		free(target);
	}
	else {
		param.shrinking = 0;
		param.eps = 1e-5;
		if (VERSION) {
		param.nu = PARAM_NU;
		}
		else {
		param.C = PARAM_C;
		param.eps = 1e-2;
		}
		printf("\nChecking svm parameters....\n");
		svm_check_parameter(&all_prob, &param);
		printf(".... SVM parameters found okay\n");
	
		printf("Start training...\n");
		struct svm_model *model = svm_train(&prob, &param);
	
		target = (double *) malloc (all_prob.l * sizeof(double));
		for (int i=0; i<(all_prob.l); i++) {
			target[i]= svm_predict(model, all_prob.x[i]);
		}

		evaluate_model(all_prob.l,
				&(all_prob.y), &target,
				&sensitivity, &specificity, 
				&precision, &recall, 
				&f1, &accuracy);
		printf("sensitivity=%f,specificity=%f,precision=%f,recall=%f,f1=%f\n", sensitivity, specificity, precision, recall, f1);
		printf("accuracy=%f\n", accuracy);
		free(target);

		test_target = (double *) malloc (test_prob.l * sizeof(double));
		for (int i=0; i<(test_prob.l); i++) {
			test_target[i]= svm_predict(model, test_prob.x[i]);
		}

		evaluate_model(test_prob.l,
				&(test_prob.y), &test_target,
				&sensitivity, &specificity, 
				&precision, &recall, 
				&f1, &accuracy);
		printf("sensitivity=%f,specificity=%f,precision=%f,recall=%f,f1=%f\n", sensitivity, specificity, precision, recall, f1);
		printf("accuracy=%f\n", accuracy);
		free(test_target);
	}

	free_svm_prob(&prob);
	free_svm_prob(&all_prob);
	free_svm_prob(&test_prob);
}
