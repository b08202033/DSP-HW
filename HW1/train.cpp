# include "hmm.h"

// from the slides, the sequence number in training set is 10000
#ifndef MAX_SEQ_NUM
#	define MAX_SEQ_NUM	10000
#endif

// from the data file, the number of observations per sequence in the dataset is 50
#ifndef MAX_OBS_NUM
#	define MAX_OBS_NUM	51
#endif

// integrated function interface for forward/backward algorithms
void fb(HMM* model, double alpha[MAX_OBS_NUM][MAX_STATE], double beta[MAX_OBS_NUM][MAX_STATE], char *obs, char mode){
	if(mode=='f'||mode=='F'){
		// initialize the alpha table
		for(int state=0; state<model->state_num; state++){
			alpha[0][state] = model->initial[state] * model->observation[(obs[0] - 'A')][state];		
		}
		
		// Induction
		for(int t=1; t<strlen(obs); t++){
			for(int state=0; state<model->state_num; state++){
				alpha[t][state] = 0;
				for(int j=0; j<model->state_num;j++){
					alpha[t][state] += alpha[t-1][j] * model->transition[j][state];
				}
				alpha[t][state] *= model->observation[obs[t] - 'A'][state];
			}
		}	
	}
	else if(mode=='b'||mode=='B'){
		// intialization
		for(int state=0; state<model->state_num; state++){
			beta[strlen(obs)-1][state] = 1;		
		}
		
		//Induction
		for(int t=strlen(obs)-2; t>=0; t--){
			for(int i=0; i<model->state_num; i++){
				beta[t][i] = 0;
				for(int j=0; j<model->state_num; j++){
					beta[t][i] += model->transition[i][j] * model->observation[obs[t+1]-'A'][j] * beta[t+1][j];
				}
			}
		}
	}
	else{
		printf("Invalid mode, please check the code.\n");
	}
} 



// the integrated function of calculating gamma and epsilon
void r3(HMM *model, int number_of_sequences, double alpha[MAX_OBS_NUM][MAX_STATE], double beta[MAX_OBS_NUM][MAX_STATE],  double epsilon[MAX_OBS_NUM][MAX_STATE][MAX_STATE], 
		double gamma1[MAX_OBS_NUM][MAX_STATE], double gamma2[MAX_OBS_NUM][MAX_STATE][MAX_OBSERV], char *obs, char mode){

	switch(mode){
		case 'e': case 'E':
			for(int t=0; t<strlen(obs) - 1; t++){ //only up to strlen(obs)-1, since the boundary condition will be violated due to t+1
				double d = 0;
				for(int i=0; i<model->state_num; i++){
					for(int j=0; j<model->state_num; j++){
						d += alpha[t][i] * model->transition[i][j] * model->observation[obs[t+1] - 'A'][j] * beta[t+1][j];
					}
				}
			
				for(int i=0; i<model->state_num; i++){
					for(int j=0; j<model->state_num; j++){
						epsilon[t][i][j] += ((alpha[t][i] * model->transition[i][j] * model->observation[obs[t+1] - 'A'][j] * beta[t+1][j]) / d) / number_of_sequences;
					}
				}
			}
			break;
		case 'r': case 'R':
			for(int t=0; t<strlen(obs); t++){
				double d = 0;
				for(int i=0; i<model->state_num; i++){
					d += alpha[t][i] * beta[t][i];
				}
			
				for(int i=0; i<model->state_num; i++){
					double a = ((alpha[t][i] * beta[t][i]) / d) / number_of_sequences;
					gamma1[t][i] += a;
					gamma2[t][i][obs[t] - 'A'] += a;
				}
			}
			break;
		default:
			printf("Invalid mode.");
			break;
	}
	
}

void prepare(HMM *model, int number_of_sequences, double alpha[MAX_OBS_NUM][MAX_STATE], double beta[MAX_OBS_NUM][MAX_STATE],  double epsilon[MAX_OBS_NUM][MAX_STATE][MAX_STATE], 
		double gamma1[MAX_OBS_NUM][MAX_STATE], double gamma2[MAX_OBS_NUM][MAX_STATE][MAX_OBSERV], char *obs)
{
	// forward
	fb(model, alpha, beta, obs, 'f');
	// backward
	fb(model, alpha, beta, obs, 'b');
	// gamma
	r3(model, number_of_sequences, alpha, beta, epsilon, gamma1, gamma2, obs, 'r');
	// epsilon
	r3(model, number_of_sequences, alpha, beta, epsilon, gamma1, gamma2, obs, 'e');
				
}

void trainer(HMM *model , char *file , int iteration){
	// the trainer of HMM
	char observation_sequence[MAX_SEQ_NUM][MAX_OBS_NUM] = {{0}};
	
	// read the training data
	int number_of_sequences = 0; 
	FILE *training_set = fopen(file, "r");
	while (fscanf(training_set , "%s", observation_sequence[number_of_sequences]) != EOF){
		number_of_sequences++; // count the number of input sequence 
	}
	fclose(training_set);
	
	int len = strlen(observation_sequence[number_of_sequences-1]);
	
	int iter = 0;
	double alpha[MAX_OBS_NUM][MAX_STATE];
	double beta[MAX_OBS_NUM][MAX_STATE];
		
	while(iter < iteration){
		// Initialize all the arrays	
		double epsilon[MAX_OBS_NUM][MAX_STATE][MAX_STATE] = {{{0}}};
		double gamma1[MAX_OBS_NUM][MAX_STATE] = {{0}};
		double gamma2[MAX_OBS_NUM][MAX_STATE][MAX_OBSERV] = {{{0}}};
		
		// update the arrays
		for(int j = 0; j < number_of_sequences; j++){
			prepare(model, number_of_sequences, alpha, beta, epsilon, gamma1, gamma2, observation_sequence[j]);
		}
		
		// re-estimate
		// Pi
		for(int i = 0; i <model->state_num; i++){
			model->initial[i] = gamma1[0][i];
		}
		// A
		for(int i=0; i<model->state_num; i++){
			for(int j=0; j<model->state_num; j++){
				double sum_e = 0;
				double sum_r = 0;
				for(int t=0; t<len-1; t++){
					sum_e += epsilon[t][i][j];
					sum_r += gamma1[t][i];
				}
				model->transition[i][j] = sum_e / sum_r;
			}
		}
		// B
		for(int k=0; k<model->observ_num; k++){
			for(int j=0; j<model->state_num; j++){
				double numerator = 0;
				double denominator = 0;
				for(int t=0; t<len; t++){
					denominator += gamma1[t][j];
					numerator += gamma2[t][j][k];
				}
				
				model->observation[k][j] = numerator / denominator;
			}
		}
		
		iter++;
	}
}

int main(int argc , char **argv){
	// Initialize the model
	HMM model;
	
	// load the model
	loadHMM(&model , argv[2]);
	
	// train the model
	trainer(&model, argv[3], atoi(argv[1]));
	
	// dump
	FILE *output_file = fopen(argv[4] , "w");
	dumpHMM(output_file , &model);
	fclose(output_file);
	
	return 0;
}








