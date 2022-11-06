# include "hmm.h"
 
# define NUMBER_OF_MODEL 5 // according to the spec
# define MAX_SEQ_NUM 2500 // according to the spec
# define MAX_OBS_NUM 51 // according to the spec
using namespace std;

int main(int argc , char **argv)
{
	HMM models[NUMBER_OF_MODEL];
	load_models(argv[1] , models, NUMBER_OF_MODEL);
	int seqs_num = 0;
	int l = 0;
	char obs_seqs[MAX_SEQ_NUM][MAX_OBS_NUM] = {{0}};
	
	FILE *test_data = fopen(argv[2], "r");
	while(fscanf(test_data, "%s", obs_seqs[seqs_num])!=EOF){
		seqs_num++;
	}
	fclose(test_data);
	
	FILE *output_file = fopen(argv[3], "w");

	for(int i=0; i<seqs_num; i++){
		double P_max = 0;
		int which_model = 0;
		char *obs = obs_seqs[i];
		int l = strlen(obs);
		for(int m=0; m<NUMBER_OF_MODEL; m++){
			// Initialize the delta table
			double delta[MAX_OBS_NUM][MAX_STATE] = {{0}};
			for(int i=0; i<models[m].state_num; i++){
				delta[0][i] = models[m].initial[i] * models[m].observation[obs[0]-'A'][i];
			}
		
			for(int t=1; t<l; t++){
				for(int j=0; j<models[m].state_num; j++){
					double v = 0;
					for(int i = 0; i<models[m].state_num; i++){
						double now = delta[t-1][i] * models[m].transition[i][j];
						if(now > v){
							v = now;
						}
					}
					delta[t][j] = v * models[m].observation[obs[t]-'A'][j];
				}
			}

			double P = 0;
			for(int i = 0; i<models[m].state_num; i++){
				if(delta[l-1][i] > P){
					P = delta[l-1][i];
				}
			}
			if(P_max < P){
				P_max = P;
				which_model = m;
			}	 
		}
		
		fprintf(output_file, "%s %e\n", models[which_model].model_name, P_max);
	}
	fclose(output_file);
	return 0;
}
