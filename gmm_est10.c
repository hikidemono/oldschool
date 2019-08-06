/**************************************************
gmm_est10.c


./GMMest01 [SPECTRUM] [OUTOUT]
**************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hkdmn_util.h"

#define NAME_MAX 200
#define NO_VOC -5000
#define FMAX 8400
#define FMIN 3600
#define FTOF 20
#define F0NUM 240
#define DBSMAX 10
#define INIT_VARI 100
#define MIN_LIMIT 0.0000000001
#define THRE 0.0001

//// prototype ////
double **dbs_read(char *filename, int size);
int dbs_size_check(char *filename);
int **dbf_read(char *filename, int size, int cand_num);
double **dbp_read(char *filename, int size, int cand_num);

//// main function ////
int main(int argc, char *argv[]){
	int i,j,k;
	int size;
	double **dbs;
	int mix_num;
	int m_num;
	int **db_f;
	double **db_p;
	double *gmm;
	FILE *fp_out;
	FILE *fp_main;
	// パラメータ //
	double *weight;
	double *mean;
	double *vari;
	double *new_weight;
	double *new_mean;
	double *new_vari;
	double temp_weight;
	double temp_score1;
	double likelihood;
	double old_likelihood;
	double **gamma;
	double **pre_gamma;
	double *N_k;
	int loop_n;
	char *filename;
	double i_vari[5]={10.0, 10.0, 10.0, 10.0, 10.0};
	double i_mean[5]={20.0, 70.0, 120.0, 170.0, 220.0};
	double i_weight[5]={0.2,0.2,0.2,0.2,0.2};


	// 引数チェック //
	mix_num = 5;
	printf("mixnum %d : ", mix_num);

	// FILEからスペクトルread //
	size = dbs_size_check(argv[1]);
	printf("// size : %d\n", size);
	dbs = dbs_read(argv[1], size);
	fp_out=myfopen("hoge_spec.txt", "w");
	for(i=0; i<F0NUM; i++){
		fprintf(fp_out,"%d %lf\n", i, dbs[0][i]);
	}
	fclose(fp_out);
	// 下準備 //
	gmm = doublemalloc(F0NUM, "gmm");
	fp_main = myfopen(argv[2], "w");
	weight = doublemalloc(mix_num, "w");
	mean = doublemalloc(mix_num, "m");
	vari = doublemalloc(mix_num, "v");
	new_weight = doublemalloc(mix_num, "w");
	new_mean = doublemalloc(mix_num, "m");
	new_vari = doublemalloc(mix_num, "v");
	gamma = double2malloc(mix_num, F0NUM, "v");
	pre_gamma = double2malloc(mix_num, F0NUM, "v");
	N_k = doublemalloc(mix_num, "NK");
	filename = charmalloc(NAME_MAX, "hoge");

	// スペクトル近似 //
	// ここから1フレームづつ //
	for(i=0; i<size; i++){
		for(j=0; j<mix_num;j++){
			weight[j] = i_weight[j];
			mean[j] = i_mean[j];
			vari[j] = i_vari[j];
		}
		// 初期パラメータ登録 //
		m_num=mix_num;j=0;
		//printf("m_num : %d\n", m_num);
		temp_weight = 0.0;
		// 初期モデルの作成
		for(k=0; k<F0NUM; k++){
			gmm[k] =  0.0;
		}
		for(j=0;j<m_num;j++){
			for(k=0; k<F0NUM; k++){
				gmm[k] +=  weight[j] * Gaussian(k, mean[j], vari[j]);
			}
		}
		// 確認 //
		fp_out = myfopen("hoge_init.txt", "w");
		for(k=0; k<F0NUM; k++){
			fprintf(fp_out, "%d %.10lf\n", k, gmm[k]);
		}
		fclose(fp_out);
		return 0;

		//// ループ ///
		loop_n=0;
		while(1){
		//for(loop_n = 0; loop_n < 5; loop_n++){
			// Log-Likelihoodの計算
			likelihood = 0.0;
			for(k=0; k<F0NUM; k++){
				temp_score1 = 0.0;
				for(j=0;j<m_num;j++){
					pre_gamma[j][k] = weight[j] * Gaussian(dbs[i][k], mean[j], vari[j]);
					temp_score1 += pre_gamma[j][k];
					//printf("%d %.10lf %.10lf\n", k, temp_score1, pre_gamma[j][k]);
				}
				likelihood += log(temp_score1);
			}

			//printf("%2dth likelihood : %.10lf\n", loop_n, likelihood);
			for(j=0; j<m_num;j++){
				//printf("   %d / w %lf / m %lf / v %lf\n", j, weight[j], mean[j], vari[j]);
			}
			if((loop_n != 0) && (abs(likelihood -old_likelihood) < THRE)){
				//printf("iteration end : %02d-loops\n", loop_n);
				break;
			}

			//// E-STEP更新 ////
			// 負担率の計算 ///
			for(k=0; k<F0NUM; k++){
				temp_score1 = 0.0;
				for(j=0;j<m_num;j++){
					temp_score1 += weight[j] * Gaussian(dbs[i][k], mean[j], vari[j]);
				}
				for(j=0; j<m_num; j++){
					//if(temp_score1 > MIN_LIMIT){
						gamma[j][k] = weight[j] * Gaussian(dbs[i][k], mean[j], vari[j]) / temp_score1;
						//printf("%d %d : %lf %lf %lf\n", j, k, gamma[j][k], pre_gamma[j][k], temp_score1);
					//}
					//else{
					//	gamma[j][k] = 0.0;
					//}
				}
			}

				/*for(j=0; j<m_num; j++){
						printf("%d : %lf\n", j,gamma[j][82]);
				}*/

			//return 0;

			for(k=0; k<F0NUM;k++){
				temp_score1 = 0.0;
				for(j=0; j<m_num; j++){
					temp_score1+=gamma[j][k];
				}
				//printf("temp : %lf\n", temp_score1);
			}

			/// M-STEP ////// パラメータの更新 //
			// N_k
			for(j=0;j<m_num;j++){
				N_k[j] = 0.0;
				for(k=0; k<F0NUM; k++){
					N_k[j] += gamma[j][k];
					//printf("%lf %lf\n", N_k[j],gamma[j][k]);
				}
			}
			/*printf("   N_k :");
			for(j=0;j<m_num;j++){
				printf(" %lf", N_k[j]);
			}
			printf("\n");*/

			// means
			for(j=0;j<m_num;j++){
				new_mean[j] = 0.0;
				for(k=0; k<F0NUM; k++){
					new_mean[j] += gamma[j][k] * k;
					//printf("%lf %lf %lf\n", new_mean[j], gamma[j][k], dbs[i][k]);
				}
				new_mean[j] /= N_k[j];
			}





			// vari
			for(j=0;j<m_num;j++){
				new_vari[j] = 0.0;
				for(k=0; k<F0NUM; k++){
					temp_score1 = dbs[i][k] - new_mean[j];
					new_vari[j] += gamma[j][k] * temp_score1 * temp_score1;
				}
				new_vari[j] /= N_k[j];
			}
			// weight
			for(j=0;j<m_num;j++){
				new_weight[j] = N_k[j] / F0NUM;
			}

			/// new->oldに代入
			for(j=0;j<m_num;j++){
				weight[j] = new_weight[j];
				mean[j] = new_mean[j];
				vari[j] = new_vari[j];
			}
			old_likelihood = likelihood;
			/*// 中間出力 //
			for(k=0; k<F0NUM; k++){
				gmm[k] =  0.0;
			}
			for(j=0;j<m_num;j++){
				for(k=0; k<F0NUM; k++){
					gmm[k] +=  weight[j] * Gaussian(k, mean[j], vari[j]);
				}
			}
			// 確認 //
			sprintf(filename, "hoge_%02dth.txt", loop_n+1);
			fp_out = myfopen(filename, "w");
			for(k=0; k<F0NUM; k++){
				fprintf(fp_out, "%d %.10lf\n", k, gmm[k]);
			}
			fclose(fp_out);*/

			loop_n++;
		} // ループ終わり

		//ファイルに出力//
		for(k=0; k<F0NUM; k++){
			gmm[k] =  0.0;
		}
		for(j=0;j<m_num;j++){
			for(k=0; k<F0NUM; k++){
				gmm[k] +=  weight[j] * Gaussian(k, mean[j], vari[j]);
			}
		}
		for(k=0; k<F0NUM; k++){
			fprintf(fp_main, "%d %d %.10lf\n", i, k, gmm[k]);
		}


	}

	return 0;


	// 後始末 //
	free(fp_main);
	//int2free(db_f,size);
	//double2free(db_p,size);
	double2free(dbs,size);
	return 0;
}



// DBS_PF.Gのサイズを確認する //
int dbs_size_check(char *filename){
	int length=0;
	FILE *fp;
	fp = myfopen(filename, "r");
	while(fscanf(fp, "%*s %*s %*s") != EOF){
		length++;
	}
	length /= F0NUM;
	fclose(fp);
	return length;
}
// DBSをconfigから読み取る //
double **dbs_read(char *filename, int size){
	int i,j;
	FILE *fp;
	double **dbs;
	fp = myfopen(filename, "r");
	dbs=double2malloc(size, F0NUM, "dbs");
	for(i=0; i<size; i++){
		for(j=0; j<F0NUM; j++){
			fscanf(fp,"%*d %*d %lf", &dbs[i][j]);
		}
	}
	fclose(fp);
	return dbs;
}
int **dbf_read(char *filename, int size, int cand_num){
	int i,j;
	FILE *fp;
	int **dbs;
	fp = myfopen(filename, "r");
	dbs = int2malloc(size,cand_num, "dbs_f");
	// ここから読み取り本番
	for(i=0; i<size; i++){
		fscanf(fp, "%*d");
		for(j=0; j<cand_num;j++){
			fscanf(fp, "%d", &dbs[i][j]);
		}
		for(j=cand_num; j<DBSMAX;j++){
			fscanf(fp, "%*d");
		}
		fscanf(fp, "%*d");
		for(j=0; j<DBSMAX;j++){
			fscanf(fp, "%*f");
		}
	}
	fclose(fp);
	return dbs;
}
double **dbp_read(char *filename, int size, int cand_num){
	int i,j;
	FILE *fp;
	double **dbs;
	fp = myfopen(filename, "r");
	dbs = double2malloc(size,cand_num, "dbs_f");
	// ここから読み取り本番
	for(i=0; i<size; i++){
		fscanf(fp, "%*d");
		for(j=0; j<DBSMAX;j++){
			fscanf(fp, "%*d");
		}
		fscanf(fp, "%*d");
		for(j=0; j<cand_num;j++){
			fscanf(fp, "%lf", &dbs[i][j]);
		}
		for(j=cand_num; j<DBSMAX;j++){
			fscanf(fp, "%*f");
		}
	}
	fclose(fp);
	return dbs;
}
