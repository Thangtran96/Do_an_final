
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <utility>

using namespace std;
// output
ofstream fo("Ans.txt");

// Cac bien hang so
const int ARRAY_SIZE = 12005;
const int ARRAY_BYTES_INT = ARRAY_SIZE * sizeof(int);

//cac bien chinh
int l = 9, d = 3;
char cDataInp[ARRAY_SIZE];
int h_dataMotif[12005];
string sDataInp[20];

//input tu file
void File_Input()
{
	l = 9; d = 2;
	FILE * pFile;
	pFile = fopen("datacu.txt", "r");
	if (pFile == NULL)
		perror("Error opening file");
	else
	{
		if (fgets(cDataInp, ARRAY_SIZE, pFile) != NULL)
			cout << "nhap du lieu thanh cong!\n";
		fclose(pFile);
	}

	for (int i = 0; i < strlen(cDataInp); ++i) {
		//A=0 C=1 G=2 T=3
		switch (cDataInp[i])
		{
		case 'A': { h_dataMotif[i] = 0; break; }
		case 'C': { h_dataMotif[i] = 1; break; }
		case 'G': { h_dataMotif[i] = 2; break; }
		case 'T': { h_dataMotif[i] = 3; break; }
		default: cout << "error chuyen sang int";
			break;
		}
	}
	int k = 0;
	string temp = cDataInp;
	cout << temp << endl;
	for (int i = 0; i < temp.size(); i += 600) {
		sDataInp[k] = temp.substr(i, 600);
		cout << k << ". " << sDataInp[k] << endl;
		k++;
	}
}

__global__ void patternBarching(const int* d_datainp, const int l, const int d, int *ans) {
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index < 600 - l) {
		//khai bao bien
		int motif_temp[40];
		int temp_val;
		int temp_dis;
		int best_dis = 99999;
		int motif_bN[40];
		int score_motif;
		int s1[40];

		//lay chuoi can duyet
		for (int i = 0; i < l; ++i) {
			motif_temp[i] = d_datainp[i + index];
			motif_bN[i] = motif_temp[i];
			s1[i] = motif_temp[i];
		}
		//ham dis_hamming
		int ans_Ham = 0;
		int temp, tempRow;
		for (int i = 0; i < 20; ++i)
		{
			tempRow = 999;
			for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
			{
				temp = 0;
				for (int k = 0; k < l; k++) {
					if (s1[k] != d_datainp[k + j]) temp++;
				}
				if (temp < tempRow) tempRow = temp;
			}
			ans_Ham += tempRow;
		}
		int sorce_Hamming = ans_Ham;
		//End ham hamming

		//chay ham patternBarching
		for (int k = 0; k < d; ++k) {
			//kiem tra chuoi tot
			//printf("\n 2 \n");
			if (best_dis < score_motif) {
				score_motif = best_dis;
				for (int i = 0; i < l; ++i) {
					motif_temp[i] = motif_bN[i];
				}
			}
			//ham bestNeighbor
			//printf("\nbestNeighbor\n");
			for (int i = 0; i < l; ++i) {
				//printf("\n 3 \n");
				//trg hop 0
				if (motif_temp[i] != 0) {
					temp_val = motif_temp[i];
					motif_temp[i] = 0;
					//temp_dis = dis_haming(d_datainp, motif_temp, l);
					//lay best neighbor
					if (temp_dis < best_dis)
					{
						best_dis = temp_dis;
						for (int j = 0; j < l; ++j) {
							motif_bN[j] = motif_temp[j];
						}
					}
					motif_temp[i] = temp_val;
				}
				//trg hop 1
				if (motif_temp[i] != 1) {
					temp_val = motif_temp[i];
					motif_temp[i] = 1;
					//temp_dis = dis_haming(d_datainp, motif_temp, l);
					//lay best neighbor
					if (temp_dis < best_dis)
					{
						best_dis = temp_dis;
						for (int j = 0; j < l; ++j) {
							motif_bN[j] = motif_temp[j];
						}
					}
					motif_temp[i] = temp_val;
				}
				//trg hop 2
				if (motif_temp[i] != 2) {
					temp_val = motif_temp[i];
					motif_temp[i] = 2;
					//temp_dis = dis_haming(d_datainp, motif_temp, l);
					//lay best neighbor
					if (temp_dis < best_dis)
					{
						best_dis = temp_dis;
						for (int j = 0; j < l; ++j) {
							motif_bN[j] = motif_temp[j];
						}
					}
					motif_temp[i] = temp_val;
				}
				//trg hop 3
				if (motif_temp[i] != 3) {
					temp_val = motif_temp[i];
					motif_temp[i] = 3;
					//temp_dis = dis_haming(d_datainp, motif_temp, l);
					//lay best neighbor
					if (temp_dis < best_dis)
					{
						best_dis = temp_dis;
						for (int j = 0; j < l; ++j) {
							motif_bN[j] = motif_temp[j];
						}
					}
					motif_temp[i] = temp_val;
				}
			}
			// END ham bestNeighbor
		}
		//printf("\n 4 \n");
		//du lieu tra lai
		//printf("\n gan du lieu vao d_motif \n");
		//End ham
	}
}

int main()
{
	File_Input();
	return 0;
}