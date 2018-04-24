
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <utility>
#include <ctime>
#include <time.h>

using namespace std;
// output
ofstream fo("Ans.txt");

// Cac bien hang so
const int ARRAY_SIZE_INP = 12005;
const int ARRAY_BYTES_INP = ARRAY_SIZE_INP * sizeof(int);
const int ARRAY_SIZE_OUT = 605;
const int ARRAY_BYTES_OUT = ARRAY_SIZE_OUT * sizeof(int);

//cac bien chinh
int l = 9, d = 2;
char cDataInp[ARRAY_SIZE_INP];
int h_dataMotif[ARRAY_SIZE_INP];
string sDataInp[20];

struct Motif_Ans
{
	int dis;
	string motif;
	int adress[20];
};

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
		if (fgets(cDataInp, ARRAY_SIZE_INP, pFile) != NULL)
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
	//cout << temp << endl;
	for (int i = 0; i < temp.size(); i += 600) {
		sDataInp[k] = temp.substr(i, 600);
		//cout << k << ". " << sDataInp[k] << endl;
		k++;
	}
}

int score_ham(string s1, string s2)
{
	int res = 0;
	for (int i = 0; i<s1.size(); ++i) if (s1[i] != s2[i]) res++;
	return res;
}

Motif_Ans dis_hamming(string s)
{
	Motif_Ans res;
	res.motif = s;
	int res_Sum = 0, temp_score = 999, temp_Adress;
	for (int i = 0; i<20; ++i)
	{
		string s1 = sDataInp[i];
		temp_score = 999;
		for (int j = 0; j < s1.size() - l + 1; ++j)
		{
			string temp_str = s1.substr(j, l);
			int score_s = score_ham(s, temp_str);
			if (score_s < temp_score)
			{
				temp_score = score_s;
				temp_Adress = j + 1;
			}
		}
		res_Sum += temp_score;
		res.adress[i] = temp_Adress;
	}
	res.dis = res_Sum;
	return res;
}

__global__ void patternBarching(const int* d_datainp, const int l, const int d, int *ans) {
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if (index < 600 - l) {
		//printf("\n %d", index);
		int ansMotif_sorce = 999;// motif tra ra
		int ansMotif_string[40];//motif tra ra

		int motif_NeSorce = 999;//kq tra ve ham NE
		int motif_NeString[40];//kq tra ve ham NE
		int temp_Sorce = 999;
		int temp_Str[40];

		//cat chuoi motif
		for (int i = 0; i < l; ++i) {
			ansMotif_string[i] = d_datainp[i + index];
			motif_NeString[i] = ansMotif_string[i];
		}
		//begin tinh hamming
		int tempRow, tempSubRow;
		for (int i = 0; i < 20; ++i)
		{
			tempRow = 999;
			for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
			{
				tempSubRow = 0;
				for (int k = 0; k < l; k++) {
					if (ansMotif_string[k] != d_datainp[k + j]) tempSubRow++;
				}
				if (tempSubRow < tempRow) tempRow = tempSubRow;
			}
			ansMotif_sorce += tempRow;
		}
		//end tinh hamming cho chuoi vao

		//begin tinh pattern branching
		for (int a = 0; a <= d; a++) {
			//kiem tra motif dis
			if (motif_NeSorce < ansMotif_sorce) {
				ansMotif_sorce = motif_NeSorce;
				for (int i = 0; i < l; ++i) {
					ansMotif_string[i] = motif_NeString[i];
					temp_Str[i] = motif_NeString[i];
				}
			}
			else
			{//gan bien Ham Ne
				for (int i = 0; i < l; ++i) {
					temp_Str[i] = ansMotif_string[i];
				}
			}//end kiem tra motif

			//begin ham bestNeighbor
			int change = -1;
			for (int b = 0; b < l; ++b) {
				//trg hop 0 A
				if (temp_Str[b] != 0) {
					change = temp_Str[b];
					temp_Str[b] = 0;
					temp_Sorce = 0;//diem dis
					//begin tinh hamming
					for (int i = 0; i < 20; ++i)
					{
						tempRow = 999;
						for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
						{
							tempSubRow = 0;
							for (int k = 0; k < l; k++) {
								if (temp_Str[k] != d_datainp[k + j]) tempSubRow++;
							}
							if (tempSubRow < tempRow) tempRow = tempSubRow;
						}
						temp_Sorce += tempRow;
					}
					//end tinh hamming cho chuoi vao
					//kiem tra dis motif Ne
					if (temp_Sorce < motif_NeSorce) {
						motif_NeSorce = temp_Sorce;
						for (int c = 0; c < l; ++c) {
							motif_NeString[c] = temp_Str[c];
						}
					}
					temp_Str[b] = change;//tra lai gia tri ban dau
				}
				//trg hop 1 C
				if (temp_Str[b] != 1) {
					change = temp_Str[b];
					temp_Str[b] = 1;
					temp_Sorce = 0;//diem dis
					//begin tinh hamming
					for (int i = 0; i < 20; ++i)
					{
						tempRow = 999;
						for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
						{
							tempSubRow = 0;
							for (int k = 0; k < l; k++) {
								if (temp_Str[k] != d_datainp[k + j]) tempSubRow++;
							}
							if (tempSubRow < tempRow) tempRow = tempSubRow;
						}
						temp_Sorce += tempRow;
					}
					//end tinh hamming cho chuoi vao
					//kiem tra dis motif Ne
					if (temp_Sorce < motif_NeSorce) {
						motif_NeSorce = temp_Sorce;
						for (int c = 0; c < l; ++c) {
							motif_NeString[c] = temp_Str[c];
						}
					}
					temp_Str[b] = change;
				}
				//trg hop 2 G
				if (temp_Str[b] != 2) {
					change = temp_Str[b];
					temp_Str[b] = 2;
					temp_Sorce = 0;//diem dis
								   //begin tinh hamming
					for (int i = 0; i < 20; ++i)
					{
						tempRow = 999;
						for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
						{
							tempSubRow = 0;
							for (int k = 0; k < l; k++) {
								if (temp_Str[k] != d_datainp[k + j]) tempSubRow++;
							}
							if (tempSubRow < tempRow) tempRow = tempSubRow;
						}
						temp_Sorce += tempRow;
					}
					//end tinh hamming cho chuoi vao
					//kiem tra dis motif Ne
					if (temp_Sorce < motif_NeSorce) {
						motif_NeSorce = temp_Sorce;
						for (int c = 0; c < l; ++c) {
							motif_NeString[c] = temp_Str[c];
						}
					}
					temp_Str[b] = change;
				}
				//trg hop 3 T
				if (temp_Str[b] != 3) {
					change = temp_Str[b];
					temp_Str[b] = 3;
					temp_Sorce = 0;//diem dis
								   //begin tinh hamming
					for (int i = 0; i < 20; ++i)
					{
						tempRow = 999;
						for (int j = i * 600; j < (i + 1) * 600 - l; ++j)
						{
							tempSubRow = 0;
							for (int k = 0; k < l; k++) {
								if (temp_Str[k] != d_datainp[k + j]) tempSubRow++;
							}
							if (tempSubRow < tempRow) tempRow = tempSubRow;
						}
						temp_Sorce += tempRow;
					}
					//end tinh hamming cho chuoi vao
					//kiem tra dis motif Ne
					if (temp_Sorce < motif_NeSorce) {
						motif_NeSorce = temp_Sorce;
						for (int c = 0; c < l; ++c) {
							motif_NeString[c] = temp_Str[c];
						}
					}
					temp_Str[b] = change;
				}
			}
		}//end Ne
		//end tinh

		int dem = 0;
		int res = 0;
		for (int i = 0; i < l; ++i) {
			res = res | (ansMotif_string[i] << dem);
			dem += 2;
			if (index == 574) printf("%d ", ansMotif_string[i]);
		}
		ans[index] = res;
	}
}


int main()
{
	File_Input();

	//test
	/*string test = "GTTCGGCGT";
	Motif_Ans testMoitf = dis_hamming(test);
	fo << testMoitf.dis << endl;
	cout<<sDataInp[0].substr(574, l) << endl;
	cout << h_dataMotif[574] << endl;*/
	//end test
	int h_dataOut[ARRAY_SIZE_OUT];
	for (int i = 0; i < 600; ++i) {
		h_dataOut[i] = -1;
	}
	//GPU khoi tao bien va bo nho
	int *d_dataMotif;
	if (cudaMalloc(&d_dataMotif, ARRAY_BYTES_INP) != cudaSuccess) {
		cout << "error allocating memory!" << endl;
		return 0;
	}
	int *d_dataOut;
	if (cudaMalloc(&d_dataOut, ARRAY_BYTES_OUT) != cudaSuccess) {
		cout << "error allocating memory!" << endl;
		cudaFree(d_dataMotif);
		return 0;
	}
	if (cudaMemcpy(d_dataMotif, h_dataMotif, ARRAY_BYTES_INP, cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "error copying memory!" << endl;
		cudaFree(d_dataMotif);
		cudaFree(d_dataOut);
		return 0;
	}
	if (cudaMemcpy(d_dataOut, h_dataOut, ARRAY_BYTES_OUT, cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "error copying memory!" << endl;
		cudaFree(d_dataMotif);
		cudaFree(d_dataOut);
		return 0;
	}

	cout << "dang chay ...." << endl;

	//khoi tao chay cuda
	int threadsPerBlock = 256;
	int blocksPerGrid = (600 + threadsPerBlock - 1) / threadsPerBlock;
	patternBarching <<<blocksPerGrid, threadsPerBlock >>> (d_dataMotif, l, d, d_dataOut);

	fo << "\nTime " << clock() / (double)1000 << " Sec" << endl;

	//copy data tro ve
	if (cudaMemcpy(h_dataOut, d_dataOut, ARRAY_BYTES_OUT, cudaMemcpyDeviceToHost) != cudaSuccess) {
		cout << "error copying memory!" << endl;
		cudaFree(d_dataMotif);
		cudaFree(d_dataOut);
		return 0;
	}
	//lay best motif
	cout << "\n du lieu tra ve" << endl;
	Motif_Ans best_motif,temp_motif_return;
	best_motif.dis = 999;
	for (int i = 0; i < 600; i++)
	{
		int chuyenStr = h_dataOut[i];
		int k = 0;
		string res = "";
		//cout << chuyenStr << endl;
		if (chuyenStr != -1) {
			//chuyen kieu in sang string
			for (int j = 0; j < l; ++j) {
				int temp = (chuyenStr >> k) & 3;
				//cout << temp << ' ';
				switch (temp)
				{
				case 0:
				{
					res += 'A'; break;
				}
				case 1:
				{
					res += 'C'; break;
				}
				case 2:
				{
					res += 'G'; break;
				}
				case 3:
				{
					res += 'T'; break;
				}
				}
				k += 2;
			}
			if (i == 574) fo << res << endl;
			//ket thuc chuyen
			//kiem tra do dai va tra vi tri
			temp_motif_return = dis_hamming(res);
			if (temp_motif_return.dis < best_motif.dis) {
				cout << "thay doi best" << endl;
				best_motif.dis = temp_motif_return.dis;
				best_motif.motif = temp_motif_return.motif;
				for (int z = 0; z < 20; ++z) {
					best_motif.adress[z] = temp_motif_return.adress[z];
				}
			}
			//end kiem tra
			cout << "------------" << endl;
			cout << temp_motif_return.motif << endl;
			cout << temp_motif_return.dis << endl;
			cout << best_motif.motif << endl;
			cout << best_motif.dis << endl;
			cout << "+++++++++++++" << endl;
		}
	}
	fo << "Best motif: " << best_motif.motif << endl << "Motif location: " << endl;
	for (int z = 0; z < 20; ++z) {
		fo << best_motif.adress[z] << ' ';
	}
	cout << "xong" << endl;

	cudaFree(d_dataMotif);
	cudaFree(d_dataOut);
	return 0;
}