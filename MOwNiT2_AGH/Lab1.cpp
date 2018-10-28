#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

// Ax = b
int PRECISION;
bool PFLAG = true;
int PARAM_K = 7, PARAM_M = 4;


template<typename T, typename T1>
T1 euklidesNorm(T X, int n, T1 _) {
	T1 sum = 0;

	for (int i = 0; i < n; i++){
		sum += (pow(X[i], 2));
	}

	return (sqrt(sum));
}

template<typename T, typename T1>
T1 maxNorm(T A, int n, T1 _) {
	T1 max = 0;

	for (int i = 0; i < n; i++){
			if (abs(A[i]) > max) max = A[i];
	}

	return max;
}

template<typename T>
void ludist(int n, T A)
{
	int i, j, k;

	for (k = 0; k < n - 1; k++)
	{
		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}
}

template<typename T, typename T1>
void lusolve(int k, int n, T A, T X, T1 _)
{
	int    i, j;
	T1 s;

	for (i = 1; i < n; i++)
	{
		s = 0;

		for (j = 0; j < i; j++) s += A[i][j] * X[j][k];

		X[i][k] -= s;
	}

	X[n - 1][k] /= A[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--)
	{
		s = 0;

		for (j = i + 1; j < n; j++) s += A[i][j] * X[j][k];

		X[i][k] = (X[i][k] - s) / A[i][i];
	}
}

template<typename T, typename T1>
double uwarunkowanie(T Pom, T A, int n, T1 _) {
	double maxA = 0;
	double maxPom = 0;
	double sumA, sumPom;

	for (size_t i = 0; i < n; i++){
		sumA = 0;
		sumPom = 0;
		for (size_t j = 0; j < n; j++){
			sumA += abs(A[i][j]);
			sumPom += abs(Pom[i][j]);
		}
		if (sumA > maxA) maxA = sumA;
		if (sumPom > maxPom) maxPom = sumA;
	}

	return maxA * maxPom;
}

//====================================================================
//====================================================================

// wypelnienie macierzy A
template<typename T>
void matrixA1(T A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0)
				A[i][j] = 1;
			else
				A[i][j] = ((1.0 / (i + j + 2.0 - 1.0)));
		}
	}
}

template<typename T>
void matrixA2(T A, int n) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j >= i)
				A[i][j] = ((2 * i + 1.0) / (j + 1.0));
			else
				A[i][j] = A[j][i];
		}
	}
}

template<typename T>
void matrixA3a(T A, int n, int k, int m) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) A[i][j] = k;
			else if (j == i + 1) A[i][j] = (1.0 / (i + m));
			else if (j == i - 1) A[i][j] = (k / (i + m + 1.0));
			else A[i][j] = 0;
		}
	}
}

// funkcja do generowania elementow wektora X ze zbioru {-1, 1}
void matrixX(int* X, int n) {
	int pom[] = { 1, -1 };
	
	for (int i = 0; i < n; i++)
		X[i] = pom[rand() % 2];
}

//====================================================================
//====================================================================

// mnozenie macierzy A i wektora X
template<typename T, typename T1, typename T2>
void matrixMul(T A, int* X, T1 B, int n, T2 _) {
	T2 sum;

	for (int i = 0; i < n; i++) {
		sum = 0;
		for (int j = 0; j < n; j++) {
			sum += (A[i][j] * X[j]);
		}
		B[i] = sum;
	}
}


template<typename T, typename T1, typename T2>
void gaus(T A, T1 B, T1 X, int n, T2 _) {	//n * 2*n/2 * n/4 + n * n/2 = n^3/4 + n^2/2 = O(n^3)
	T2 m;

	// eliminacja wspoł
	for (int i = 0; i < n-1; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			m = -A[j][i] / A[i][i];
			for (int k = i+1; k < n; k++)
			{
				A[j][k] += m * A[i][k];
			}
			B[j] += m * B[i];
		}
	}

	T2 s;
	for (int i = n-1; i >= 0; i--)
	{
		s = B[i];
		for (int j = n - 1; j >= i + 1; j--)
		{
			s -= (A[i][j] * X[j]);
		}
		X[i] = (s / A[i][i]);
	}
}

template<typename T, typename T1, typename T2>  //t1 = float* | t = float**
void trojDiagonalnaSolver(T A, T1 B, T1 X, int n, T2 _) {
	//a_i = A[i][i]; b_i = A[i][i-1]; c_i = A[i][i+1]
	
	T2* l = new T2[n];			//1...n-1
	T2* u = new T2[n];			//0...n-1

	//pierwzsa
	u[0] = A[0][0];
	for (int i = 1; i < n; i++){
		l[i] = (A[i][i - 1] / u[i - 1]);
		u[i] = A[i][i] - l[i] * A[i - 1][i];
	}

	T2* y = new T2[n];

	//rozwiazywanie macierzy
	y[0] = B[0];
	for (int i = 1; i < n; i++){
		y[i] = B[i] - l[i]*y[i-1];
	}

	X[n - 1] = y[n - 1] / u[n - 1];
	for (int i = n-2; i >= 0; i--){
		X[i] = (y[i] - A[i][i + 1] * X[i + 1]) / u[i];
	}
}


//=============
//TRESCI ZADAN
//=============

template<typename T, typename T1, typename T2>
void zad1_2(T A, int* X, T1 X_po, T1 B, int n, T2 pomocnicza) {
	matrixMul(A, X, B, n, pomocnicza);

	auto start = chrono::high_resolution_clock::now();
	gaus(A, B, X_po, n, pomocnicza);
	auto elapsed = chrono::high_resolution_clock::now() - start;

	long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

	cout << "\nNorma Euklidesowska: " << euklidesNorm(X_po, n, pomocnicza) << "   " << euklidesNorm(X, n, pomocnicza) << endl;
	cout << "Norma Max: " << maxNorm(X_po, n, pomocnicza) << endl;
	cout << "Czas :" << microseconds << "microseconds" << endl;

}

template<typename T, typename T1, typename T2>
void zad3a(T A, int* X, T1 X_po, T1 B, int n, T2 pomocnicza) {
	zad1_2(A, X, X_po, B, n, pomocnicza);

	cout << "\n================\nMacierz trojdiagonalna\n";
	matrixA3a(A, n, PARAM_K, PARAM_M);

	auto start = chrono::high_resolution_clock::now();
	trojDiagonalnaSolver(A, B, X_po, n, pomocnicza);
	auto elapsed = chrono::high_resolution_clock::now() - start;

	long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

	cout << "\nNorma Euklidesowska: " << euklidesNorm(X_po, n, pomocnicza) << "   " << euklidesNorm(X, n, pomocnicza) << endl;
	cout << "Norma Max: " << maxNorm(X_po, n, pomocnicza) << endl;
	cout << "Czas :" << microseconds << "microseconds" << endl;
}

//==========================================================

template<typename T>
void wypelnienie_odwrotnej(T X, int n) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) X[i][j] = 0;
		X[i][i] = 1;
	}
}

template<typename T, typename T1>
void zwolnij_miejsce(T A, T Pom, T1 X_po, T1 B, int n) {
	for (size_t i = 0; i < n; i++){
		delete[] A[i];
		delete[] Pom[i];
	}

	delete[] A;
	delete[] Pom;
	delete[] X_po;
	delete[] B;
}

template<typename T, typename T1, typename T2>
void zadania(T A, T Pom, T1 X_po, int* X, T1 B, int n, T2 pomocnicza) {
	cout << "\n=======================================\nZADANIE NR 1";
	matrixA1(A, n);
	zad1_2(A, X, X_po, B, n, pomocnicza);

	matrixA1(A, n);
	ludist(n, A);
	wypelnienie_odwrotnej(Pom, n);
	for (size_t k = 0; k < n; k++)
		lusolve(k, n, A, Pom, pomocnicza);
	cout << "\nUwarunkowanie: " << uwarunkowanie(Pom, A, n, pomocnicza) << fixed << endl;

	cout << "\n=======================================\nZADANIE NR 2";
	matrixA2(A, n);
	zad1_2(A, X, X_po, B, n, pomocnicza);

	matrixA2(A, n);
	ludist(n, A);
	wypelnienie_odwrotnej(Pom, n);
	for (size_t k = 0; k < n; k++)
		lusolve(k, n, A, Pom, pomocnicza);
	cout << "\nUwarunkowanie: " << uwarunkowanie(Pom, A, n, pomocnicza) << fixed << endl;

	cout << "\n=======================================\nZADANIE NR 3a";
	matrixA3a(A, n, PARAM_K, PARAM_M);
	zad3a(A, X, X_po, B, n, pomocnicza);

	zwolnij_miejsce(A, Pom, X_po, B, n);
}

int main()
{
    freopen(".\\Outputs\\out1.txt","w",stdout);

	srand(time(NULL));
	int n;

	int N[] = { 2, 4, 8, 10, 20, 50, 100, 250, 1000 };

	for (size_t i = 0; i < 9; i++) {
		n = N[i];
		int* X = new int[n];
		matrixX(X, n);

		cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\nROZMIAR MACIERZY TO: " << n << endl;
		for (size_t j = 0; j < 2; j++)
		{
			cout << "\n##################################\nPRECYZJA: ";
			if (j == 0) {
				cout << "Double" << endl;

				double* X_po;
				double** A;
				double** Pom;
				double* B;
				double pomocnicza = 0;

				A = new double *[n];
				Pom = new double *[n];
				for (int i = 0; i < n; i++) {
					A[i] = new double[n];
					Pom[i] = new double[n];
				}

				B = new double[n];

				X_po = new double[n];
				for (int i = 0; i < n; i++)
					X_po[i] = 0;

				zadania(A, Pom, X_po, X, B, n, pomocnicza);
			}
			else {
				cout << "Float" << endl;
				float* X_po;
				float** A;
				float** Pom;
				float* B;
				float pomocnicza = 0;

				A = new float *[n];
				Pom = new float *[n];
				for (int i = 0; i < n; i++) {
					A[i] = new float[n];
					Pom[i] = new float[n];
				}

				B = new float[n];

				X_po = new float[n];
				for (int i = 0; i < n; i++)
					X_po[i] = 0;

				zadania(A, Pom, X_po, X, B, n, pomocnicza);
			}			
		}
	}
}

