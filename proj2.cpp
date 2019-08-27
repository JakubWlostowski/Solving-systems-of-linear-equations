#define A1A 5 + 5
#define A2A -1
#define A3A -1
#define A1C 3
#define A2C -1
#define A3C -1


#define N 959

#include<math.h>
#include<stdio.h>
#include<time.h>
#include "proj2.h"


	Matrix::Matrix(int height, int width) {
		this->h = height;
		this->w = width;
		this->numberOfIter = 0;
		matrix = new double*[h]; 
		for (int i = 0; i < h; i++)
		{
			matrix[i] = new double[w];
		}

		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
			{
				matrix[i][j] = 0;
			}
		}
	}

	Matrix::Matrix() {
		this->h = 1;
		this->w = 1;

		this->matrix = new double*[h];
		for (int i = 0; i < h; i++)
		{
			this->matrix[i] = new double[w];
		}
		this->matrix[0][0] = 0;

	}

	Matrix::~Matrix() {
		for (int i = 0; i < h; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;
	}

	int Matrix::getHeight() {
		return h;
	}

	int Matrix::getWidth() {
		return w;
	}

	void Matrix::setHeight(int newHeight) {
		this->h = newHeight;
	}

	void Matrix::setWidth(int newWidth) {
		this->w = newWidth;
	}

	void Matrix::initA(int a1, int a2, int a3) {
		for (int i = 0; i < h - 2; i++)
		{
			for (int j = 0; j < w - 2; j++)
			{
				if (i == j)
				{
					matrix[i][j] = a1;
					matrix[i + 1][j] = a2;
					matrix[i][j + 1] = a2;
					matrix[i + 2][j] = a3;
					matrix[i][j + 2] = a3;
				}
			}
		}

		matrix[h - 2][w - 2] = a1;
		matrix[h - 1][w - 1] = a1;

		matrix[h - 1][w - 2] = a2;
		matrix[h - 2][w - 1] = a2;
	}

	Matrix::Matrix(int height) {
		this->h = height;
		this->w = 1;
		this->matrix = new double*[h];
		for (int i = 0; i < h; i++)
		{
			this->matrix[i] = new double[w];
		}

		for (int i = 1; i <= h; i++) {
			matrix[i][0] = sin(i * (1 + 1));
		}
	}

	Matrix::Matrix(const Matrix& second_matrix) {
		this->h = second_matrix.h;
		this->w = second_matrix.w;
		this->numberOfIter = second_matrix.numberOfIter;

		this->matrix = new double*[h];
		for (int i = 0; i < h; i++)
		{
			this->matrix[i] = new double[w];
		}
		for (int i = 0; i < h; i++) {
			for (int j = 0; j < w; j++) {
				this->matrix[i][j] = second_matrix.matrix[i][j];
			}
		}
	}



	double Matrix::norm(Matrix another_matrix) {
		double normValue = 0;

		for (int i = 0; i < another_matrix.h; i++) {
			normValue += pow(another_matrix.matrix[i][0], 2);
		}
		return sqrt(normValue);
	}

	Matrix Matrix::diag() {
		Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.h; i++) {
			for (int j = 0; j < third_matrix.w; j++) {
				if (i != j) {
					third_matrix.matrix[i][j] = 0;
				}
			}
		}
		return third_matrix;
	}

	Matrix Matrix::tril() {
		Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.h; i++) {
			for (int j = 0; j < third_matrix.w; j++) {
				if (i <= j) {
					third_matrix.matrix[i][j] = 0;
				}
			}
		}
		return third_matrix;
	}

	Matrix Matrix::triu() {
		Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.h; i++) {
			for (int j = 0; j < third_matrix.w; j++) {
				if (i >= j) {
					third_matrix.matrix[i][j] = 0;
				}
			}
		}
		return third_matrix;
	}

	void Matrix::ones() {
		for (int i = 0; i < this->h; i++) {
			for (int j = 0; j < this->w; j++) {
					matrix[i][j] = 1;
			}
		}
	}

	Matrix Matrix::operator*(double number) {
		Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.h; i++) {
			for (int j = 0; j < third_matrix.w; j++) {
				third_matrix.matrix[i][j] = matrix[i][j] * number;
			}
		}
		return third_matrix;
	}

	Matrix Matrix::operator-(Matrix second_matrix) {
		Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.getHeight(); i++) {
			for (int j = 0; j < third_matrix.getWidth(); j++) {
				third_matrix.matrix[i][j] = matrix[i][j] - second_matrix.matrix[i][j];
			}
		}
		return third_matrix;
	}

	

	 Matrix Matrix::operator+(Matrix second_matrix) {
		 Matrix third_matrix(*this);
		for (int i = 0; i < third_matrix.getHeight(); i++) {
			for (int j = 0; j < third_matrix.getWidth(); j++) {
				third_matrix.matrix[i][j] = matrix[i][j] + second_matrix.matrix[i][j];
			}
		}
		return third_matrix;
	}

	 Matrix Matrix::operator*(Matrix second_matrix) {
		 Matrix third_matrix(this->getHeight(), second_matrix.getWidth());

		 for (int i = 0; i < third_matrix.getHeight(); i++) {
			 for (int j = 0; j < third_matrix.getWidth(); j++) {
				 for (int k = 0; k < this->getWidth(); k++) {
					 third_matrix.matrix[i][j] += matrix[i][k] * second_matrix.matrix[k][j];
				 }
				 
			 }
		 }
		 return third_matrix;
	 }

	 Matrix Matrix::operator=(Matrix second_matrix) {
		 
		 this->~Matrix();
		 this->h = second_matrix.h;
		 this->w = second_matrix.w;
		 this->numberOfIter = second_matrix.numberOfIter;
		 this->matrix = new double*[h];
		 for (int i = 0; i < this->h; i++)
		 {
			 this->matrix[i] = new double[w];
		 }


		 for (int i = 0; i < this->h; i++) {
			 for (int j = 0; j < this->w; j++) {
				this->matrix[i][j] = second_matrix.matrix[i][j];
			}
		}
		 return *this;
	 }

	 Matrix Matrix::operator/(Matrix second_matrix) {

		 Matrix third_matrix(this->getHeight(), 1);
		 third_matrix.matrix[0][0] = second_matrix.matrix[0][0] / this->matrix[0][0];

		 double sum = 0;
		 for (int i = 0; i < this->getHeight(); i++) {
			 sum = 0;
			 for (int j = 0; j < i; j++) {
				 sum += this->matrix[i][j] * third_matrix.matrix[j][0];

			 }
			 third_matrix.matrix[i][0] = (second_matrix.matrix[i][0] - sum) / this->matrix[i][i];
		 }

		 return third_matrix; 

	 }

	 double Matrix::Jacobi() {
		 Matrix D(h, w);
		 D = this->diag(); // diagonala

		 Matrix L(h, w);
		 L = this->tril(); // macierz trojkatna dolna

		 Matrix U(h, w);
		 U = this->triu(); // macierz trojkatna gorna

		 Matrix r(h, 1);
		 r.ones(); // macierz jednostkowa do wyniku

		 Matrix b(h,1); // wektor b

		 

		 for (int i = 0; i < b.h; i++) {
			 b.matrix[i][0] = sin(i * (1 + 1));
		 }
		 clock_t time;
		 time = clock();

		 while (Matrix::norm((*this * r) - b) > pow(10, -9)) {
			
			 r = (D * -1.0) / ((L + U) * r) + D / b;
			 this->numberOfIter++;
		 }
		 time = clock() - time;
		 return ((double)time / CLOCKS_PER_SEC);
	}

	 double Matrix::Gauss() {
		  Matrix D(h, w);
		  D = this->diag(); // diagonala

		  Matrix L(h, w);
		  L = this->tril(); // macierz trojkatna dolna

		  Matrix U(h, w);
		  U = this->triu(); // macierz trojkatna gorna

		  Matrix r(h, 1);
		  r.ones(); // macierz jednostkowa do wyniku

		  Matrix b(h, 1); // wektor b

		  for (int i = 0; i < b.h; i++) {
			  b.matrix[i][0] = sin(i * (1 + 1));
		  }
		  clock_t time;
		  time = clock();

		  while (Matrix::norm((*this * r) - b) > pow(10, -9)) {


			  r = ((D + L) * -1) / (U*r) + (D + L) / b;
			  this->numberOfIter++;
		  }
		  time = clock() - time;
		  return ((double)time / CLOCKS_PER_SEC);
	 }

	 double Matrix::LU() {
		 Matrix L(h, w);
		 Matrix U(h, w);
		 clock_t time;

		 time = clock();
		 for (int i = 0; i < h; i++) {
			 L.matrix[i][i] = 1;
		 }

		 for (int i = 0; i < h; i++) {
			 for (int j = i; j < w; j++) {
				 U.matrix[i][j] += this->matrix[i][j];
				 for (int k = 0; k < i; k++) {
					 U.matrix[i][j] -= L.matrix[i][k] * U.matrix[k][j];
				 }
			 }


			 for (int j = i + 1; j < h; j++) {
				 L.matrix[j][i] += this->matrix[i][j];
				 for (int k = 0; k < i; k++) {
					 L.matrix[j][i] -= L.matrix[j][k] * U.matrix[k][i];
				 }
				 L.matrix[j][i] /= U.matrix[i][i];
			 }

		 }

		 
		 // podstawienie w przód Ly = b

		 Matrix y(h, 1);

		 Matrix b(h, 1); // wektor b

		 for (int i = 0; i < h; i++) {
			 b.matrix[i][0] = sin(i * (1 + 1));
		 }
		 
		 y = L / b;

		 // podstawienie wstecz Ux = y

		 Matrix x(h, 1);
		 double sum = 0;
		 for (int i = h - 1; i >= 0; i--) {
			 sum = 0;
			 for (int j = h - 1; j >= i; j--) {
				 sum += U.matrix[i][j] * x.matrix[j][0];
			 }
			 x.matrix[i][0] = (y.matrix[i][0] - sum) / U.matrix[i][i];
		 }


		 time = clock() - time;
		 return ((double)time / CLOCKS_PER_SEC);

	 }
	 
	

int main() {
	// PODPUNKT A

	Matrix j(N, N);
	j.initA(A1A, A2A, A3A);

	// PODPUNKT B

	double time;
	Matrix g(N, N);
	g.initA(A1A, A2A, A3A);
	time = j.Jacobi();
	printf("Czas rozwiazania metoda Jacobiego: %f oraz liczba iteracji: %d \n", time, j.numberOfIter);
	time = g.Gauss();
	printf("Czas rozwiazania metoda Gaussa-Seidla: %f oraz liczba iteracji: %d \n", time, g.numberOfIter);*/

	// PODPUNKT C

	Matrix j2(N, N);
	Matrix g2(N, N);
	j2.initA(A1C, A2C, A3C);
	g2.initA(A1C, A2C, A3C);

	// PODPUNKT D

	double norm;
	Matrix metodaLU(N, N);
	metodaLU.initA(A1C, A2C, A3C);
	norm = metodaLU.LU();
	printf("Norma z residuum w przypadku rozwiazania metoda LU wynosi: %E \n", norm);

	// PODPUNKT E

	int tabN[5] = { 100, 500, 1000, 2000, 3000 };
	double timeNJ[5];
	double timeNG[5];
	double timeNLU[5];


	for (int i = 0; i < 5; i++) {
		Matrix metLU(tabN[i], tabN[i]);
		metLU.initA(A1A, A2A, A3A);
		timeNLU[i] = metLU.LU();
	}

	for (int i = 0; i < 5; i++) {
		printf("%F \n", timeNLU[i]);
	}
	
	getchar();
	return 0;
}