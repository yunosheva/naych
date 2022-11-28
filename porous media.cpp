#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
using namespace std;


int h = 1;
int dt = 1;
size_t n = 100;
size_t m = 30;
double teta = 1. / 3.;
int tau = 1;
size_t M = 300;
double full_rho;
double full_rho_ux;
double full_rho_uy;
double g = 0;


/* равновесные функции распределения, sp - скалярное произведение, u2 - вектор в квадрате */
double F_e(double sp, double u2, double w, double rho) {
	return w * rho* (1 + sp / teta + sp * sp / (2. * teta * teta) - u2 / 2. / teta);
}

/* решеточное уравнение Больцмана */
double F(double f, double f_eq, double f_eq1) {
	return f + (f_eq - f) / tau + f_eq1 - f_eq;
}


/* задаем начальную скорость и изменение скорости */
vector<vector<double>> ux(n + 2, vector<double>(m + 2));
vector<vector<double>> uy(n + 2, vector<double>(m + 2));


/* задаем начальную плотность */
vector<vector<double>> rho(n + 2, vector<double>(m + 2));
vector<vector<double>> mask(n + 2, vector<double>(m + 2));


void SaveVTKFile(int tStep)
{
	stringstream fname;
	fname << "VTK/adv_";
	if (tStep < 10) fname << "0";
	if (tStep < 100) fname << "0";
	if (tStep < 1000) fname << "0";
	if (tStep < 10000) fname << "0";
	if (tStep < 100000) fname << "0";
	if (tStep < 1000000) fname << "0";
	if (tStep < 10000000) fname << "0";
	fname << tStep << ".vtk";
	ofstream vtk_file(fname.str().c_str());
	vtk_file << "# vtk DataFile Version 3.0\n";
	vtk_file << "Immiscible displacement\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << n << " " << m << " 1\n";
	vtk_file << "X_COORDINATES " << n << " double\n";
	for (int i = 0; i < n; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Y_COORDINATES " << m << " double\n";
	for (int i = 0; i < m; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Z_COORDINATES 1 double\n0\n";
	vtk_file << "POINT_DATA " << n * m << endl;
	vtk_file << "SCALARS rho double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++) vtk_file << rho[i][j] << " ";
	vtk_file << endl;
	vtk_file << "VECTORS uflow double\n";
	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++) vtk_file << ux[i][j] + g / 2 << "  " << uy[i][j] << "  0.0" << " ";
	vtk_file << endl;

	vtk_file << "SCALARS mask double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++) vtk_file << mask[i][j] << " ";
	vtk_file << endl;

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}
double sum = 0;
int f1() {
	system("mkdir VTK");
	/* задаем начальную плотность */
	for (int j = 1; j < m + 1; j++) {
		for (int i = 1; i < n + 1; i++) {
			rho[i][j] = 1.0;
		};
	};

	for (int j = 0; j < m + 2; j++) {
		rho[n + 1][j] = 0.9;
	};

	for (int j = 0; j < m + 2; j++) {
		rho[0][j] = 1.1;
	};

	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < m + 1; j++) {
			sum += rho[i][j];
		};
	};

	cout << " Summa = " << sum << endl;

	/* задаем вектор возможных скоростей частиц */
	vector<vector<double>> c(9, vector<double>(2));

	c[0][0] = 0;
	c[0][1] = 0;

	c[1][0] = 1;
	c[1][1] = 0;

	c[2][0] = 0;
	c[2][1] = 1;

	c[3][0] = -1;
	c[3][1] = 0;

	c[4][0] = 0;
	c[4][1] = -1;

	c[5][0] = 1;
	c[5][1] = 1;

	c[6][0] = -1;
	c[6][1] = 1;

	c[7][0] = -1;
	c[7][1] = -1;

	c[8][0] = 1;
	c[8][1] = -1;

	/* задаем коэффициенты w_k*/
	vector<double> w(9);
	w[0] = 4. / 9.;
	for (int i = 1; i < 5; i++) {
		w[i] = 1. / 9.;
	}

	for (int i = 5; i < 9; i++) {
		w[i] = 1. / 36.;
	}

	/* маска для стенок */
	for (int j = 0; j <= m + 1; j++) {
		for (int i = 0; i <= n + 1; i++) {
			mask[i][j] = 0.0;
		};
	};
	mask[49][14] = 1;
	mask[50][14] = 1;
	for (int i = 48; i <= 51; i++) {
		mask[i][15] = 1;
		mask[i][16] = 1;
	}
	mask[49][17] = 1;
	mask[50][17] = 1;

	/* задаем начальные одночастичные функции f, равновесные функции распределения f_eq */
	vector<vector<vector<double>>> f(9, vector<vector<double>>(n + 2, vector<double>(m + 2)));
	vector<vector<vector<double>>> buf(9, vector<vector<double>>(n + 2, vector<double>(m + 2)));


	for (size_t i = 1; i < n + 1; i++) {
		for (size_t j = 1; j < m + 1; j++) {
			for (int k = 0; k < 9; k++) {
				f[k][i][j] = F_e(c[k][0] * ux[i][j] + c[k][1] * uy[i][j], ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j], w[k], rho[i][j]);
			};
		};
	};


	vector <int> dx = { 0, -1, 0, 1, 0, -1, 1, 1, -1 };
	vector <int> dy = { 0, 0, -1, 0, 1, -1, -1, 1, 1 };

	/* запускаем цикл по времени */
	for (int t = 0; t < M; t++) {

		/* учтем движение частиц */

		buf = f;

		/* твердые стенки */
		for (size_t j = 1; j < n + 1; j++) {
			buf[2][j][0] = buf[4][j][1];
			buf[4][j][m + 1] = buf[2][j][m];
			buf[6][j + 1][0] = buf[8][j][1];
			buf[5][j - 1][0] = buf[7][j][1];
			buf[8][j - 1][m + 1] = buf[6][j][m];
			buf[7][j + 1][m + 1] = buf[5][j][m];
		}

		/* учет задания градиента давлений */
		for (size_t j = 1; j < m + 1; j++) {
			for (size_t k = 1; k < 9; k++) {
				buf[k][0][j] = buf[k][1][j] + F_e(0, 0, w[k], rho[0][j]) - F_e(0, 0, w[k], rho[1][j]);
				buf[k][n + 1][j] = buf[k][n][j] + F_e(0, 0, w[k], rho[n + 1][j]) - F_e(0, 0, w[k], rho[n][j]);
			}
		}

		/* отражение от квадратика */
		/*for (size_t i = 5; i < 7; i++) {
			//нижняя стенка
			buf[4][i][5] = buf[2][i][4];
			buf[7][i][5] = buf[5][i - 1][4];
			buf[8][i][5] = buf[6][i + 1][4];

			// правая стенка
			buf[1][6][i] = buf[3][7][i];
			buf[5][6][i] = buf[7][7][i + 1];
			buf[8][6][i] = buf[6][7][i - 1];

			// верхняя стенка
			buf[2][i][6] = buf[4][i][7];
			buf[5][i][6] = buf[7][i + 1][7];
			buf[6][i][6] = buf[8][i - 1][7];

			// левая стенка 
			buf[3][5][i] = buf[1][4][i];
			buf[7][5][i] = buf[5][4][i - 1];
			buf[6][5][i] = buf[8][4][i + 1];
		}*/
		

		/* нижний квадратик */
		// правая стенка
		buf[1][50][14] = buf[3][51][14];
		buf[8][50][15] = buf[6][51][14];

		// левая стенка 
		buf[3][49][14] = buf[1][48][14];
		buf[7][49][15] = buf[5][48][14];
	
		for (size_t i = 49; i < 51; i++) {
			//нижняя стенка
			buf[4][i][14] = buf[2][i][13];
			buf[7][i][14] = buf[5][i - 1][13];
			buf[8][i][14] = buf[6][i + 1][13];
		}

		/* правый квадратик */
		// верхняя стенка
		buf[2][51][16] = buf[4][51][17];
		buf[5][50][16] = buf[7][51][17];

		// нижняя стенка 
		buf[4][51][15] = buf[2][51][14];
		buf[8][50][15] = buf[6][51][14];

		for (size_t i = 15; i < 17; i++) {
			//левая стенка
			buf[1][51][i] = buf[3][52][i];
			buf[5][51][i] = buf[7][52][i + 1];
			buf[8][51][i] = buf[6][52][i - 1];
		}

		/* верхний квадратик */
		// правая стенка
		buf[1][50][17] = buf[3][51][17];
		buf[5][50][16] = buf[7][51][17];

		// левая стенка 
		buf[3][49][16] = buf[1][48][16];
		buf[6][49][16] = buf[8][48][17];

		for (size_t i = 49; i < 51; i++) {
			// верхняя стенка
			buf[2][i][17] = buf[4][i][18];
			buf[5][i][17] = buf[7][i + 1][18];
			buf[6][i][17] = buf[8][i - 1][18];
		}

		/* левый квадратик */
		// верхняя стенка
		buf[2][48][16] = buf[4][48][17];

		// нижняя стенка 
		buf[4][48][15] = buf[2][48][14];

		for (size_t i = 15; i < 17; i++) {
			//левая стенка
			buf[3][48][i] = buf[1][47][i];
			buf[7][48][i] = buf[5][47][i - 1];
			buf[6][48][i] = buf[8][47][i + 1];
		}
		


		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1 ; j++) {
				for (size_t k = 1; k < 9; k++) {
					if (mask[i][j] == 0) {
						f[k][i][j] = buf[k][i + dx[k]][j + dy[k]];
					}
				}
			}
		};

		/*for (int k = 1; k < 9; k++) {
			
			f[k][5][5] = buff[k][0];
			f[k][5][6] = buff[k][1];
			f[k][6][5] = buff[k][2];
			f[k][6][6] = buff[k][3];
		};

		/* посчитаем новую плотность */
		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				rho[i][j] = f[0][i][j];
				for (int k = 1; k < 9; k++) {
					rho[i][j] += f[k][i][j];
				};
			};
		};

		/* проверяем закон сохранения массы */
		full_rho = 0;
		for (int i = 1; i < n + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				full_rho += rho[i][j];
			};
		};
		cout << " Density for " << t << " step = " << full_rho << endl;


		/* посчитаем новую скорость вещества в узле */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				ux[i][j] = f[1][i][j] * c[1][0] / rho[i][j];
				uy[i][j] = f[1][i][j] * c[1][1] / rho[i][j];
				for (int k = 2; k < 9; k++) {
					ux[i][j] += f[k][i][j] * c[k][0] / rho[i][j];
					uy[i][j] += f[k][i][j] * c[k][1] / rho[i][j];
				}
			};
		};

		/* проверяем закон сохранения импульса */
		full_rho_ux = 0;
		full_rho_uy = 0;
		for (int i = 1; i < n + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				full_rho_ux += rho[i][j] * ux[i][j];
				full_rho_uy += rho[i][j] * uy[i][j];
			};
		};
		cout << " Impulse for " << t << " step (x) = " << full_rho_ux << endl;
		cout << " Impulse for " << t << " step (y) = " << full_rho_uy << endl;


		/* сделаем одну итерацию */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (int k = 0; k < 9; k++) {
					f[k][i][j] = F(f[k][i][j],
						F_e(c[k][0] * ux[i][j] + c[k][1] * uy[i][j], ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j], w[k], rho[i][j]),
						F_e(c[k][0] * (ux[i][j] + g) + c[k][1] * uy[i][j], (ux[i][j] + g) * (ux[i][j] + g) + uy[i][j] * uy[i][j], w[k], rho[i][j]));
				};
			};
		};

		SaveVTKFile(t);

	}


	return 0;
}

int f2();
int main() {
	int choice;
	std::cin >> choice;

	switch (choice) {
	case 1: f1(); break;
	case 2: f2(); break;
	}
}