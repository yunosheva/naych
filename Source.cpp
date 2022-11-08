#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
using namespace std;


int h = 1;
int dt = 1;
size_t n = 40;
size_t m = 10;
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

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}
double sum = 0;
int main() {
	system("mkdir VTK");
	/* задаем начальную плотность */
	for (int j = 1; j < m + 1; j++) {
		for (int i = 1; i < n + 1; i++) {
			rho[i][j] = 1.0;
		};
	};

	for (int j = 0 ; j < m + 2; j++) {
			rho[n + 1][j] = 0.95;
	};

	for (int j = 0; j < m + 2; j++) {
		rho[0][j] = 1.05;
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

		for (size_t j = 1; j < n + 1; j++) {
			buf[2][j][0] = buf[4][j][1];
			buf[4][j][m + 1] = buf[2][j][m];
			buf[6][j + 1][0] = buf[8][j][1];
			buf[5][j - 1][0] = buf[7][j][1];
			buf[8][j - 1][m + 1] = buf[6][j][m];
			buf[7][j + 1][m + 1] = buf[5][j][m];
		}

		for (size_t j = 1; j < m + 1; j++) {
			for (size_t k = 1; k < 9; k++) {
				buf[k][0][j] = buf[k][1][j] + F_e(0, 0, w[k], rho[0][j]) - F_e(0, 0, w[k], rho[1][j]);
				buf[k][n + 1][j] = buf[k][n][j] + F_e(0, 0, w[k], rho[n + 1][j]) - F_e(0, 0, w[k], rho[n][j]);
			}
		}


		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (size_t k = 1; k < 9; k++) {
					f[k][i][j] = buf[k][i + dx[k]][j + dy[k]];
				}
			}
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
			for (int j = 1; j < m + 1 ; j++) {
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

