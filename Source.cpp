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
double g = 1e-5;


/* ����������� ������� ������������� */
double F_e(double sp, double u2, double w, double rho) {
	return w * rho* (1 + sp / teta + sp * sp / (2. * teta * teta) - u2 / 2. / teta);
}

/* ���������� ��������� ��������� */
double F(double f, double f_eq, double f_eq1) {
	return f + (f_eq - f) / tau + f_eq1 - f_eq;
}

/* ������ ��������� */

double Par(int y) {
	return -3*g/2/(tau - 1./2.)*(y - 1./2.)*(y - 1./2. - m);
}

/* ������ ��������� �������� � ��������� �������� */
vector<vector<double>> ux(n, vector<double>(m));
vector<vector<double>> uy(n, vector<double>(m));




/* ������ ��������� ��������� */
vector<vector<double>> rho(n, vector<double>(m));


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
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++) vtk_file << rho[i][j] << " ";
	vtk_file << endl;
	vtk_file << "VECTORS uflow double\n";
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++) vtk_file << ux[i][j] + g / 2 << "  " << uy[i][j] << "  0.0" << " ";
	vtk_file << endl;
	vtk_file << "VECTORS ux double\n";
		for (int i = 0; i < m; i++) 
			for (int j = 0; j < n; j++)vtk_file << Par(dt*(i+1)) <<  "  0.0  " << "  0.0" << " ";
	vtk_file << endl;

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}
double sum = 0;
int main() {
	system("mkdir VTK");
	/* ������ ��������� ��������� */
	for (int j = 0; j < m ; j++) {
		for (int i = 0; i < n; i++) {
			rho[i][j] = 1.0;
		};
	};
	/*for (int j = m / 2; j < m; j++) {
		for (int i = 0; i < n; i++) {
			rho[i][j] = 2.0;
		};
	};*/

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			sum += rho[i][j];
		};
	};

	cout << " Summa = " << sum << endl;

	/* ������ ������ ��������� ��������� ������ */
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

	/* ������ ������������ w_k*/
	vector<double> w(9);
	w[0] = 4. / 9.;
	for (int i = 1; i < 5; i++) {
		w[i] = 1. / 9.;
	}

	for (int i = 5; i < 9; i++) {
		w[i] = 1. / 36.;
	}

	/* ������ ��������� ������������� ������� f, ����������� ������� ������������� f_eq */
	vector<vector<vector<double>>> f(9, vector<vector<double>>(n + 2, vector<double>(m + 2)));
	vector<vector<vector<double>>> buf(9, vector<vector<double>>(n + 2, vector<double>(m + 2)));

	for (size_t i = 1; i < n + 1; i++) {
		for (size_t j = 1; j < m + 1; j++) {
			for (int k = 0; k < 9; k++) {
				f[k][i][j] = F_e(c[k][0] * ux[i - 1][j - 1] + c[k][1] * uy[i - 1][j - 1], ux[i - 1][j - 1] * ux[i - 1][j - 1] + uy[i - 1][j - 1] * uy[i - 1][j - 1], w[k], rho[i - 1][j - 1]);
			};
		};
	};

	vector <int> dx = { 0, -1, 0, 1, 0, -1, 1, 1, -1 };
	vector <int> dy = { 0, 0, -1, 0, 1, -1, -1, 1, 1 };

	/* ��������� ���� �� ������� */
	for (int t = 0; t < M; t++) {

		/* ����� �������� ������ */

		buf = f;

		for (int j = 1; j < m + 1; j++) {
			buf[1][0][j] = buf[1][n][j];
			buf[3][n + 1][j] = buf[3][1][j];
		}
		for (int j = 1; j < n + 1; j++) {
			buf[2][j][0] = buf[4][j][1];
			buf[4][j][m + 1] = buf[2][j][m];
			buf[6][j + 1][0] = buf[8][j][1];
			buf[5][j - 1][0] = buf[7][j][1];
			buf[8][j - 1][m + 1] = buf[6][j][m];
			buf[7][j + 1][m + 1] = buf[5][j][m];
		}

		for (size_t j = 1; j < m; j++) {
			buf[5][0][j] = buf[5][n][j];
			buf[6][n + 1][j] = buf[6][1][j];
		}

		for (size_t j = 2; j < m + 1; j++) {
			buf[7][n + 1][j] = buf[7][1][j];
			buf[8][0][j] = buf[8][n][j];
		}


		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (size_t k = 1; k < 9; k++) {
					f[k][i][j] = buf[k][i + dx[k]][j + dy[k]];
				}
			}
		};


		/* ��������� ����� ��������� */
		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				rho[i - 1][j - 1] = f[0][i][j];
				for (int k = 1; k < 9; k++) {
					rho[i - 1][j - 1] += f[k][i][j];
				};
			};
		};

		/* ��������� ����� ���������� ����� */
		full_rho = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				full_rho += rho[i][j];
			};
		};
		cout << " Density for " << t << " step = " << full_rho << endl;


		/* ��������� ����� �������� �������� � ���� */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				ux[i - 1][j - 1] = f[1][i][j] * c[1][0] / rho[i - 1][j - 1];
				uy[i - 1][j - 1] = f[1][i][j] * c[1][1] / rho[i - 1][j - 1];
				for (int k = 2; k < 9; k++) {
					ux[i - 1][j - 1] += f[k][i][j] * c[k][0] / rho[i - 1][j - 1];
					uy[i - 1][j - 1] += f[k][i][j] * c[k][1] / rho[i - 1][j - 1];
				}
			};
		};

		/* ��������� ����� ���������� �������� */
		full_rho_ux = 0;
		full_rho_uy = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				full_rho_ux += rho[i][j] * ux[i][j];
				full_rho_uy += rho[i][j] * uy[i][j];
			};
		};
		cout << " Impulse for " << t << " step (x) = " << full_rho_ux << endl;
		cout << " Impulse for " << t << " step (y) = " << full_rho_uy << endl;


		/* ������� ���� �������� */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (int k = 0; k < 9; k++) {
					f[k][i][j] = F(f[k][i][j],
						F_e(c[k][0] * ux[i - 1][j - 1] + c[k][1] * uy[i - 1][j - 1], ux[i - 1][j - 1] * ux[i - 1][j - 1] + uy[i - 1][j - 1] * uy[i - 1][j - 1], w[k], rho[i - 1][j - 1]),
						F_e(c[k][0] * (ux[i - 1][j - 1] + g) + c[k][1] * uy[i - 1][j - 1], (ux[i - 1][j - 1] + g) * (ux[i - 1][j - 1] + g) + uy[i - 1][j - 1] * uy[i - 1][j - 1], w[k], rho[i - 1][j - 1]));
				};
			};
		};

		/*if (t % 10 == 0)
		{
			SaveVTKFile(t);
		}*/
		SaveVTKFile(t);
		/*for (size_t i = 0; i < n ; i++) {
			for (size_t j = 0; j < m ; j++) {
				cout << uy[i][j] << " ";
			}
			cout << "\n";
		}*/
	}
	

	return 0;
}

