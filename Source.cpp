#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
using namespace std;


size_t n = 50;
size_t m = 20;
size_t k = 20;
double teta = 1. / 3.;
int tau = 1;
size_t M = 40;
double full_rho;
double full_rho_ux;
double full_rho_uy;
double full_rho_uz;
double g = 1e-2;



/* равновесные функции распределения, sp - скалярное произведение, u2 - вектор скорости в квадрате */
double F_e(double sp, double u2, double w, double rho) {
	return w * rho* (1 + sp / teta + sp * sp / (2. * teta * teta) - u2 / 2. / teta);
}

/* решеточное уравнение Больцмана */
double F(double f, double f_eq, double f_eq1) {
	return f + (f_eq - f) / tau + f_eq1 - f_eq;
}


/* задаем начальную скорость*/
vector<vector<vector<double>>> ux(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> uy(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> uz(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));



/* задаем начальную плотность и маску*/
vector<vector<vector<double>>> rho(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> mask(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));



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
	vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << n << " " << m << " " << k;
	vtk_file << "X_COORDINATES " << n << " double\n";
	for (int i = 0; i < n; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Y_COORDINATES " << m << " double\n";
	for (int i = 0; i < m; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "Z_COORDINATES " << k << " double\n";
	for (int i = 0; i < k; i++) vtk_file << i << " ";
	vtk_file << endl;
	vtk_file << "POINT_DATA " << n * m * k << endl;

	vtk_file << "SCALARS rho double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < k + 1; l++)
		for (int j = 1; j < m + 1; j++)
			for (int i = 1; i < n + 1; i++)
				vtk_file << rho[i][j][l] << " ";
	vtk_file << endl;

	vtk_file << "VECTORS uflow double\n";
	for (int l = 1; l < k + 1; l++)
		for (int j = 1; j < m + 1; j++)
			for (int i = 1; i < n + 1; i++)
				vtk_file << ux[i][j][l] + g/2 << " " << uy[i][j][l] << " " << uz[i][j][l] << " ";
	vtk_file << endl;

	/*vtk_file << "SCALARS mask double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++)
			for (int l = 1; l < k + 1; l++) vtk_file << mask[i][j][l] << " ";
	vtk_file << endl;*/

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}
double sum = 0;
int main() {
	system("mkdir VTK");
	/* задаем начальную плотность */
	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++)
			for (int l = 1; l < k / 2; l++)  rho[i][j][l] = 1;

	for (int j = 1; j < m + 1; j++)
		for (int i = 1; i < n + 1; i++)
			for (int l = k / 2; l < k + 1; l++)  rho[i][j][l] = 1;

	for (int j = 0; j < m + 2; j++)
		for (int i = 0; i < n + 2; i++)
			for (int l = 0; l < k + 2; l++) {
				sum += rho[i][j][l];
				ux[i][j][l] = uy[i][j][l] = uz[i][j][l] = 0.0;
			}

	cout << " Summa = " << sum << endl;

	/* задаем вектор возможных скоростей частиц */
	vector<vector<int>> c(19, vector<int>(3));
	vector <int> dx = { 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0};
	vector <int> dy = { 0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, -1, 1};
	vector <int> dz = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1,-1, -1, 1, 1, -1, -1};

	for (int j = 1; j < 19; j++) {
		c[j][0] = dx[j];
		c[j][1] = dy[j];
		c[j][2] = dz[j];
	}

	/* задаем коэффициенты w_k*/
	vector<double> w(19);

	w[0] = 1. / 3.;

	for (int i = 1; i < 7; i++) {
		w[i] = 1. / 18.;
	}

	for (int i = 7; i < 19; i++) {
		w[i] = 1. / 36.;
	}


	/* задаем начальные одночастичные функции f, равновесные функции распределения f_eq */
	vector<vector<vector<vector<double>>>> f(19, vector<vector<vector<double>>>(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2))));
	vector<vector<vector<vector<double>>>> buf(19, vector<vector<vector<double>>>(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2))));


	for (size_t i = 1; i < n + 1; i++) {
		for (size_t j = 1; j < m + 1; j++) {
			for (size_t l = 1; l < k + 1; l++) {
				for (size_t s = 0; s < 19; s++) {
					f[s][i][j][l] = F_e(c[s][0] * ux[i][j][l] + c[s][1] * uy[i][j][l] + c[s][2] * uz[i][j][l]
						, ux[i][j][l] * ux[i][j][l] + uy[i][j][l] * uy[i][j][l] + uz[i][j][l] * uz[i][j][l], w[s], rho[i][j][l]);
				}
			};
		};
	};


	/* запускаем цикл по времени */
	for (int t = 0; t < M; t++) {

		/* учтем движение частиц */

		buf = f;

		for (size_t i = 1; i < k + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				buf[1][0][j][i] = buf[1][n][j][i];
				buf[2][n + 1][j][i] = buf[2][1][j][i];
			};
		};
		/* wall 
		for (size_t i = 1; i < k + 1; i++) {
			for (size_t j = 1; j < n + 1; j++) {
				buf[3][j][0][i] = buf[4][j][1][i];
				buf[4][j][m + 1][i] = buf[3][j][m][i];
			};
		};*/

		for (size_t i = 1; i < k + 1; i++) {
			for (size_t j = 1; j < n + 1; j++) {
				buf[3][j][0][i] = buf[3][j][m][i];
				buf[4][j][m + 1][i] = buf[4][j][1][i];
			};
		}
		
		/* wall */
		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				buf[5][i][j][0] = buf[6][i][j][1];
				buf[6][i][j][k + 1] = buf[5][i][j][k];
			};
		};

		/*
		for (size_t i = 0; i < k + 1; i++) {
			buf[7][0][0][i] = buf[7][n][m][i];
		}*/
		
		for (size_t i = 0; i < k + 1; i++) {
			buf[7][0][0][i] = buf[7][n][m][i];
		}
		for (size_t j = 1; j < m; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[7][0][j][l] = buf[7][n][j][l];
			}
		}
		/* wall 
		for (size_t j = 0; j < n ; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[7][j][0][l] = buf[9][j + 1][1][l];
			}
		}*/

		for (size_t j = 1; j < n; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[7][j][0][l] = buf[7][j][m][l];
			}
		}

		
		for (size_t l = 0; l < k + 1; l++) {
			buf[8][n + 1][0][l] = buf[8][1][m][l];
		}
		for (size_t j = 1; j < m; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[8][n + 1][j][l] = buf[8][1][j][l];
			}
		}
		/* wall 
		for (size_t j = 2; j < n + 2; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[8][j][0][l] = buf[10][j - 1][1][l];
			}
		}*/
		for (size_t j = 2; j < n + 1; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[8][j][0][l] = buf[8][j][m][l];
			}
		}

		
		for (size_t l = 0; l < k + 1; l++) { 
			buf[9][n + 1][m + 1][l] = buf[9][1][1][l]; 
		}
		for (size_t j = 2; j < m + 1; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[9][n + 1][j][l] = buf[9][1][j][l];
			}
		}
		/* wall 
		for (size_t j = 2; j < n + 2; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[9][j][m + 1][l] = buf[7][j - 1][m][l];
			}
		}*/

		for (size_t j = 2; j < n + 1; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[9][j][m + 1][l] = buf[9][j][1][l];
			}
		}

		
		for (size_t l = 0; l < k + 1; l++) {
			buf[10][0][m + 1][l] = buf[10][n][1][l];
		}
		for (size_t j = 2; j < m + 1; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[10][0][j][l] = buf[10][n][j][l];
			}
		}
		/* wall 
		for (size_t j = 0; j < n; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[10][j][m + 1][l] = buf[8][j + 1][m][l];
			}
		}*/

		for (size_t j = 1; j < n; j++) {
			for (size_t l = 0; l < k + 1; l++) {
				buf[10][j][m + 1][l] = buf[10][j][1][l];
			}
		}

		/* не забыть убрать
		for (size_t i = 0; i < m + 1; i++) {
			buf[11][0][i][0] = buf[11][n][i][k];
		}*/
		
		for (size_t j = 1; j < k; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[11][0][l][j] = buf[11][n][l][j];
			}
		}
		/* wall */
		for (size_t j = 0; j < n; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[11][j][l][0] = buf[13][j + 1][l][1];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < m + 1; l++) {
			buf[12][n + 1][l][0] = buf[12][1][l][k];
		}*/
		for (size_t j = 1; j < k; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[12][n + 1][l][j] = buf[12][1][l][j];
			}
		}
		/* wall */
		for (size_t j = 2; j < n + 2; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[12][j][l][0] = buf[14][j - 1][l][1];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < m + 1; l++) {
			buf[13][n + 1][l][k + 1] = buf[13][1][l][1];
		}*/
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[13][n + 1][l][j] = buf[13][1][l][j];
			}
		}
		/* wall */
		for (size_t j = 2; j < n + 2; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[13][j][l][k + 1] = buf[11][j - 1][l][k];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < m + 1; l++) {
			buf[14][0][l][k + 1] = buf[14][n][l][1];
		}*/
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[14][0][l][j] = buf[14][n][l][j];
			}
		}
		/* wall */
		for (size_t j = 0; j < n; j++) {
			for (size_t l = 0; l < m + 1; l++) {
				buf[14][j][l][k + 1] = buf[12][j + 1][l][k];
			}
		}

		/* не забыть убрать 
		for (size_t i = 0; i < n + 1; i++) {
			buf[15][i][0][0] = buf[15][i][m][k];
		}*/
		for (size_t j = 1; j < k; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[15][l][0][j] = buf[15][l][m][j];
			}
		}
		/* wall */
		for (size_t j = 0; j < m; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[15][l][j][0] = buf[17][l][j + 1][1];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < n + 1; l++) {
			buf[16][l][m + 1][0] = buf[16][l][1][k];
		}*/
		for (size_t j = 1; j < k; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[16][l][m + 1][j] = buf[16][l][1][j];
			}
		}
		/* wall */
		for (size_t j = 2; j < m + 2; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[16][l][j][0] = buf[18][l][j - 1][1];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < n + 1; l++) {
			buf[17][l][m + 1][k + 1] = buf[17][l][1][1];
		}*/
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[17][l][m + 1][j] = buf[17][l][1][j];
			}
		}
		/* wall */
		for (size_t j = 2; j < m + 2; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[17][l][j][k + 1] = buf[15][l][j - 1][k];
			}
		}

		/* не забыть убрать 
		for (size_t l = 0; l < n + 1; l++) {
			buf[18][l][0][k + 1] = buf[18][l][m][1];
		}*/
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[18][l][0][j] = buf[18][l][m][j];
			}
		}
		/* wall */
		for (size_t j = 0; j < m; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[18][l][j][k + 1] = buf[16][l][j + 1][k];
			}
		}

		/*
		for (size_t j = 1; j < k; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[16][l][m + 1][j] = buf[16][l][1][j];
			}
		}
		for (size_t j = 2; j < m + 2; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[16][l][j][0] = buf[18][l][j - 1][1];
			}
		}
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[17][l][m + 1][j] = buf[17][l][1][j];
			}
		}
		for (size_t j = 2; j < m + 2; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[17][l][j][k + 1] = buf[15][l][j - 1][k];
			}
		}
		for (size_t j = 2; j < k + 1; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[18][l][0][j] = buf[18][l][m][j];
			}
		}
		for (size_t j = 0; j < m; j++) {
			for (size_t l = 0; l < n + 1; l++) {
				buf[18][l][j][k + 1] = buf[16][l][j + 1][k];
			}
		}*/

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (size_t l = 1; l < k + 1; l++) {
					for (size_t s = 1; s < 19 ; s++) {
						f[s][i][j][l] = buf[s][i - dx[s]][j - dy[s]][l - dz[s]];
					}
				}
			}
		};

		/* посчитаем новую плотность */

		for (int i = 1; i < n + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				for (int l = 1; l < k + 1; l++) {
					rho[i][j][l] = f[0][i][j][l];
					for (int s = 1; s < 19; s++) {
						rho[i][j][l] += f[s][i][j][l];
					};
				}
			};
		};

		/* проверяем закон сохранения массы */
		full_rho = 0;
		for (int i = 1; i < n + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				for (int l = 1; l < k + 1; l++) {
					full_rho += rho[i][j][l];
				}
			};
		};
		cout << " Density for " << t << " step = " << full_rho << endl;


		/* посчитаем новую скорость вещества в узле */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (size_t l = 1; l < k + 1; l++){
					ux[i][j][l] = f[1][i][j][l] * c[1][0] / rho[i][j][l];
					uy[i][j][l] = f[1][i][j][l] * c[1][1] / rho[i][j][l];
					uz[i][j][l] = f[1][i][j][l] * c[1][2] / rho[i][j][l];
					for (int s = 2; s < 19; s++) {
						ux[i][j][l] += f[s][i][j][l] * c[s][0] / rho[i][j][l];
						uy[i][j][l] += f[s][i][j][l] * c[s][1] / rho[i][j][l];
						uz[i][j][l] += f[s][i][j][l] * c[s][2] / rho[i][j][l];
					}
				}
			};
		};

		/* проверяем закон сохранения импульса */
		full_rho_ux = 0;
		full_rho_uy = 0;
		full_rho_uz = 0;
		for (int i = 1; i < n + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				for (int l = 1; l < k + 1; l++) {
					full_rho_ux += rho[i][j][l] * ux[i][j][l];
					full_rho_uy += rho[i][j][l] * uy[i][j][l];
					full_rho_uz += rho[i][j][l] * uz[i][j][l];
				}
			};
		};
		cout << " Impulse for " << t << " step (x) = " << full_rho_ux << endl;
		cout << " Impulse for " << t << " step (y) = " << full_rho_uy << endl;
		cout << " Impulse for " << t << " step (z) = " << full_rho_uz << endl;


		/* сделаем одну итерацию */

		for (size_t i = 1; i < n + 1; i++) {
			for (size_t j = 1; j < m + 1; j++) {
				for (size_t l = 1; l < k + 1; l++) {
					for (size_t s = 0; s < 19; s++) {
						f[s][i][j][l] = F(f[s][i][j][l],
						F_e(c[s][0] * ux[i][j][l] + c[s][1] * uy[i][j][l] + c[s][2] * uz[i][j][l]
							, ux[i][j][l] * ux[i][j][l] + uy[i][j][l] * uy[i][j][l] + uz[i][j][l] * uz[i][j][l], w[s], rho[i][j][l]),
						F_e(c[s][0] * (ux[i][j][l] + g) + c[s][1] * uy[i][j][l] + c[s][2] * uz[i][j][l]
							, (ux[i][j][l] + g) * (ux[i][j][l] + g) + uy[i][j][l] * uy[i][j][l] + uz[i][j][l] * uz[i][j][l], w[s], rho[i][j][l]));
					}
				}
			}
		};

		SaveVTKFile(t);

	}


	return 0;
}
