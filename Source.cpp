#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>
using namespace std;


size_t n = 30;
size_t m = 20;
size_t k = 20;
double teta = 1. / 3.;
int tau = 1;
size_t M = 2000;
double full_rho;
double full_rho_ux;
double full_rho_uy;
double full_rho_uz;
double g = 1e-3; /* добавка к скорости */
double k_koeff = 0.01;
double A = -0.0456;
double T = 0.8; /* температура */
double omega = 0.040; /* ацентрический фактор из википедии для azota */
double sum = 0;
double sumG = 0;



/* равновесные функции распределения, sp - скалярное произведение, u2 - вектор скорости в квадрате */
double F_e(double sp, double u2, double w, double rho) {
	return w * rho* (1 + sp / teta + sp * sp / (2. * teta * teta) - u2 / 2. / teta);
}

/* решеточное уравнение Больцмана */
double F(double f, double f_eq, double f_eq1, double f_eq2) {
	return f + (f_eq - f) / tau + (f_eq1 - f_eq) + (f_eq2 - f_eq);
}

/* уравнение состояния Пенга-Робинсона */
double a(const double temperature, double omega){
	double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
	double a = pow((1 + m * (1 - sqrt(temperature))), 2);
	return a;
}

double PressurePengRobinson(double rho, const double temperature, double omega){
	double pressure = 1 / 0.307 * (temperature / (1. / rho - 0.253) -
		1.487 * a(temperature, omega) / (1. / rho / rho + 2 * 0.253 / rho - 0.253 * 0.253));
	return pressure;
}


/* задаем начальную скорость и изменение скорости*/
vector<vector<vector<double>>> ux(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> uy(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> uz(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> dux(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> duy(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> duz(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));



/* задаем начальную плотность, "эффективную" плотность и маску*/
vector<vector<vector<double>>> rho(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> mask(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));
vector<vector<vector<double>>> Fi(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2)));


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
				vtk_file << ux[i][j][l] + g / 2 << " " << uy[i][j][l] << " " << uz[i][j][l] << " ";
	vtk_file << endl;

	vtk_file << "SCALARS mask double 1\n";
	vtk_file << "LOOKUP_TABLE default\n";
	for (int l = 1; l < k + 1; l++)
		for (int j = 1; j < m + 1; j++)
			for (int i = 1; i < n + 1; i++)
				vtk_file << mask[i][j][l] << " ";
	vtk_file << endl;

	vtk_file.close();

	cout << endl << "File " << fname.str() << " written" << endl << endl;
}


int main() {
	system("mkdir VTK");
	  
		/* задаем начальную плотность */

		for (int j = 1; j < m + 1; j++)
			for (int l = 1; l < k + 1; l++)
				for (int i = n / 2 - 5; i < n / 2 + 5; i++)
					rho[i][j][l] = 2.5335;


		for (int j = 1; j < m + 1; j++)
			for (int l = 1; l < k + 1; l++)
				for (int i = 1; i < n / 2 - 8; i++)
					rho[i][j][l] = 0.09638;

		for (int j = 1; j < m + 1; j++)
			for (int l = 1; l < k + 1; l++)
				for (int i = n / 2 + 8; i < n + 1; i++)
					rho[i][j][l] = 0.09638;

		for (int j = 1; j < m + 1; j++)
			for (int l = 1; l < k + 1; l++)
				for (int i = n / 2 + 5; i < n / 2 + 8; i++)
					rho[i][j][l] = 0.09638 + i * 0.03;

		for (int j = 1; j < m + 1; j++)
			for (int l = 1; l < k + 1; l++)
				for (int i = n / 2 - 8; i < n / 2 - 5; i++)
					rho[i][j][l] = 0.09638 + i * 0.03;

		for (int j = 0; j < m + 2; j++) {
			for (int i = 0; i < n + 2; i++) {
				for (int l = 0; l < k + 2; l++) {
					sum += rho[i][j][l];
					ux[i][j][l] = uy[i][j][l] = uz[i][j][l] = dux[i][j][l] = duy[i][j][l] = duz[i][j][l] = 0.0;
				}
			}
		}
		cout << " Summa = " << sum << endl;

		/* маска */
		for (int j = 1; j < m + 1; j++) {
			for (int i = 1; i < n + 1; i++) {
				for (int l = 1; l < k + 1; l++) {
					if (sqrt((i - 20)* (i - 20) + (j - 10)* (j - 10) + (l - 10)* (l - 10)) <= 5)
						mask[i][j][l] = 1.0;
					else mask[i][j][l] = 0.0;
				}
			}
		}


		/* задаем вектор возможных скоростей частиц */
		vector<vector<int>> c(19, vector<int>(3));
		vector <int> dx = { 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0 };
		vector <int> dy = { 0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, -1, 1 };
		vector <int> dz = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1,-1, -1, 1, 1, -1, -1 };

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

		/* задаем коэффициенты G_k*/
		vector<double> G(19);

		for (int i = 0; i < 7; i++) {
			G[i] = 1.;
		}

		for (int i = 7; i < 19; i++) {
			G[i] = 1. / 2.;
		}



		/* задаем начальные одночастичные функции f, равновесные функции распределения f_eq */
		vector<vector<vector<vector<double>>>> f(19, vector<vector<vector<double>>>(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2))));
		vector<vector<vector<vector<double>>>> buf(19, vector<vector<vector<double>>>(n + 2, vector<vector<double>>(m + 2, vector<double>(k + 2))));

#pragma omp parallel for
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
		for (int t = 0; t <= M ; t++) {

			/* учтем движение частиц */

			buf = f;

			for (size_t i = 1; i < k + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					buf[1][0][j][i] = buf[1][n][j][i];
					buf[2][n + 1][j][i] = buf[2][1][j][i];
				};
			};
			/* wall */
			for (size_t i = 1; i < k + 1; i++) {
				for (size_t j = 1; j < n + 1; j++) {
					buf[3][j][0][i] = buf[4][j][1][i];
					buf[4][j][m + 1][i] = buf[3][j][m][i];
				};
			};

			/* wall */
			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					buf[5][i][j][0] = buf[6][i][j][1];
					buf[6][i][j][k + 1] = buf[5][i][j][k];
				};
			};


			/* wall */
			for (size_t j = 0; j < n; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[7][j][0][l] = buf[9][j + 1][1][l];
				}
			}

			for (size_t j = 1; j < m; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[7][0][j][l] = buf[7][n][j][l];
				}
			}


			/* wall */
			for (size_t j = 2; j < n + 2; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[8][j][0][l] = buf[10][j - 1][1][l];
				}
			}
			for (size_t j = 1; j < m; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[8][n + 1][j][l] = buf[8][1][j][l];
				}
			}


			/* wall */
			for (size_t j = 2; j < n + 2; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[9][j][m + 1][l] = buf[7][j - 1][m][l];
				}
			}

			for (size_t j = 2; j < m + 1; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[9][n + 1][j][l] = buf[9][1][j][l];
				}
			}


			/* wall */
			for (size_t j = 0; j < n; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[10][j][m + 1][l] = buf[8][j + 1][m][l];
				}
			}
			for (size_t j = 2; j < m + 1; j++) {
				for (size_t l = 0; l < k + 1; l++) {
					buf[10][0][j][l] = buf[10][n][j][l];
				}
			}

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

			/* wall */
			for (size_t j = 0; j < k; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[15][l][0][j] = buf[17][l][1][j + 1];
				}
			}
			/* wall */
			for (size_t j = 0; j < m; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[15][l][j][0] = buf[17][l][j + 1][1];
				}
			}


			/* wall */
			for (size_t j = 0; j < k; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[16][l][m + 1][j] = buf[18][l][m][j + 1];
				}
			}
			/* wall */
			for (size_t j = 2; j < m + 2; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[16][l][j][0] = buf[18][l][j - 1][1];
				}
			}


			/* wall */
			for (size_t j = 2; j < k + 2; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[17][l][m + 1][j] = buf[15][l][m][j - 1];
				}
			}
			/* wall */
			for (size_t j = 2; j < m + 2; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[17][l][j][k + 1] = buf[15][l][j - 1][k];
				}
			}


			/* wall */
			for (size_t j = 2; j < k + 2; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[18][l][0][j] = buf[16][l][1][j - 1];
				}
			}
			/* wall */
			for (size_t j = 0; j < m; j++) {
				for (size_t l = 0; l < n + 1; l++) {
					buf[18][l][j][k + 1] = buf[16][l][j + 1][k];
				}
			}

			/* препятствия */
			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					for (size_t l = 1; l < k + 1; l++) {
						if (mask[i][j][l] == 1.0) {

							if (mask[i + 1][j][l] != 1.0)
								buf[1][i][j][l] = buf[2][i + 1][j][l];

							if (mask[i + 1][j + 1][l] != 1.0)
								buf[7][i][j][l] = buf[9][i + 1][j + 1][l];

							if (mask[i][j + 1][l] != 1.0)
								buf[3][i][j][l] = buf[4][i][j + 1][l];

							if (mask[i - 1][j + 1][l] != 1.0)
								buf[8][i][j][l] = buf[10][i - 1][j + 1][l];

							if (mask[i - 1][j][l] != 1.0)
								buf[2][i][j][l] = buf[1][i - 1][j][l];

							if (mask[i - 1][j - 1][l] != 1.0)
								buf[9][i][j][l] = buf[7][i - 1][j - 1][l];

							if (mask[i][j - 1][l] != 1.0)
								buf[4][i][j][l] = buf[3][i][j - 1][l];

							if (mask[i + 1][j - 1][l] != 1.0)
								buf[10][i][j][l] = buf[8][i + 1][j - 1][l];

							if (mask[i + 1][j][l + 1] != 1.0)
								buf[11][i][j][l] = buf[13][i + 1][j][l + 1];

							if (mask[i][j][l + 1] != 1.0)
								buf[5][i][j][l] = buf[6][i][j][l + 1];

							if (mask[i - 1][j][l + 1] != 1.0)
								buf[12][i][j][l] = buf[14][i - 1][j][l + 1];

							if (mask[i - 1][j][l - 1] != 1.0)
								buf[13][i][j][l] = buf[11][i - 1][j][l - 1];

							if (mask[i][j][l - 1] != 1.0)
								buf[6][i][j][l] = buf[5][i][j][l - 1];

							if (mask[i + 1][j][l - 1] != 1.0)
								buf[14][i][j][l] = buf[12][i + 1][j][l - 1];

							if (mask[i][j + 1][l + 1] != 1.0)
								buf[15][i][j][l] = buf[17][i][j + 1][l + 1];

							if (mask[i][j - 1][l + 1] != 1.0)
								buf[16][i][j][l] = buf[18][i][j - 1][l + 1];

							if (mask[i][j - 1][l -1] != 1.0)
								buf[17][i][j][l] = buf[15][i][j - 1][l - 1];

							if (mask[i][j + 1][l - 1] != 1.0)
								buf[18][i][j][l] = buf[16][i][j + 1][l - 1];
						}
					}
				}
			}
			
	#pragma omp parallel for
			for (size_t i = 1; i < n + 1; i++) {
					for (size_t j = 1; j < m + 1; j++) {
						for (size_t l = 1; l < k + 1; l++) {
							for (size_t s = 1; s < 19; s++) {
								if(mask[i][j][l] == 0.0)
								f[s][i][j][l] = buf[s][i - dx[s]][j - dy[s]][l - dz[s]];
							}
						}
					}
				}; 

			/* посчитаем новую плотность */
			double rho_min = 100;
			double rho_max = -100;

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
			for (int i = 1; i < n + 1; i++) {
				for (int j = 1; j < m + 1; j++) {
					for (int l = 1; l < k + 1; l++) {
						if (rho_min > rho[i][j][l]) rho_min = rho[i][j][l];
						if (rho_max < rho[i][j][l]) rho_max = rho[i][j][l];
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
			/*cout << " Density for " << t << " step = " << full_rho << endl;
			cout << " rho_min = " << rho_min << "     , rho_max = " << rho_max << endl;

			/* посчитаем новую скорость вещества в узле */
#pragma omp parallel for
			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					for (size_t l = 1; l < k + 1; l++) {
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
			/*cout << " Impulse for " << t << " step (x) = " << full_rho_ux << endl;
			cout << " Impulse for " << t << " step (y) = " << full_rho_uy << endl;
			cout << " Impulse for " << t << " step (z) = " << full_rho_uz << endl;


			/* посчитаем эффективную плотность */
#pragma omp parallel for
			for (size_t i = 0; i < n + 2; i++) {
				for (size_t j = 0; j < m + 2; j++) {
					for (size_t l = 0; l < k + 2; l++) {
						Fi[i][j][l] = sqrt(rho[i][j][l] * teta - k_koeff * PressurePengRobinson(rho[i][j][l], T, omega));
					}
				}
			};

			/* задаем эффективную плотность на границе 
			в препятствиях на границе сделать сумму фишек по жидкости умноженные на G и разделить на сумму весов
			дождаться установления макс и мин плотность взять за начальные полоску жидкости с плавным переходом в пар или наоборот*/

			Fi[0][0][0] = Fi[1][1][1];
			Fi[n + 1][0][0] = Fi[n][1][1];
			Fi[0][m + 1][0] = Fi[1][m][1];
			Fi[0][0][k + 1] = Fi[1][1][k];
			Fi[n + 1][m + 1][0] = Fi[n][m][1];
			Fi[n + 1][0][k + 1] = Fi[n][1][k];
			Fi[0][m + 1][k + 1] = Fi[1][m][k];
			Fi[n + 1][m + 1][k + 1] = Fi[n][m][k];

			for (int i = 1; i < n + 1; i++) {
				Fi[i][0][0] = Fi[i][1][1];
				Fi[i][0][k + 1] = Fi[i][1][k];
				Fi[i][m + 1][0] = Fi[i][m][1];
				Fi[i][m + 1][k + 1] = Fi[i][m][k];
			}

			for (int i = 1; i < m + 1; i++) {
				Fi[0][i][0] = Fi[1][i][1];
				Fi[n + 1][i][k + 1] = Fi[n][i][k];
				Fi[n + 1][i][0] = Fi[n][i][1];
				Fi[0][i][k + 1] = Fi[1][i][k];
			}

			for (int i = 1; i < k + 1; i++) {
				Fi[0][0][i] = Fi[1][1][i];
				Fi[n + 1][m + 1][i] = Fi[n][m][i];
				Fi[n + 1][0][i] = Fi[n][1][i];
				Fi[0][m + 1][i] = Fi[1][m][i];
			}

			for (int i = 1; i < n + 1; i++) {
				for (int j = 1; j < m + 1; j++) {
					Fi[i][j][k + 1] = Fi[i][j][k];
					Fi[i][j][0] = Fi[i][j][1];
				}
			}

			for (int i = 1; i < n + 1; i++) {
				for (int j = 1; j < k + 1; j++) {
					Fi[i][m + 1][j] = Fi[i][m][j];
					Fi[i][0][j] = Fi[i][1][j];
				}
			}

			for (int i = 1; i < m + 1; i++) {
				for (int j = 1; j < k + 1; j++) {
					Fi[n + 1][i][j] = Fi[1][i][j];
					Fi[0][i][j] = Fi[n][i][j];
				}
			}

			/* препятствия для эффективной плотности */

			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					for (size_t l = 1; l < k + 1; l++) {
						if (mask[i][j][l] == 1.0) {
							sumG = 0;
							Fi[i][j][l] = 0;
							for (int s = 1; s < 19; s++) {
								if (mask[i + dx[s]][j + dy[s]][l + dz[s]] != 1.0) {
									Fi[i][j][l] += Fi[i + dx[s]][j + dy[s]][l + dz[s]] * G[s];
									sumG += G[s];
								}
							}
							/*
							if (mask[i + 1][j][l] != 1.0) {
								Fi[i][j][l] += Fi[i + 1][j][l] * G[1];
								sumG += G[1];
							}

							if (mask[i + 1][j + 1][l] != 1.0)
								Fi[i][j][l] += Fi[i + 1][j + 1][l] * G[7];

							if (mask[i][j + 1][l] != 1.0)
								Fi[i][j][l] += Fi[i][j + 1][l] * G[3];

							if (mask[i - 1][j + 1][l] != 1.0)
								Fi[i][j][l] += Fi[i - 1][j + 1][l] * G[8];

							if (mask[i - 1][j][l] != 1.0)
								Fi[i][j][l] += Fi[i - 1][j][l] * G[2];

							if (mask[i - 1][j - 1][l] != 1.0)
								Fi[i][j][l] += Fi[i - 1][j - 1][l] * G[9];

							if (mask[i][j - 1][l] != 1.0)
								Fi[i][j][l] += Fi[i][j - 1][l] * G[4];

							if (mask[i + 1][j - 1][l] != 1.0)
								Fi[i][j][l] += Fi[i + 1][j - 1][l] * G[10];

							if (mask[i + 1][j][l + 1] != 1.0)
								Fi[i][j][l] += Fi[i + 1][j][l + 1] * G[11];

							if (mask[i][j][l + 1] != 1.0)
								Fi[i][j][l] += Fi[i][j][l + 1] * G[5];

							if (mask[i - 1][j][l + 1] != 1.0)
								Fi[i][j][l] += Fi[i - 1][j][l + 1] * G[12];

							if (mask[i - 1][j][l - 1] != 1.0)
								Fi[i][j][l] += Fi[i - 1][j][l - 1] * G[13];

							if (mask[i][j][l - 1] != 1.0)
								Fi[i][j][l] += Fi[i][j][l - 1] * G[6];

							if (mask[i + 1][j][l - 1] != 1.0)
								Fi[i][j][l] += Fi[i + 1][j][l - 1] * G[14];

							if (mask[i][j + 1][l + 1] != 1.0)
								Fi[i][j][l] += Fi[i][j + 1][l + 1] * G[15];

							if (mask[i][j - 1][l + 1] != 1.0)
								Fi[i][j][l] += Fi[i][j - 1][l + 1] * G[16];

							if (mask[i][j - 1][l - 1] != 1.0)
								Fi[i][j][l] += Fi[i][j - 1][l - 1] * G[17];

							if (mask[i][j + 1][l - 1] != 1.0)
								Fi[i][j][l] += Fi[i][j + 1][l - 1] * G[18];*/
							if (sumG != 0)
								Fi[i][j][l] = Fi[i][j][l] / sumG;

						}
					}
				}
			}

			/* посчитаем изменение скорости */
	#pragma omp parallel for
			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					for (size_t l = 1; l < k + 1; l++) {
						dux[i][j][l] = 0.;
						duy[i][j][l] = 0.;
						duz[i][j][l] = 0.;
						if (mask[i][j][l] == 0.0) {
							for (size_t s = 1; s < 19; s++) {
								dux[i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i + dx[s]][j][l] * dx[s] +
									A * G[s] * Fi[i + dx[s]][j][l] * Fi[i + dx[s]][j][l] * dx[s]) / rho[i][j][l];
								duy[i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i][j + dy[s]][l] * dy[s] +
									A * G[s] * Fi[i][j + dy[s]][l] * Fi[i][j + dy[s]][l] * dy[s]) / rho[i][j][l];
								duz[i][j][l] += 1. / 3. * ((1 - 2 * A) * Fi[i][j][l] * G[s] * Fi[i][j][l + dz[s]] * dz[s] +
									A * G[s] * Fi[i][j][l + dz[s]] * Fi[i][j][l + dz[s]] * dz[s]) / rho[i][j][l];
							}
						}
					}
				}
			};

			/* сделаем одну итерацию */
		#pragma omp parallel for
			for (size_t i = 1; i < n + 1; i++) {
				for (size_t j = 1; j < m + 1; j++) {
					for (size_t l = 1; l < k + 1; l++) {
						for (size_t s = 0; s < 19; s++) {
							f[s][i][j][l] = F(f[s][i][j][l],
								F_e(c[s][0] * ux[i][j][l] + c[s][1] * uy[i][j][l] + c[s][2] * uz[i][j][l]
									, ux[i][j][l] * ux[i][j][l] + uy[i][j][l] * uy[i][j][l] + uz[i][j][l] * uz[i][j][l], w[s], rho[i][j][l]),
								F_e(c[s][0] * (ux[i][j][l] + g) + c[s][1] * uy[i][j][l] + c[s][2] * uz[i][j][l]
									, (ux[i][j][l] + g) * (ux[i][j][l] + g) + uy[i][j][l] * uy[i][j][l] + uz[i][j][l] * uz[i][j][l], w[s], rho[i][j][l]),
								F_e(c[s][0] * (ux[i][j][l] + dux[i][j][l]) + c[s][1] * (uy[i][j][l] + duy[i][j][l]) + c[s][2] * (uz[i][j][l] + duz[i][j][l])
									, (ux[i][j][l] + dux[i][j][l]) * (ux[i][j][l] + dux[i][j][l]) + (uy[i][j][l] + duy[i][j][l]) *
									(uy[i][j][l] + duy[i][j][l]) + (uz[i][j][l] + duz[i][j][l]) * (uz[i][j][l] + duz[i][j][l]), w[s], rho[i][j][l]));
						}
					}
				}
			};

			if (t % 50 == 0) 
			{ 
				SaveVTKFile(t); 
				cout << " Density for " << t << " step = " << full_rho << endl;
				cout << " rho_min = " << rho_min << "     , rho_max = " << rho_max << endl;
				cout << " Impulse for " << t << " step (x) = " << full_rho_ux << endl;
				cout << " Impulse for " << t << " step (y) = " << full_rho_uy << endl;
				cout << " Impulse for " << t << " step (z) = " << full_rho_uz << endl;
			}

		}


	return 0;
}
