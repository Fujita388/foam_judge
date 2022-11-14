#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include <cmath>
#include "split.h"


using namespace std;


int index_of_density(int mx, int my, int mz) {
	int index = mx + my * 21 + mz * 21 * 21;  // 密度データの中でのインデックス
	return index;
}


// 気泡破壊を判定する関数
void gas_volume(double d, double thresh) {
	ifstream ifile("rescale.lammpstrj");  // 読み込むファイルのパスを指定
	ofstream ofile("check.dat");  // 書き出すファイルのパスを指定

	const int N = 120000;
	int num_atoms;  // 粒子数
	double L;  // ボックスサイズ
	double V = pow(d, 3.0);  // セルの体積
	int Lx, Ly, Lz;  // 各次元のセルの個数
	int i_dump = 0;  // dumpfileの各行の通し番号
	int i_step = 0;  // step内でのインデックス 
	string str_dump;
	static double pos_data[N][3];  // 1stepあたりの粒子の座標データの配列
	vector<double> density;  // 密度データの配列

	while (getline(ifile, str_dump)) {  // ifileを1行ずつstr_dumpに読み込む
		if (i_dump == 3) {
			num_atoms = stoi(str_dump);  
		}
		if (i_dump == 5) {
			L = stod(split(str_dump, ' ')[1]);  // 立方体のシミュレーションボックス
			Lx = Ly = Lz = (int)(L / d);  // 一次元のセルの個数
			density.resize(Lx * Ly * Lz);  // densityをリサイズ
		}
		if (split(str_dump, ' ').size() == 8) {  // 座標データの行の場合
			i_step += 1;

			// pos_data[][]にstep内の座標データを保存
			pos_data[i_step-1][0] = stod(split(str_dump, ' ')[2]);
			pos_data[i_step-1][1] = stod(split(str_dump, ' ')[3]);  
			pos_data[i_step-1][2] = stod(split(str_dump, ' ')[4]);

			/////座標データから密度計算/////
			if (i_step == num_atoms) {
				fill(density.begin(), density.end(), 0);  // densityを0で初期化
				int gas = 0;  // 閾値以下のセルを気泡としてカウント
				for (int i = 0; i < num_atoms; i++) {
					int mx = int(pos_data[i][0] / d);
					int my = int(pos_data[i][1] / d);
					int mz = int(pos_data[i][2] / d);
					if (mx >= Lx) mx -= Lx;  // はみだした値に周期境界条件を適用
					if (my >= Ly) my -= Ly;
					if (mz >= Lz) mz -= Lz;
					int i_density = mx + my * Lx + mz * Lx * Ly;  // 密度データの中でのインデックス
					density[i_density] += 1.0 / V;
				}

				// 上下(z=0とz=21)のセルの10%以上が気泡となった時、気泡破壊と判定する
				int check_up = 0;  // z=0のセルを気泡かどうかを判定
				int check_bottom = 0;  // z=21のセルを気泡かどうかを判定
				for (int i = 0; i < 21; i++) {
					for (int j = 0; j < 21; j++) {
						if (density[i+j*21+0*21*21] <= 0.1) check_up += 1;
						if (density[i+j*21+21*21*21] <= 0.1) check_bottom += 1;
					}
				}

				if (check_up / (21 * 21) >= 0.1 && check_bottom / (21 * 21) >= 0.1) ofile << "Yes rupture" << '\n';
				else ofile << "Not rupture" << '\n';

				i_step = 0;  // i_stepを初期化
			}
		}
		i_dump += 1;
	}
}


int main() {
	gas_volume(3.0, 0.1);
	return 0;
}
