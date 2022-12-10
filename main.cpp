#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

const int N = 2000000;
using namespace std;

// 空白でsplitする関数
std::vector<std::string> mysplit(std::string &str) {
  static std::vector<std::string> r;
  r.clear();
  int s = 0;
  for (int i = 0; i < str.length(); i++) {
    if (str[i] == ' ') {
      r.push_back(str.substr(s, i - s));
      s = i + 1;
      i++;
    }
  }
  r.push_back(str.substr(s, str.length() - s));
  return r;
}

void save_vtk(int frame, int Lx, int Ly, int Lz, std::vector<double> &density) {
  char filename[1000];
  sprintf(filename, "test%03d.vtk", frame);
  std::cout << filename << std::endl;
  std::ofstream ofs(filename);
  ofs << "# vtk DataFile Version 1.0" << std::endl;
  ofs << "test" << std::endl;
  ofs << "ASCII" << std::endl;
  ofs << "DATASET STRUCTURED_POINTS" << std::endl;
  ofs << "DIMENSIONS " << Lx << " " << Ly << " " << Lz << std::endl;
  ofs << "ORIGIN 0.0 0.0 0.0" << std::endl;
  ofs << "SPACING 1.0 1.0 1.0" << std::endl;
  ofs << std::endl;
  ofs << "POINT_DATA"
      << " " << density.size() << std::endl;
  ofs << std::endl;
  ofs << "SCALARS density float" << std::endl;
  ofs << "LOOKUP_TABLE default" << std::endl;
  for (auto &d : density) {
    ofs << d << std::endl;
  }
}

//密度の計算
vector<double> calc_density(double pos_data[N][3], int num_atoms, double V, int Lx, int Ly, int Lz, double d) {
  vector<double> density(Lx * Ly * Lz, 0);
  for (int i = 0; i < num_atoms; i++) {
    int mx = int(pos_data[i][0] / d);
    int my = int(pos_data[i][1] / d);
    int mz = int(pos_data[i][2] / d);
    if (mx >= Lx) mx -= Lx; // はみだした値に周期境界条件を適用
    if (my >= Ly) my -= Ly;
    if (mz >= Lz) mz -= Lz;
    int i_density = mx + my * Lx + mz * Lx * Ly; // 密度データの中でのインデックス
    density[i_density] += 1.0 / V;
  }
  return density;
}

int index2pos(int x, int y, int z, int Lx, int Ly) {
  return x + y * Lx + z * Lx * Ly;
}

//破壊をチェックする関数
bool check_rapture(int Lx, int Ly, int Lz, std::vector<double> &density, double thresh) {
  // 上下(x,y,z=0とx,y,z=84)のセルの10%以上が気泡となった時、気泡破壊と判定する
  int x_up = 0; // x=0のセルを気泡かどうかを判定
  int y_up = 0;
  int z_up = 0;
  int x_bottom = 0; // x=84のセルを気泡かどうかを判定
  int y_bottom = 0;
  int z_bottom = 0;
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      if (density[index2pos(0, i, j, Lx, Ly)] <= thresh) x_up += 1;
      if (density[index2pos(i, 0, j, Lx, Ly)] <= thresh) y_up += 1;
      if (density[index2pos(i, j, 0, Lx, Ly)] <= thresh) z_up += 1;
      if (density[index2pos(Lx - 0, i, j, Lx, Ly)] <= thresh) x_bottom += 1;
      if (density[index2pos(i, Ly - 1, j, Lx, Ly)] <= thresh) y_bottom += 1;
      if (density[index2pos(i, j, Lz - 1, Lx, Ly)] <= thresh) z_bottom += 1;
    }
  }

  double x_up_density = static_cast<double>(x_up) / Ly * Lz;
  double x_bottom_density = static_cast<double>(x_bottom) / Ly * Ly;
  double y_up_density = static_cast<double>(y_up) / Lz * Lx;
  double y_bottom_density = static_cast<double>(y_bottom) / Lz * Lx;
  double z_up_density = static_cast<double>(z_up) / Lx * Ly;
  double z_bottom_density = static_cast<double>(z_bottom) / Lx * Ly;

  if (x_up_density >= 0.1 && x_bottom_density >= 0.1) return true;
  if (y_up_density >= 0.1 && y_bottom_density >= 0.1) return true;
  if (z_up_density >= 0.1 && z_bottom_density >= 0.1) return true;

  return false;
}

// 気泡破壊を判定する関数
void gas_volume(double d, double thresh, std::string input_file, std::string output_file) {
  ifstream ifile(input_file);  // 読み込むファイルのパスを指定
  ofstream ofile(output_file); // 書き出すファイルのパスを指定

  std::cout << "Reading: " << input_file << std::endl;
  int num_atoms;        // 粒子数
  double L;             // ボックスサイズ
  double V = d * d * d; // セルの体積
  int Lx, Ly, Lz;       // 各次元のセルの個数
  int i_dump = 0;       // dumpfileの各行の通し番号
  int i_step = 0;       // step内でのインデックス
  string str_dump;
  static double pos_data[N][3]; // 1stepあたりの粒子の座標データの配列
  int frame = 0;

  // ifileを1行ずつstr_dumpに読み込む
  while (getline(ifile, str_dump) && frame < 10) { // とりあえず10フレームだけ処理する(全部処理したい場合は消すこと)
    if (i_dump == 3) {
      num_atoms = stoi(str_dump);
    }
    if (i_dump == 5) {
      L = stod(mysplit(str_dump)[1]); // 立方体のシミュレーションボックス
      Lx = Ly = Lz = (int)(L / d);    // 一次元のセルの個数
    }
    auto list = mysplit(str_dump);
    if (list.size() == 8) { // 座標データの行の場合
      i_step += 1;

      // pos_data[][]にstep内の座標データを保存
      pos_data[i_step - 1][0] = stod(list[2]);
      pos_data[i_step - 1][1] = stod(list[3]);
      pos_data[i_step - 1][2] = stod(list[4]);

      // 全ての原子を読み込んだ
      if (i_step == num_atoms) {
        std::cout << "Frame: " << frame << std::endl;
        auto density = calc_density(pos_data, num_atoms, V, Lx, Ly, Lz, d);
        if (check_rapture(Lx, Ly, Lz, density, thresh)) {
          ofile << "Yes Rapture" << std::endl;
        } else {
          ofile << "Not Rapture" << std::endl;
        }
        save_vtk(frame, Lx, Ly, Lz, density);
        i_step = 0; // i_stepを初期化
        frame++;
      }
    }
    i_dump += 1;
  }
  std::cout << "Generated: " << output_file << std::endl;
}

int main() {
  gas_volume(1.4875, 0.1, "rescale.lammpstrj", "f_j_no_surf01.dat");
}
