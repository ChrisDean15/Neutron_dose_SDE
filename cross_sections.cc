#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#ifndef CS
#define CS

struct CS_1d {

  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
  }
  CS_1d(const std::string filename, const double cuttoff) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    double tmp_val, tmp_val_old = 0, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (tmp_val - cuttoff) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
        }
      }
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          rate.push_back(tmp_val * lin_inter_val +
                         (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
    }
    file.close();
  }

  double hydrogen_cm_to_lab(double ang, const double E) {
    ang = M_PI - ang;
    double mp = 938.346;
    double p = sqrt(E * (E + 2 * mp));
    double u = p / (E + 2 * mp);
    double g = 1 / sqrt(1 - u * u);
    double e = E + mp;
    double v_ratio = u * (e - u * p) / (p - u * e);
    double out;
    if (fabs(g * (cos(ang) + v_ratio)) == 0) {
      out = M_PI / 2;
    } else {
      out = atan(sin(ang) / (g * (cos(ang) + v_ratio)));
    }
    if (out < 0) {
      out += M_PI;
    }
    return out;
  }

  CS_1d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
        }
      }
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      rate.push_back(total_rate);
    }
    file.close();
  }

  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}

  CS_1d() : energy(), rate() {}

  double evaluate(const double e) const {
    double ret = 0;
    int r;
    double tol = 1e-7;
    if (energy.size() > 0) {
      if (e <= energy[0]) {
        ret = rate[0];
      } else if (e >= energy.back()) {
        ret = rate.back();
      } else {
        r = std::distance(energy.begin(),
                          std::lower_bound(energy.begin(), energy.end(), e));
        if (energy[r] - energy[r - 1] > tol) {
          ret =
              ((energy[r] - e) * rate[r - 1] + (e - energy[r - 1]) * rate[r]) /
              (energy[r] - energy[r - 1]);
        } else {
          ret = energy[r];
        }
      }
    }
    return ret;
  }

  std::vector<double> energy, rate;
};

struct CS_3d {

  CS_3d(const std::vector<std::string> filenames) : energy(), yield(), exit_energy(), cdf(), rvalue() {
    for (size_t i=0; i!=filenames.size(); i++) {
      std::vector<double> tmp_energy, tmp_yield;
      std::vector<std::vector<double>> tmp_exit_energy, tmp_cdf, tmp_rvalue;
      std::string filename = filenames[i];
      std::ifstream file;
      file.open(filename);
      std::string line, token;
      getline(file, line);
      std::stringstream iss;
      iss << line;
      while (getline(iss, token, ' ')) {
        tmp_energy.push_back(atof(token.c_str()));
      }
      getline(file, line);
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_yield.push_back(atof(token.c_str()));
      }
      std::vector<double> tmp_vec;
      while (getline(file, line)) {
        tmp_vec.clear();
        std::stringstream iss3;
        iss3 << line;
        while (getline(iss3, token, ' ')) {
          tmp_vec.push_back(atof(token.c_str()));
        }
        tmp_exit_energy.push_back(tmp_vec);
        getline(file, line);
        tmp_vec.clear();
        std::stringstream iss4;
        iss4 << line;
        while (getline(iss4, token, ' ')) {
          tmp_vec.push_back(atof(token.c_str()));
        }
        tmp_cdf.push_back(tmp_vec);
        getline(file, line);
        tmp_vec.clear();
        std::stringstream iss5;
        iss5 << line;
        while (getline(iss5, token, ' ')) {
          tmp_vec.push_back(atof(token.c_str()));
        }
        tmp_rvalue.push_back(tmp_vec);
      }
        file.close();
        energy.push_back(tmp_energy);
        exit_energy.push_back(tmp_exit_energy);
        cdf.push_back(tmp_cdf);
        rvalue.push_back(tmp_rvalue);
        yield.push_back(tmp_yield);
    }
  }

  CS_3d(const CS_3d &other)
      : energy(other.energy), yield(other.yield), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  CS_3d() : energy(), yield(), exit_energy(), cdf(), rvalue() {}
  void sample_secondary_particles(const double e, gsl_rng *gen, std::vector<std::vector<std::vector<double>>> &secondary_particles) const {
    double tmp_rvalue, tmp_energy_cm;
    for(size_t i=0;i!= energy.size();i++){
      std::vector<std::vector<double>> tmp_particle_type;
      size_t energy_index = std::distance(
        energy[i].begin(), std::lower_bound(energy[i].begin(), energy[i].end(), e));
      double yield_temp = linear_interpolate_from_energy_index(e,energy_index,i);
      int particles_to_generate = floor(yield_temp+gsl_rng_uniform(gen));
      for(int j = 0; j!=particles_to_generate; j++){
      double u = gsl_rng_uniform(gen);
      if (energy_index == 0) {
      sample_from_energy_index(0, u,
                                tmp_energy_cm,tmp_rvalue,i);
    } else if (energy_index == energy[i].size()) {
      sample_from_energy_index(energy_index - 1, u,
                                tmp_energy_cm,tmp_rvalue,i);
    } else {
      sample_from_energy_index(energy_index, u,
                                tmp_energy_cm,tmp_rvalue,i);
      }
      tmp_particle_type.push_back({tmp_energy_cm, tmp_rvalue});

    }
    secondary_particles.push_back(tmp_particle_type);

  }
  return;
}
  double linear_interpolate_from_energy_index(double e, size_t energy_index, int i )const {
    if (energy_index == 0) {
      return yield[i][0];
    } else if (energy_index == yield[i].size()){
      return yield[i][yield[i].size() -1];
    } else {
      double tmp_interp = ((e-energy[i][energy_index-1])/(energy[i][energy_index]-energy[i][energy_index-1]));
      return yield[i][energy_index]*tmp_interp + yield[i][energy_index-1]*(1-tmp_interp);
    }
  }
  void sample_from_energy_index(const int energy_index, const double u,
                                double &out_energy_cm,
                                double &out_rvalue, int i) const {
    double diff = 0;
    double tol = 1e-7;
    int density_index =
        std::distance(cdf[i][energy_index].begin(),
                      std::lower_bound(cdf[i][energy_index].begin(),
                                       cdf[i][energy_index].end(), u));
    if (density_index == 0) {
      out_energy_cm = exit_energy[i][energy_index][0];
      out_rvalue = rvalue[i][energy_index][0];
    } else if (density_index == int(cdf[i][energy_index].size())) {
      out_energy_cm = exit_energy[i][energy_index].back();
      out_rvalue = rvalue[i][energy_index].back();
    } else {
      if (cdf[i][energy_index][density_index] -
              cdf[i][energy_index][density_index - 1] >
          tol) {
        diff = (u - cdf[i][energy_index][density_index - 1]) /
               (cdf[i][energy_index][density_index] -
                cdf[i][energy_index][density_index - 1]);
        out_energy_cm =
            exit_energy[i][energy_index][density_index] * diff +
            exit_energy[i][energy_index][density_index - 1] * (1 - diff);
        out_rvalue = rvalue[i][energy_index][density_index] * diff +
                     rvalue[i][energy_index][density_index - 1] * (1 - diff);
      } else {
        out_energy_cm = exit_energy[i][energy_index][density_index];
        out_rvalue = rvalue[i][energy_index][density_index];
      }
    }
    return;
  }
  std::vector<std::vector<double>> energy, yield;
  std::vector<std::vector<std::vector<double>>> exit_energy, cdf, rvalue;
};

struct CS_2d {

  CS_2d(const std::string filename, const double cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double tmp_val, tmp_val_old = 0, total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      tmp_vec.clear();
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      total_rate = tmp_vec.back();
      for (double &i : tmp_vec) {
        i /= total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }

  double hydrogen_cm_to_lab(double ang, const double E) {
    ang = M_PI - ang;
    double mp = 938.346;
    double p = sqrt(E * (E + 2 * mp));
    double u = p / (E + 2 * mp);
    double g = 1 / sqrt(1 - u * u);
    double e = E + mp;
    double v_ratio = u * (e - u * p) / (p - u * e);
    double out;
    if (fabs(g * (cos(ang) + v_ratio)) == 0) {
      out = M_PI / 2;
    } else {
      out = atan(sin(ang) / (g * (cos(ang) + v_ratio)));
    }
    if (out < 0) {
      out += M_PI;
    }
    return out;
  }

  CS_2d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_vec.clear();
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      for (double &i : tmp_vec) {
        i = (i - bottom_rate) / total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }

  CS_2d(const CS_2d &other)
      : energy(other.energy), exit_angle(other.exit_angle), cdf(other.cdf) {}

  CS_2d() : energy(), exit_angle(), cdf() {}

  double sample_from_energy_index(const size_t energy_index,
                                  const double u) const {
    double ret = 0;
    double diff = 0;
    double tol = 1e-7;
    int density_index =
        std::distance(cdf[energy_index].begin(),
                      std::lower_bound(cdf[energy_index].begin(),
                                       cdf[energy_index].end(), u));
    if (density_index == 0) {
      ret = exit_angle[energy_index][0];
    } else if (density_index == int(cdf[energy_index].size())) {
      ret = exit_angle[energy_index].back();
    } else {
      if (cdf[energy_index][density_index] -
              cdf[energy_index][density_index - 1] >
          tol) {
        diff = (u - cdf[energy_index][density_index - 1]) /
               (cdf[energy_index][density_index] -
                cdf[energy_index][density_index - 1]);
        ret = exit_angle[energy_index][density_index] * diff +
              exit_angle[energy_index][density_index - 1] * (1 - diff);
      } else {
        ret = exit_angle[energy_index][density_index];
      }
    }
    return ret;
  }

  double sample(const double e, gsl_rng *gen) const {
    int energy_index = std::distance(
        energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    double u = gsl_rng_uniform(gen);
    double out_angle_cm;
    double out_angle_cm_2;
    double diff;
    double tol = 1e-7;
    if (energy_index == 0) {
      out_angle_cm = sample_from_energy_index(0, u);
    } else if (energy_index == int(energy.size())) {
      out_angle_cm = sample_from_energy_index(energy_index - 1, u);
    } else {
      out_angle_cm = sample_from_energy_index(energy_index, u);
      if (energy[energy_index] - energy[energy_index - 1] > tol) {
        out_angle_cm_2 = sample_from_energy_index(energy_index - 1, u);
        diff = (e - energy[energy_index - 1]) /
               (energy[energy_index] - energy[energy_index - 1]);
        out_angle_cm = out_angle_cm * diff + out_angle_cm_2 * (1 - diff);
      }
    }
    return out_angle_cm;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_angle, cdf;
};

#endif
