#include "cross_sections.cc"
#include <cstdlib>
#include <vector>
#include <cmath>

#ifndef MAT
#define MAT

struct Atom {
  // General constructor with nucleus of atom equal to incident particle
  Atom(const int a_atom0, const int z_atom0, const int a_inc0, const int z_inc0, const std::vector<std::vector<int>> charge_mass_numbers0, const double I_in0, const std::vector<double> I_out0, const double c_in0, const std::vector<double> c_out0, const std::string ne_r,
       const std::string el_ruth_cs, const std::vector<std::string> ne_ea,
       const double cutoff, const double back_cutoff)
      : a_atom(a_atom0), z_atom(z_atom0), a_inc(a_inc0), z_inc(z_inc0), charge_mass_numbers(charge_mass_numbers0), I_in(I_in0), I_out(I_out0), c_in(c_in0),c_out(c_out0), el_ruth_rate(el_ruth_cs, cutoff,back_cutoff), ne_rate(ne_r),
        el_ruth_angle_cdf(el_ruth_cs, cutoff,back_cutoff), ne_energy_angle(ne_ea) {}
  // General constructor with nucleus of atom not equal to incident particle
   Atom(const int a_atom0, const int z_atom0, const int a_inc0, const int z_inc0, const std::vector<std::vector<int>> charge_mass_numbers0, const double I_in0, const std::vector<double> I_out0, const double c_in0, const std::vector<double> c_out0, const std::string ne_r,
       const std::string el_ruth_cs, const std::vector<std::string> ne_ea,
       const double cutoff)
      : a_atom(a_atom0), z_atom(z_atom0), a_inc(a_inc0), z_inc(z_inc0), charge_mass_numbers(charge_mass_numbers0), I_in(I_in0), I_out(I_out0), c_in(c_in0), c_out(c_out0), el_ruth_rate(el_ruth_cs, cutoff), ne_rate(ne_r),
        el_ruth_angle_cdf(el_ruth_cs, cutoff), ne_energy_angle(ne_ea) {}
  // Hydrogen constructor when incident particle is not a proton (zero non-elastic rate)
   Atom(const int a_atom0, const int z_atom0, const int a_inc0, const int z_inc0,
       const std::string el_ruth_cs, const double cutoff)
      : a_atom(a_atom0), z_atom(z_atom0), a_inc(a_inc0), z_inc(z_inc0), charge_mass_numbers(), I_in(), I_out(), c_in(),c_out(), el_ruth_rate(el_ruth_cs, cutoff), ne_rate(),
        el_ruth_angle_cdf(el_ruth_cs, cutoff), ne_energy_angle() {}
  // Hydrogen constructor when incident particle is a proton (zero non-elastic rate)
   Atom(const int a_atom0, const int z_atom0, const int a_inc0, const int z_inc0, 
       const std::string el_ruth_cs, const double cutoff, const double back_cutoff)
      : a_atom(a_atom0), z_atom(z_atom0), a_inc(a_inc0), z_inc(z_inc0), charge_mass_numbers(), I_in(), I_out(), c_in(),  c_out(), el_ruth_rate(el_ruth_cs, cutoff, back_cutoff), ne_rate(),
        el_ruth_angle_cdf(el_ruth_cs, cutoff, back_cutoff), ne_energy_angle() {}
  // Copy constructor
  Atom(const Atom &other)
      : a_atom(other.a_atom), z_atom(other.z_atom), a_inc(other.a_inc), z_inc(other.z_inc), charge_mass_numbers(other.charge_mass_numbers), I_in(other.I_in), I_out(other.I_out), c_in(other.c_in),c_out(other.c_out), el_ruth_rate(other.el_ruth_rate),
        ne_rate(other.ne_rate), el_ruth_angle_cdf(other.el_ruth_angle_cdf),
        ne_energy_angle(other.ne_energy_angle) {}

  double s(const double a_a, const double a_c, const double z_a, const double z_c, const double I) const {
    double n_c = a_c - z_c;
    double n_a = a_a - z_a;
    double ret = 15.68 * (a_c - a_a) -
                 28.07 * (pow(n_c - z_c, 2) / a_c - pow(n_a - z_a, 2) / a_a) -
                 18.56 * (pow(a_c, 2.0 / 3.0) - pow(a_a, 2.0 / 3.0)) +
                 33.22 * (pow(n_c - z_c, 2) / pow(a_c, 4.0 / 3) -
                          pow(n_a - z_a, 2) / pow(a_a, 4.0 / 3)) -
                 0.717 * (z_c * z_c / pow(a_c, 1.0 / 3) -
                          z_a * z_a / pow(a_a, 1.0 / 3)) +
                 1.211 * (z_c * z_c / a_c - z_a * z_a / a_a)-I;
    return ret;
  }
  std::vector<double> compute_new_angle(std::vector<double> ang, const double alpha,
                         const double beta) const {
    double omega_new1 =
        sin(ang[0]) * cos(ang[1]) * cos(alpha) +
        (cos(ang[0]) * cos(ang[1]) * sin(beta) - sin(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new2 =
        sin(ang[0]) * sin(ang[1]) * cos(alpha) +
        (cos(ang[0]) * sin(ang[1]) * sin(beta) + cos(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new3 =
        cos(ang[0]) * cos(alpha) - sin(ang[0]) * sin(beta) * sin(alpha);
    double magnitude =
        std::sqrt(omega_new1 * omega_new1 + omega_new2 * omega_new2 +
                  omega_new3 * omega_new3);
    omega_new1 /= magnitude;
    omega_new2 /= magnitude;
    omega_new3 /= magnitude;
    return {acos(omega_new3), atan2(omega_new2, omega_new1)};
  }
  void sample_nonelastic_collision(double &e, std::vector<double> ang,
                                   gsl_rng *gen, std::vector<std::vector<std::vector<double>>> &particle_list) const {
    ne_energy_angle.sample_secondary_particles(e, gen, particle_list); //should be an output of all possible particles to handle
    for(size_t i = 0; i!=particle_list.size(); i++) {
      for (size_t j = 0; j!=particle_list[i].size(); j++) {
        if (i==1) {
        //std::cout << particle_list[i][j][0] << "outgoing energy value" << std::endl;
        }
        double eps_a = a_atom * e / (a_atom + a_inc);
        double eps_b = (a_atom + a_inc) * particle_list[i][j][0] / (a_atom + a_inc - charge_mass_numbers[i][1]);
        double e_a = eps_a + s(a_atom,a_atom+a_inc,z_atom,z_atom+z_inc,I_in);
        double e_b = eps_b + s(a_atom+a_inc-charge_mass_numbers[i][1],a_atom+a_inc,z_atom+z_inc-charge_mass_numbers[i][0],z_atom+z_inc,I_out[i]);
        double x1 = fmin(e_a, 130) * e_b / e_a;
        double x3 = fmin(e_a, 41) * e_b / e_a;
        double aval = 0.04 * x1 + 1.8 * 1e-6 * pow(x1, 3) + c_in*c_out[i]*6.7 * 1e-7 * pow(x3, 4);
        double cdfc2 = particle_list[i][j][1] * cosh(aval) - sinh(aval);
        double cdfc1 = 2 * sinh(aval);
        double u2 = gsl_rng_uniform(gen);
        double z1 = cdfc1 * u2 + cdfc2;
        double z2 =
            (z1 + sqrt(pow(z1, 2) - pow(particle_list[i][j][1], 2) + 1)) / (particle_list[i][j][1] + 1);
        particle_list[i][j][1] = log(z2) / aval;
        double out_energy_lab =
            particle_list[i][j][0] + e*a_inc*charge_mass_numbers[i][1] / pow(a_inc+a_atom, 2) +
            2 * sqrt(particle_list[i][j][0] * e *a_inc*charge_mass_numbers[i][1]) * particle_list[i][j][1] / (a_atom + a_inc);
        particle_list[i][j][1] = sqrt(particle_list[i][j][0] / out_energy_lab) * particle_list[i][j][1] +
                              sqrt(e*a_inc*charge_mass_numbers[i][1] / out_energy_lab) / (a_atom + a_inc);
        particle_list[i][j][0] = out_energy_lab;
        double beta = 2 * M_PI * gsl_rng_uniform(gen);
        std::vector<double> ang_tmp = compute_new_angle(ang, acos(particle_list[i][j][1]), beta);
        particle_list[i][j][1] = ang_tmp[0];
        particle_list[i][j].push_back(ang_tmp[1]);
      }
    }
    e = 0;
  }
  const int a_atom, z_atom, a_inc, z_inc;
  const std::vector<std::vector<int>> charge_mass_numbers;
  const double I_in;
  const std::vector<double> I_out;
  const double c_in;
  const std::vector<double> c_out;
  CS_1d el_ruth_rate, ne_rate;
  CS_2d el_ruth_angle_cdf;
  CS_3d ne_energy_angle;
};

struct Material {

  Material(std::vector<Atom> &atoms, const std::vector<int> &id,
           const std::vector<double> &x0, const double d0, const double I0)
      : density(d0), I(I0 / 1e6), x(x0), at() {
    for (unsigned int i = 0; i < id.size(); i++) {
      at.push_back(atoms[id[i]]);
    }
  }

  Material() : density(), I(), x(), at() {}

  void read_material(const std::string filename, std::vector<Atom> &atoms) {
    std::ifstream file;
    file.open(filename + ".txt");
    std::string line, token;
    getline(file, line);
    density = atof(line.c_str());
    getline(file, line);
    I = atof(line.c_str()) / 1e6;
    while (getline(file, line)) {
      std::stringstream iss;
      iss << line;
      getline(iss, token, ' ');
      at.push_back(atoms[atoi(token.c_str())]);
      getline(iss, token, ' ');
      x.push_back(atof(token.c_str()));
    }
    file.close();
    return;
  }

  double bethe_bloch(const double e) const {
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mprocsq = 938.346; // mass of proton * speed of light squared, MeV
    double mpcsq = mprocsq*at[0].a_inc; //approx mass of particle * speed of light squared, MeV
    double z_part = at[0].z_inc;
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      ret += x[i] * pow(z_part,2)* 0.3072 * at[i].z_atom * density *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq * at[i].a_atom); // MeV / cm
    }
    return ret;
  }

  double multiple_scattering_sd(const double e, const double dt) const {
    double mprocsq = 938.346; // mass of proton * speed of light squared, MeV
    double mpcsq = mprocsq*at[0].a_inc; //approx mass of particle * speed of light squared, MeV 
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double c = 29979245800; // speed of light
    double vel = sqrt(betasq) * c;
    double p = pv / sqrt(betasq); // momentum in MeV / c.
    // effective chi_c_sq is just the sum of individual elements
    double chi_c_sq = 0;
    std::vector<double> chi_a_sq_vec(at.size());
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_c_sq += pow(at[i].z_inc,2)*x[i] * at[i].z_atom * (at[i].z_atom + 1.0) / at[i].a_atom;
      chi_a_sq_vec[i] = 2.007e-5 * pow(at[i].z_atom, 2.0 / 3) *
                        (1 + 3.34 * pow(at[i].z_inc*at[i].z_atom / (137 * vel), 2)) / (p * p);
    }
    chi_c_sq *= 0.157 * dt * density / (pv * pv);
    // effective chi_a_sq is a weighted average on the log-scale
    double chi_a_sq = 0;
    double denom = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_a_sq +=
          x[i] * at[i].z_atom * (at[i].z_atom + 1) * log(chi_a_sq_vec[i]) / at[i].a_atom;
      denom += x[i] * at[i].z_atom * (at[i].z_atom + 1.0) / at[i].a_atom;
    }
    chi_a_sq = exp(chi_a_sq / denom);
    double omega = chi_c_sq / chi_a_sq;
    double F = 0.98;
    double v = omega / (2 * (1 - F));
    double ret = sqrt(chi_c_sq * ((1 + v) * log(1 + v) / v - 1) / (1 + F * F));
    return ret;
  }

  double energy_straggling_sd(const double e) const {
    double alpha = 1 / 137.0;
    double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
    double log_c = log(29979245800);                              // cm / s
    double log_avogadro = log(6) + 23 * log(10);
    double mprocsq = 938.346; // mass of proton * speed of light squared, MeV
    double mpcsq = mprocsq*at[0].a_inc; //approx mass of particle * speed of light squared, MeV 
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double z = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      z += x[i] * at[i].z_atom / at[i].a_atom; // electrons per average molecule
    }
    double log_molecule_density =
        log(density) + log_avogadro; // molecules / cm^3
    double ret =
        4 * pow(at[0].z_inc,2)* M_PI * z * (1 - betasq / 2) / sqrt(1 - betasq) *
        exp(2 * (log(alpha) + log_hbar + log_c) + log_molecule_density);
    return sqrt(ret);
  }

  double nonelastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10);
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a_atom; // average molar mass
      ret += x[i] * at[i].ne_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm
  }

  double rutherford_and_elastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10);
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a_atom; // average molar mass
      ret += x[i] * at[i].el_ruth_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm
  }

  void compute_new_angle(std::vector<double> &ang, const double alpha,
                         const double beta) const {
    double omega_new1 =
        sin(ang[0]) * cos(ang[1]) * cos(alpha) +
        (cos(ang[0]) * cos(ang[1]) * sin(beta) - sin(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new2 =
        sin(ang[0]) * sin(ang[1]) * cos(alpha) +
        (cos(ang[0]) * sin(ang[1]) * sin(beta) + cos(ang[1]) * cos(beta)) *
            sin(alpha);
    double omega_new3 =
        cos(ang[0]) * cos(alpha) - sin(ang[0]) * sin(beta) * sin(alpha);
    double magnitude =
        std::sqrt(omega_new1 * omega_new1 + omega_new2 * omega_new2 +
                  omega_new3 * omega_new3);
    omega_new1 /= magnitude;
    omega_new2 /= magnitude;
    omega_new3 /= magnitude;
    ang[0] = acos(omega_new3);
    ang[1] = atan2(omega_new2, omega_new1);
    return;
  }

  void nonelastic_scatter(std::vector<double> &ang, double &e,
                          gsl_rng *gen, std::vector<std::vector<std::vector<double>>> &secondary_particles) const {
    double rate = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += x[i] * at[i].ne_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    int ind = 0;
    double tmp = x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    }
    // ENDF non-elastic scattering, both energy + angle from CM to LAB
    at[ind].sample_nonelastic_collision(e, ang, gen, secondary_particles);
    return;
  }

  void rutherford_elastic_scatter(std::vector<double> &ang, double &e,
                                  gsl_rng *gen) const {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double rate = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += x[i] * at[i].el_ruth_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    int ind = 0;
    double tmp = x[ind] * at[ind].el_ruth_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += x[ind] * at[ind].el_ruth_rate.evaluate(e) / rate;
    }
    double alpha;
    alpha = at[ind].el_ruth_angle_cdf.sample(e, gen);
    compute_new_angle(ang, alpha, beta);
    return;
  }

  double density; // density, g / cm^3
  double I;       // mean excitation energy, MeV
  std::vector<double> x;
  std::vector<Atom> at;
};

#endif
