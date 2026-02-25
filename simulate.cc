#include "cross_sections.cc"
#include "grid.cc"
#include "material.cc"
#include "particle_beam.cc"
#include "particle_bank.cc"
#include <cfloat>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Call " << argv[0] << " config-path" << std::endl;
    return 1;
  }
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  double rutherford_cutoff, backscatter_cutoff;
  libconfig::Config cfg;
  cfg.readFile(argv[1]);
  cfg.lookupValue("rutherford_cutoff", rutherford_cutoff);
  cfg.lookupValue("backscatter_cutoff", backscatter_cutoff);
  std::ifstream file;
  std::string line, token;
  std::vector<std::string> name_vec;
  std::vector<std::vector<int>> charge_mass_numbers;
  std::vector<double> c_in_vec, c_out_vec, I_vec;
  std::vector<Atom> atoms;
  std::string name;
  file.open("./Materials/particles.txt");
  while (getline(file,line)) {
    std::stringstream iss;
    std::vector<int> tmp_charge_mass;
    iss << line;
    getline(iss, token, ' ');
    name_vec.push_back(token);
    getline(iss, token, ' ');
    tmp_charge_mass.push_back(atoi(token.c_str())); //z
    getline(iss, token, ' ');
    tmp_charge_mass.push_back(atoi(token.c_str()));//a
    charge_mass_numbers.push_back(tmp_charge_mass);
    getline(iss, token, ' ');
    c_in_vec.push_back(atof(token.c_str())); //M_a
    getline(iss, token, ' ');
    c_out_vec.push_back(atof(token.c_str())); //m_b
    getline(iss, token, ' ');
    I_vec.push_back(atof(token.c_str())); //I_a
  }
  file.close();
  std::vector<std::vector<Material>> materials_each_secondary;
  for (size_t i = 0; i!=name_vec.size(); i++){
      std::vector<Atom> atoms;
      std::string name;
      double a;
      int z;
      file.open("./Materials/atoms.txt");
      while (getline(file, line)) {
        std::vector<std::string> file_vec_ne;
        std::stringstream iss;
        iss << line;
        getline(iss, token, ' ');
        name = token.c_str();
        for(size_t j = 0; j!=name_vec.size(); j++){
          file_vec_ne.push_back("./Splines/" + name + "_ne_energyangle_cdf_"+name_vec[i]+"_"+name_vec[j]+".txt");
        }
        getline(iss, token, ' ');
        z = atoi(token.c_str());
        getline(iss, token, ' ');
        a = atoi(token.c_str());
        if (name == "hydrogen" && name_vec[i] == "proton") {
          Atom tmp(a, z, charge_mass_numbers[i][1],charge_mass_numbers[i][0], "./Splines/" + name + "_el_ruth_cross_sec_"+name_vec[i]+".txt",
                  rutherford_cutoff, backscatter_cutoff); // special case for proton-hydrogen, no non-elastic, special elastic
          atoms.push_back(tmp);
        } else if (name == "hydrogen") {
              // special case for hydrogen target, no non-elastic
          Atom tmp(a, z, charge_mass_numbers[i][1],charge_mass_numbers[i][0], "./Splines/" + name + "_el_ruth_cross_sec_"+name_vec[i]+".txt",
                  rutherford_cutoff);
          atoms.push_back(tmp);
        } else if (name == name_vec[i]) { // special case for self-interaction
          Atom tmp(a, z, charge_mass_numbers[i][1],charge_mass_numbers[i][0], charge_mass_numbers, I_vec[i], I_vec, c_in_vec[i], c_out_vec,
            "./Splines/" + name + "_ne_rate_"+name_vec[i]+".txt",
            "./Splines/" + name + "_el_ruth_cross_sec_"+name_vec[i]+".txt",
            file_vec_ne, rutherford_cutoff, backscatter_cutoff);
          atoms.push_back(tmp);
        } else {
          Atom tmp(a, z, charge_mass_numbers[i][1],charge_mass_numbers[i][0], charge_mass_numbers, I_vec[i], I_vec, c_in_vec[i], c_out_vec,
            "./Splines/" + name + "_ne_rate_"+name_vec[i]+".txt",
            "./Splines/" + name + "_el_ruth_cross_sec_"+name_vec[i]+".txt",
            file_vec_ne, rutherford_cutoff);
          atoms.push_back(tmp);
        }
      }
      file.close();
      std::vector<std::string> material_names;
      file.open("./Materials/materials.txt");
      while (getline(file, line)) {
        material_names.push_back(line);
      }
      file.close();

      std::vector<Material> materials(material_names.size());
      for (unsigned int i = 0; i < material_names.size(); i++) {
        materials[i].read_material("./Materials/" + material_names[i], atoms);
      }
      materials_each_secondary.push_back(materials);
  }
  double nozzle_radius, initial_x_sd;
  std::vector<double> x(3, 0), w(2, 0);
  w[0] = M_PI / 2;
  cfg.lookupValue("nozzle_radius", nozzle_radius);
  cfg.lookupValue("initial_x_sd", initial_x_sd);

  double initial_e_mean, initial_e_sd;
  cfg.lookupValue("initial_e_mean", initial_e_mean);
  cfg.lookupValue("initial_e_sd", initial_e_sd);

  double dt, absorption_e;
  int nrep;
  cfg.lookupValue("step_size", dt);
  cfg.lookupValue("absorption_energy", absorption_e);
  cfg.lookupValue("replicates", nrep);

  const libconfig::Setting &root = cfg.getRoot();
  std::vector<double> change_points_x(root["change_points_x"].getLength() + 2);
  change_points_x[0] = -DBL_MAX;
  for (unsigned int i = 0; i < change_points_x.size() - 2; i++) {
    change_points_x[i + 1] = root["change_points_x"][i];
  }
  change_points_x[change_points_x.size() - 1] = DBL_MAX;
  std::vector<double> change_points_y(root["change_points_y"].getLength());
  for (unsigned int i = 0; i < change_points_y.size(); i++) {
    change_points_y[i] = root["change_points_y"][i];
  }
  std::vector<int> material_tmp(2, 0);
  std::vector<std::vector<int>> interval_materials(
      root["interval_materials_hi"].getLength(), material_tmp);
  for (unsigned int i = 0; i < interval_materials.size(); i++) {
    interval_materials[i][0] = root["interval_materials_lo"][i];
    interval_materials[i][1] = root["interval_materials_hi"][i];
  }
  double grid_dx;
  cfg.lookupValue("grid_dx", grid_dx);
  std::string out_path;
  cfg.lookupValue("out_path", out_path);
  particle_bank p_bank(initial_e_mean + 3 * initial_e_sd, dt, absorption_e, materials_each_secondary, charge_mass_numbers,
                  change_points_x, change_points_y, interval_materials);
  //int n = p.energy.size();
  //Grid grid(n * grid_dx / dt, grid_dx);

for (int i = 0; i < nrep; i++) {
    p_bank.add_particle(initial_e_mean,x,w,0);
    if (i % 500 == 490){
      p_bank.evaluate_particles_proton_neutron_storage(gen);
    }
  }
  p_bank.evaluate_particles_proton_neutron_storage(gen);
  nlohmann::json json_output = nlohmann::json::array();
  for (size_t j = 0; j!=p_bank.particles[1].size(); j++){
      nlohmann::json neutron_info;
      neutron_info["id"]=j;
      neutron_info["energy_MeV"]=p_bank.particles[1][j].energy;
      neutron_info["position_cm_x"]=p_bank.particles[1][j].position[0];
      neutron_info["position_cm_y"]=p_bank.particles[1][j].position[1];
      neutron_info["position_cm_z"]=p_bank.particles[1][j].position[2];
      neutron_info["angle_x"]=sin(p_bank.particles[1][j].angle[0])*cos(p_bank.particles[1][j].angle[1]);
      neutron_info["angle_y"]=sin(p_bank.particles[1][j].angle[0])*sin(p_bank.particles[1][j].angle[1]);
      neutron_info["angle_z"]=cos(p_bank.particles[1][j].angle[0]);
      json_output.push_back(neutron_info);
  }
  std::ofstream outputfile("Output/output_neutrons.json");
  outputfile << json_output.dump(8);
  outputfile.close();
  gsl_rng_free(gen);
  return 1;
}
