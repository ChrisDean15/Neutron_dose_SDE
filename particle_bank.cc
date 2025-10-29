#include "material.cc"
#include "particle_beam.cc"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <vector>
#include <map>

#ifndef PBA
#define PBA
struct initial_condition {
    initial_condition(double energy0, std::vector<double> position0, std::vector<double> angle0) :  
    energy(energy0), position(position0), angle(angle0){
    }
    double energy;
    std::vector<double> position;
    std::vector<double> angle;  
};
struct particle_bank {
    particle_bank(const double e0_m_3sd, const double dt0, const double absorption_e0, std::vector<std::vector<Material>> materials0, std::vector<std::vector<int>> charge_mass_numbers0, std::vector<double> change_points0, std::vector<int> interval_materials0)
        : dt(dt0),absorption_e(absorption_e0), change_points(change_points0), interval_materials(interval_materials0), particles(materials0.size()), materials(materials0), charge_mass_numbers(charge_mass_numbers0),
          particle_traject(charge_mass_numbers0[0][0], charge_mass_numbers0[0][1], e0_m_3sd, dt0, absorption_e0, change_points0, interval_materials0, materials0[0]) {}
    void add_particle(double &energy0, std::vector<double> &position0, std::vector<double> &angle0, int label) 
    {
        particles[label].push_back(initial_condition(energy0,position0,angle0));
    }
    void evaluate_particles_proton_neutron_storage(gsl_rng *gen){
        while (!particles[0].empty()) {
            particle_traject.reset(charge_mass_numbers[0][1], charge_mass_numbers[0][0], particles[0].back().energy, 
            particles[0].back().position, particles[0].back().angle);
            particles[0].pop_back();
            std::vector<std::vector<std::vector<double>>> secondary_particles;
            particle_traject.simulate(dt, absorption_e, change_points, interval_materials,
                     materials[0], gen, secondary_particles);
            for(unsigned int i=0;i!=secondary_particles.size();i++){
                for(unsigned int j=0;j!=secondary_particles[i].size();j++){
                    double sec_energy=secondary_particles[i][j][0];
                    std::vector<double> sec_position={particle_traject.x[particle_traject.ix-1][0], particle_traject.x[particle_traject.ix-1][1], particle_traject.x[particle_traject.ix-1][2]};
                    std::vector<double> sec_angle={secondary_particles[i][j][1], secondary_particles[i][j][2]};
                    add_particle(sec_energy, sec_position, sec_angle, i);
                    if (i==1) {
                        std::cout << "Neutron energy: " << sec_energy << " MeV at position (" << sec_position[0] << ", " << sec_position[1] << ", " << sec_position[2] << ") cm" << std::endl;
                        std::cout << "Neutron angle: " << sec_angle[0] << " rad polar, " << sec_angle[1] << " rad azimuthal" << std::endl;
                    }
                }
            }
        }
    }
    double dt, absorption_e;
    std::vector<double> change_points;
    std::vector<int> interval_materials;
    std::vector<std::vector<initial_condition>> particles;
    std::vector<std::vector<Material>> materials;
    std::vector<std::vector<int>> charge_mass_numbers;
    particle_path particle_traject;
};
#endif