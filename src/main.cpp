#include <cmath>
#include <iostream>
#include <chrono>
#include <sstream>
#include <fstream>

#include "TScheme.h"
#include "omp.h"
#include <mpi.h>

class MPIController : public TScheme::TSender {

    void send_and_receive_density(
            double *const global_density, const double *const density, const size_t local_size, const size_t frame_size,
            const size_t total_computer_count, size_t *const sizes
    ) const override {
        if (global_density != nullptr) {
            std::copy_n(density, local_size, global_density);
            size_t shift = local_size;
            for (size_t i = 1; i < total_computer_count; i++) {
                MPI_Recv(
                        global_density + shift, sizes[i] * frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                );
                shift += sizes[i] * frame_size;
            }
        } else {
            MPI_Send(density, local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    void send_and_receive_forces(
            const double *const global_ax, const double *const global_ay, const double *const global_az,
            const double *const ax, const double *const ay, const double *const az,
            const double *const next_ax, const double *const next_ay, const double *const next_az,
            const double *const prev_ax, const double *const prev_ay, const double *const prev_az,
            size_t local_size, size_t frame_size,
            size_t total_computer_count, size_t *const sizes
    ) const override {
        if (global_ax != nullptr) {
            std::copy_n(global_ax, local_size, ax);
            std::copy_n(global_ay, local_size, ay);
            std::copy_n(global_az, local_size, az);

            std::copy_n(global_ax + local_size, frame_size, next_ax);
            std::copy_n(global_ay + local_size, frame_size, next_ay);
            std::copy_n(global_az + local_size, frame_size, next_az);

            size_t shift = local_size;
            for (size_t i = 1; i < total_computer_count; i++) {
                MPI_Send(global_ax + shift, sizes[i] * frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(global_ay + shift, sizes[i] * frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(global_az + shift, sizes[i] * frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

                MPI_Send(global_ax + shift - frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(global_ay + shift - frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(global_az + shift - frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

                if (i + 1 == total_computer_count) {
                    MPI_Send(global_ax, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(global_ay, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(global_az, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                } else {
                    MPI_Send(global_ax + shift + sizes[i] * frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(global_ay + shift + sizes[i] * frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                    MPI_Send(global_az + shift + sizes[i] * frame_size, frame_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                }

                shift += sizes[i] * frame_size;

            }
            std::copy_n(global_ax + shift - frame_size, frame_size, prev_ax);
            std::copy_n(global_ay + shift - frame_size, frame_size, prev_ay);
            std::copy_n(global_az + shift - frame_size, frame_size, prev_az);
        } else {
            MPI_Recv(ax, local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(ay, local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(az, local_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Recv(next_ax, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(next_ay, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(next_az, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Recv(prev_ax, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(prev_ay, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(prev_az, frame_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    void send_previous_x(double *const f, size_t size, size_t computer_index) const override {
        MPI_Send(f, size, MPI_DOUBLE, computer_index, 0, MPI_COMM_WORLD);
    }

    void send_next_x(double *const f, size_t size, size_t computer_index) const override {
        MPI_Send(f, size, MPI_DOUBLE, computer_index, 0, MPI_COMM_WORLD);
    }

    void receive_next_x(double *const f, size_t size, size_t computer_index) const override {
        MPI_Recv(f, size, MPI_DOUBLE, computer_index, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    void receive_previous_x(double *const f, size_t size, size_t computer_index) const override {
        MPI_Recv(f, size, MPI_DOUBLE, computer_index, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    ~MPIController() override {

    }
};
//

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int myRank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    //#######################################//
    //## Physical parameters of the system ##//
    //#######################################//

    ///////////////////////////////////////////
    //ion temperature  (kelvin)
    const double Ti = 300;
    //electron temperature (kelvin)
    const double Te = 9000;
    // time between ion-neutral collisions in seconds
    const double wc = 361745;
    //equilibrium concentaration of ions [1/cm^3]
    const double ni = 3. * pow(10., 7.);
    //particle charge in electron units
    const double q = 10000;
    //electricity field in cgs units
    const double Eext = 0.257e-2;
    //ion mass
    const double mi = 6.6464 * pow(10., -23.);
    //In debay radious units
    //x-dimension of the computational box in Debye radious
    const double lx = 15;
    //x-dimension of the computational box in Debye radious
    const double ly = 6;
    //x-dimension of the computational box in Debye radious
    const double lz = 6;
    //coordinate-x of the dust particle
    const double r0x = 3;
    //coordinate-y of the dust particle
    const double r0y = 3;
    //coordinate-z of the dust particle
    const double r0z = 3;
    ///////////////////////////////////////////

    //########################################//
    //## Numerical parameters of the system ##//
    //########################################//

    ///////////////////////////////////////////
    //number of cells in cordinate space x
    const size_t nx = 10;
    //number of cells in cordinate space y
    const size_t ny = 10;
    //number of ells in cordinate space z
    const size_t nz = 10;
    //integration step in velocity space x
    const double hvx = 0.2;
    //integration step in velocity space y
    const double hvy = 0.2;
    //integration step in velocity space z
    const double hvz = 0.2;
    //maximum velocity yz
    const double vmaxyz = 1;
    //minimum velocity yz
    const double vminyz = -1;
    //integration step for calculation initial conditions
    const double hxi = 0.002;
    //integration limit for calculation initial conditions
    const double ximax = 5;
    ///////////////////////////////////////////


    //#######################################//
    //##     Define output frequency       ##//
    //#######################################//

    ///////////////////////////////////////////
    size_t ITmax = 3;
    size_t T_output = 20;
    ///////////////////////////////////////////


    //#######################################//
    //##         Inizialization            ##//
    //#######################################//

//    omp_set_num_threads(44);


    const auto start = std::chrono::high_resolution_clock::now();
    const time_t start_time_t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Start scheme crating at " << std::ctime(&start_time_t) << std::endl;

    TScheme::TWorkController scheme = TScheme::TBuilder()
            .set_electron_temperature(Te)
            .set_ion_temperature(Ti)
            .set_ion_concentration(ni)
            .set_dust_particle_charge(q)
            .set_external_electricity_field(Eext)
            .set_ion_neutral_collisions_frequency(wc)
            .set_ion_mass(mi)
            .set_box_split_step_count(nx, ny, nz)
            .set_box_size(lx, ly, lz)
            .set_particle_position(r0x, r0y, r0z)
            .set_velocity_step(hvx, hvy, hvz)
            .set_velocity_bound_for_yz(vminyz, vmaxyz)
            .set_integration_parameters(ximax, hxi)
            .set_computer_index(myRank)
            .set_max_computer_index(numprocs)
            .set_sender(new MPIController())
            .build_work_controller();

    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken by inizialisation: "
              << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0
              << " seconds" << std::endl;

    //####################################//
    //##         Computation            ##//
    //####################################//
    std::cout << "##############################################\n"
                 "##******* Comutations start ****************##\n"
                 "##############################################\n";
    const auto total_start = std::chrono::high_resolution_clock::now();
    const time_t total_start_time_t = std::chrono::system_clock::to_time_t(total_start);
    std::cout << "Start time: " << std::ctime(&total_start_time_t);
    for (size_t t = 0; t < ITmax; ++t) {
        std::cout << "# Step " << t << std::endl;
        const auto step_start = std::chrono::high_resolution_clock::now();
        scheme.next_step();

        if (myRank == 0) {
            std::cout << "Current dimmensionless time: " << t * scheme.get_dt() << std::endl;
            std::cout << "Current full charge in the computational box: " << scheme.get_full_charge() << std::endl;
            if (t % T_output == 0) {
                std::ofstream file;

                std::stringstream filename;
                filename << "./density_t" << t << ".dat";
                file.open(filename.str().c_str());
                scheme.write_density(file);
                file.close();

                std::stringstream filename2;
                filename2 << "./potential_t" << t << ".dat";
                file.open(filename2.str().c_str());
                scheme.write_potential(file);
                file.close();
            }
        }

        const auto step_stop = std::chrono::high_resolution_clock::now();
        std::cout << "Time taken by one step: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(step_stop - step_start).count() / 1000000.0
                  << " seconds" << std::endl;
    }
    const auto total_stop = std::chrono::high_resolution_clock::now();
    const time_t total_stop_time_t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Stop time: " << std::ctime(&total_stop_time_t);
    std::cout << "Time taken by scheme: "
              << std::chrono::duration_cast<std::chrono::microseconds>(total_stop - total_start).count() / 1000000.0
              << " seconds" << std::endl;

    MPI_Finalize();
    return 0;
}
