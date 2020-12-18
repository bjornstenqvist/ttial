#pragma once
#include "auxiliary.h"
#include "nodes.h"
#include "visual.h"

/**
 * @brief Solves dc/dt = nabla dot j
 */
void getFlux(mat &lambda, mat &mesh_hor, mat &mesh_ver, double dx, double dy, double dt, mat &j_ver_up, mat &j_ver_down, mat &j_hor_left, mat &j_hor_right) {
    int R = lambda.rows();
    int C = lambda.cols();
    assert(lambda.rows() == R);
    assert(mesh_hor.rows() == R-2);
    assert(mesh_ver.rows() == R-1);
    assert(lambda.cols() == C);
    assert(mesh_hor.cols() == C);
    assert(mesh_ver.cols() == C);
    //mat delta_conc = mat::Zero(R,C);

    for(int c = 0; c < C; c++)
        j_ver_down(0,c) -= ( lambda(0,c) - lambda(1,c) ) / dy * mesh_ver(0,c);

    for(int r = 1; r < R-1; r++) {
        j_ver_up(r,0) += ( lambda(r,0) - lambda(r-1,0) ) / dy * mesh_ver(r-1,0);
        j_ver_down(r,0) -= ( lambda(r,0) - lambda(r+1,0) ) / dy * mesh_ver(r,0);
        j_hor_right(r,0) += ( lambda(r,0) - lambda(r,C-1) ) / dx * mesh_hor(r-1,C-1);
        j_hor_left(r,0) -= ( lambda(r,0) - lambda(r,1) ) / dx * mesh_hor(r-1,0);

        for(int c = 1; c < C-1; c++) {
            j_ver_up(r,c) += ( lambda(r,c) - lambda(r-1,c) ) / dy * mesh_ver(r-1,c);
            j_ver_down(r,c) -= ( lambda(r,c) - lambda(r+1,c) ) / dy * mesh_ver(r,c);
            j_hor_right(r,c) += ( lambda(r,c) - lambda(r,c-1) ) / dx * mesh_hor(r-1,c-1);
            j_hor_left(r,c) -= ( lambda(r,c) - lambda(r,c+1) ) / dx * mesh_hor(r-1,c);
        }

        j_ver_up(r,C-1) += ( lambda(r,C-1) - lambda(r-1,C-1) ) / dy * mesh_ver(r-1,C-1);
        j_ver_down(r,C-1) -= ( lambda(r,C-1) - lambda(r+1,C-1) ) / dy * mesh_ver(r,C-1);
        j_hor_right(r,C-1) += ( lambda(r,C-1) - lambda(r,C-2) ) / dx * mesh_hor(r-1,C-2);
        j_hor_left(r,C-1) -= ( lambda(r,C-1) - lambda(r,0) ) / dx * mesh_hor(r-1,C-1);
    }

    for(int c = 0; c < C; c++)
        j_ver_up(R-1,c) += ( lambda(R-1,c) - lambda(R-2,c) ) / dy * mesh_ver(R-2,c);

    //delta_conc *= dt;
    j_ver_up *= dt;
    j_ver_down *= dt;
    j_hor_left *= dt;
    j_hor_right *= dt;

    // node 'r' surrounded by mesh_ver 'r-1' and mesh_ver 'r'
    //return delta_conc;
}

// Give warning if negative
void correctNegative(mat &a) {
    for(int r = 0; r < a.rows(); r++)
        for(int c = 0; c < a.cols(); c++)
            if(a(r,c) < 0.0)
                a(r,c) = 0.0;
}

double evaporateVolume(double v_0, double vt1, double t, double t1, int k=1) {
    double A = std::pow(t1+1.0,double(k)) * (v_0 - vt1) / ( std::pow(t1+1.0,double(k)) - 1.0  );
    double B = v_0 - A;
    return ( A/std::pow(t+1.0,double(k)) + B );
}

void nonsteadystate(InputData ipd, mat varpi_hor, mat varpi_ver, mat s_nodes) {
    int R = s_nodes.rows();
    mat conc_t = 0.0*s_nodes;
    conc_t.row(0).setConstant(ipd.c_out);
    mat lambda_t = conc_t.cwiseProduct(s_nodes.cwiseInverse());

    vec mass_out = vec::Zero(ipd.time_steps);
    vec volume_change = vec::Zero(ipd.time_steps);

    double time_since_added = 0.0;
    double dx = ipd.width/double(lambda_t.cols());
    double dy = ipd.height/double(lambda_t.rows()-1);

    int cnt = 1;
    double dV = 1.0;
    double dV_above_init = dV;
    double dV_above_current = dV_above_init;
    writeMatrixToFile(ipd.output_folder+"conc_0.txt",conc_t);
    writeMatrixToFile(ipd.output_folder+"lambda_0.txt",lambda_t);

    ProgressBar pb(ipd.time_steps,ipd.display,70); // 70 is with of output
    for(int n = 0; n < ipd.time_steps; n++) {

        mat j_ver_up = mat::Zero(lambda_t.rows(),lambda_t.cols());
        mat j_ver_down = mat::Zero(lambda_t.rows(),lambda_t.cols());
        mat j_hor_left = mat::Zero(lambda_t.rows(),lambda_t.cols());
        mat j_hor_right = mat::Zero(lambda_t.rows(),lambda_t.cols());

        getFlux(lambda_t,varpi_hor,varpi_ver,dx,dy,ipd.dt,j_ver_up,j_ver_down,j_hor_left,j_hor_right);
        //conc_t += delta_conc;
        conc_t -= ( j_ver_up - j_ver_down ) / dy;
        conc_t -= ( j_hor_right - j_hor_left ) / dx;
        correctNegative(conc_t);
        mass_out(n) = conc_t.row(R-1).sum()*dV; // save flow out of system
        conc_t.row(R-1).setZero(); // set zero concentration at bottom boundary (i.e. sink)
        lambda_t = conc_t.cwiseProduct(s_nodes.cwiseInverse()); // update standard activity
        if( ( n + 1 ) % ipd.sample == 0) {
            writeMatrixToFile(ipd.output_folder+"conc_"+std::to_string(cnt)+".txt",conc_t);
            writeMatrixToFile(ipd.output_folder+"lambda_"+std::to_string(cnt)+".txt",lambda_t);

            mat j_ver = 0.5*(j_ver_up + j_ver_down);
            j_ver.row(0) = j_ver_down.row(0);
            j_ver.row(j_ver.rows()-1) = j_ver_up.row(j_ver.rows()-1);
            mat j_hor = 0.5*(j_hor_left + j_hor_right);
            writeMatrixToFile(ipd.output_folder+"j_ver_"+std::to_string(cnt)+".txt",j_ver);
            writeMatrixToFile(ipd.output_folder+"j_hor_"+std::to_string(cnt)+".txt",j_hor);
            cnt++;
        }

        time_since_added += ipd.dt;
        if(time_since_added > ipd.time_periodic) {
            time_since_added = 0.0;
            conc_t.row(0).setConstant(ipd.c_out); // set 'c_out' concentration at top boundary (i.e. reservoir)
            lambda_t = conc_t.cwiseProduct(s_nodes.cwiseInverse()); // update standard activity
            dV_above_current = dV_above_init;
        }

        if(ipd.evaporate) {
            conc_t.row(0) *= dV_above_current; // calculate mass (which does not evaporate)
            dV_above_current = evaporateVolume(dV_above_init,dV_above_init*ipd.evap_left,time_since_added,ipd.evap_time); // evaporate part of volume
            conc_t.row(0) /= dV_above_current; // calculate concentration based on conserved mass and new volume
        }
        volume_change(n) = dV_above_current;
        ++pb;
        if( n % 100 )
            pb.display();
    }
    pb.end_of_process();
    writeMatrixToFile(ipd.output_folder+"volume_change.txt",volume_change);
    writeMatrixToFile(ipd.output_folder+"mass_out.txt",mass_out);

}
