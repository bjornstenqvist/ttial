#pragma once
#include "auxiliary.h"
#include "nodes.h"

mat updateConc(mat &lambda, mat &mesh_hor, mat &mesh_ver, double dx, double dy, double dt) {
    int R = lambda.rows();
    int C = lambda.cols();
    assert(lambda.rows() == R);
    assert(mesh_hor.rows() == R-2);
    assert(mesh_ver.rows() == R-1);
    assert(lambda.cols() == C);
    assert(mesh_hor.cols() == C);
    assert(mesh_ver.cols() == C);
    mat delta_conc = mat::Zero(R,C);

    for(int c = 0; c < C; c++)
        delta_conc(0,c) -= ( lambda(0,c) - lambda(1,c) ) / dy * mesh_ver(0,c);

    for(int r = 1; r < R-1; r++) {
        delta_conc(r,0) -= ( lambda(r,0) - lambda(r-1,0) ) / dy * mesh_ver(r-1,0);
        delta_conc(r,0) -= ( lambda(r,0) - lambda(r+1,0) ) / dy * mesh_ver(r,0);
        delta_conc(r,0) -= ( lambda(r,0) - lambda(r,C-1) ) / dx * mesh_hor(r-1,C-1);
        delta_conc(r,0) -= ( lambda(r,0) - lambda(r,1) ) / dx * mesh_hor(r-1,0);

        for(int c = 1; c < C-1; c++) {
            delta_conc(r,c) -= ( lambda(r,c) - lambda(r-1,c) ) / dy * mesh_ver(r-1,c);
            delta_conc(r,c) -= ( lambda(r,c) - lambda(r+1,c) ) / dy * mesh_ver(r,c);
            delta_conc(r,c) -= ( lambda(r,c) - lambda(r,c-1) ) / dx * mesh_hor(r-1,c-1);
            delta_conc(r,c) -= ( lambda(r,c) - lambda(r,c+1) ) / dx * mesh_hor(r-1,c);
        }

        delta_conc(r,C-1) -= ( lambda(r,C-1) - lambda(r-1,C-1) ) / dy * mesh_ver(r-1,C-1);
        delta_conc(r,C-1) -= ( lambda(r,C-1) - lambda(r+1,C-1) ) / dy * mesh_ver(r,C-1);
        delta_conc(r,C-1) -= ( lambda(r,C-1) - lambda(r,C-2) ) / dx * mesh_hor(r-1,C-2);
        delta_conc(r,C-1) -= ( lambda(r,C-1) - lambda(r,0) ) / dx * mesh_hor(r-1,C-1);
    }

    for(int c = 0; c < C; c++)
        delta_conc(R-1,c) -= ( lambda(R-1,c) - lambda(R-2,c) ) / dy * mesh_ver(R-2,c);

    delta_conc *= dt;

    // node 'r' surrounded by mesh_ver 'r-1' and mesh_ver 'r'
    return delta_conc;
}

void correctNegative(mat &a) {
    for(int r = 0; r < a.rows(); r++)
        for(int c = 0; c < a.cols(); c++)
            if(a(r,c) < 0.0)
                a(r,c) = 0.0;
}

double evaporateVolume(double v_0, double v_inf, double t) {
    int k = 3;
    double volume = ( v_0 - v_inf ) / std::pow(t+1.0,double(k)) + v_inf;
    return volume;
}

void nonsteadystate(InputData ipd, mat s_nodes, mat varpi_hor, mat varpi_ver) {
    int R = s_nodes.rows();
    mat conc_t = 0.0*s_nodes;
    conc_t.row(0).setConstant(ipd.c_out);
    mat lambda_t = conc_t.cwiseProduct(s_nodes.cwiseInverse());

    vec mass_out = vec::Zero(ipd.time_steps);
    vec volume_change = vec::Zero(ipd.time_steps);

    double time_since_added = 0.0;
    double dx = ipd.width/double(lambda_t.cols());
    double dy = ipd.height/double(lambda_t.rows()-1);

    int cnt = 0;
    double dV = 1.0;
    double dV_above_init = dV;
    double dV_above_steady_state = 0.5*dV_above_init;
    double dV_above_current = dV_above_init;
    for(int n = 0; n < ipd.time_steps; n++) {

        mat delta_conc = updateConc(lambda_t,varpi_hor,varpi_ver,dx,dy,ipd.dt);
        conc_t += delta_conc;
        mass_out(n) = conc_t.row(R-1).sum()*dV; // save flow out of system
        conc_t.row(R-1).setZero(); // set zero concentration at bottom boundary (i.e. sink)
        correctNegative(conc_t);
        lambda_t = conc_t.cwiseProduct(s_nodes.cwiseInverse()); // update standard activity
        if (n % ipd.sample == 0) {
            writeMatrixToFile("Output/conc_"+std::to_string(cnt)+".txt",conc_t);
            writeMatrixToFile("Output/lambda_"+std::to_string(cnt)+".txt",lambda_t);
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
            dV_above_current = evaporateVolume(dV_above_init,dV_above_steady_state,time_since_added); // evaporate part of volume
            conc_t.row(0) /= dV_above_current; // calculate concentration based on conserved mass and new volume
        }
        volume_change(n) = dV_above_current;
    }
    writeMatrixToFile("Output/volume_change.txt",volume_change);
    writeMatrixToFile("Output/mass_out.txt",mass_out);

}
