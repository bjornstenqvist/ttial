#pragma once
#include "auxiliary.h"
#include "nodes.h"

/**
 * @brief Calculate the effective resistance for a mesh of vertical and horizontal resistors. Using V = R*I where V is the potential difference 
 * over the resistance R, where the current I flows.
 * @param res_hor mesh of horizontal resistors
 * @param res_ver mesh of vertical resistors
 * @param DeltaV potential difference over the entire mesh of vertical resistors
 * @param output_file name of output-file
 */
void calcNodeAbsActivities(mat &res_hor, mat &res_ver, double DeltaV, double dy, std::string output_folder, std::string init_guess, std::string output_file, int max_iterations, double max_error, bool display, int sample) {
    // initialize nodes
    std::vector<Node> nodes;
    int nbr_of_nodes = res_hor.rows()*res_hor.cols();
    nodes.resize(nbr_of_nodes);
    int cnt = 0;
    for(unsigned int r = 1; r <= res_hor.rows(); r++)
        for(unsigned int c = 1; c <= res_hor.cols(); c++) {
            Node node_tmp;
            node_tmp.initNode(res_hor,res_ver,c,r,res_hor.cols(),res_hor.rows());
            nodes.at(cnt++) = node_tmp;
        }

    mat Rinv = mat::Zero(nbr_of_nodes,nbr_of_nodes); // initialize Rinv matrix
    vec I = vec::Zero(nbr_of_nodes);                 // initialize I matrix
    fillMatricesFromNodes(nodes,res_hor.cols(),res_hor.rows(),Rinv,I,DeltaV); // fill matrices and vector

    // calculates potential - heart of program
    vec V;
    if(max_iterations > 0) {
        double error = 2.0*max_error;
        vec V_guess = loadVector(init_guess);
        if(V_guess.size() == 0) {
            V_guess = I*0.0;
            appendDataToFile(output_file,"initial guess in Conjugate Gradient Method is a zero-vector\n");
        } else {
            appendDataToFile(output_file,"initial guess in Conjugate Gradient Method has been loaded from file: "+init_guess+"\n");
        }
        V = conjugateGradientMethod(Rinv,V_guess,I,error,max_error,max_iterations,display,sample,output_file,output_folder);
        writeMatrixToFile(output_folder+"V_vec_last_iteration.txt",V);
        appendDataToFile(output_file,"method conjugate gradient (numerical), error "+to_string_precision((Rinv*V - I).array().square().sum()/double(Rinv.rows()))+"\n");
    } else {
        //mat R = Rinv.inverse();                             // more effective for small systems (I think...)
        mat R = splitInverse(Rinv);                           // more effective for large systems
        V = R*I; // gets potential throughout the mesh
        appendDataToFile(output_file,"method inverse (exact), error "+to_string_precision((Rinv*V - I).array().square().sum()/double(Rinv.rows()))+"\n");
    }

    mat V_nodes = convertVector2Matrix(V,res_hor.rows(),res_hor.cols()); // convert potential in nodes from vector- to matrix-shape
    V_nodes = addBoundariesToMatrix(V_nodes,DeltaV,0.0);                  // add boundaries to potential matrix
    writeMatrixToFile(output_folder+"V_matrix.txt",V_nodes);                  // write potential to file

    // assert steady-state, and prepares for calculation of effective resistance
    vec I_top = vec::Zero(res_ver.cols());   // current going into top node from below
    vec I_bot = vec::Zero(res_ver.cols());   // current going into bottom node from above
    for(unsigned int c = 0; c < res_ver.cols(); c++) {
        I_top(c) = ( V_nodes(0,c) - V_nodes(1,c) ) / dy / nodes.at(c).upperRes;
        I_bot(c) = ( V_nodes(res_ver.rows()-1,c) - V_nodes(res_ver.rows(),c) ) / dy / nodes.at( (res_hor.rows()-1)*res_hor.cols() + c).lowerRes;
    }
    double I_top_sum = I_top.sum();  // current going into system
    double I_bot_sum = I_bot.sum();  // current going out of system
    appendDataToFile(output_file,"I_bot_sum "+to_string_precision(I_bot_sum)+"\n");
    appendDataToFile(output_file,"I_top_sum "+to_string_precision(I_top_sum)+"\n");
    if( std::fabs(1.0-I_bot_sum/I_top_sum) < 1e-6 )  // makes sure flow out of mesh equals that into mesh, i.e. steady-state
        std::cerr << "flow into system is different from that out of system. Check 'I_bot_sum' and 'I_top_sum' in output-file for more information.\n";

    // calculate effective resistance
    double Reff = ( DeltaV - 0.0 ) / I_top_sum;
    appendDataToFile(output_file,"R_eff "+to_string_precision(Reff)+"\n"); // generate output to file
}

void steadystate(InputData ipd, mat varpi_hor, mat varpi_ver, mat s_hor, mat s_ver) {
    mat res_hor = varpi_hor.cwiseInverse(); // set horizontal resistance
    mat res_ver = varpi_ver.cwiseInverse(); // set vertical resistance
    writeMatrixToFile(ipd.output_folder+"res_hor.txt",res_hor);
    writeMatrixToFile(ipd.output_folder+"res_ver.txt",res_ver);

    double dy = ipd.height/double(ipd.R+1);
    calcNodeAbsActivities(res_hor,res_ver,ipd.DLambda,dy,ipd.output_folder,ipd.init_guess,ipd.output_file,ipd.max_iterations,ipd.max_error,ipd.display,ipd.sample); // perform main calculations, i.e. get potential in nodes
    // generate output and end program
    mat Lambda_mat = loadMatrix(ipd.output_folder+"V_matrix.txt"); // load output from 'calcNodeAbsActivities' which is written to file
    calcProp(varpi_hor,varpi_ver,s_hor,s_ver,Lambda_mat,ipd.height,ipd.width,ipd.c_out,ipd.S_out,ipd.output_folder,ipd.output_file); // calculate concentrations and fluxes
}
