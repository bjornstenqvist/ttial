#include <chrono>
#include "src/io.h" 
#include "src/geometry.h" 
#include "src/auxiliary.h"

/**
 * @brief Calculate the effective resistance for a mesh of vertical and horizontal resistors. Using V = R*I where V is the potential difference 
 * over the resistance R, where the current I flows.
 * @param res_hor mesh of horizontal resistors
 * @param res_ver mesh of vertical resistors
 * @param DeltaV potential difference over the entire mesh of vertical resistors
 * @param output_file name of output-file
 */
void calcNodeAbsActivities(mat &res_hor, mat &res_ver, double DeltaV, double dy, std::string output_file) {

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

    // calculates inverse of Rinv - heart of program
    //mat R = Rinv.inverse();                             // more effective for small systems (I think...)
    mat R = splitInverse(Rinv);                           // more effective for large systems

    // gets potential throughout the mesh
    vec V = R*I;
    mat V_nodes = convertVector2Matrix(V,res_hor.rows(),res_hor.cols()); // convert potential in nodes from vector- to matrix-shape
    V_nodes = addBoundariesToMatrix(V_nodes,DeltaV,0.0);                  // add boundaries to potential matrix
    writeMatrixToFile("V_matrix.txt",V_nodes);                  // write potential to file

    // assert steady-state, and prepares for calculation of effective resistance
    vec I_top = vec::Zero(res_ver.cols());   // current going into top node from below
    vec I_bot = vec::Zero(res_ver.cols());   // current going into bottom node from above
    for(unsigned int c = 0; c < res_ver.cols(); c++) {
        I_top(c) = ( V_nodes(0,c) - V_nodes(1,c) ) / dy / nodes.at(c).upperRes;
        I_bot(c) = ( V_nodes(res_ver.rows()-1,c) - V_nodes(res_ver.rows(),c) ) / dy / nodes.at( (res_hor.rows()-1)*res_hor.cols() + c).lowerRes;
    }
    double I_top_sum = I_top.sum();  // current going into system
    double I_bot_sum = I_bot.sum();  // current going out of system
    assert(("flow into system is different than that out of system" && std::fabs(1.0-I_bot_sum/I_top_sum) < 1e-6));  // makes sure flow out of mesh equals that into mesh, i.e. steady-state

    // calculate effective resistance
    double Reff = ( DeltaV - 0.0 ) / I_top_sum;

    appendDataToFile(output_file,"I_bot_sum "+to_string_precision(I_bot_sum)+"\n");
    appendDataToFile(output_file,"I_top_sum "+to_string_precision(I_top_sum)+"\n");
    appendDataToFile(output_file,"Reff "+to_string_precision(Reff)+"\n"); // generate output to file
}

int main() {
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now(); // start timer
    std::string output_file = "output.txt"; // set name of output-file
    std::remove("output.txt"); // remove any eventual old output-file

    // initialize program
    InputData ipd("input.txt"); // load data in input-file
    srand(ipd.seed); // set seed

    // get geometry
    mat s_hor, s_ver, D_hor, D_ver, varpi_hor, varpi_ver, res_hor, res_ver; // initiate matrices
    double height, width; // initiate height and width of system
    Geometry geo;
    geo.getGeometry(s_hor,s_ver,D_hor,D_ver,height,width,ipd,output_file); // get geometry
    varpi_hor = s_hor.cwiseProduct(D_hor); // set horizontal momentum permittivity
    varpi_ver = s_ver.cwiseProduct(D_ver); // set vertical momentum permittivity
    res_hor = varpi_hor.cwiseInverse(); // set horizontal resistance
    res_ver = varpi_ver.cwiseInverse(); // set vertical resistance

    double dy = ipd.height/double(ipd.R);
    calcNodeAbsActivities(res_hor,res_ver,ipd.DLambda,dy,output_file); // perform main calculations, i.e. get potential in nodes

    // generate output and end program
    mat Lambda_mat = loadMatrix("V_matrix.txt"); // load output from 'calcNodeAbsActivities' which is written to file
    calcProp(varpi_hor,varpi_ver,s_hor,s_ver,Lambda_mat,height,width,output_file); // calculate concentrations and fluxes

    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now(); // stop timer
    double timed = double(std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count())/1000.0; // calculate elapsed time
    appendDataToFile(output_file,getTimeString(timed)); // write elapsed time to file

    std::cout << "Calculation successfully terminated!" << std::endl;
    return 0;
}