#include <chrono>
#include "src/io.h" 
#include "src/geometry.h" 
#include "src/auxiliary.h"
#include "src/steadystate.h"
#include "src/nonsteadystate.h"

int main() {
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now(); // start timer
    std::string output_file = "output.txt"; // set name of output-file
    std::remove("output.txt"); // remove any eventual old output-file

    // initialize program
    InputData ipd("input.txt"); // load data in input-file
    srand(ipd.seed); // set seed

    // get geometry
    mat s_hor, s_ver, s_nodes, D_hor, D_ver, varpi_hor, varpi_ver; // initiate matrices
    Geometry geo;
    geo.getGeometry(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd,output_file); // get geometry
    varpi_hor = s_hor.cwiseProduct(D_hor); // set horizontal momentum permittivity
    varpi_ver = s_ver.cwiseProduct(D_ver); // set vertical momentum permittivity

    mat res_hor = varpi_hor.cwiseInverse(); // set horizontal resistance
    mat res_ver = varpi_ver.cwiseInverse(); // set vertical resistance
    steadystate(ipd,res_hor,res_ver,varpi_hor,varpi_ver,s_hor,s_ver,height,width,output_file);

    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now(); // stop timer
    double timed = double(std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count())/1000.0; // calculate elapsed time
    appendDataToFile(output_file,getTimeString(timed)); // write elapsed time to file

    std::cout << "TTIAL successfully terminated!" << std::endl;
    return 0;
}