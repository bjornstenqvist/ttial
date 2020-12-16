#include <chrono>
#include "src/io.h" 
#include "src/geometry.h" 
#include "src/auxiliary.h"
#include "src/steadystate.h"

int main(int argc, const char **argv) {
    std::string input_file = "";
    if (argc > 1)
        input_file = argv[1];
    if(input_file.length() == 0) {
        std::cerr << "No input-file given. Exit program..." << std::endl;
        return 1;
    }

    // initialize program
    InputData ipd(input_file); // load data in input-file
    srand(ipd.seed); // set seed

    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now(); // start timer
    std::remove(ipd.output_file.c_str()); // remove any eventual old output-file

    // get geometry
    mat s_hor, s_ver, s_nodes, D_hor, D_ver, varpi_hor, varpi_ver; // initiate matrices
    Geometry geo;
    geo.getGeometry(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd); // get geometry
    varpi_hor = s_hor.cwiseProduct(D_hor); // set horizontal momentum permittivity
    varpi_ver = s_ver.cwiseProduct(D_ver); // set vertical momentum permittivity

    steadystate(ipd,varpi_hor,varpi_ver,s_hor,s_ver);

    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now(); // stop timer
    double timed = double(std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count())/1000.0; // calculate elapsed time
    appendDataToFile(ipd.output_file,getTimeString(timed)); // write elapsed time to file

    std::cout << "TTIAL successfully terminated!" << std::endl;
    return 0;
}