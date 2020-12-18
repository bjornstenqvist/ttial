#pragma once
#include "auxiliary.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <vector>
#include <sstream>
#include <cmath>

/**
 * @brief A class for input data
 */
class InputData {

    public:
        int R, C, N, seed, model_nbr;
        double d, s, g, t, omega, S_mv, S_mh, S_bv, S_bh, S_mtv, S_mth, D_mtv, D_mth, D_mv, D_mh, D_bv, D_bh, DLambda, z_break, height, width, c_out, S_out, evap_time, evap_left; // DLambda superfluous?
        bool load_external;
        std::string input_folder, output_folder, output_file;
        const double inf = std::numeric_limits<double>::infinity();

        int time_steps, sample;
        double time_periodic;
        double dt;
        bool evaporate, display;

        /**
         * @brief Load input data from file
         * @param filename Name of input data file
         */
        void loadInputData(std::string filename) {
            std::string line;
            std::ifstream myfile (filename);
            R = -1;
            C = -1;
            N = -1;
            seed = -1; // random seed
            model_nbr = -1;

            d = -0.1;
            s = -0.1;
            g = -0.1;
            t = -0.1;
            omega = -0.1; // random off-set
            S_mv = -0.1;
            S_mh = -0.1;
            S_bv = -0.1;
            S_bh = -0.1;
            S_mtv = -0.1;
            S_mth = -0.1;
            D_mtv = -0.1;
            D_mth = -0.1;
            D_mv = -0.1;
            D_mh = -0.1;
            D_bv = -0.1;
            D_bh = -0.1;
            DLambda = -0.1;
            z_break = -0.1;
            height = -0.1;
            width = -0.1;
            evap_time = -0.1;
            evap_left = -0.1;

            S_out = -0.1;
            c_out = -0.1;

            load_external = false;

            input_folder = "";
            output_folder = "";
            output_file = "output.txt";

            time_steps = -1;
            sample = -1;

            time_periodic = -0.1;
            dt = -0.1;

            evaporate = false;
            display = true;

            if (myfile.is_open())
            {
                while ( getline (myfile,line) )
                {
                    std::vector<std::string> result;
                    std::istringstream iss(line);
                    for(std::string line; iss >> line; )
                        result.push_back(line);

                    if (result.at(0).compare("Nr") == 0)
                        R = std::stoi(result.at(1));
                    if (result.at(0).compare("Nc") == 0)
                        C = std::stoi(result.at(1));
                    if (result.at(0).compare("N") == 0)
                        N = std::stoi(result.at(1));
                    if (result.at(0).compare("seed") == 0)
                        seed = std::stoi(result.at(1));
                    if (result.at(0).compare("model_nbr") == 0)
                        model_nbr = std::stoi(result.at(1));
                    if (result.at(0).compare("d") == 0)
                        d = std::stod(result.at(1));
                    if (result.at(0).compare("s") == 0)
                        s = std::stod(result.at(1));
                    if (result.at(0).compare("g") == 0)
                        g = std::stod(result.at(1));
                    if (result.at(0).compare("t") == 0)
                        t = std::stod(result.at(1));
                    if (result.at(0).compare("omega") == 0)
                        omega = std::stod(result.at(1));
                    if (result.at(0).compare("S_mv") == 0)
                        S_mv = std::stod(result.at(1));
                    if (result.at(0).compare("S_mh") == 0)
                        S_mh = std::stod(result.at(1));
                    if (result.at(0).compare("S_bv") == 0)
                        S_bv = std::stod(result.at(1));
                    if (result.at(0).compare("S_bh") == 0)
                        S_bh = std::stod(result.at(1));
                    if (result.at(0).compare("S_mtv") == 0)
                        S_mtv = std::stod(result.at(1));
                    if (result.at(0).compare("S_mth") == 0)
                        S_mth = std::stod(result.at(1));
                    if (result.at(0).compare("D_mtv") == 0)
                        D_mtv = std::stod(result.at(1));
                    if (result.at(0).compare("D_mth") == 0)
                        D_mth = std::stod(result.at(1));
                    if (result.at(0).compare("D_mv") == 0)
                        D_mv = std::stod(result.at(1));
                    if (result.at(0).compare("D_mh") == 0)
                        D_mh = std::stod(result.at(1));
                    if (result.at(0).compare("D_bv") == 0)
                        D_bv = std::stod(result.at(1));
                    if (result.at(0).compare("D_bh") == 0)
                        D_bh = std::stod(result.at(1));
                    if (result.at(0).compare("z_break") == 0)
                        z_break = std::stod(result.at(1));
                    if (result.at(0).compare("height") == 0)
                        height = std::stod(result.at(1));
                    if (result.at(0).compare("width") == 0)
                        width = std::stod(result.at(1));
                    if (result.at(0).compare("S_out") == 0)
                        S_out = std::stod(result.at(1));
                    if (result.at(0).compare("c_out") == 0)
                        c_out = std::stod(result.at(1));
                    if (result.at(0).compare("load_external") == 0)
                        if(!result.at(1).compare("true"))
                            load_external = true;
                    if (result.at(0).compare("input_folder") == 0)
                        input_folder = result.at(1)+"/";
                    if (result.at(0).compare("output_folder") == 0)
                        output_folder = result.at(1)+"/";
                    if (result.at(0).compare("output_file") == 0)
                        output_file = result.at(1);

                    if (result.at(0).compare("time_steps") == 0)
                        time_steps = std::stoi(result.at(1));
                    if (result.at(0).compare("sample") == 0)
                        sample = std::stoi(result.at(1));
                    if (result.at(0).compare("time_periodic") == 0)
                        time_periodic = std::stod(result.at(1));
                    if (result.at(0).compare("dt") == 0)
                        dt = std::stod(result.at(1));
                    if (result.at(0).compare("evaporate") == 0)
                        if(!result.at(1).compare("true"))
                            evaporate = true;
                    if (result.at(0).compare("display") == 0)
                        if(!result.at(1).compare("false"))
                            display = false;
                    if (result.at(0).compare("evap_time") == 0)
                        evap_time = std::stod(result.at(1));
                    if (result.at(0).compare("evap_left") == 0)
                        evap_left = std::stod(result.at(1));
                }
                myfile.close();
            }
            else {
                std::cerr << "Unable to open input-file\n";
                exit(EXIT_FAILURE);
            }

            if(evaporate && ( evap_time < 0.0 || evap_left < 0.0 ) ) {
                std::cerr << "Could not load 'evap_time' or 'evap_left'\n";
                exit(EXIT_FAILURE);
            }

            if(S_out < 0.0) {
                std::cerr << "Could not load solubility outside membrane 'S_out'\n";
                exit(EXIT_FAILURE);
            }

            if(c_out < 0.0) {
                std::cerr << "Could not load concentration outside membrane 'c_out'\n";
                exit(EXIT_FAILURE);
            }

            if(seed < 0)
                seed = time(NULL);
            DLambda = c_out / S_out;
        }

        /**
         * @brief Constructor
         */
        InputData(std::string filename) { loadInputData(filename); }
};
