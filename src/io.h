#pragma once
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <vector>
#include <sstream>
#include <cmath>
#include "auxiliary.h"

/**
 * @brief A class for input data
 */
class InputData {

    public:
        int R, C, N, seed, model_nbr, max_iterations;
        double d, d_comp, s, g, t, omega, S_mv, S_mh, S_bv, S_bh, S_comp_bv, S_comp_bh, S_mtv, S_mth, D_mtv, D_mth, D_mv, D_mh, D_bv, D_bh, D_comp_bv, D_comp_bh, DLambda, z_break, height, width, c_out, S_out, evap_time, evap_left, max_error; // DLambda superfluous?
        bool load_external; //!< Load external input?
        std::string input_folder, output_folder, output_file;
        const double inf = std::numeric_limits<double>::infinity();

        std::vector<double> recX, recY, recW, recH, recVD, recVS, recHD, recHS;

        int time_steps, sample;
        double time_periodic; //!< Time between applying updated boundary conditions
        double dt; //!< Time-step
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
            max_iterations = -1;

            d = -0.1;
            d_comp = -0.1;
            s = -0.1;
            g = -0.1;
            t = -0.1;
            omega = -0.1; // random off-set
            S_mv = -0.1;
            S_mh = -0.1;
            S_bv = -0.1;
            S_bh = -0.1;
            S_comp_bv = -0.1;
            S_comp_bh = -0.1;
            S_mtv = -0.1;
            S_mth = -0.1;
            D_mtv = -0.1;
            D_mth = -0.1;
            D_mv = -0.1;
            D_mh = -0.1;
            D_bv = -0.1;
            D_bh = -0.1;
            D_comp_bv = -0.1;
            D_comp_bh = -0.1;
            DLambda = -0.1;
            z_break = -0.1;
            height = -0.1;
            width = -0.1;
            evap_time = -0.1;
            evap_left = -0.1;
            max_error = -0.1;

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

            recX.resize(0);
            recY.resize(0);
            recW.resize(0);
            recH.resize(0);
            recVD.resize(0);
            recVS.resize(0);
            recHD.resize(0);
            recHS.resize(0);


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
                    if (result.at(0).compare("max_iterations") == 0)
                        max_iterations = std::stoi(result.at(1));
                    if (result.at(0).compare("d") == 0)
                        d = std::stod(result.at(1));
                    if (result.at(0).compare("d_comp") == 0)
                        d_comp = std::stod(result.at(1));
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
                    if (result.at(0).compare("S_comp_bv") == 0)
                        S_comp_bv = std::stod(result.at(1));
                    if (result.at(0).compare("S_comp_bh") == 0)
                        S_comp_bh = std::stod(result.at(1));
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
                    if (result.at(0).compare("D_comp_bv") == 0)
                        D_comp_bv = std::stod(result.at(1));
                    if (result.at(0).compare("D_comp_bh") == 0)
                        D_comp_bh = std::stod(result.at(1));
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
                    if (result.at(0).compare("max_error") == 0)
                        max_error = std::stod(result.at(1));

                    if (result.at(0).compare("recX") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recX.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recY") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recY.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recW") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recW.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recH") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recH.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recVS") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recVS.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recVD") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recVD.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recHS") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recHS.push_back(std::stod(result.at(k)));
                    }
                    if (result.at(0).compare("recHD") == 0) {
                        for(unsigned int k = 1; k < result.size(); k++)
                            recHD.push_back(std::stod(result.at(k)));
                    }




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

            writeLoadedData();
        }

        void writeVectorToLoadedValues(std::string filename, std::vector<double> vec_out, std::string name) {
            for(unsigned int k = 0; k < vec_out.size(); k++)
                name += " "+to_string_precision(vec_out.at(k));
            if (vec_out.size() > 0) appendDataToFile(filename,name+"\n");
        }

        void writeLoadedData() {
            std::string filename = "loaded_and_available_data.txt";
            std::remove(filename.c_str()); // remove any eventual old output-file
            if (R >= 0) appendDataToFile(filename,"Nr "+to_string_precision(R)+"\n");
            if (C >= 0) appendDataToFile(filename,"Nc "+to_string_precision(C)+"\n");
            if (N >= 0) appendDataToFile(filename,"N "+to_string_precision(N)+"\n");
            if (seed >= 0) appendDataToFile(filename,"seed "+to_string_precision(seed)+"\n");
            if (model_nbr >= 0) appendDataToFile(filename,"model_nbr "+to_string_precision(model_nbr)+"\n");
            if (max_iterations >= 0) appendDataToFile(filename,"max_iterations "+to_string_precision(max_iterations)+"\n");
            if (d >= -0.05) appendDataToFile(filename,"d "+to_string_precision(d)+"\n");
            if (d_comp >= -0.05) appendDataToFile(filename,"d_comp "+to_string_precision(d_comp)+"\n");
            if (s >= -0.05) appendDataToFile(filename,"s "+to_string_precision(s)+"\n");
            if (g >= -0.05) appendDataToFile(filename,"g "+to_string_precision(g)+"\n");
            if (t >= -0.05) appendDataToFile(filename,"t "+to_string_precision(t)+"\n");
            if (omega >= -0.05) appendDataToFile(filename,"omega "+to_string_precision(omega)+"\n");
            if (S_mv >= -0.05) appendDataToFile(filename,"S_mv "+to_string_precision(S_mv)+"\n");
            if (S_mh >= -0.05) appendDataToFile(filename,"S_mh "+to_string_precision(S_mh)+"\n");
            if (S_bv >= -0.05) appendDataToFile(filename,"S_bv "+to_string_precision(S_bv)+"\n");
            if (S_bh >= -0.05) appendDataToFile(filename,"S_bh "+to_string_precision(S_bh)+"\n");
            if (S_comp_bv >= -0.05) appendDataToFile(filename,"S_comp_bv "+to_string_precision(S_comp_bv)+"\n");
            if (S_comp_bh >= -0.05) appendDataToFile(filename,"S_comp_bh "+to_string_precision(S_comp_bh)+"\n");
            if (S_mtv >= -0.05) appendDataToFile(filename,"S_mtv "+to_string_precision(S_mtv)+"\n");
            if (S_mth >= -0.05) appendDataToFile(filename,"S_mth "+to_string_precision(S_mth)+"\n");
            if (D_mtv >= -0.05) appendDataToFile(filename,"D_mtv "+to_string_precision(D_mtv)+"\n");
            if (D_mth >= -0.05) appendDataToFile(filename,"D_mth "+to_string_precision(D_mth)+"\n");
            if (D_mv >= -0.05) appendDataToFile(filename,"D_mv "+to_string_precision(D_mv)+"\n");
            if (D_mh >= -0.05) appendDataToFile(filename,"D_mh "+to_string_precision(D_mh)+"\n");
            if (D_bv >= -0.05) appendDataToFile(filename,"D_bv "+to_string_precision(D_bv)+"\n");
            if (D_bh >= -0.05) appendDataToFile(filename,"D_bh "+to_string_precision(D_bh)+"\n");
            if (D_comp_bv >= -0.05) appendDataToFile(filename,"D_comp_bv "+to_string_precision(D_comp_bv)+"\n");
            if (D_comp_bh >= -0.05) appendDataToFile(filename,"D_comp_bh "+to_string_precision(D_comp_bh)+"\n");
            if (z_break >= -0.05) appendDataToFile(filename,"z_break "+to_string_precision(z_break)+"\n");
            if (height >= -0.05) appendDataToFile(filename,"height "+to_string_precision(height)+"\n");
            if (width >= -0.05) appendDataToFile(filename,"width "+to_string_precision(width)+"\n");
            if (c_out >= -0.05) appendDataToFile(filename,"c_out "+to_string_precision(c_out)+"\n");
            if (S_out >= -0.05) appendDataToFile(filename,"S_out "+to_string_precision(S_out)+"\n");
            if (evap_time >= -0.05) appendDataToFile(filename,"evap_time "+to_string_precision(evap_time)+"\n");
            if (evap_left >= -0.05) appendDataToFile(filename,"evap_left "+to_string_precision(evap_left)+"\n");
            if (max_error >= -0.05) appendDataToFile(filename,"max_error "+to_string_precision(max_error)+"\n");

            if (time_steps >= 0) appendDataToFile(filename,"time_steps "+to_string_precision(time_steps)+"\n");
            if (sample >= 0) appendDataToFile(filename,"sample "+to_string_precision(sample)+"\n");
            if (time_periodic >= -0.05) appendDataToFile(filename,"time_periodic "+to_string_precision(time_periodic)+"\n");
            if (dt >= -0.05) appendDataToFile(filename,"dt "+to_string_precision(dt)+"\n");

            writeVectorToLoadedValues(filename,recX,"recX");
            writeVectorToLoadedValues(filename,recY,"recY");
            writeVectorToLoadedValues(filename,recW,"recW");
            writeVectorToLoadedValues(filename,recH,"recH");
            writeVectorToLoadedValues(filename,recVD,"recVD");
            writeVectorToLoadedValues(filename,recVS,"recVS");
            writeVectorToLoadedValues(filename,recHD,"recHD");
            writeVectorToLoadedValues(filename,recHS,"recHS");

            if (load_external)
                appendDataToFile(filename,"load_external true\n");
            else
                appendDataToFile(filename,"load_external false\n");

            if (evaporate)
                appendDataToFile(filename,"evaporate true\n");
            else
                appendDataToFile(filename,"evaporate false\n");

            if (display)
                appendDataToFile(filename,"display true\n");
            else
                appendDataToFile(filename,"display false\n");

            appendDataToFile(filename,"input_folder "+input_folder+"\n");
            appendDataToFile(filename,"output_folder "+output_folder+"\n");
            appendDataToFile(filename,"output_file "+output_file+"\n");
        }

        /**
         * @brief Constructor
         */
        InputData(std::string filename) { loadInputData(filename); }
};
