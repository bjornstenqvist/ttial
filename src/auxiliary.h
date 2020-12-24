#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <vector>

typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

void removeRow(mat &inm, unsigned int index) {
    unsigned int rows = inm.rows()-1;
    unsigned int cols = inm.cols();
    if( index < rows )
        inm.block(index,0,rows-index,cols) = inm.block(index+1,0,rows-index,cols);
    inm.conservativeResize(rows,cols);
}

void removeColumn(mat &inm, unsigned int index) {
    unsigned int rows = inm.rows();
    unsigned int cols = inm.cols()-1;
    if( index < cols )
        inm.block(0,index,rows,cols-index) = inm.block(0,index+1,rows,cols-index);
    inm.conservativeResize(rows,cols);
}

void insertRow(mat &inm, unsigned int index, double value=1.0) {
    unsigned int rows = inm.rows();
    inm.conservativeResize(rows+1, inm.cols());
    mat tmp = inm.bottomRows(rows-index);
    inm.bottomRows(rows-index -1) = tmp.topRows(tmp.rows()-1);
    //inm.row(index) *= 0.0;

    inm.row(index) = vec::Ones(inm.row(index).cols());
    inm.row(index) *= value;
}

void insertColumn(mat &inm, unsigned int index, double value=1.0) {
    unsigned int cols = inm.cols();
    inm.conservativeResize(inm.rows(), cols+1);
    mat tmp = inm.rightCols(cols-index);
    inm.rightCols(cols-index -1) = tmp.leftCols(tmp.cols()-1);

    inm.col(index) = vec::Ones(inm.col(index).rows());
    inm.col(index) *= value;
}

void expandMatrix(mat &inm, std::vector<bool> empty_vec) {
    assert(inm.rows() == inm.cols());
    mat ones = inm.cwiseProduct(inm.transpose().cwiseInverse());
    assert( 1e-9 > 1.0 - ones.minCoeff()/ones.maxCoeff());
    for(unsigned int k = 0; k < empty_vec.size(); k++) {
        if(empty_vec.at(k)) {
            insertColumn(inm,k);
            insertRow(inm,k);
        }
    }
    assert(inm.rows() == inm.cols());
}

std::vector<bool> reduceMatrix(mat &inm, vec inv) {
    assert(inm == inm.transpose());
    std::vector<bool> empty_vec(0);
    for(unsigned int r = 0; r < inm.rows(); r++) {
        double row_sum = inm.row(r).cwiseAbs().sum();
        if(row_sum < 1e-9 && inv(r) < 1e-9) {
            empty_vec.push_back(true);
        } else {
            empty_vec.push_back(false);
        }
    }

    for(int k = empty_vec.size()-1; k > -1; k--) {
        if(empty_vec.at(k)) {
            removeRow(inm,k);
            removeColumn(inm,k);
        }
    }
    assert(inm == inm.transpose());
    return empty_vec;
}

/**
 * @brief Gives shortest distance using PBC in one dimension.
 * @param x coordinate
 * @param Lx periodic length
 */
double PBC_1D(const double x, const double Lx) {
    if(x > Lx/2.0)
        return (x - Lx);
    if(x < -Lx/2.0)
        return (x + Lx);
    return x;
}

std::string to_string_precision(const double val, const int precision = 15)
{
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << val;
    return out.str();
}

/**
 * @brief Check such that matrices have compatible sizes for the mesh-approach. Assumes vertical matrices have one more row than horizontal ones, 
 * and all have identical number of columns.
 * @param Av vertical matrix of type A
 * @param Ah horizontal matrix of type A
 * @param Bv vertical matrix of type B
 * @param Bh horizontal matrix of type B
 */
bool compatibleSizes(const mat &Av, const mat &Ah, const mat &Bv, const mat &Bh) {
    if(Av.cols() != Ah.cols())
        return false;
    if(Ah.cols() != Bv.cols())
        return false;
    if(Bv.cols() != Bh.cols())
        return false;
    if(Av.rows() != Bv.rows())
        return false;
    if(Ah.rows() != Bh.rows())
        return false;
    if(Av.rows()-1 != Ah.rows())
        return false;
    return true;
}

/**
 * @brief Load matrix from file.
 * @note Assumes input is exactly of the form N x M, i.e. each colum has same size, and each row has same size.
 */
mat loadMatrix(const std::string filename) {
    std::string line;
    std::ifstream myfile (filename);

    std::vector<std::vector<double>> mat_v(0);
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            std::vector<std::string> result;
            std::istringstream iss(line);
            for(std::string line; iss >> line; )
                result.push_back(line);

            std::vector<double> result_d(0);
            for(unsigned int i = 0; i < result.size(); i++)
                result_d.push_back(std::stod(result.at(i)));
            mat_v.push_back(result_d);
        }
        myfile.close();
    } 

    mat mat0 = mat::Zero(mat_v.size(),mat_v.at(0).size());
    for(unsigned int i = 0; i < mat_v.size(); i++)
        for(unsigned int j = 0; j < mat_v.at(0).size(); j++)
            mat0(i,j) = mat_v.at(i).at(j);

    return mat0;
}

/**
 * @brief Generates a random number which is uniformely distributed between the input parameters.
 * @param minv lower bound of random value
 * @param maxv upper bound of random value
 */
double randomUnit(const double minv, const double maxv)
{
    double f = (double)rand() / RAND_MAX;
    return minv + f * (maxv - minv);
}


/**
 * @brief Appends a string of text to a file.
 * @param output_file name of file
 * @param text string to be appended
 */
void appendDataToFile(const std::string output_file, const std::string text) {
    std::ofstream outfile;
    outfile.open(output_file, std::ios_base::app);
    outfile << text; 
}



/**
 * @brief Converts a vector to matrix.
 * @param R rows in matrix
 * @param C columns in matrix
 * @note Assumes length of vector is R*C.
 * @note Reads first C values in vector as first row in matrix, and so forth.
 */
mat convertVector2Matrix(const vec &in, const int R, const int C) {
    assert(in.size() == R*C);
    mat out = mat::Zero(R,C);

    int cnt = 0;
    for(int r = 0; r < R; r++) {
        for(int c = 0; c < C; c++) {
            out(r,c) = in(cnt);
            cnt++;
        }
    }
    return out;
}

/**
 * @brief A more effective way o calculate the inverse of a matrix.
 * @note Does only work for square matrices
 * @todo Add reference to method
 */
mat splitInverse(const mat &mat0) {

    int R = mat0.rows();
    int C = mat0.cols();
    assert(R == C);

    if(R < 50)
        return mat0.inverse();

    int sr = int(double(R)/2.0);
    int sc = int(double(C)/2.0);


    mat Ei = mat0.topLeftCorner(sr,sc);
    mat F = mat0.topRightCorner(sr,C-sc);
    mat G = mat0.bottomLeftCorner(R-sr,sc);
    mat H = mat0.bottomRightCorner(R-sr,C-sc);


    //Ei = Ei.inverse();
    Ei = splitInverse(Ei);
    G = G*Ei;
    mat Si = H - G*F;
    //Si = Si.inverse();
    Si = splitInverse(Si);
    //Hi = Hi.inverse();
    //Hi = splitInverse(Hi);

    F = Ei*F*Si;
    mat mat_inverse = mat::Zero(R,C);
    mat_inverse.block(0,0,sr,sc) = Ei + F*G;
    mat_inverse.block(0,sc,sr,C-sc) = -F;
    mat_inverse.block(sr,0,R-sr,sc) = -Si*G;
    mat_inverse.block(sr,sc,R-sr,C-sc) = Si;

    return mat_inverse;
}

void writeMatrixToFile(const std::string filename, const mat &in) {
    std::ofstream file(filename);
    if (file.is_open())
    {
        for(int r = 0; r < in.rows(); r++) {
            for(int c = 0; c < in.cols(); c++)
                file << in(r,c) << ' ';
            file << std::endl;
        }
    }
}

std::vector<bool> findEmptyRowsInVector(const vec &vec0, const double limit=1e-12) {

    std::vector<bool> empty_vec;
    empty_vec.resize(vec0.size());
    std::fill(empty_vec.begin(), empty_vec.end(), true);
    for(int r = 0; r < vec0.size(); r++)
        if(std::fabs(vec0(r)) > limit)
            empty_vec.at(r) = false;
    return empty_vec;
}

std::vector<bool> findEmptyRowsInMatrix(const mat &mat0, const double limit=1e-12) {
    assert(mat0 == mat0.transpose());

    std::vector<bool> empty_vec;
    empty_vec.resize(mat0.rows());
    std::fill(empty_vec.begin(), empty_vec.end(), true);
    for(int r = 0; r < mat0.rows(); r++)
        for(int c = 0; c < mat0.cols(); c++)
            if(std::fabs(mat0(r,c)) > limit) {
                empty_vec.at(r) = false;
                break;
            }
    return empty_vec;
}

std::vector<bool> findEmptyRows(const mat &Rinv, const vec &Iex, const double limit=1e-12) {
    std::vector<bool> empty_mat = findEmptyRowsInMatrix(Rinv,limit);
    std::vector<bool> empty_vec = findEmptyRowsInVector(Iex,limit);
    assert(empty_mat.size() == empty_vec.size());

    std::vector<bool> empty_both;
    empty_both.resize(empty_mat.size());
    std::fill(empty_both.begin(), empty_both.end(), false);

    for(unsigned int r = 0; r < empty_both.size(); r++)
        if(empty_mat.at(r) && empty_vec.at(r))
            empty_both.at(r) = true;
    return empty_both;
}

std::string getTimeString(double timed) {
    int days = int(timed/86400.0);
    timed -= double(days)*86400.0;
    int hours = int(timed/3600.0);
    timed -= double(hours)*3600.0;
    int minutes = int(timed/60.0);
    timed -= double(minutes)*60.0;
    int seconds = int(timed);
    std::string str = "Time "+std::to_string(days)+"D_"+std::to_string(hours)+"H_"+std::to_string(minutes)+"M_"+std::to_string(seconds)+"S\n";
    return str;
}

/**
 * @brief Adds one row to the input matrix from above and one from below, with given values.
 * @param top value for each column in the top row of the new matrix
 * @param bottom value for each column in the bottom row of the new matrix
 */
mat addBoundariesToMatrix(const mat &in, const double top, const double bottom) {
    int R = in.rows();
    int C = in.cols();
    mat out = mat::Zero(R+2,C); // create an equal but two rows larger matrix
    out.row(0) = vec::Ones(C)*top; // add top-value to first row
    out.block(1,0,R,C) = in; // add matrix elements in the middle
    out.row(R+1) = vec::Ones(C)*bottom; // add bottom-value to last row
    return out;
}

void calcProp(mat &varpi_hor, mat &varpi_ver, mat &s_hor, mat &s_ver, mat &V_nodes, double height, double width, double c_out, double S_out, std::string output_folder, std::string output_file) {
    assert(compatibleSizes(varpi_ver,varpi_hor,s_ver,s_hor));
    assert(V_nodes.cols() == s_hor.cols());
    assert(V_nodes.rows() == s_hor.rows()+2);

    double dy = height/double(varpi_ver.rows());
    double dx = width/double(varpi_hor.cols());
    mat concV = mat::Zero(varpi_ver.rows(),varpi_ver.cols());
    mat concH = mat::Zero(varpi_hor.rows(),varpi_hor.cols());
    mat jv = mat::Zero(varpi_ver.rows(),varpi_ver.cols());
    mat jh = mat::Zero(varpi_hor.rows(),varpi_hor.cols());

    for(int r = 0; r < varpi_ver.rows(); r++)
        for(int c = 0; c < varpi_ver.cols(); c++) {
            concV(r,c) = 0.5 * ( V_nodes(r,c) + V_nodes(r+1,c) ) * s_ver(r,c);
            jv(r,c) = - varpi_ver(r,c) * ( V_nodes(r,c) - V_nodes(r+1,c) ) / dy;
        }
    for(int r = 0; r < varpi_hor.rows(); r++) {
        for(int c = 0; c < varpi_hor.cols()-1; c++) {
            concH(r,c) = 0.5 * ( V_nodes(r,c) + V_nodes(r,c+1) ) * s_hor(r,c);
            jh(r,c) = - varpi_hor(r,c) * ( V_nodes(r,c) - V_nodes(r,c+1) ) / dx;
        }
        concH(r,varpi_hor.cols()-1) = 0.5 * ( V_nodes(r,varpi_hor.cols()-1) + V_nodes(r,0) ) * s_hor(r,varpi_hor.cols()-1);
        jh(r,varpi_hor.cols()-1) = - varpi_hor(r,varpi_hor.cols()-1) * ( V_nodes(r,varpi_hor.cols()-1) - V_nodes(r,0) ) / dx;
    }

    mat jva, jha;
    jva.resize(jv.rows(),2);
    jha.resize(jh.rows(),2);
    jva.col(0) = vec::LinSpaced(jva.rows(),0.5*dy,height-0.5*dy);
    jha.col(0) = vec::LinSpaced(jha.rows(),dy,height-dy);
    jva.col(1) = jv.rowwise().mean();
    jha.col(1) = jh.cwiseAbs().rowwise().mean();

    double c_in = c_out * s_ver.row(0).mean()/S_out;
    appendDataToFile(output_file,"P_eff "+to_string_precision(-jva.col(1).mean()/c_out)+"\n"); // generate output to file
    appendDataToFile(output_file,"D_eff "+to_string_precision(-jva.col(1).mean()/c_in*height)+"\n"); // generate output to file
    appendDataToFile(output_file,"j_ver "+to_string_precision(jva.col(1).mean())+"\n"); // generate output to file
    appendDataToFile(output_file,"j_hor "+to_string_precision(jha.col(1).mean())+"\n"); // generate output to file
    writeMatrixToFile(output_folder+"j_ver_vector.txt", jva);
    writeMatrixToFile(output_folder+"j_hor_vector.txt", jha);
    writeMatrixToFile(output_folder+"j_ver_matrix.txt", jv);
    writeMatrixToFile(output_folder+"j_hor_matrix.txt", jh);
    writeMatrixToFile(output_folder+"conc_ver_matrix.txt", concV);
    writeMatrixToFile(output_folder+"conc_hor_matrix.txt", concH);
}
