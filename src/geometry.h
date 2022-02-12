#pragma once
#include "auxiliary.h"
#include "nodes.h"

/**
 * @brief A class for geometry 
 */
class Geometry
{
    public:
        void fillGeometryWithTunnels(mat &A_mat, mat &B_mat, double A_sca, double B_sca, int x_start);
        void fillGeometryWithCirclesHor(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height);
        void fillGeometryWithCirclesVer(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height);
        void fillGeometryWithCirclesNodes(mat &A_mat, std::vector<double> A_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height);
        void fillGeometryWithRectangles(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double dx, double dy, double L, double offset);
        void fillGeometryWithRectanglesNodes(mat &A_mat, std::vector<double> A_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height);
        void fillGeometryWithRectanglesHor(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height);
        void fillGeometryWithRectanglesVer(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height);
        void brickAndMortar(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void brickAndMortarEnglishBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void brickAndMortarFlemishBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void brickAndMortarMonkBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void brickAndMortarSussexBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void custom_made(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void brickAndMortarDualLinear(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name);
        void getGeometry(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, InputData &ipd);
        void checkInputForBrickAndMortar(InputData ipd);
        void getBricks(InputData ipd, std::vector<double> &xv, std::vector<double> &yv, std::vector<double> &dv, std::vector<double> &tv, std::vector<double> &Sv, std::vector<double> &Sh, std::vector<double> &Dv, std::vector<double> &Dh, int model, int K);
};

/**
 * @brief Fills geometry with tunnels
 * @param A_mat Matrix A
 * @param B_mat Matrix B
 * @param A_sca Value in tunnel in matrix A
 * @param B_sca Value in tunnel in matrix B
 * @param x_start Starting-column in matrix for tunnel
 */
void Geometry::fillGeometryWithTunnels(mat &A_mat, mat &B_mat, double A_sca, double B_sca, int x_start) {
    A_mat(0,x_start) = A_sca;
    B_mat(0,x_start) = B_sca;
    int dir = 0;

    for(int r = 1; r < A_mat.rows(); r++) {
        A_mat(r,x_start) = A_sca;
        B_mat(r,x_start) = B_sca;

        if(dir == 0) {
            double rand = randomUnit(-1.5,1.5);
            if(rand < -0.5) {
                dir = -1;
            } else if(rand > 0.5) {
                dir = 1;
            }
        } else if(dir == -1) {
            if(randomUnit(-1.0,1.0) > 0.0)
                dir = 0;
        } else if(dir == 1) {
            if(randomUnit(-1.0,1.0) < 0.0)
                dir = 0;
        } else
            std::cout << "Something is wrong!" << std::endl;

        x_start += dir;
        if(x_start < 0) x_start = A_mat.cols()-1;
        if(x_start >= A_mat.cols()) x_start = 0;

        A_mat(r,x_start) = A_sca;
        B_mat(r,x_start) = B_sca;
    }
    A_mat(A_mat.rows()-1,x_start) = A_sca;
    B_mat(B_mat.rows()-1,x_start) = B_sca;
}

/**
 * @brief Fills geometry with circles
 * @param A_mat Matrix A
 * @param B_mat Matrix B
 * @param A_sca Value inside circles in matrix A
 * @param B_sca Value inside circles in matrix B
 * @param xv Vector of x-values for circle centers
 * @param yv Vector of y-values for circle centers
 * @param rcv Vector of radius-values for circles
 * @param dx bin-width
 * @param dy bin-height
 * @param L Lateral period
 * @param offset For vertical components 'offset=0.0', while for horizontal components 'offset=1.0'
 */
void Geometry::fillGeometryWithCirclesHor(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height) {
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows()+1);
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = ( double(c) + 0.5 ) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 1.0 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = std::fabs(PBC_1D(x-xv.at(i),width));
                double yn = std::fabs(y - yv.at(i));

                if( xn*xn + yn*yn < rcv.at(i)*rcv.at(i)) {
                    A_mat(r,c) = A_sca.at(i);
                    B_mat(r,c) = B_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::fillGeometryWithCirclesVer(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height) {
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows());
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = double(c) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 0.5 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = std::fabs(PBC_1D(x-xv.at(i),width));
                double yn = std::fabs(y - yv.at(i));

                if( xn*xn + yn*yn < rcv.at(i)*rcv.at(i)) {
                    A_mat(r,c) = A_sca.at(i);
                    B_mat(r,c) = B_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::fillGeometryWithCirclesNodes(mat &A_mat, std::vector<double> A_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double width, double height) {
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows());
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = double(c) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 1.0 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = std::fabs(PBC_1D(x-xv.at(i),width));
                double yn = std::fabs(y - yv.at(i));

                if( xn*xn + yn*yn < rcv.at(i)*rcv.at(i)) {
                    A_mat(r,c) = A_sca.at(i);
                    break;
                }
            }
        }
}

/**
 * @brief Fills geometry with rectangles
 * @param A_mat Matrix A
 * @param B_mat Matrix B
 * @param A_sca Value inside circles in matrix A
 * @param B_sca Value inside circles in matrix B
 * @param xv Vector of x-values for rectangle centers
 * @param yv Vector of y-values for rectangle centers
 * @param d Rectangle width
 * @param t Rectangle height
 * @param dx bin-width
 * @param dy bin-height
 * @param L Lateral period
 * @param offset For vertical components 'offset=0.0', while for horizontal components 'offset=1.0'
 */
void Geometry::fillGeometryWithRectangles(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double dx, double dy, double L, double offset) {
    assert(xv.size() == yv.size());
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = ( double(c) + 0.5*offset ) * dx; // center of resistance
            double y = ( double(2*r) + 0.5 + 0.5*offset ) * dy; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = (PBC_1D(x-xv.at(i),L));
                double yn = (y - yv.at(i));

                if(xn >= -dv.at(i)/2.0 && xn < dv.at(i)/2.0 && yn >= -tv.at(i)/2.0 && yn < tv.at(i)/2.0) {
                    A_mat(r,c) = A_sca.at(i);
                    B_mat(r,c) = B_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::fillGeometryWithRectanglesNodes(mat &A_mat, std::vector<double> A_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height) {
    assert(xv.size() == yv.size());
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows()+1);
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = double(c) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 1.0 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = (PBC_1D(x-xv.at(i),width));
                double yn = (y - yv.at(i));

                if(xn >= -dv.at(i)/2.0 && xn < dv.at(i)/2.0 && yn >= -tv.at(i)/2.0 && yn < tv.at(i)/2.0) {
                    A_mat(r,c) = A_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::fillGeometryWithRectanglesHor(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height) {
    assert(xv.size() == yv.size());
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows()+1);
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = ( double(c) + 0.5 ) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 1.0 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = (PBC_1D(x-xv.at(i),width));
                double yn = (y - yv.at(i));

                if(xn >= -dv.at(i)/2.0 && xn < dv.at(i)/2.0 && yn >= -tv.at(i)/2.0 && yn < tv.at(i)/2.0) {
                    A_mat(r,c) = A_sca.at(i);
                    B_mat(r,c) = B_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::fillGeometryWithRectanglesVer(mat &A_mat, mat &B_mat, std::vector<double> A_sca, std::vector<double> B_sca, std::vector<double> xv, std::vector<double> yv, std::vector<double> dv, std::vector<double> tv, double width, double height) {
    assert(xv.size() == yv.size());
    double dx = width/double(A_mat.cols());
    double dy = height/double(A_mat.rows());
    for(unsigned int r = 0; r < A_mat.rows(); r++)
        for(unsigned int c = 0; c < A_mat.cols(); c++) {
            double x = double(c) * dx + dx*1e-9; // center of resistance
            double y = ( double(r) + 0.5 ) * dy + dy*1e-9; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = (PBC_1D(x-xv.at(i),width));
                double yn = (y - yv.at(i));

                if(xn >= -dv.at(i)/2.0 && xn < dv.at(i)/2.0 && yn >= -tv.at(i)/2.0 && yn < tv.at(i)/2.0) {
                    A_mat(r,c) = A_sca.at(i);
                    B_mat(r,c) = B_sca.at(i);
                    break;
                }
            }
        }
}

void Geometry::checkInputForBrickAndMortar(InputData ipd) {
    assert(("brick thickness cannot be negatve" && ipd.t >= 0));
    assert(("brick width cannot be negatve" && ipd.d >= 0));
    assert(("vertical spacing between bricks cannot be negatve" && ipd.g >= 0));
    assert(("lateral spacing between bricks cannot be negatve" && ipd.s >= 0));
    assert(("number of layers of bricks cannot be negatve" && ipd.N >= 0));
    assert(("vertical resistance in mortar cannot be negative" && ipd.S_mv >= 0));
    assert(("horizontal resistance in mortar cannot be negative" && ipd.S_mh >= 0));
    assert(("vertical resistance in bricks cannot be negative" && ipd.S_bv >= 0));
    assert(("horizontal resistance in bricks cannot be negative" && ipd.S_bh >= 0));
    assert(("vertical mobility in mortar cannot be negative" && ipd.D_mv >= 0));
    assert(("horizontal mobility in mortar cannot be negative" && ipd.D_mh >= 0));
    assert(("vertical mobility in bricks cannot be negative" && ipd.D_bv >= 0));
    assert(("horizontal mobility in bricks cannot be negative" && ipd.D_bh >= 0));
    assert(("number of columns in mesh needs to be positive" && ipd.C > 0));
    assert(("number of rows in mesh cannot be negative" && ipd.R >= 0));
}

void Geometry::getBricks(InputData ipd, std::vector<double> &xv, std::vector<double> &yv, std::vector<double> &dv, std::vector<double> &tv, std::vector<double> &Sv, std::vector<double> &Sh, std::vector<double> &Dv, std::vector<double> &Dh, int model, int K = 2) {
    xv.resize(0);
    yv.resize(0);
    dv.resize(0);
    tv.resize(0);
    Sv.resize(0);
    Sh.resize(0);
    Dv.resize(0);
    Dh.resize(0);

    if(model == 0) { // Stretcher bond (Brick and Mortar original)
        double width = ipd.d + ipd.s;
        for(int n = 0; n < ipd.N; n += 2) {
            for(int i = 0; i < 2; i++) {
                dv.push_back(ipd.d);
                tv.push_back(ipd.t);
                Sv.push_back(ipd.S_bv);
                Sh.push_back(ipd.S_bh);
                Dv.push_back(ipd.D_bv);
                Dh.push_back(ipd.D_bh);
            }
            yv.push_back(ipd.t/2.0 + n*(ipd.t+ipd.g));
            yv.push_back(ipd.t/2.0 + (n+1)*(ipd.t+ipd.g));

            if(ipd.omega < 0.0) {
                xv.push_back(randomUnit(0.0,1.0)*width);
                xv.push_back(randomUnit(0.0,1.0)*width);
            } else {
                xv.push_back(ipd.d/2.0);
                xv.push_back(ipd.d/2.0 + ipd.omega*width/(1.0+ipd.omega));
            }
        }
    } else if(model == 2) { // English Bond (Blockförband)
        double width = ipd.d + ipd.s;
        for(int n = 0; n < ipd.N; n += 2) {

            tv.push_back(ipd.t);
            Sv.push_back(ipd.S_bv);
            Sh.push_back(ipd.S_bh);
            Dv.push_back(ipd.D_bv);
            Dh.push_back(ipd.D_bh);

            tv.push_back(ipd.t);
            Sv.push_back(ipd.S_comp_bv);
            Sh.push_back(ipd.S_comp_bh);
            Dv.push_back(ipd.D_comp_bv);
            Dh.push_back(ipd.D_comp_bh);

            tv.push_back(ipd.t);
            Sv.push_back(ipd.S_comp_bv);
            Sh.push_back(ipd.S_comp_bh);
            Dv.push_back(ipd.D_comp_bv);
            Dh.push_back(ipd.D_comp_bh);

            dv.push_back(ipd.d);
            dv.push_back((ipd.d-ipd.s)/2.0);
            dv.push_back((ipd.d-ipd.s)/2.0);

            yv.push_back(ipd.t/2.0 + n*(ipd.t+ipd.g));
            yv.push_back(ipd.t/2.0 + (n+1)*(ipd.t+ipd.g));
            yv.push_back(ipd.t/2.0 + (n+1)*(ipd.t+ipd.g));

            if(ipd.omega < 0.0) {
                xv.push_back(randomUnit(0.0,1.0)*width);
                xv.push_back(randomUnit(0.0,1.0)*width);
            } else {
                xv.push_back(ipd.d/2.0);
                xv.push_back(ipd.d/2.0 + ipd.omega*width/(1.0+ipd.omega));
            }
            xv.push_back(xv.back()+dv.back()+ipd.s);
        }
    } else if(model == 3) { // K=2 => Flemish Bond (Vendiskt (eller Götiskt) förband), K=3 => Monk bond (Munkförband), K = 4 => Sussex bond
        //int K = 3; // total number of bricks in each row (including complmentary one)
        double width = ipd.d*double(K-1) + double(K)*ipd.s + ipd.d_comp;
        for(int n = 0; n < ipd.N; n += 2) {
            for(int m = 0; m < 2; m++) { // for line 'n' and line 'n+1'
                for(int k = 0; k < K-1; k++) { // bricks
                    dv.push_back(ipd.d);
                    tv.push_back(ipd.t);
                    Sv.push_back(ipd.S_bv);
                    Sh.push_back(ipd.S_bh);
                    Dv.push_back(ipd.D_bv);
                    Dh.push_back(ipd.D_bh);
                    yv.push_back(ipd.t/2.0 + (n+m)*(ipd.t+ipd.g));
                    if(ipd.omega < 0.0 && k == 0) {
                        xv.push_back(randomUnit(0.0,1.0)*width);
                    } else if(ipd.omega >= 0.0 && k == 0 && m == 0) {
                        xv.push_back(0.0);
                    } else if(ipd.omega >= 0.0 && k == 0 && m == 1) {
                        //xv.push_back(0.0 + (width - (ipd.d/2.0 + ipd.s)*double(K-1) - ipd.d_comp/2.0) * 2.0*ipd.omega/(1.0+ipd.omega));
                        xv.push_back(0.0 + (width - 2.0*ipd.d_comp - 2.0*ipd.s - 2.0*double(K-1)*(ipd.d+ipd.s))/2.0 * 2.0*ipd.omega/(1.0+ipd.omega));
                    } else {
                        xv.push_back(xv.back()+ipd.d+ipd.s);
                    }
                }
                // complementary brick
                dv.push_back(ipd.d_comp);
                tv.push_back(ipd.t);
                Sv.push_back(ipd.S_comp_bv);
                Sh.push_back(ipd.S_comp_bh);
                Dv.push_back(ipd.D_comp_bv);
                Dh.push_back(ipd.D_comp_bh);
                yv.push_back(ipd.t/2.0 + (n+m)*(ipd.t+ipd.g));
                xv.push_back(xv.back()+ipd.d/2.0+ipd.d_comp/2.0+ipd.s);
            }
        }
    }
}

void Geometry::custom_made(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    if(ipd.recX.size() != ipd.recY.size()) { std::cerr << "'recX' and 'recY' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recW.size()) { std::cerr << "'recX' and 'recW' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recH.size()) { std::cerr << "'recX' and 'recH' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recVS.size()) { std::cerr << "'recX' and 'recVS' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recVD.size()) { std::cerr << "'recX' and 'recVD' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recHS.size()) { std::cerr << "'recX' and 'recHS' is not same size\n"; exit(EXIT_FAILURE); }
    if(ipd.recX.size() != ipd.recHD.size()) { std::cerr << "'recX' and 'recHD' is not same size\n"; exit(EXIT_FAILURE); }

    name = "custom_made";
    width = ipd.width;
    height = ipd.height;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv = ipd.recX;
    std::vector<double> yv = ipd.recY;
    std::vector<double> dv = ipd.recW;
    std::vector<double> tv = ipd.recH;
    std::vector<double> VS = ipd.recVS;
    std::vector<double> VD = ipd.recVD;
    std::vector<double> HS = ipd.recHS;
    std::vector<double> HD = ipd.recHD;

    fillGeometryWithRectanglesVer(s_ver,D_ver,VS,VD,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,HS,HD,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,HS,xv,yv,dv,tv,width,height);
    //fillGeometryWithCirclesHor(s_hor,D_hor,HS,HD,xv,yv,dv,width,height);
    //fillGeometryWithCirclesVer(s_ver,D_ver,VS,VD,xv,yv,dv,width,height);
    //fillGeometryWithCirclesNodes(s_ver,VS,xv,yv,dv,width,height);

    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}





/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar system
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 */
void Geometry::brickAndMortar(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    Geometry::checkInputForBrickAndMortar(ipd);

    name = "Brick and Mortar";
    width = ipd.d + ipd.s;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh;
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,0);

    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,width,height);
    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar (English bond style) system
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 */
void Geometry::brickAndMortarEnglishBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    Geometry::checkInputForBrickAndMortar(ipd);

    name = "Brick and Mortar English Bond";
    width = ipd.d + ipd.s;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh;
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,2);

    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,width,height);
    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar (Flemish bond style) system
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 */
void Geometry::brickAndMortarFlemishBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    Geometry::checkInputForBrickAndMortar(ipd);

    name = "Brick and Mortar Flemish Bond";
    int K = 2;
    width = double(K-1)*ipd.d + double(K)*ipd.s + ipd.d_comp;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh;
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,3,K);

    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,width,height);
    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar (Monk bond style) system
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 */
void Geometry::brickAndMortarMonkBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    Geometry::checkInputForBrickAndMortar(ipd);

    name = "Brick and Mortar Monk Bond";
    int K = 3;
    width = double(K-1)*ipd.d + double(K)*ipd.s + ipd.d_comp;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh;
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,3,K);

    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,width,height);
    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar (Sussex bond style) system
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 */
void Geometry::brickAndMortarSussexBond(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    Geometry::checkInputForBrickAndMortar(ipd);

    name = "Brick and Mortar Sussex Bond";
    int K = 4;
    width = double(K-1)*ipd.d + double(K)*ipd.s + ipd.d_comp;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh;
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,3,K);

    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height);
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,width,height);
    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get solubility and diffusion coefficient matrices for a brick and mortar system where the solubility is a linear function of the depth.
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param height Height of system (will be set)
 * @param width Width of system (will be set)
 * @param ipd Input data for the system
 * @param name Name of system (will be set)
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'. The solubility 'S' is described as a function of the depth 'z' as
 * @f[
 *     S(z) = \left[S_{mt} + \left(\frac{S_{m} - S_{mt}}{z_{break}} \right)\cdot z\right]\theta(z_{break}-z) + S_m\cdot\theta(z-z_{break})
 * @f]
 * where @f$ \theta(z) @f$ is the Heaviside step function.
 */
void Geometry::brickAndMortarDualLinear(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
    assert(("brick thickness cannot be negatve" && ipd.t >= 0));
    assert(("brick width cannot be negatve" && ipd.d >= 0));
    assert(("vertical spacing between bricks cannot be negatve" && ipd.g >= 0));
    assert(("lateral spacing between bricks cannot be negatve" && ipd.s >= 0));
    assert(("number of layers of bricks cannot be negatve" && ipd.N >= 0));
    assert(("vertical resistance in mortar cannot be negative" && ipd.S_mv >= 0));
    assert(("horizontal resistance in mortar cannot be negative" && ipd.S_mh >= 0));
    assert(("vertical resistance in bricks cannot be negative" && ipd.S_bv >= 0));
    assert(("horizontal resistance in bricks cannot be negative" && ipd.S_bh >= 0));
    assert(("vertical mobility in mortar cannot be negative" && ipd.D_mv >= 0));
    assert(("horizontal mobility in mortar cannot be negative" && ipd.D_mh >= 0));
    assert(("vertical mobility in bricks cannot be negative" && ipd.D_bv >= 0));
    assert(("horizontal mobility in bricks cannot be negative" && ipd.D_bh >= 0));
    assert(("number of columns in mesh needs to be positive" && ipd.C > 0));
    assert(("number of rows in mesh cannot be negative" && ipd.R >= 0));

    name = "Brick and Mortar Gradient";
    width = ipd.d + ipd.s;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;
    //double res_width = width / double( ipd.C );
    double res_hight = height / double( 2*ipd.R + 1 );

    assert(("vertical top resistance in mortar cannot be negative" && ipd.S_mtv >= 0));
    assert(("horizontal top resistance in mortar cannot be negative" && ipd.S_mth >= 0));
    assert(("vertical top resistance in mortar cannot be negative" && ipd.D_mtv >= 0));
    assert(("horizontal top resistance in mortar cannot be negative" && ipd.D_mth >= 0));
    assert(("breaking-point between non-zero and zero resistance gradient in mortar cannot be negative" && ipd.z_break >= 0));
    assert(("breaking-point between non-zero and zero resistance gradient in mortar cannot lie outside of system" && ipd.z_break <= height));

    // fill geometry with mortar ("background")
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_nodes = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;

    double k = (ipd.S_mtv - ipd.S_mv)/(0.0 - ipd.z_break);
    double m = ipd.S_mtv;
    for(int r = 0; r < s_ver.rows(); r++) {
        double y = ( double(2*r) + 0.5 ) * res_hight; // the 'two' counts one vertical layer and one horizontal layer
        if( y < ipd.z_break)
            s_ver.row(r) = vec::Ones(s_ver.cols())*( k*y + m );
    }

    k = (ipd.S_mth - ipd.S_mh)/(0.0 - ipd.z_break);
    m = ipd.S_mth;
    for(int r = 0; r < s_hor.rows(); r++) {
        double y = ( double(2*r) + 0.5 ) * res_hight; // the 'two' counts one vertical layer and one horizontal layer
        if( y < ipd.z_break) {
            s_hor.row(r) = vec::Ones(s_hor.cols())*( k*y + m );
            s_nodes.row(r) = vec::Ones(s_nodes.cols())*( k*y + m );
        }
    }

    k = (ipd.D_mtv - ipd.D_mv)/(0.0 - ipd.z_break);
    m = ipd.D_mtv;
    for(int r = 0; r < D_ver.rows(); r++) {
        double y = ( double(2*r) + 0.5 ) * res_hight; // the 'two' counts one vertical layer and one horizontal layer
        if( y < ipd.z_break)
            D_ver.row(r) = vec::Ones(D_ver.cols())*( k*y + m );
    }

    k = (ipd.D_mth - ipd.D_mh)/(0.0 - ipd.z_break);
    m = ipd.D_mth;
    for(int r = 0; r < D_hor.rows(); r++) {
        double y = ( double(2*r) + 0.5 ) * res_hight; // the 'two' counts one vertical layer and one horizontal layer
        if( y < ipd.z_break)
            D_hor.row(r) = vec::Ones(D_hor.cols())*( k*y + m );
    }

    // calculate brick centers
    /* OLD
       std::vector<double> xv(0);
       std::vector<double> yv(0);
       std::vector<double> dv(0);
       std::vector<double> tv(0);
       std::vector<double> Sv(0);
       std::vector<double> Sh(0);
       std::vector<double> Dv(0);
       std::vector<double> Dh(0);
       for(int n = 0; n < ipd.N; n += 2) {
       Sv.push_back(ipd.S_bv);
       Sh.push_back(ipd.S_bh);
       Dv.push_back(ipd.D_bv);
       Dh.push_back(ipd.D_bh);
       dv.push_back(ipd.d);
       tv.push_back(ipd.t);
       yv.push_back(ipd.t/2.0 + n*(ipd.t+ipd.g));
       yv.push_back(ipd.t/2.0 + (n+1)*(ipd.t+ipd.g));

       if(ipd.omega < 0.0) {
       xv.push_back(randomUnit(0.0,1.0)*width);
       xv.push_back(randomUnit(0.0,1.0)*width);
       } else {
       xv.push_back(ipd.d/2.0);
       xv.push_back(ipd.d/2.0 + ipd.omega*(ipd.d+ipd.s)/(1.0+ipd.omega));
       }
       }
     */ // OLD

    // calculate brick centers
    std::vector<double> xv, yv, dv, tv, Sv, Sh, Dv, Dh; // UNTESTED (HERE)
    Geometry::getBricks(ipd,xv,yv,dv,tv,Sv,Sh,Dv,Dh,0); // UNTESTED (HERE)

    // fill geometry with bricks ("foreground")
    //fillGeometryWithRectangles(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,res_width,res_hight,width,0.0); // OLD
    //fillGeometryWithRectangles(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,res_width,res_hight,width,1.0); // OLD
    fillGeometryWithRectanglesVer(s_ver,D_ver,Sv,Dv,xv,yv,dv,tv,width,height); // UNTESTED (HERE)
    fillGeometryWithRectanglesHor(s_hor,D_hor,Sh,Dh,xv,yv,dv,tv,width,height); // UNTESTED (HERE)
    fillGeometryWithRectanglesNodes(s_nodes,Sh,xv,yv,dv,tv,ipd.width,ipd.height);

    insertRow(s_nodes,0,ipd.S_out);
    s_nodes.conservativeResize(s_nodes.rows()+1, s_nodes.cols());
    s_nodes.row(s_nodes.rows()-1) = ipd.S_out*vec::Ones(s_nodes.row(s_nodes.rows()-1).cols());
}

/**
 * @brief Get geometry and properties of membrane
 * @param s_hor Matrix of horizontal solubilities (will be set)
 * @param s_ver Matrix of vertical solubilities (will be set)
 * @param D_hor Matrix of horizontal diffusion coefficients (will be set)
 * @param D_ver Matrix of vertical diffusion coefficients (will be set)
 * @param ipd Input data for the system
 * @note All input except 'ipd' is set using this function. All input data is supplied in 'ipd'.
 * @todo Check 's_nodes' size for load_external.
 */
void Geometry::getGeometry(mat &s_hor, mat &s_ver, mat &s_nodes, mat &D_hor, mat &D_ver, InputData &ipd) {
    if(ipd.load_external) {
        s_ver = loadMatrix(ipd.input_folder+"sv_matrix.txt");
        s_hor = loadMatrix(ipd.input_folder+"sh_matrix.txt");
        s_nodes = loadMatrix(ipd.input_folder+"sn_matrix.txt");
        D_ver = loadMatrix(ipd.input_folder+"Dv_matrix.txt");
        D_hor = loadMatrix(ipd.input_folder+"Dh_matrix.txt");
        assert(("externally loaded matrices do not have compatible sizes" && compatibleSizes(s_ver,s_hor,D_ver,D_hor)));
        assert(("height of system must be positive" && ipd.height > 0.0));
        assert(("width of system must be positive" && ipd.width > 0.0));
        appendDataToFile(ipd.output_file,"model external\n");
    } else {
        std::string name = "";
        if(ipd.model_nbr == 6) {
            custom_made(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
        } else if(ipd.model_nbr == 5) {
            brickAndMortarSussexBond(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else if(ipd.model_nbr == 4) {
            brickAndMortarMonkBond(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else if(ipd.model_nbr == 3) {
            brickAndMortarFlemishBond(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else if(ipd.model_nbr == 2) {
            brickAndMortarEnglishBond(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else if(ipd.model_nbr == 1) {
            brickAndMortarDualLinear(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else if(ipd.model_nbr == 0) {
            brickAndMortar(s_hor,s_ver,s_nodes,D_hor,D_ver,ipd.height,ipd.width,ipd,name);
            appendDataToFile(ipd.output_file,"K_{M_ver/out} "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{B_ver/out} "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_ver/B_ver} "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
            appendDataToFile(ipd.output_file,"K_{M_hor/B_hor} "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        } else {
            std::cerr << "model not found\n";
            exit(EXIT_FAILURE);
        }
        appendDataToFile(ipd.output_file,"model "+name+"\n");
        writeMatrixToFile(ipd.output_folder+"sh_matrix.txt",s_hor);
        writeMatrixToFile(ipd.output_folder+"sv_matrix.txt",s_ver);
        writeMatrixToFile(ipd.output_folder+"sn_matrix.txt",s_nodes);
        writeMatrixToFile(ipd.output_folder+"Dh_matrix.txt",D_hor);
        writeMatrixToFile(ipd.output_folder+"Dv_matrix.txt",D_ver);
    }
    appendDataToFile(ipd.output_file,"height "+to_string_precision(ipd.height) +"\n");
    appendDataToFile(ipd.output_file,"width "+to_string_precision(ipd.width) +"\n");
}
