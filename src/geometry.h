#pragma once
#include "auxiliary.h"
#include "nodes.h"

void fillGeometryWithTunnels(mat &Am_v, mat &Am_h, mat &Bm_v, mat &Bm_h, double As, double Bs, int xc) {
    Am_v(0,xc) = As;
    Bm_v(0,xc) = Bs;
    Am_h(0,xc) = As;
    Bm_h(0,xc) = Bs;
    int dir = 0;

    for(int r = 1; r < Am_h.rows(); r++) {
        Am_v(r,xc) = As;
        Bm_v(r,xc) = Bs;
        Am_h(r,xc) = As;
        Bm_h(r,xc) = Bs;

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

        xc += dir;
        if(xc < 0) xc = Am_v.cols()-1;
        if(xc >= Am_v.cols()) xc = 0;

        Am_v(r,xc) = As;
        Bm_v(r,xc) = Bs;
        Am_h(r,xc) = As;
        Bm_h(r,xc) = Bs;
    }
    Am_v(Am_v.rows()-1,xc) = As;
    Bm_v(Am_v.rows()-1,xc) = Bs;
}

void fillGeometryWithCircles(mat &Am, mat &Bm, double As, double Bs, std::vector<double> xv, std::vector<double> yv, std::vector<double> rcv, double dx, double dy, double L, double offset) {
    for(unsigned int r = 0; r < Am.rows(); r++)
        for(unsigned int c = 0; c < Am.cols(); c++) {
            double x = ( double(c) + 0.5 ) * dx; // center of resistance
            double y = ( double(2*r) + 0.5 + offset ) * dy; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = std::fabs(PBC_1D(x-xv.at(i),L));
                double yn = std::fabs(y - yv.at(i));

                if( xn*xn + yn*yn < rcv.at(i)*rcv.at(i)) {
                    Am(r,c) = As;
                    Bm(r,c) = Bs;
                    break;
                }
            }
        }
}

void fillGeometryWithRectangles(mat &Am, mat &Bm, double As, double Bs, std::vector<double> xv, std::vector<double> yv, double d, double t, double dx, double dy, double L, double offset) {
    assert(xv.size() == yv.size());
    for(unsigned int r = 0; r < Am.rows(); r++)
        for(unsigned int c = 0; c < Am.cols(); c++) {
            double x = ( double(c) + 0.5 ) * dx; // center of resistance
            double y = ( double(2*r) + 0.5 + offset ) * dy; // the 'two' counts one vertical layer and one horizontal layer

            for(unsigned int i = 0; i < xv.size(); i++) {
                double xn = std::fabs(PBC_1D(x-xv.at(i),L));
                double yn = std::fabs(y - yv.at(i));

                if(xn < d/2.0 && yn < t/2.0) {
                    Am(r,c) = As;
                    Bm(r,c) = Bs;
                    break;
                }
            }
        }
}

void brickAndMortar(mat &s_hor, mat &s_ver, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
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

    name = "Brick and Mortar";
    width = ipd.d + ipd.s;
    height = double( ipd.N )*ipd.t + double( ipd.N-1 )*ipd.g;
    double res_width = width / double( ipd.C );
    double res_hight = height / double( 2*ipd.R + 1 );

    // fill geometry with mortar ("background")
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    D_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.D_mv;
    D_hor = mat::Ones(ipd.R,ipd.C)*ipd.D_mh;

    // calculate brick centers
    std::vector<double> xv(0);
    std::vector<double> yv(0);
    std::vector<double> rcv(0);
    for(int n = 0; n < ipd.N; n += 2) {
        yv.push_back(ipd.t/2.0 + n*(ipd.t+ipd.g));
        yv.push_back(ipd.t/2.0 + (n+1)*(ipd.t+ipd.g));
        rcv.push_back(ipd.t/2.0);
        rcv.push_back(ipd.t/2.0);

        if(ipd.omega < 0.0) {
            xv.push_back(randomUnit(0.0,1.0)*width);
            xv.push_back(randomUnit(0.0,1.0)*width);
        } else {
            xv.push_back(ipd.d/2.0);
            xv.push_back(ipd.d/2.0 + ipd.omega*(ipd.d+ipd.s)/(1.0+ipd.omega));
        }
    }

    double y = ipd.g/2.0;
    std::vector<double> xlv(0);
    std::vector<double> ylv(0);
    while(y < height) {
        xlv.push_back(0.0);
        ylv.push_back(y);
        y += ipd.g;
    }

    fillGeometryWithRectangles(s_ver,D_ver,ipd.S_bv,ipd.D_bv,xv,yv,ipd.d,ipd.t,res_width,res_hight,width,0.0);
    fillGeometryWithRectangles(s_hor,D_hor,ipd.S_bh,ipd.D_bh,xv,yv,ipd.d,ipd.t,res_width,res_hight,width,1.0);
}

void brickAndMortarDualLinear(mat &s_hor, mat &s_ver, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string &name) {
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
    double res_width = width / double( ipd.C );
    double res_hight = height / double( 2*ipd.R + 1 );

    assert(("vertical top resistance in mortar cannot be negative" && ipd.S_mtv >= 0));
    assert(("horizontal top resistance in mortar cannot be negative" && ipd.S_mth >= 0));
    assert(("breaking-point between non-zero and zero resistance gradient in mortar cannot be negative" && ipd.z_break >= 0));
    assert(("breaking-point between non-zero and zero resistance gradient in mortar cannot lie outside of system" && ipd.z_break <= height));

    // fill geometry with mortar ("background")
    s_hor = mat::Ones(ipd.R,ipd.C)*ipd.S_mh;
    s_ver = mat::Ones(ipd.R+1,ipd.C)*ipd.S_mv;
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
        if( y < ipd.z_break)
            s_hor.row(r) = vec::Ones(s_hor.cols())*( k*y + m );
    }

    // calculate brick centers
    std::vector<double> xv(0);
    std::vector<double> yv(0);
    for(int n = 0; n < ipd.N; n += 2) {
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

    // fill geometry with bricks ("foreground")
    fillGeometryWithRectangles(s_ver,D_ver,ipd.S_bv,ipd.D_bv,xv,yv,ipd.d,ipd.t,res_width,res_hight,width,0.0);
    fillGeometryWithRectangles(s_hor,D_hor,ipd.S_bh,ipd.D_bh,xv,yv,ipd.d,ipd.t,res_width,res_hight,width,1.0);
}

void getGeometry(mat &s_hor, mat &s_ver, mat &D_hor, mat &D_ver, double &height, double &width, InputData ipd, std::string output_file) {
    if(ipd.load_external) {
        s_ver = loadMatrix("sv_matrix.txt");
        s_hor = loadMatrix("sh_matrix.txt");
        D_ver = loadMatrix("Dv_matrix.txt");
        D_hor = loadMatrix("Dh_matrix.txt");
        assert(("externally loaded matrices do not have compatible sizes" && compatibleSizes(s_ver,s_hor,D_ver,D_hor)));
        assert(("height of system must be positive" && ipd.height > 0.0));
        assert(("width of system must be positive" && ipd.width > 0.0));
        appendDataToFile(output_file,"model  external\n");
    } else {
        std::string name = "";
        if(ipd.model_nbr == 1) {
            brickAndMortarDualLinear(s_hor,s_ver,D_hor,D_ver,height,width,ipd,name);
        } else if(ipd.model_nbr == 0) {
            brickAndMortar(s_hor,s_ver,D_hor,D_ver,height,width,ipd,name);
        } else {
            std::cerr << "model not found\n";
            exit(EXIT_FAILURE);
        }
        appendDataToFile(output_file,"model  "+name+"\n");
        appendDataToFile(output_file,"K_{M_ver/out}  "+to_string_precision(ipd.S_mv/ipd.S_out) +"\n");
        appendDataToFile(output_file,"K_{B_ver/out}  "+to_string_precision(ipd.S_bv/ipd.S_out) +"\n");
        appendDataToFile(output_file,"K_{M_ver/B_ver}  "+to_string_precision(ipd.S_mv/ipd.S_bv) +"\n");
        appendDataToFile(output_file,"K_{M_hor/B_hor}  "+to_string_precision(ipd.S_mh/ipd.S_bh) +"\n");
        writeMatrixToFile("sh_matrix.txt",s_hor);
        writeMatrixToFile("sv_matrix.txt",s_ver);
        writeMatrixToFile("Dh_matrix.txt",D_hor);
        writeMatrixToFile("Dv_matrix.txt",D_ver);
    }
}