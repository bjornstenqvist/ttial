#pragma once
#include "auxiliary.h"

/**
 * @brief Get ID of left node
 */
int getLeftID(const int c, const int r, const int C) {
    int ID = (r-1)*C + c - 1;
    if(c == 1) ID += C;
    return ID;
}

/**
 * @brief Get ID of right node
 */
int getRightID(const int c, const int r, const int C) {
    int ID = (r-1)*C + c + 1;
    if(c == C) ID -= C;
    return ID;
}

/**
 * @brief Get ID of upper node
 */
int getUpperID(const int c, const int r, const int C) {
    if(r == 1) return 0;
    return ( (r-1)*C + c - C );
}

/**
 * @brief Get ID of lower node
 */
int getLowerID(const int c, const int r, const int C, const int R) {
    if(r == R) return -1;
    return ( (r-1)*C + c + C );
}

/**
 * @brief Get resistance left of node
 */
double getLeftRes(const mat res_hor, const int c, const int r,const int C) {
    if(c == 1) return res_hor(r-1,C-1);
    return res_hor(r-1,c-2);
}

struct node {
    int ID, leftID, rightID, upperID, lowerID;
    double leftRes, rightRes, upperRes, lowerRes;

    /**
     * @brief Initialize node with ID's of surrounding nodes and resistances
     */
    void initNode(const mat res_hor, const mat res_ver, const int c, const int r, const int C, const int R) {
        ID = (r-1)*C + c;
        leftID = getLeftID(c,r,C);
        rightID = getRightID(c,r,C);
        upperID = getUpperID(c,r,C);
        lowerID = getLowerID(c,r,C,R);
        leftRes = getLeftRes(res_hor,c,r,C);
        rightRes = res_hor(r-1,c-1);
        upperRes = res_ver(r-1,c-1);
        lowerRes = res_ver(r,c-1);
    }
};

/**
 * @brief Constructs matrices such that mass (or rather mass flow) is conserved in each node. 
 */
void fillMatricesFromNodes(const std::vector<node> nodes, const int C, const int R, mat &Rinv, vec &Iex, const double Vtop) {
    // First row
    for(int k = 1; k <= C; k++) {
        node node_tmp = nodes.at(k-1);
        int ID = node_tmp.ID;
        Rinv(ID-1,ID-1) += 1/node_tmp.leftRes + 1.0/node_tmp.rightRes + 1.0/node_tmp.upperRes + 1.0/node_tmp.lowerRes;
        Rinv(ID-1,node_tmp.leftID-1) -= 1.0/node_tmp.leftRes;
        Rinv(ID-1,node_tmp.rightID-1) -= 1.0/node_tmp.rightRes;
        if(R > 2)
            Rinv(ID-1,node_tmp.lowerID-1) -= 1.0/node_tmp.lowerRes;
        Iex(ID-1) = Iex(ID-1) + Vtop/node_tmp.upperRes;
    }
    // Middle rows
    for(int k = C+1; k <= ((R-1)*C); k++) {
        node node_tmp = nodes.at(k-1);
        int ID = node_tmp.ID;
        Rinv(ID-1,ID-1) += 1/node_tmp.leftRes + 1.0/node_tmp.rightRes + 1.0/node_tmp.upperRes + 1.0/node_tmp.lowerRes;
        Rinv(ID-1,node_tmp.leftID-1) -= 1.0/node_tmp.leftRes;
        Rinv(ID-1,node_tmp.rightID-1) -= 1.0/node_tmp.rightRes;
        Rinv(ID-1,node_tmp.upperID-1) -= 1.0/node_tmp.upperRes;
        Rinv(ID-1,node_tmp.lowerID-1) -= 1.0/node_tmp.lowerRes;
    }
    // Last row
    if(R > 1) {
        for( int k = (R-1)*C+1; k <= (R*C); k++) {
            node node_tmp = nodes.at(k-1);
            int ID = node_tmp.ID;
            Rinv(ID-1,ID-1) += 1.0/node_tmp.leftRes + 1.0/node_tmp.rightRes + 1.0/node_tmp.upperRes + 1.0/node_tmp.lowerRes;
            Rinv(ID-1,node_tmp.leftID-1) -= 1.0/node_tmp.leftRes;
            Rinv(ID-1,node_tmp.rightID-1) -= 1.0/node_tmp.rightRes;
            if(R > 2)
                Rinv(ID-1,node_tmp.upperID-1) -= 1.0/node_tmp.upperRes;
            //Iex(ID-1) = Iex(ID-1) + 0.0/node_tmp.lowerRes;
        }
    }
}
