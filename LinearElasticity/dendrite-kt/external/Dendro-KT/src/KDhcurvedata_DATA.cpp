/**
 * @author Masado Ishii
 * @brief The actual table data generated by another script.
 */
#include "../include/KDhcurvedata_decl.h"

template <>
const char HilbertData<2>::m_rotations[] = {
    0,1,3,2, 0,1,3,2,
    0,2,3,1, 0,3,1,2,
    3,2,0,1, 2,3,1,0,
    3,1,0,2, 2,1,3,0,
};

template <>
const int HilbertData<2>::m_HILBERT_TABLE[] = {
    1,0,3,0,
    0,2,1,1,
    2,1,2,3,
    3,3,0,2,
};

template <>
const char HilbertData<3>::m_rotations[] = {
    0,1,3,2,6,7,5,4, 0,1,3,2,7,6,4,5,
    0,4,6,2,3,7,5,1, 0,7,3,4,1,6,2,5,
    0,4,5,1,3,7,6,2, 0,3,7,4,1,2,6,5,
    0,2,3,1,5,7,6,4, 0,3,1,2,7,4,6,5,
    0,2,6,4,5,7,3,1, 0,7,1,6,3,4,2,5,
    0,1,5,4,6,7,3,2, 0,1,7,6,3,2,4,5,
    5,1,0,4,6,2,3,7, 2,1,5,6,3,0,4,7,
    5,7,6,4,0,2,3,1, 4,7,5,6,3,0,2,1,
    5,7,3,1,0,2,6,4, 4,3,5,2,7,0,6,1,
    5,4,0,1,3,2,6,7, 2,3,5,4,1,0,6,7,
    5,4,6,7,3,2,0,1, 6,7,5,4,1,0,2,3,
    5,1,3,7,6,2,0,4, 6,1,5,2,7,0,4,3,
    3,1,5,7,6,4,0,2, 6,1,7,0,5,2,4,3,
    3,2,6,7,5,4,0,1, 6,7,1,0,5,4,2,3,
    3,2,0,1,5,4,6,7, 2,3,1,0,5,4,6,7,
    3,7,5,1,0,4,6,2, 4,3,7,0,5,2,6,1,
    3,7,6,2,0,4,5,1, 4,7,3,0,5,6,2,1,
    3,1,0,2,6,4,5,7, 2,1,3,0,5,6,4,7,
    6,2,0,4,5,1,3,7, 2,5,1,6,3,4,0,7,
    6,7,5,4,0,1,3,2, 4,5,7,6,3,2,0,1,
    6,7,3,2,0,1,5,4, 4,5,3,2,7,6,0,1,
    6,4,0,2,3,1,5,7, 2,5,3,4,1,6,0,7,
    6,4,5,7,3,1,0,2, 6,5,7,4,1,2,0,3,
    6,2,3,7,5,1,0,4, 6,5,1,2,7,4,0,3,
};

template <>
const int HilbertData<3>::m_HILBERT_TABLE[] = {
    1,5,17,2,11,20,17,23,
    0,10,21,21,2,16,5,13,
    3,9,22,9,1,4,15,12,
    2,14,4,1,23,14,8,11,
    5,13,3,7,18,18,0,10,
    4,0,12,19,6,3,6,22,
    8,11,21,18,5,7,5,17,
    19,16,1,4,19,6,11,8,
    15,15,0,10,20,9,3,7,
    7,2,17,2,10,8,14,21,
    13,1,16,22,9,11,6,22,
    23,6,20,9,0,10,12,12,
    22,17,5,13,19,14,11,11,
    10,4,14,12,7,23,17,23,
    16,3,13,15,6,3,9,18,
    8,8,19,14,5,13,2,16,
    20,7,20,17,4,1,12,15,
    15,12,0,16,18,21,0,6,
    20,9,23,6,4,4,19,14,
    7,2,15,5,7,23,18,20,
    16,3,16,22,8,0,21,19,
    19,14,1,1,22,17,20,9,
    12,15,2,10,21,18,23,10,
    11,8,18,21,3,13,22,13,
};

template <>
const char HilbertData<4>::m_rotations[] = {
    0,1,3,2,6,7,5,4,12,13,15,14,10,11,9,8, 0,1,3,2,7,6,4,5,15,14,12,13,8,9,11,10,
    0,8,12,4,6,14,10,2,3,11,15,7,5,13,9,1, 0,15,7,8,3,12,4,11,1,14,6,9,2,13,5,10,
    0,8,9,1,3,11,10,2,6,14,15,7,5,13,12,4, 0,3,7,4,15,12,8,11,1,2,6,5,14,13,9,10,
    0,4,6,2,3,7,5,1,9,13,15,11,10,14,12,8, 0,7,3,4,1,6,2,5,15,8,12,11,14,9,13,10,
    0,4,12,8,9,13,5,1,3,7,15,11,10,14,6,2, 0,7,15,8,1,6,14,9,3,4,12,11,2,5,13,10,
    0,2,3,1,9,11,10,8,12,14,15,13,5,7,6,4, 0,3,1,2,15,12,14,13,7,4,6,5,8,11,9,10,
    0,2,6,4,12,14,10,8,9,11,15,13,5,7,3,1, 0,15,1,14,3,12,2,13,7,8,6,9,4,11,5,10,
    0,1,9,8,12,13,5,4,6,7,15,14,10,11,3,2, 0,1,15,14,7,6,8,9,3,2,12,13,4,5,11,10,
    0,2,6,4,5,7,3,1,9,11,15,13,12,14,10,8, 0,7,1,6,3,4,2,5,15,8,14,9,12,11,13,10,
    0,8,9,1,5,13,12,4,6,14,15,7,3,11,10,2, 0,3,15,12,7,4,8,11,1,2,14,13,6,5,9,10,
    0,8,10,2,6,14,12,4,5,13,15,7,3,11,9,1, 0,15,3,12,7,8,4,11,1,14,2,13,6,9,5,10,
    0,1,5,4,6,7,3,2,10,11,15,14,12,13,9,8, 0,1,7,6,3,2,4,5,15,14,8,9,12,13,11,10,
    0,1,9,8,10,11,3,2,6,7,15,14,12,13,5,4, 0,1,7,6,15,14,8,9,3,2,4,5,12,13,11,10,
    0,4,6,2,10,14,12,8,9,13,15,11,3,7,5,1, 0,15,3,12,1,14,2,13,7,8,4,11,6,9,5,10,
    0,4,5,1,9,13,12,8,10,14,15,11,3,7,6,2, 0,3,15,12,1,2,14,13,7,4,8,11,6,5,9,10,
    0,2,10,8,9,11,3,1,5,7,15,13,12,14,6,4, 0,7,1,6,15,8,14,9,3,4,2,5,12,11,13,10,
    0,4,5,1,3,7,6,2,10,14,15,11,9,13,12,8, 0,3,7,4,1,2,6,5,15,12,8,11,14,13,9,10,
    0,8,10,2,3,11,9,1,5,13,15,7,6,14,12,4, 0,7,3,4,15,8,12,11,1,6,2,5,14,9,13,10,
    0,8,12,4,5,13,9,1,3,11,15,7,6,14,10,2, 0,7,15,8,3,4,12,11,1,6,14,9,2,5,13,10,
    0,2,3,1,5,7,6,4,12,14,15,13,9,11,10,8, 0,3,1,2,7,4,6,5,15,12,14,13,8,11,9,10,
    0,2,10,8,12,14,6,4,5,7,15,13,9,11,3,1, 0,15,1,14,7,8,6,9,3,12,2,13,4,11,5,10,
    0,1,5,4,12,13,9,8,10,11,15,14,6,7,3,2, 0,1,15,14,3,2,12,13,7,6,8,9,4,5,11,10,
    0,1,3,2,10,11,9,8,12,13,15,14,6,7,5,4, 0,1,3,2,15,14,12,13,7,6,4,5,8,9,11,10,
    0,4,12,8,10,14,6,2,3,7,15,11,9,13,5,1, 0,15,7,8,1,14,6,9,3,12,4,11,2,13,5,10,
    12,4,0,8,9,1,5,13,15,7,3,11,10,2,6,14, 2,5,13,10,1,6,14,9,3,4,12,11,0,7,15,8,
    12,14,15,13,9,11,10,8,0,2,3,1,5,7,6,4, 8,11,9,10,15,12,14,13,7,4,6,5,0,3,1,2,
    12,14,6,4,0,2,10,8,9,11,3,1,5,7,15,13, 4,11,5,10,3,12,2,13,7,8,6,9,0,15,1,14,
    12,13,9,8,0,1,5,4,6,7,3,2,10,11,15,14, 4,5,11,10,7,6,8,9,3,2,12,13,0,1,15,14,
    12,13,15,14,6,7,5,4,0,1,3,2,10,11,9,8, 8,9,11,10,7,6,4,5,15,14,12,13,0,1,3,2,
    12,8,0,4,6,2,10,14,15,11,3,7,5,1,9,13, 2,13,5,10,3,12,4,11,1,14,6,9,0,15,7,8,
    12,8,9,13,15,11,10,14,6,2,3,7,5,1,0,4, 14,13,9,10,15,12,8,11,1,2,6,5,0,3,7,4,
    12,4,6,14,15,7,5,13,9,1,3,11,10,2,0,8, 14,9,13,10,1,6,2,5,15,8,12,11,0,7,3,4,
    12,8,9,13,5,1,0,4,6,2,3,7,15,11,10,14, 6,5,9,10,7,4,8,11,1,2,14,13,0,3,15,12,
    12,14,6,4,5,7,15,13,9,11,3,1,0,2,10,8, 12,11,13,10,3,4,2,5,15,8,14,9,0,7,1,6,
    12,14,10,8,9,11,15,13,5,7,3,1,0,2,6,4, 12,11,13,10,15,8,14,9,3,4,2,5,0,7,1,6,
    12,4,5,13,9,1,0,8,10,2,3,11,15,7,6,14, 6,5,9,10,1,2,14,13,7,4,8,11,0,3,15,12,
    12,4,6,14,10,2,0,8,9,1,3,11,15,7,5,13, 6,9,5,10,1,14,2,13,7,8,4,11,0,15,3,12,
    12,13,9,8,10,11,15,14,6,7,3,2,0,1,5,4, 12,13,11,10,15,14,8,9,3,2,4,5,0,1,7,6,
    12,13,5,4,6,7,15,14,10,11,3,2,0,1,9,8, 12,13,11,10,3,2,4,5,15,14,8,9,0,1,7,6,
    12,8,10,14,6,2,0,4,5,1,3,7,15,11,9,13, 6,9,5,10,7,8,4,11,1,14,2,13,0,15,3,12,
    12,13,5,4,0,1,9,8,10,11,3,2,6,7,15,14, 4,5,11,10,3,2,12,13,7,6,8,9,0,1,15,14,
    12,14,10,8,0,2,6,4,5,7,3,1,9,11,15,13, 4,11,5,10,7,8,6,9,3,12,2,13,0,15,1,14,
    12,14,15,13,5,7,6,4,0,2,3,1,9,11,10,8, 8,11,9,10,7,4,6,5,15,12,14,13,0,3,1,2,
    12,8,0,4,5,1,9,13,15,11,3,7,6,2,10,14, 2,5,13,10,3,4,12,11,1,6,14,9,0,7,15,8,
    12,8,10,14,15,11,9,13,5,1,3,7,6,2,0,4, 14,9,13,10,15,8,12,11,1,6,2,5,0,7,3,4,
    12,4,5,13,15,7,6,14,10,2,3,11,9,1,0,8, 14,13,9,10,1,2,6,5,15,12,8,11,0,3,7,4,
    12,4,0,8,10,2,6,14,15,7,3,11,9,1,5,13, 2,13,5,10,1,14,6,9,3,12,4,11,0,15,7,8,
    12,13,15,14,10,11,9,8,0,1,3,2,6,7,5,4, 8,9,11,10,15,14,12,13,7,6,4,5,0,1,3,2,
    15,13,12,14,6,4,5,7,3,1,0,2,10,8,9,11, 10,9,11,8,5,6,4,7,13,14,12,15,2,1,3,0,
    15,11,3,7,6,2,10,14,12,8,0,4,5,1,9,13, 10,13,5,2,11,12,4,3,9,14,6,1,8,15,7,0,
    15,11,9,13,12,8,10,14,6,2,0,4,5,1,3,7, 10,13,9,14,11,12,8,15,5,2,6,1,4,3,7,0,
    15,7,6,14,12,4,5,13,9,1,0,8,10,2,3,11, 10,9,13,14,5,6,2,1,11,8,12,15,4,7,3,0,
    15,7,3,11,9,1,5,13,12,4,0,8,10,2,6,14, 10,5,13,2,9,6,14,1,11,4,12,3,8,7,15,0,
    15,14,12,13,9,8,10,11,3,2,0,1,5,4,6,7, 10,11,9,8,13,12,14,15,5,4,6,7,2,3,1,0,
    15,14,6,7,3,2,10,11,9,8,0,1,5,4,12,13, 10,11,5,4,13,12,2,3,9,8,6,7,14,15,1,0,
    15,13,9,11,3,1,5,7,6,4,0,2,10,8,12,14, 10,5,11,4,9,6,8,7,13,2,12,3,14,1,15,0,
    15,14,6,7,5,4,12,13,9,8,0,1,3,2,10,11, 10,11,13,12,5,4,2,3,9,8,14,15,6,7,1,0,
    15,11,9,13,5,1,3,7,6,2,0,4,12,8,10,14, 10,5,9,6,11,4,8,7,13,2,14,1,12,3,15,0,
    15,11,10,14,6,2,3,7,5,1,0,4,12,8,9,13, 10,9,5,6,11,8,4,7,13,14,2,1,12,15,3,0,
    15,13,5,7,6,4,12,14,10,8,0,2,3,1,9,11, 10,13,11,12,5,2,4,3,9,14,8,15,6,1,7,0,
    15,13,9,11,10,8,12,14,6,4,0,2,3,1,5,7, 10,13,11,12,9,14,8,15,5,2,4,3,6,1,7,0,
    15,7,6,14,10,2,3,11,9,1,0,8,12,4,5,13, 10,9,5,6,13,14,2,1,11,8,4,7,12,15,3,0,
    15,7,5,13,9,1,3,11,10,2,0,8,12,4,6,14, 10,5,9,6,13,2,14,1,11,4,8,7,12,3,15,0,
    15,14,10,11,9,8,12,13,5,4,0,1,3,2,6,7, 10,11,13,12,9,8,14,15,5,4,2,3,6,7,1,0,
    15,7,5,13,12,4,6,14,10,2,0,8,9,1,3,11, 10,13,9,14,5,2,6,1,11,12,8,15,4,3,7,0,
    15,11,10,14,12,8,9,13,5,1,0,4,6,2,3,7, 10,9,13,14,11,8,12,15,5,6,2,1,4,7,3,0,
    15,11,3,7,5,1,9,13,12,8,0,4,6,2,10,14, 10,5,13,2,11,4,12,3,9,6,14,1,8,7,15,0,
    15,14,12,13,5,4,6,7,3,2,0,1,9,8,10,11, 10,11,9,8,5,4,6,7,13,12,14,15,2,3,1,0,
    15,14,10,11,3,2,6,7,5,4,0,1,9,8,12,13, 10,11,5,4,9,8,6,7,13,12,2,3,14,15,1,0,
    15,13,5,7,3,1,9,11,10,8,0,2,6,4,12,14, 10,5,11,4,13,2,12,3,9,6,8,7,14,1,15,0,
    15,13,12,14,10,8,9,11,3,1,0,2,6,4,5,7, 10,9,11,8,13,14,12,15,5,6,4,7,2,1,3,0,
    15,7,3,11,10,2,6,14,12,4,0,8,9,1,5,13, 10,13,5,2,9,14,6,1,11,12,4,3,8,15,7,0,
    3,7,15,11,9,13,5,1,0,4,12,8,10,14,6,2, 8,7,15,0,9,6,14,1,11,4,12,3,10,5,13,2,
    3,2,0,1,9,8,10,11,15,14,12,13,5,4,6,7, 2,3,1,0,13,12,14,15,5,4,6,7,10,11,9,8,
    3,2,6,7,15,14,10,11,9,8,12,13,5,4,0,1, 14,15,1,0,13,12,2,3,9,8,6,7,10,11,5,4,
    3,1,9,11,15,13,5,7,6,4,12,14,10,8,0,2, 14,1,15,0,9,6,8,7,13,2,12,3,10,5,11,4,
    3,1,0,2,6,4,5,7,15,13,12,14,10,8,9,11, 2,1,3,0,5,6,4,7,13,14,12,15,10,9,11,8,
    3,11,15,7,6,14,10,2,0,8,12,4,5,13,9,1, 8,15,7,0,11,12,4,3,9,14,6,1,10,13,5,2,
    3,11,9,1,0,8,10,2,6,14,12,4,5,13,15,7, 4,3,7,0,11,12,8,15,5,2,6,1,10,13,9,14,
    3,7,6,2,0,4,5,1,9,13,12,8,10,14,15,11, 4,7,3,0,5,6,2,1,11,8,12,15,10,9,13,14,
    3,11,9,1,5,13,15,7,6,14,12,4,0,8,10,2, 12,3,15,0,11,4,8,7,13,2,14,1,10,5,9,6,
    3,2,6,7,5,4,0,1,9,8,12,13,15,14,10,11, 6,7,1,0,5,4,2,3,9,8,14,15,10,11,13,12,
    3,2,10,11,9,8,0,1,5,4,12,13,15,14,6,7, 6,7,1,0,9,8,14,15,5,4,2,3,10,11,13,12,
    3,7,5,1,9,13,15,11,10,14,12,8,0,4,6,2, 12,3,15,0,13,2,14,1,11,4,8,7,10,5,9,6,
    3,7,6,2,10,14,15,11,9,13,12,8,0,4,5,1, 12,15,3,0,13,14,2,1,11,8,4,7,10,9,5,6,
    3,1,9,11,10,8,0,2,6,4,12,14,15,13,5,7, 6,1,7,0,9,14,8,15,5,2,4,3,10,13,11,12,
    3,1,5,7,6,4,0,2,10,8,12,14,15,13,9,11, 6,1,7,0,5,2,4,3,9,14,8,15,10,13,11,12,
    3,11,10,2,6,14,15,7,5,13,12,4,0,8,9,1, 12,15,3,0,11,8,4,7,13,14,2,1,10,9,5,6,
    3,1,5,7,15,13,9,11,10,8,12,14,6,4,0,2, 14,1,15,0,13,2,12,3,9,6,8,7,10,5,11,4,
    3,2,10,11,15,14,6,7,5,4,12,13,9,8,0,1, 14,15,1,0,9,8,6,7,13,12,2,3,10,11,5,4,
    3,2,0,1,5,4,6,7,15,14,12,13,9,8,10,11, 2,3,1,0,5,4,6,7,13,12,14,15,10,11,9,8,
    3,11,15,7,5,13,9,1,0,8,12,4,6,14,10,2, 8,7,15,0,11,4,12,3,9,6,14,1,10,5,13,2,
    3,11,10,2,0,8,9,1,5,13,12,4,6,14,15,7, 4,7,3,0,11,8,12,15,5,6,2,1,10,9,13,14,
    3,7,5,1,0,4,6,2,10,14,12,8,9,13,15,11, 4,3,7,0,5,2,6,1,11,12,8,15,10,13,9,14,
    3,7,15,11,10,14,6,2,0,4,12,8,9,13,5,1, 8,15,7,0,9,14,6,1,11,12,4,3,10,13,5,2,
    3,1,0,2,10,8,9,11,15,13,12,14,6,4,5,7, 2,1,3,0,13,14,12,15,5,6,4,7,10,9,11,8,
    10,2,0,8,12,4,6,14,15,7,5,13,9,1,3,11, 2,13,1,14,5,10,6,9,3,12,0,15,4,11,7,8,
    10,11,15,14,12,13,9,8,0,1,5,4,6,7,3,2, 8,9,15,14,11,10,12,13,7,6,0,1,4,5,3,2,
    10,11,3,2,0,1,9,8,12,13,5,4,6,7,15,14, 4,5,3,2,11,10,12,13,7,6,0,1,8,9,15,14,
    10,14,12,8,0,4,6,2,3,7,5,1,9,13,15,11, 4,11,7,8,5,10,6,9,3,12,0,15,2,13,1,14,
    10,14,15,11,3,7,6,2,0,4,5,1,9,13,12,8, 8,11,7,4,9,10,6,5,15,12,0,3,14,13,1,2,
    10,8,0,2,3,1,9,11,15,13,5,7,6,4,12,14, 2,5,3,4,13,10,12,11,1,6,0,7,14,9,15,8,
    10,8,12,14,15,13,9,11,3,1,5,7,6,4,0,2, 14,9,15,8,13,10,12,11,1,6,0,7,2,5,3,4,
    10,2,3,11,15,7,6,14,12,4,5,13,9,1,0,8, 14,13,1,2,9,10,6,5,15,12,0,3,8,11,7,4,
    10,8,12,14,6,4,0,2,3,1,5,7,15,13,9,11, 6,9,7,8,5,10,4,11,1,14,0,15,2,13,3,12,
    10,11,3,2,6,7,15,14,12,13,5,4,0,1,9,8, 12,13,3,2,11,10,4,5,15,14,0,1,8,9,7,6,
    10,11,9,8,12,13,15,14,6,7,5,4,0,1,3,2, 12,13,15,14,11,10,8,9,3,2,0,1,4,5,7,6,
    10,2,6,14,12,4,0,8,9,1,5,13,15,7,3,11, 6,9,1,14,5,10,2,13,7,8,0,15,4,11,3,12,
    10,2,3,11,9,1,0,8,12,4,5,13,15,7,6,14, 6,5,1,2,9,10,14,13,7,4,0,3,8,11,15,12,
    10,14,12,8,9,13,15,11,3,7,5,1,0,4,6,2, 12,11,15,8,13,10,14,9,3,4,0,7,2,5,1,6,
    10,14,6,2,3,7,15,11,9,13,5,1,0,4,12,8, 12,11,3,4,13,10,2,5,15,8,0,7,14,9,1,6,
    10,8,9,11,3,1,0,2,6,4,5,7,15,13,12,14, 6,5,7,4,9,10,8,11,1,2,0,3,14,13,15,12,
    10,14,6,2,0,4,12,8,9,13,5,1,3,7,15,11, 4,11,3,12,5,10,2,13,7,8,0,15,6,9,1,14,
    10,11,9,8,0,1,3,2,6,7,5,4,12,13,15,14, 4,5,7,6,11,10,8,9,3,2,0,1,12,13,15,14,
    10,11,15,14,6,7,3,2,0,1,5,4,12,13,9,8, 8,9,7,6,11,10,4,5,15,14,0,1,12,13,3,2,
    10,8,0,2,6,4,12,14,15,13,5,7,3,1,9,11, 2,13,3,12,5,10,4,11,1,14,0,15,6,9,7,8,
    10,8,9,11,15,13,12,14,6,4,5,7,3,1,0,2, 14,13,15,12,9,10,8,11,1,2,0,3,6,5,7,4,
    10,2,6,14,15,7,3,11,9,1,5,13,12,4,0,8, 14,9,1,6,13,10,2,5,15,8,0,7,12,11,3,4,
    10,2,0,8,9,1,3,11,15,7,5,13,12,4,6,14, 2,5,1,6,13,10,14,9,3,4,0,7,12,11,15,8,
    10,14,15,11,9,13,12,8,0,4,5,1,3,7,6,2, 8,11,15,12,9,10,14,13,7,4,0,3,6,5,1,2,
    9,11,15,13,5,7,3,1,0,2,6,4,12,14,10,8, 8,7,9,6,11,4,10,5,15,0,14,1,12,3,13,2,
    9,8,0,1,5,4,12,13,15,14,6,7,3,2,10,11, 2,3,13,12,5,4,10,11,1,0,14,15,6,7,9,8,
    9,8,10,11,15,14,12,13,5,4,6,7,3,2,0,1, 14,15,13,12,9,8,10,11,1,0,2,3,6,7,5,4,
    9,1,5,13,15,7,3,11,10,2,6,14,12,4,0,8, 14,1,9,6,13,2,10,5,15,0,8,7,12,3,11,4,
    9,1,0,8,10,2,3,11,15,7,6,14,12,4,5,13, 2,1,5,6,13,14,10,9,3,0,4,7,12,15,11,8,
    9,13,15,11,10,14,12,8,0,4,6,2,3,7,5,1, 8,15,11,12,9,14,10,13,7,0,4,3,6,1,5,2,
    9,13,5,1,0,4,12,8,10,14,6,2,3,7,15,11, 4,3,11,12,5,2,10,13,7,0,8,15,6,1,9,14,
    9,11,10,8,0,2,3,1,5,7,6,4,12,14,15,13, 4,7,5,6,11,8,10,9,3,0,2,1,12,15,13,14,
    9,13,5,1,3,7,15,11,10,14,6,2,0,4,12,8, 12,3,11,4,13,2,10,5,15,0,8,7,14,1,9,6,
    9,8,10,11,3,2,0,1,5,4,6,7,15,14,12,13, 6,7,5,4,9,8,10,11,1,0,2,3,14,15,13,12,
    9,8,12,13,5,4,0,1,3,2,6,7,15,14,10,11, 6,7,9,8,5,4,10,11,1,0,14,15,2,3,13,12,
    9,11,3,1,5,7,15,13,12,14,6,4,0,2,10,8, 12,3,13,2,11,4,10,5,15,0,14,1,8,7,9,6,
    9,11,10,8,12,14,15,13,5,7,6,4,0,2,3,1, 12,15,13,14,11,8,10,9,3,0,2,1,4,7,5,6,
    9,1,5,13,12,4,0,8,10,2,6,14,15,7,3,11, 6,1,9,14,5,2,10,13,7,0,8,15,4,3,11,12,
    9,1,3,11,10,2,0,8,12,4,6,14,15,7,5,13, 6,1,5,2,9,14,10,13,7,0,4,3,8,15,11,12,
    9,13,12,8,10,14,15,11,3,7,6,2,0,4,5,1, 12,15,11,8,13,14,10,9,3,0,4,7,2,1,5,6,
    9,1,3,11,15,7,5,13,12,4,6,14,10,2,0,8, 14,1,13,2,9,6,10,5,15,0,12,3,8,7,11,4,
    9,8,12,13,15,14,10,11,3,2,6,7,5,4,0,1, 14,15,9,8,13,12,10,11,1,0,6,7,2,3,5,4,
    9,8,0,1,3,2,10,11,15,14,6,7,5,4,12,13, 2,3,5,4,13,12,10,11,1,0,6,7,14,15,9,8,
    9,13,15,11,3,7,5,1,0,4,6,2,10,14,12,8, 8,7,11,4,9,6,10,5,15,0,12,3,14,1,13,2,
    9,13,12,8,0,4,5,1,3,7,6,2,10,14,15,11, 4,7,11,8,5,6,10,9,3,0,12,15,2,1,13,14,
    9,11,3,1,0,2,10,8,12,14,6,4,5,7,15,13, 4,3,5,2,11,12,10,13,7,0,6,1,8,15,9,14,
    9,11,15,13,12,14,10,8,0,2,6,4,5,7,3,1, 8,15,9,14,11,12,10,13,7,0,6,1,4,3,5,2,
    9,1,0,8,12,4,5,13,15,7,6,14,10,2,3,11, 2,1,13,14,5,6,10,9,3,0,12,15,4,7,11,8,
    5,1,0,4,12,8,9,13,15,11,10,14,6,2,3,7, 2,1,13,14,3,0,12,15,5,6,10,9,4,7,11,8,
    5,7,15,13,12,14,6,4,0,2,10,8,9,11,3,1, 8,15,9,14,7,0,6,1,11,12,10,13,4,3,5,2,
    5,7,3,1,0,2,6,4,12,14,10,8,9,11,15,13, 4,3,5,2,7,0,6,1,11,12,10,13,8,15,9,14,
    5,13,12,4,0,8,9,1,3,11,10,2,6,14,15,7, 4,7,11,8,3,0,12,15,5,6,10,9,2,1,13,14,
    5,13,15,7,3,11,9,1,0,8,10,2,6,14,12,4, 8,7,11,4,15,0,12,3,9,6,10,5,14,1,13,2,
    5,4,0,1,3,2,6,7,15,14,10,11,9,8,12,13, 2,3,5,4,1,0,6,7,13,12,10,11,14,15,9,8,
    5,4,12,13,15,14,6,7,3,2,10,11,9,8,0,1, 14,15,9,8,1,0,6,7,13,12,10,11,2,3,5,4,
    5,1,3,7,15,11,9,13,12,8,10,14,6,2,0,4, 14,1,13,2,15,0,12,3,9,6,10,5,8,7,11,4,
    5,4,12,13,9,8,0,1,3,2,10,11,15,14,6,7, 6,7,9,8,1,0,14,15,5,4,10,11,2,3,13,12,
    5,7,3,1,9,11,15,13,12,14,10,8,0,2,6,4, 12,3,13,2,15,0,14,1,11,4,10,5,8,7,9,6,
    5,7,6,4,12,14,15,13,9,11,10,8,0,2,3,1, 12,15,13,14,3,0,2,1,11,8,10,9,4,7,5,6,
    5,1,9,13,12,8,0,4,6,2,10,14,15,11,3,7, 6,1,9,14,7,0,8,15,5,2,10,13,4,3,11,12,
    5,1,3,7,6,2,0,4,12,8,10,14,15,11,9,13, 6,1,5,2,7,0,4,3,9,14,10,13,8,15,11,12,
    5,13,12,4,6,14,15,7,3,11,10,2,0,8,9,1, 12,15,11,8,3,0,4,7,13,14,10,9,2,1,5,6,
    5,13,9,1,3,11,15,7,6,14,10,2,0,8,12,4, 12,3,11,4,15,0,8,7,13,2,10,5,14,1,9,6,
    5,4,6,7,3,2,0,1,9,8,10,11,15,14,12,13, 6,7,5,4,1,0,2,3,9,8,10,11,14,15,13,12,
    5,13,9,1,0,8,12,4,6,14,10,2,3,11,15,7, 4,3,11,12,7,0,8,15,5,2,10,13,6,1,9,14,
    5,7,6,4,0,2,3,1,9,11,10,8,12,14,15,13, 4,7,5,6,3,0,2,1,11,8,10,9,12,15,13,14,
    5,7,15,13,9,11,3,1,0,2,10,8,12,14,6,4, 8,7,9,6,15,0,14,1,11,4,10,5,12,3,13,2,
    5,4,0,1,9,8,12,13,15,14,10,11,3,2,6,7, 2,3,13,12,1,0,14,15,5,4,10,11,6,7,9,8,
    5,4,6,7,15,14,12,13,9,8,10,11,3,2,0,1, 14,15,13,12,1,0,2,3,9,8,10,11,6,7,5,4,
    5,1,9,13,15,11,3,7,6,2,10,14,12,8,0,4, 14,1,9,6,15,0,8,7,13,2,10,5,12,3,11,4,
    5,1,0,4,6,2,3,7,15,11,10,14,12,8,9,13, 2,1,5,6,3,0,4,7,13,14,10,9,12,15,11,8,
    5,13,15,7,6,14,12,4,0,8,10,2,3,11,9,1, 8,15,11,12,7,0,4,3,9,14,10,13,6,1,5,2,
    6,7,15,14,10,11,3,2,0,1,9,8,12,13,5,4, 8,9,7,6,15,14,0,1,11,10,4,5,12,13,3,2,
    6,4,0,2,10,8,12,14,15,13,9,11,3,1,5,7, 2,13,3,12,1,14,0,15,5,10,4,11,6,9,7,8,
    6,4,5,7,15,13,12,14,10,8,9,11,3,1,0,2, 14,13,15,12,1,2,0,3,9,10,8,11,6,5,7,4,
    6,2,10,14,15,11,3,7,5,1,9,13,12,8,0,4, 14,9,1,6,15,8,0,7,13,10,2,5,12,11,3,4,
    6,2,0,4,5,1,3,7,15,11,9,13,12,8,10,14, 2,5,1,6,3,4,0,7,13,10,14,9,12,11,15,8,
    6,14,15,7,5,13,12,4,0,8,9,1,3,11,10,2, 8,11,15,12,7,4,0,3,9,10,14,13,6,5,1,2,
    6,14,10,2,0,8,12,4,5,13,9,1,3,11,15,7, 4,11,3,12,7,8,0,15,5,10,2,13,6,9,1,14,
    6,7,5,4,0,1,3,2,10,11,9,8,12,13,15,14, 4,5,7,6,3,2,0,1,11,10,8,9,12,13,15,14,
    6,14,10,2,3,11,15,7,5,13,9,1,0,8,12,4, 12,11,3,4,15,8,0,7,13,10,2,5,14,9,1,6,
    6,4,5,7,3,1,0,2,10,8,9,11,15,13,12,14, 6,5,7,4,1,2,0,3,9,10,8,11,14,13,15,12,
    6,4,12,14,10,8,0,2,3,1,9,11,15,13,5,7, 6,9,7,8,1,14,0,15,5,10,4,11,2,13,3,12,
    6,7,3,2,10,11,15,14,12,13,9,8,0,1,5,4, 12,13,3,2,15,14,0,1,11,10,4,5,8,9,7,6,
    6,7,5,4,12,13,15,14,10,11,9,8,0,1,3,2, 12,13,15,14,3,2,0,1,11,10,8,9,4,5,7,6,
    6,2,10,14,12,8,0,4,5,1,9,13,15,11,3,7, 6,9,1,14,7,8,0,15,5,10,2,13,4,11,3,12,
    6,2,3,7,5,1,0,4,12,8,9,13,15,11,10,14, 6,5,1,2,7,4,0,3,9,10,14,13,8,11,15,12,
    6,14,12,4,5,13,15,7,3,11,9,1,0,8,10,2, 12,11,15,8,3,4,0,7,13,10,14,9,2,5,1,6,
    6,2,3,7,15,11,10,14,12,8,9,13,5,1,0,4, 14,13,1,2,15,12,0,3,9,10,6,5,8,11,7,4,
    6,4,12,14,15,13,5,7,3,1,9,11,10,8,0,2, 14,9,15,8,1,6,0,7,13,10,12,11,2,5,3,4,
    6,4,0,2,3,1,5,7,15,13,9,11,10,8,12,14, 2,5,3,4,1,6,0,7,13,10,12,11,14,9,15,8,
    6,14,15,7,3,11,10,2,0,8,9,1,5,13,12,4, 8,11,7,4,15,12,0,3,9,10,6,5,14,13,1,2,
    6,14,12,4,0,8,10,2,3,11,9,1,5,13,15,7, 4,11,7,8,3,12,0,15,5,10,6,9,2,13,1,14,
    6,7,3,2,0,1,5,4,12,13,9,8,10,11,15,14, 4,5,3,2,7,6,0,1,11,10,12,13,8,9,15,14,
    6,7,15,14,12,13,5,4,0,1,9,8,10,11,3,2, 8,9,15,14,7,6,0,1,11,10,12,13,4,5,3,2,
    6,2,0,4,12,8,10,14,15,11,9,13,5,1,3,7, 2,13,1,14,3,12,0,15,5,10,6,9,4,11,7,8,
};

template <>
const int HilbertData<4>::m_HILBERT_TABLE[] = {
    1,7,95,18,156,187,78,184,123,105,48,117,156,35,64,32,
    0,122,104,104,43,49,27,68,2,157,190,82,22,164,185,85,
    3,121,111,140,44,54,111,58,1,13,72,75,158,151,181,178,
    2,159,169,188,4,77,9,74,45,159,55,62,110,133,103,130,
    5,152,170,152,3,142,109,86,46,26,52,69,19,137,116,81,
    4,90,6,23,171,53,153,165,108,79,139,136,108,65,39,36,
    7,89,5,154,172,50,12,150,107,107,28,129,189,63,31,134,
    6,0,88,106,155,25,155,177,124,8,51,102,141,30,59,182,
    9,81,15,167,191,174,2,164,103,81,33,134,57,52,45,129,
    8,138,102,56,32,127,32,48,10,6,80,88,162,165,178,181,
    11,137,101,59,39,39,113,67,9,87,21,74,168,152,171,155,
    10,12,86,173,144,17,160,170,136,38,86,108,58,31,71,111,
    13,11,85,119,151,37,85,182,143,3,126,116,61,44,49,177,
    12,150,186,60,14,84,7,89,36,36,175,53,114,130,117,133,
    15,149,185,63,13,20,83,75,35,161,35,70,120,123,104,107,
    14,82,8,125,184,82,34,156,96,112,16,122,62,66,30,159,
    17,163,177,147,23,10,91,88,31,68,177,61,128,136,107,104,
    16,129,115,99,30,129,69,57,18,94,14,89,176,155,184,152,
    19,130,116,130,29,41,66,55,17,145,183,85,5,150,170,82,
    18,73,20,1,182,92,148,151,117,67,131,123,182,51,36,39,
    21,74,19,132,181,181,47,159,118,64,11,137,98,56,44,156,
    20,22,75,180,166,15,65,185,133,42,133,111,146,45,60,108,
    23,21,76,4,165,179,70,171,134,100,93,103,134,32,50,35,
    22,164,178,178,16,135,97,81,24,71,40,54,0,122,102,86,
    42,138,111,56,31,141,118,59,23,6,72,88,25,163,177,163,
    119,92,10,13,176,73,162,158,119,51,136,139,24,67,26,46,
    190,82,3,125,183,78,37,149,112,112,0,122,27,68,25,161,
    142,2,86,173,160,5,160,170,135,33,79,101,26,28,69,113,
    167,14,93,9,167,184,50,187,128,114,76,110,29,27,70,43,
    47,159,186,60,18,77,7,89,30,166,189,63,28,129,115,115,
    165,148,174,169,17,74,116,87,29,36,52,55,31,130,116,143,
    117,126,100,121,24,49,32,54,16,164,75,83,30,164,178,191,
    153,158,169,174,9,132,9,76,39,26,57,69,33,137,101,81,
    80,72,16,122,188,181,30,159,100,56,8,125,32,56,34,156,
    83,91,2,164,187,63,15,167,99,107,45,129,35,63,33,134,
    131,128,115,112,36,41,62,55,14,154,14,95,34,150,186,82,
    105,121,110,126,35,61,27,68,13,13,180,73,37,149,185,85,
    84,17,77,170,148,12,60,173,140,31,133,111,36,38,60,108,
    87,3,94,116,147,44,155,177,139,11,59,119,39,37,59,182,
    179,163,176,160,10,10,106,90,32,58,40,54,38,138,102,86,
    145,16,85,119,157,34,92,186,126,19,126,116,41,47,55,175,
    97,81,17,167,174,174,22,164,109,93,38,138,40,54,42,127,
    173,79,13,10,173,65,151,148,110,90,120,128,43,53,41,29,
    25,149,177,63,1,20,91,75,44,146,172,60,42,121,111,121,
    171,160,187,163,2,122,88,80,43,71,35,68,45,122,104,96,
    123,139,112,115,46,39,66,69,3,89,170,84,44,152,170,144,
    28,129,101,59,45,124,98,56,4,94,21,74,47,159,169,169,
    125,9,78,14,158,168,95,176,125,103,64,100,46,40,48,24,
    80,83,16,157,188,191,30,157,112,104,0,126,66,55,47,49,
    37,161,189,70,20,1,75,91,34,156,186,50,105,127,105,48,
    179,147,176,155,10,6,106,92,27,61,24,49,136,120,106,51,
    131,128,99,107,26,29,57,52,14,154,7,93,184,154,168,50,
    33,141,113,67,38,138,108,51,21,74,4,94,153,153,175,53,
    84,17,87,109,160,5,152,174,140,31,143,109,71,42,54,52,
    135,15,79,185,146,2,60,173,132,110,76,110,159,41,53,55,
    183,78,11,137,180,73,158,158,98,56,3,125,111,48,40,54,
    132,4,93,9,161,171,50,187,135,33,81,99,166,33,63,57,
    176,160,179,163,0,122,80,80,40,54,32,58,8,142,98,56,
    158,153,174,169,19,87,116,87,26,39,69,57,11,141,97,59,
    106,79,23,10,175,65,165,148,109,86,38,140,172,60,38,58,
    180,92,1,13,183,85,37,147,113,51,123,139,118,59,37,61,
    128,131,112,115,41,36,55,62,5,84,170,84,12,146,190,60,
    110,126,105,121,27,68,35,61,22,164,83,83,15,145,189,63,
    154,18,78,14,157,34,82,188,127,117,64,100,124,34,56,62,
    105,140,110,133,40,58,43,71,13,20,180,79,151,162,180,65,
    153,158,188,181,9,132,21,78,41,46,62,66,103,132,114,64,
    34,146,175,53,7,89,18,77,37,149,182,65,131,131,113,67,
    87,3,84,183,147,44,144,183,126,19,130,112,49,25,68,66,
    157,8,92,102,154,176,95,176,141,16,59,119,129,26,67,69,
    109,93,12,150,189,63,17,167,106,90,128,128,177,70,27,68,
    83,80,2,135,174,178,22,160,99,96,45,135,52,69,28,71,
    38,127,98,48,33,134,101,64,6,23,88,72,179,161,179,70,
    162,162,180,73,11,137,119,79,40,54,24,71,8,142,106,90,
    94,19,74,72,155,25,163,181,143,3,140,118,61,44,58,118,
    164,20,73,75,145,30,85,182,127,117,48,117,124,34,51,186,
    116,76,21,74,175,53,165,165,97,81,31,134,172,50,38,138,
    91,75,22,77,191,188,2,166,107,115,28,133,57,62,45,166,
    114,132,114,76,41,29,55,66,15,167,185,78,12,154,190,95,
    139,131,113,79,39,26,113,65,7,84,4,77,168,144,171,160,
    187,161,179,78,6,1,80,72,35,161,27,64,120,123,96,112,
    33,141,97,81,28,129,57,57,21,74,9,87,171,155,168,152,
    157,8,82,80,154,176,78,184,124,8,56,96,127,24,64,32,
    135,15,81,83,166,15,63,191,132,110,93,103,161,43,50,35,
    34,146,190,82,7,89,14,84,47,159,62,62,117,133,114,130,
    37,145,189,85,20,13,75,83,25,61,177,61,123,120,107,104,
    109,86,12,84,172,60,12,144,106,79,128,136,175,65,29,36,
    183,85,11,87,180,92,158,151,118,59,11,143,113,51,46,39,
    38,142,98,86,42,58,111,58,6,10,88,80,165,162,181,178,
    170,95,7,89,190,82,44,156,113,67,123,123,118,64,37,149,
    122,6,90,88,161,171,70,171,142,45,86,108,166,33,65,101,
    77,5,89,91,144,17,147,172,133,42,121,107,58,31,61,172,
    120,120,106,90,27,68,43,49,12,150,173,92,15,145,180,73,
    100,127,105,93,32,127,40,50,20,23,83,91,162,165,191,174,
    148,153,175,92,21,87,18,94,36,41,175,51,114,143,117,126,
    168,154,168,95,8,125,102,93,26,46,69,52,11,132,97,76,
    72,88,0,94,181,169,47,155,96,99,16,124,62,57,30,124,
    114,130,103,143,43,49,40,54,15,145,97,81,22,164,191,191,
    142,2,80,190,160,5,144,183,142,45,96,98,71,42,58,118,
    94,19,87,109,155,25,147,172,141,16,99,97,141,30,57,189,
    47,159,188,188,18,77,21,74,34,146,98,56,105,121,100,140,
    187,161,187,70,6,1,88,91,33,149,101,63,139,131,99,115,
    119,90,10,23,182,53,148,165,102,93,100,138,186,50,32,138,
    185,78,9,137,173,73,151,158,101,64,103,137,108,67,39,46,
    136,120,96,104,26,29,69,66,8,150,102,82,184,154,184,95,
    180,73,1,1,183,78,37,149,111,48,105,121,98,56,44,156,
    157,8,92,102,145,30,85,182,120,0,104,106,49,25,49,177,
    154,18,91,7,167,184,50,187,127,117,107,105,134,32,50,35,
    38,138,108,51,33,141,113,67,6,6,106,90,179,163,171,155,
    128,131,107,99,29,26,52,57,5,84,109,89,5,144,172,152,
    153,148,173,53,9,74,4,77,41,36,108,53,103,130,110,133,
    158,162,174,178,19,137,116,81,42,54,111,54,3,142,109,86,
    83,80,2,135,191,188,2,166,104,112,110,122,55,66,43,159,
    168,152,176,160,8,142,106,90,26,26,113,67,11,137,119,79,
    125,9,78,14,156,187,78,184,132,110,112,114,161,43,66,27,
    77,5,77,170,146,2,60,173,131,28,115,113,166,33,65,101,
    97,81,17,167,172,50,12,150,116,76,114,130,175,53,29,29,
    75,91,18,164,188,191,30,157,115,107,117,129,62,57,30,124,
    31,141,118,59,42,138,111,56,19,74,116,74,165,153,181,169,
    100,121,117,126,32,54,24,49,20,13,119,73,162,151,182,73,
    25,147,183,163,1,6,72,80,25,61,118,68,123,120,112,96,
    91,91,22,164,189,63,17,167,105,121,42,127,172,50,38,138,
    142,2,86,173,166,15,65,185,122,120,90,104,161,43,70,43,
    89,1,77,170,147,44,144,183,121,123,133,111,61,44,58,118,
    45,124,98,56,28,129,101,59,0,122,88,88,171,155,179,163,
    139,123,115,112,39,46,69,66,7,125,14,95,168,156,184,95,
    47,157,188,169,18,94,21,87,47,124,62,55,117,126,114,143,
    162,158,178,174,11,132,97,76,40,127,40,48,8,125,102,93,
    119,92,10,13,182,92,148,151,106,126,128,120,175,49,29,41,
    176,160,168,152,0,122,102,86,28,129,69,69,16,135,97,81,
    84,17,87,109,144,17,147,172,130,128,126,116,68,29,49,177,
    154,18,95,18,157,34,92,186,129,131,67,115,141,30,59,182,
    183,78,11,137,190,82,44,156,114,130,19,132,66,66,47,159,
    180,77,1,20,173,65,151,148,113,133,123,131,108,65,39,36,
    33,134,101,64,38,127,98,48,21,132,21,76,153,165,169,181,
    110,133,105,140,43,71,40,58,22,135,83,75,22,166,191,178,
    179,167,187,70,10,23,88,91,27,134,35,70,136,128,104,107,
    103,143,114,130,40,54,43,49,11,137,185,85,151,151,180,73,
    150,10,92,102,167,184,95,176,138,136,51,102,134,32,48,24,
    125,9,76,4,156,187,70,171,137,139,79,101,149,39,65,101,
    148,148,175,53,21,74,18,77,38,138,186,60,100,140,105,121,
    25,147,177,147,1,6,91,88,37,141,189,59,131,139,115,99,
    72,80,0,135,181,188,47,166,98,140,3,142,98,58,44,146,
    97,87,17,145,174,191,22,157,97,143,31,141,52,57,28,124,
    120,136,104,96,29,26,66,69,12,142,190,86,5,144,170,144,
    162,151,178,191,11,145,97,85,46,41,52,55,19,143,116,143,
    190,84,3,142,190,144,44,146,112,96,0,135,66,62,47,166,
    91,83,22,157,189,147,17,145,107,99,28,124,189,61,31,141,
    42,140,111,140,38,146,98,60,23,20,72,75,153,148,169,188,
    139,139,113,67,37,149,101,59,7,89,4,94,187,147,179,163,
    167,14,95,18,150,148,92,186,134,100,48,117,138,36,51,186,
    137,13,79,185,149,151,65,185,125,103,76,110,156,35,70,43,
    184,144,168,152,12,150,102,86,27,68,24,71,136,136,106,90,
    132,4,76,4,159,153,53,169,135,33,79,101,146,45,60,108,
    109,93,12,150,168,152,5,154,97,81,31,134,52,52,28,129,
    106,94,23,6,175,155,165,153,119,51,136,139,182,51,36,39,
    34,156,186,50,7,154,7,95,37,161,189,70,131,123,115,107,
    176,155,179,147,0,157,80,88,24,49,27,61,0,124,96,104,
    105,125,100,48,40,156,32,48,13,1,75,72,151,158,178,181,
    110,126,114,130,47,159,55,55,22,164,185,85,2,157,190,82,
    87,3,84,183,152,158,160,170,143,3,140,118,54,46,71,111,
    120,128,104,112,27,161,27,70,12,154,190,95,15,167,185,78,
    173,79,13,10,180,160,158,162,108,79,139,136,113,71,46,26,
    72,72,0,122,179,163,25,161,98,56,3,125,118,64,37,149,
    145,16,85,119,164,162,73,178,124,8,51,102,127,24,48,24,
    74,23,94,116,163,165,155,177,140,31,143,109,58,31,61,172,
    30,166,189,63,22,164,75,75,47,159,186,60,117,133,105,121,
    148,165,169,174,21,167,9,76,36,29,55,52,114,134,103,76,
    28,135,99,115,28,166,57,69,4,77,7,84,171,160,168,144,
    94,19,94,116,153,47,169,175,141,16,59,119,124,34,51,186,
    190,82,3,125,170,95,168,152,118,64,11,137,113,67,46,46,
    88,72,4,122,169,181,171,159,99,96,45,135,57,62,45,166,
    44,146,172,60,5,89,170,89,25,149,177,63,123,131,107,115,
    187,163,171,160,6,10,173,90,35,68,43,71,120,136,108,90,
    42,140,109,121,42,58,172,54,23,20,91,83,165,162,174,191,
    114,130,110,126,41,41,175,53,15,145,180,73,12,150,173,92,
    167,14,93,9,154,176,174,168,134,100,93,103,127,24,52,40,
    128,120,112,104,25,68,177,68,5,150,170,82,17,145,183,85,
    80,83,16,157,178,174,176,164,96,99,16,124,69,52,24,129,
    106,90,23,23,177,70,179,163,109,93,38,138,189,63,31,134,
    135,15,79,185,162,22,178,180,142,45,86,108,71,42,71,111,
    132,4,72,21,161,171,181,179,125,103,64,100,156,35,64,32,
    37,149,182,65,20,20,180,73,34,146,175,53,105,121,117,133,
    158,153,181,188,19,87,183,74,46,41,66,62,19,143,118,130,
    131,139,119,67,26,39,182,67,14,89,18,94,184,152,176,155,
    151,162,191,178,15,137,185,81,41,46,55,52,103,132,103,76,
    102,93,14,150,186,50,184,150,119,90,136,128,182,53,36,29,
    173,73,13,1,185,78,187,149,108,67,139,123,101,64,35,149,
    100,127,100,48,34,138,186,56,20,23,75,72,148,153,188,169,
    28,129,99,99,33,141,189,63,4,94,7,89,179,163,187,147,
    77,5,84,183,146,2,188,190,133,42,140,118,146,45,62,98,
    145,16,83,97,145,30,191,189,126,19,143,109,49,25,61,172,
    168,152,184,144,8,142,190,82,24,71,27,68,0,122,96,96,
};

