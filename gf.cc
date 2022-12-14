/**
 * @file gf.cc
 * @author Michal Findra, xfindr00
 * @brief Galois field operations
 * @date 12.12.2022
 */
#include "gf.hh"

#include <math.h>

BMS_GF::BMS_GF() {
}

int GF_field_exp[256] = {1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19,
                         53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34,
                         102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112,
                         144, 171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104,
                         184, 211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98,
                         166, 241, 8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220,
                         127, 129, 152, 179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16,
                         48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125,
                         135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 251, 22,
                         58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 195,
                         94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218,
                         117, 159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223,
                         122, 142, 137, 128, 155, 182, 193, 88, 232, 35, 101, 175, 234, 37,
                         111, 177, 200, 67, 197, 84, 252, 31, 33, 99, 165, 244, 7, 9, 27,
                         45, 119, 153, 176, 203, 70, 202, 69, 207, 74, 222, 121, 139, 134,
                         145, 168, 227, 62, 66, 198, 81, 243, 14, 18, 54, 90, 238, 41, 123,
                         141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 57, 75, 221, 124,
                         132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1};

int GF_field_log[256] = {-1, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223,
                         3, 100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105,
                         28, 193, 125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114,
                         154, 201, 9, 120, 101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130,
                         69, 53, 147, 218, 142, 150, 143, 219, 189, 54, 208, 206, 148, 19,
                         92, 210, 241, 64, 70, 131, 56, 102, 221, 253, 48, 191, 6, 139, 98,
                         179, 37, 226, 152, 34, 136, 145, 16, 126, 110, 72, 195, 163, 182,
                         30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 43, 121, 10, 21, 155,
                         159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 175, 88, 168,
                         80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 44,
                         215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81,
                         160, 127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164,
                         118, 123, 183, 204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161,
                         108, 170, 85, 41, 157, 151, 178, 135, 144, 97, 190, 220, 252, 188,
                         149, 207, 205, 55, 63, 91, 209, 83, 57, 132, 60, 65, 162, 109, 71,
                         20, 42, 158, 93, 86, 242, 211, 171, 68, 17, 146, 217, 35, 32, 46,
                         137, 180, 124, 184, 38, 119, 153, 227, 165, 103, 74, 237, 222, 197,
                         49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7};

int BMS_GF::gf_mul(int x, int y) {
    if (x == 0 || y == 0)
        return 0;
    return GF_field_exp[((GF_field_log[x] + GF_field_log[y]) % 255)];  // should be gf_exp[(gf_log[x]+gf_log[y])%255] if gf_exp wasn't oversized
}

int BMS_GF::gf_pow(unsigned int x, int power) {
    if ((GF_field_log[x] * power) < 0) {
        return GF_field_exp[(255 + (GF_field_log[x] * power)) % 255];
    } else
        return GF_field_exp[(GF_field_log[x] * power) % 255];
}

std::vector<int> BMS_GF::gf_poly_mul(std::vector<int> p, std::vector<int> q) {
    // Multiply two polynomials, inside Galois Field
    //  Pre-allocate the result array
    std::vector<int> r(p.size() + q.size() - 1);
    for (int j = 0; j < q.size(); ++j) {
        for (int i = 0; i < p.size(); ++i) {
            r[i + j] ^= gf_mul(p[i], q[j]);  // equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))
                                             // -- you can see it's your usual polynomial multiplication
        }
    }
    return r;
}

int BMS_GF::gf_div(int x, int y) {
    if (y == 0) {
        exit(1);
    }
    if (x == 0) {
        return 0;
    }
    return GF_field_exp[(GF_field_log[x] + 255 - GF_field_log[y]) % 255];
}

int BMS_GF::gf_sub(int x, int y) {
    return x ^ y;
}

std::vector<int> BMS_GF::gf_poly_add(std::vector<int> p, std::vector<int> q) {
    std::vector<int> r(std::max(p.size(), q.size()), 0);
    for (std::size_t i = 0; i < p.size(); ++i) {
        r[i + r.size() - p.size()] = p[i];
    }
    for (std::size_t i = 0; i < q.size(); ++i) {
        r[i + r.size() - q.size()] ^= q[i];
    }
    return r;
}

int BMS_GF::gf_inverse(int x) {
    return GF_field_exp[255 - GF_field_log[x]];
}

std::vector<int> BMS_GF::gf_poly_div(const std::vector<int>& dividend, const std::vector<int>& divisor) {
    std::vector<int> msg_out = dividend;  // Copy the dividend
    for (int i = 0; i < dividend.size() - (divisor.size() - 1); i++) {
        int coef = msg_out[i];  // precaching
        if (coef != 0) {
            for (int j = 1; j < divisor.size(); j++) {
                if (divisor[j] != 0) {
                    msg_out[i + j] ^= gf_mul(divisor[j], coef);  // equivalent to the more mathematically correct
                }
            }
        }
    }

    int separator = -(divisor.size() - 1);
    std::vector<int> b;
    for (auto it = msg_out.end() - 1; it >= msg_out.end() - abs(separator); it--) {
        b.insert(b.begin(), (*it));
    }
    return b;  // return quotient, remainder
}

std::vector<int> BMS_GF::gf_poly_scale(const std::vector<int>& p, int x) {
    std::vector<int> r(p.size());
    for (int i = 0; i < p.size(); i++) {
        r[i] = gf_mul(p[i], x);
    }
    return r;
}

int BMS_GF::gf_poly_eval(const std::vector<int>& poly, int x) {
    int y = poly[0];
    for (int i = 1; i < poly.size(); ++i) {
        y = gf_mul(y, x) ^ poly[i];
    }
    return y;
}
