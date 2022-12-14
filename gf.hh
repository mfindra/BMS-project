/**
 * @file gf.hh
 * @author Michal Findra, xfindr00
 * @brief Galois field operations
 * @date 12.12.2022
 */
#include <vector>
#ifndef GF_H
#define GF_H

class BMS_GF {
   public:
    BMS_GF();
    int gf_mul(int x, int y);
    int gf_pow(unsigned int x, int power);
    std::vector<int> gf_poly_mul(std::vector<int> p, std::vector<int> q);
    int gf_div(int x, int y);
    int gf_sub(int x, int y);
    std::vector<int> gf_poly_add(std::vector<int> p, std::vector<int> q);
    int gf_inverse(int x);
    std::vector<int> gf_poly_div(const std::vector<int>& dividend, const std::vector<int>& divisor);
    std::vector<int> gf_poly_scale(const std::vector<int>& p, int x);
    int gf_poly_eval(const std::vector<int>& poly, int x);

   private:
};

#endif