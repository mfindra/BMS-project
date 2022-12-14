/**
 * @file rs.hh
 * @author Michal Findra, xfindr00
 * @brief Reed-Salomon code operations
 * @date 12.12.2022
 */
#include <vector>

#include "gf.hh"
#ifndef RS_H
#define RS_H

#define GENERATOR 3
#define FCR 1

class BMS_RS {
   public:
    BMS_RS();
    std::vector<int> rs_correct_msg(const std::vector<int>& msg_in, int nsym);
    std::vector<int> rs_encode_msg(const std::vector<int>& msg_in, int nsym);

   private:
    // initialize space for galois fields
    int GF_field_exp[256];
    int GF_field_log[256];

    // galois fields functions
    BMS_GF bms_gf;

    // Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
    std::vector<int> rs_generator_poly(int nsym);

    std::vector<int> rs_calc_syndromes(const std::vector<int>& msg, int nsym);

    bool rs_check(const std::vector<int>& msg, int nsym);

    std::vector<int> rs_find_errata_locator(const std::vector<int>& e_pos);

    std::vector<int> rs_find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym);

    std::vector<int> rs_find_errors(std::vector<int>& err_loc, int nmess);

    std::vector<int> rs_find_error_locator(const std::vector<int>& synd, int nsym, int erase_count);

    std::vector<int> rs_forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess);

    std::vector<int> rs_correct_errata(std::vector<int>& msg_in, std::vector<int>& synd, std::vector<int>& err_pos);
};

#endif