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
    // constructor
    BMS_RS();

    /**
     * @brief Correct message if it is possible to correct it.
     *
     * @param msg_in Original message + ecc
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> rs_correct_msg(const std::vector<int>& msg_in, int nsym);

    /**
     * @brief Encoding message using reed-salomon alg. with generator 3, fcr 1
     *
     * @param msg_in Message to be encoded
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> rs_encode_msg(const std::vector<int>& msg_in, int nsym);

   private:
    // initialize space for galois fields
    int GF_field_exp[256];
    int GF_field_log[256];

    // galois fields functions
    BMS_GF bms_gf;

    /**
     * @brief Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
     *
     * @param nsym  Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> rs_generator_poly(int nsym);

    /**
     * @brief The syndromes polynomial is calculated from the given codeword and the number of error-correcting symbols.
     * This can be thought of as performing a Fourier Transform, with the Chien search being the inverse operation.
     *
     * @param msg Message to be checked
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> rs_calc_syndromes(const std::vector<int>& msg, int nsym);

    /**
     * @brief Checks message syndrome.
     *
     * @param msg Message to be checked
     * @param nsym Number of error correcting symbols
     * @return true If message is correct
     * @return false Otherwise
     */
    bool rs_check(const std::vector<int>& msg, int nsym);

    /**
     * @brief Computes locator polynomial from position of errors
     *
     * @param e_pos Error positions in code word
     * @return std::vector<int>
     */
    std::vector<int> rs_find_errata_locator(const std::vector<int>& e_pos);

    /**
     * @brief Computes error evaluator polynomial
     *
     * @param synd Syndromes
     * @param err_loc Error locations
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> rs_find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym);

    /**
     * @brief By brute-force finds the root of error polynomial
     *
     * @param err_loc Error locations
     * @param nmess Message length
     * @return std::vector<int>
     */
    std::vector<int> rs_find_errors(std::vector<int>& err_loc, int nmess);

    /**
     * @brief The Berlekamp-Massey algorithm is used to find the error/errata locator and evaluator polynomials.
     *
     * @param synd Syndromes
     * @param nsym Number of error correcting symbols
     * @param erase_count Number of erasures
     * @return std::vector<int>
     */
    std::vector<int> rs_find_error_locator(const std::vector<int>& synd, int nsym, int erase_count);

    /**
     * @brief
     *
     * @param synd Computes Forney syndromes
     * @param pos Position
     * @param nmess Message length
     * @return std::vector<int>
     */
    std::vector<int> rs_forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess);

    /**
     * @brief Calculates error magnitude polynomial
     *
     * @param msg_in Message to be checked
     * @param synd Syndromes
     * @param err_pos List of the errors positions
     * @return std::vector<int>
     */
    std::vector<int> rs_correct_errata(std::vector<int>& msg_in, std::vector<int>& synd, std::vector<int>& err_pos);
};

#endif