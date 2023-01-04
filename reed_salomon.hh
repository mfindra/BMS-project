/**
 * @file reed-salomon.hh
 * @author Michal Findra, xfindr00
 * @brief Reed-Salomon code operations
 * @date 12.12.2022
 * Resources used:
 * - Andrew Brown, Stephen Larroque - unireedsolomon
    (https://pypi.org/project/unireedsolomon/)
    MIT licensed (https://mit-license.org)
 * - Wikiversity - Reed–Solomon codes for coders
    (https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)
    licensed under CC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)
 */
#include <vector>

#include "galois_field.hh"
#ifndef RS_H
#define RS_H

#define GENERATOR 3
#define FCR 1

class REED_SOLOMON {
   public:
    // constructor
    REED_SOLOMON();

    /**
     * @brief Correct message if it is possible to correct it.
     *
     * @param msg_in Original message + ecc
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> correct_msg(const std::vector<int>& msg_in, int nsym);

    /**
     * @brief Encoding message using reed-salomon alg. with generator 3, fcr 1
     *
     * @param msg_in Message to be encoded
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> encode_msg(const std::vector<int>& msg_in, int nsym);

   private:
    // initialize space for galois fields
    int galois_field_exponential[256];
    int galois_field_logarithmical[256];

    // galois fields functions
    GALOIS_FIELD gf;

    /**
     * @brief Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
     *
     * @param nsym  Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> generator_poly(int nsym);

    /**
     * @brief The syndromes polynomial is calculated from the given codeword and the number of error-correcting symbols.
     * This can be thought of as performing a Fourier Transform, with the Chien search being the inverse operation.
     *
     * @param msg Message to be checked
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> calc_syndromes(const std::vector<int>& msg, int nsym);

    /**
     * @brief Checks message syndrome
     *
     * @param msg Message to be checked
     * @param nsym Number of error correcting symbols
     * @return true If message is correct
     * @return false Otherwise
     */
    bool check(const std::vector<int>& msg, int nsym);

    /**
     * @brief Computes locator polynomial from position of errors
     *
     * @param e_pos Error positions in code word
     * @return std::vector<int>
     */
    std::vector<int> find_errata_locator(const std::vector<int>& e_pos);

    /**
     * @brief Computes error evaluator polynomial
     *
     * @param synd Syndromes
     * @param err_loc Error locations
     * @param nsym Number of error correcting symbols
     * @return std::vector<int>
     */
    std::vector<int> find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym);

    /**
     * @brief By brute-force finds the root of error polynomial
     *
     * @param err_loc Error locations
     * @param nmess Message length
     * @return std::vector<int>
     */
    std::vector<int> find_errors(std::vector<int>& err_loc, int nmess);

    /**
     * @brief The Berlekamp-Massey algorithm is used to find the error/errata locator and evaluator polynomials.
     *
     * @param synd Syndromes
     * @param nsym Number of error correcting symbols
     * @param erase_count Number of erasures
     * @return std::vector<int>
     */
    std::vector<int> find_error_locator(const std::vector<int>& synd, int nsym, int erase_count);

    /**
     * @brief Computes Forney syndromes
     *
     * @param synd Syndromes
     * @param pos Position
     * @param nmess Message length
     * @return std::vector<int>
     */
    std::vector<int> forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess);

    /**
     * @brief Calculates error magnitude polynomial
     *
     * @param msg_in Message to be checked
     * @param synd Syndromes
     * @param err_pos List of the errors positions
     * @return std::vector<int>
     */
    std::vector<int> correct_errata(std::vector<int>& msg_in, std::vector<int>& synd, std::vector<int>& err_pos);
};

#endif