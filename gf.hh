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
    // constructor
    BMS_GF();

    /**
     * @brief Multiplication of two operands in Galois fields
     *
     * @param x Multiplicand
     * @param y Multiplier
     * @return int
     */
    int gf_mul(int x, int y);

    /**
     * @brief Raising x to power in Galois fields
     *
     * @param x Base
     * @param power Exponent
     * @return int
     */
    int gf_pow(unsigned int x, int power);

    /**
     * @brief Polynomial multiplication of two vectors in Galois fields
     *
     * @param p First vector
     * @param q Second vector
     * @return std::vector<int>
     */
    std::vector<int> gf_poly_mul(std::vector<int> p, std::vector<int> q);

    /**
     * @brief Division in Galois fields
     *
     * @param x Dividend
     * @param y Divisor
     * @return int
     */
    int gf_div(int x, int y);

    /**
     * @brief Subtraction in Galois fields
     *
     * @param x Minuend
     * @param y Subtrahend
     * @return int
     */
    int gf_sub(int x, int y);

    /**
     * @brief Polynomial addition
     *
     * @param p Augend
     * @param q Addend
     * @return std::vector<int>
     */
    std::vector<int> gf_poly_add(std::vector<int> p, std::vector<int> q);

    /**
     * @brief Element inversion In Galois fields
     *
     * @param x base
     * @return int
     */
    int gf_inverse(int x);

    /**
     * @brief Polynomial division in Galois fields
     *
     * @param x Dividend
     * @param y Divisor
     * @return std::vector<int>
     */
    std::vector<int> gf_poly_div(const std::vector<int>& x, const std::vector<int>& y);

    /**
     * @brief Polynomial scale in Galois fields
     *
     * @param p base
     * @param x scaler
     * @return std::vector<int>
     */
    std::vector<int> gf_poly_scale(const std::vector<int>& p, int x);

    /**
     * @brief Evaluate polynomial in Galois fields
     *
     * @param poly base
     * @param x evaluator
     * @return int
     */
    int gf_poly_eval(const std::vector<int>& poly, int x);

   private:
};

#endif