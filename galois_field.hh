/**
 * @file galois-field.hh
 * @author Michal Findra, xfindr00
 * @brief Galois field operations
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
#ifndef GF_H
#define GF_H

class GALOIS_FIELD {
   public:
    // constructor
    GALOIS_FIELD();

    /**
     * @brief Multiplication of two operands in Galois fields
     *
     * @param x Multiplicand
     * @param y Multiplier
     * @return int
     */
    int multiplication(int x, int y);

    /**
     * @brief Raising x to power in Galois fields
     *
     * @param x Base
     * @param power Exponent
     * @return int
     */
    int pow(unsigned int x, int power);

    /**
     * @brief Polynomial multiplication of two vectors in Galois fields
     *
     * @param p First vector
     * @param q Second vector
     * @return std::vector<int>
     */
    std::vector<int> polynomial_multiplication(std::vector<int> p, std::vector<int> q);

    /**
     * @brief Division in Galois fields
     *
     * @param x Dividend
     * @param y Divisor
     * @return int
     */
    int division(int x, int y);

    /**
     * @brief Subtraction in Galois fields
     *
     * @param x Minuend
     * @param y Subtrahend
     * @return int
     */
    int subtraction(int x, int y);

    /**
     * @brief Polynomial addition in Galois fields
     *
     * @param p Augend
     * @param q Addend
     * @return std::vector<int>
     */
    std::vector<int> polynomial_addition(std::vector<int> p, std::vector<int> q);

    /**
     * @brief Element inversion In Galois fields
     *
     * @param x base
     * @return int
     */
    int inverse(int x);

    /**
     * @brief Polynomial division in Galois fields
     *
     * @param x Dividend
     * @param y Divisor
     * @return std::vector<int>
     */
    std::vector<int> polynomial_division(const std::vector<int>& x, const std::vector<int>& y);

    /**
     * @brief Polynomial scale in Galois fields
     *
     * @param p base
     * @param x scaler
     * @return std::vector<int>
     */
    std::vector<int> polynomial_scale(const std::vector<int>& p, int x);

    /**
     * @brief Evaluate polynomial in Galois fields
     *
     * @param poly base
     * @param x evaluator
     * @return int
     */
    int polynomial_evaluation(const std::vector<int>& poly, int x);

   private:
};

#endif