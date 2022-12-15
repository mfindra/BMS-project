/**
 * @file rs.cc
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

// other libs
#include "reed_salomon.hh"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// Constructor
REED_SOLOMON::REED_SOLOMON() {
}

// Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
std::vector<int> REED_SOLOMON::generator_poly(int nsym) {
    std::vector<int> g;
    g.push_back(1);
    for (int i = 1; i <= nsym; i++) {
        std::vector<int> tmp;
        tmp.push_back(1);
        tmp.push_back(gf.pow(3, i));
        g = gf.polynomial_multiplication(g, tmp);
    }
    return g;
}

// Encoding message using reed-salomon alg. with generator 3, fcr 1
std::vector<int> REED_SOLOMON::encode_msg(const std::vector<int>& msg_in, int nsym) {
    if (msg_in.size() + nsym > 255) {
        std::cerr << "ERROR: message too long" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> gen = generator_poly(nsym);
    std::vector<int> msg_out(msg_in.size() + gen.size() - 1);

    // initialize the first part of the output message with the input message
    for (int i = 0; i < msg_in.size(); i++) {
        msg_out[i] = msg_in[i];
    }

    // Use galois filed operations to calculate rest of the encoded message
    for (int i = 0; i < msg_in.size(); i++) {
        int coef = msg_out[i];

        if (coef != 0) {
            for (int j = 1; j < gen.size(); j++) {
                msg_out[i + j] ^= gf.multiplication(gen[j], coef);
            }
        }
    }
    return msg_out;
}

// Syndrome polynomial calculation for given code word
std::vector<int> REED_SOLOMON::calc_syndromes(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd(nsym);
    // Equivalent to fourier transformation
    for (int i = 0; i < nsym; ++i) {
        synd[i] = gf.polynomial_evaluation(msg, gf.pow(3, i + 1));
    }

    // Add 0 coefficient for the the constant, error position wont be shifted.
    synd.insert(synd.begin(), 0);

    return synd;
}

// Checks message syndrome
bool REED_SOLOMON::check(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd = calc_syndromes(msg, nsym);
    // checks if all elements are 0 -> no errors
    auto it = std::max_element(synd.begin(), synd.end());
    return (*it == 0);
}

// Computes locator polynomial from position of errors
std::vector<int> REED_SOLOMON::find_errata_locator(const std::vector<int>& e_pos) {
    // init t o1 because of multiplication
    std::vector<int> e_loc{1};
    // erasures_loc = product(1 - x*alpha**i) for each in erasures_pos
    for (const int i : e_pos) {
        e_loc = gf.polynomial_multiplication(e_loc, gf.polynomial_addition({1}, {gf.pow(3, i), 0}));
    }

    return e_loc;
}

// Computes error evaluator polynomial
std::vector<int> REED_SOLOMON::find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym) {
    // Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    std::vector<int> v{1};
    std::vector<int> tmp(nsym + 1, 0);
    v.insert(v.end(), tmp.begin(), tmp.end());

    std::vector<int> remainder = gf.polynomial_division(gf.polynomial_multiplication(synd, err_loc), v);
    return remainder;
}

// By brute-force finds the root of error polynomial
std::vector<int> REED_SOLOMON::find_errors(std::vector<int>& err_loc, int nmess) {
    int errs = err_loc.size() - 1;
    std::vector<int> err_pos;
    // Chien search
    for (int i = 0; i < nmess; ++i) {
        if (gf.polynomial_evaluation(err_loc, gf.pow(3, i)) == 0) {
            err_pos.push_back(nmess - 1 - i);
        }
    }
    // number of errors = length of the errata locator polynomial
    if (err_pos.size() != errs) {
        // couldn't find error locations
        std::cerr << "ERROR: too many errors" << std::endl;
        exit(EXIT_FAILURE);
    }
    return err_pos;
}
// The Berlekamp-Massey algorithm is used to find the error/errata locator and evaluator polynomials.
std::vector<int> REED_SOLOMON::find_error_locator(const std::vector<int>& synd, int nsym, int erase_count) {
    std::vector<int> err_loc{1};
    std::vector<int> old_loc{1};
    // Berlekamp-Massey algorithm
    int synd_shift = synd.size() - nsym;
    for (int i = 0; i < nsym - erase_count; i++) {
        int K = i + synd_shift;  //(erase_loc.empty()) ? (i + synd_shift) : (erase_count + i + synd_shift);
        int delta = synd[K];
        for (int j = 1; j < err_loc.size(); j++) {
            delta ^= gf.multiplication(err_loc.at(err_loc.size() - (j + 1)), synd[K - j]);
        }

        old_loc.push_back(0);
        // Errata locator and evaluator estimation
        if (delta != 0) {
            if (old_loc.size() > err_loc.size()) {
                std::vector<int> new_loc = gf.polynomial_scale(old_loc, delta);
                old_loc = gf.polynomial_scale(err_loc, gf.inverse(delta));
                err_loc = new_loc;
            }
            err_loc = gf.polynomial_addition(err_loc, gf.polynomial_scale(old_loc, delta));
        }
    }

    // check if there are not too many errors to correct
    int count = 0;
    for (int i = 0; i < err_loc.size(); i++) {
        if (err_loc[i] == 0) {
            count++;
        } else {
            break;
        }
    }

    err_loc.erase(err_loc.begin(), err_loc.begin() + count);

    int errs = err_loc.size() - 1;
    if ((errs)*2 > nsym) {
        std::cerr << "ERROR: too many errors" << std::endl;
        exit(EXIT_FAILURE);
    }
    return err_loc;
}

// Computes Forney syndromes
std::vector<int> REED_SOLOMON::forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess) {
    std::vector<int> erase_pos_reversed;
    // Preparation of coefficient degree positions
    for (auto p : pos) {
        erase_pos_reversed.push_back(nmess - 1 - p);
    }

    // Copy
    std::vector<int> fsynd(synd.begin() + 1, synd.end());

    // Forney syndrome calculation
    for (int i = 0; i < pos.size(); i++) {
        int x = gf.pow(2, erase_pos_reversed[i]);
        for (int j = 0; j < fsynd.size() - 1; j++) {
            fsynd[j] = gf.multiplication(fsynd[j], x) ^ fsynd[j + 1];
        }
    }

    return fsynd;
}

// Calculates error magnitude polynomial
std::vector<int> REED_SOLOMON::correct_errata(std::vector<int>& msg_in, std::vector<int>& synd, std::vector<int>& err_pos) {
    // Calculate errata locator
    std::vector<int> coef_pos;
    for (int p : err_pos) {
        coef_pos.push_back(msg_in.size() - 1 - p);
    }

    // Errata evaluator polynomial calculation
    std::vector<int> err_loc = find_errata_locator(coef_pos);
    reverse(synd.begin(), synd.end());
    std::vector<int> err_eval = find_error_evaluator(synd, err_loc, err_loc.size() - 1);

    std::vector<int> X;
    for (int i = 0; i < coef_pos.size(); i++) {
        int l = 255 - coef_pos[i];
        X.push_back(gf.pow(3, -l));
    }

    // Forney alg
    std::vector<int> E(msg_in.size(), 0);
    int Xlength = X.size();
    for (int i = 0; i < X.size(); i++) {
        int Xi_inv = gf.inverse(X[i]);

        // Calculate formal derivate
        std::vector<int> err_loc_prime_tmp;
        for (int j = 0; j < Xlength; j++) {
            if (j != i) {
                err_loc_prime_tmp.push_back(gf.subtraction(1, gf.multiplication(Xi_inv, X[j])));
            }
        }
        int err_loc_prime = 1;

        // Calculate denominator of the Forney algorithm
        for (auto coef : err_loc_prime_tmp) {
            err_loc_prime = gf.multiplication(err_loc_prime, coef);
        }

        auto y = gf.polynomial_evaluation(err_eval, Xi_inv);
        y = gf.multiplication(gf.pow(X[i], 1 - 1), y);

        if (err_loc_prime == 0) {
            std::cerr << "ERROR: could not find magnitude" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Compute magnitude
        auto magnitude = gf.division(y, err_loc_prime);
        E[err_pos[i]] = magnitude;
    }

    // Apply value correction
    msg_in = gf.polynomial_addition(msg_in, E);
    return msg_in;
}

// Correct message if it is possible to correct it
std::vector<int> REED_SOLOMON::correct_msg(const std::vector<int>& msg_in, int nsym) {
    // check message leng
    if (msg_in.size() > 255) {
        std::cerr << "ERROR: message too long" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> erase_pos;
    std::vector<int> msg_out = msg_in;

    // Calculate syndrome polynomial
    auto synd = calc_syndromes(msg_out, nsym);

    // If no errors return message
    if (*std::max_element(synd.begin(), synd.end()) == 0) {
        // return {msg_out.begin(), msg_out.end() - nsym}, {msg_out.end() - nsym, msg_out.end()};
        return msg_out;
    }

    // Calculate Forney syndromes
    auto fsynd = forney_syndromes(synd, erase_pos, msg_out.size());

    // Compute error locator polynomial
    auto err_loc = find_error_locator(fsynd, nsym, erase_pos.size());

    // Locate errors in message
    reverse(err_loc.begin(), err_loc.end());
    auto err_pos = find_errors(err_loc, msg_out.size());

    if (err_pos.empty()) {
        std::cerr << "ERROR: error location failed" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Correct errors
    std::vector<int> tmp = erase_pos;
    tmp.insert(tmp.end(), err_pos.begin(), err_pos.end());
    msg_out = correct_errata(msg_out, synd, tmp);

    // Final check if message was corrected
    if (!check(msg_out, nsym)) {
        std::cerr << "ERROR: message uncorrectable" << std::endl;
        exit(EXIT_FAILURE);
    }

    return msg_out;
}
