/**
 * @file rs.cc
 * @author Michal Findra, xfindr00
 * @brief Reed-Salomon code operations
 * @date 12.12.2022
 */

// other libs
#include "rs.hh"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

BMS_RS::BMS_RS() {
}

std::vector<int> BMS_RS::rs_generator_poly(int nsym) {
    // Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
    std::vector<int> g;
    g.push_back(1);
    for (int i = 1; i <= nsym; i++) {
        std::vector<int> tmp;
        tmp.push_back(1);
        tmp.push_back(bms_gf.gf_pow(3, i));
        g = bms_gf.gf_poly_mul(g, tmp);
    }
    return g;
}

std::vector<int> BMS_RS::rs_encode_msg(const std::vector<int>& msg_in, int nsym) {
    if (msg_in.size() + nsym > 255) {
        std::cerr << "ERROR: message too long" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> gen = rs_generator_poly(nsym);
    std::vector<int> msg_out(msg_in.size() + gen.size() - 1);

    // initialize the first part of the output message with the input message
    for (int i = 0; i < msg_in.size(); i++) {
        msg_out[i] = msg_in[i];
    }

    for (int i = 0; i < msg_in.size(); i++) {
        int coef = msg_out[i];

        if (coef != 0) {
            for (int j = 1; j < gen.size(); j++) {
                msg_out[i + j] ^= bms_gf.gf_mul(gen[j], coef);
            }
        }
    }
    return msg_out;
}

std::vector<int> BMS_RS::rs_calc_syndromes(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd(nsym);
    for (int i = 0; i < nsym; ++i) {
        synd[i] = bms_gf.gf_poly_eval(msg, bms_gf.gf_pow(3, i + 1));
    }

    // Note the "[0] +" : we add a 0 coefficient for the lowest degree (the constant). This effectively shifts the syndrome, and will shift every computations depending on the syndromes (such as the errors locator polynomial, errors evaluator polynomial, etc. but not the errors positions).
    // This is not necessary, you can adapt subsequent computations to start from 0 instead of skipping the first iteration (ie, the often seen range(1, n-k+1)),
    synd.insert(synd.begin(), 0);

    return synd;
}

bool BMS_RS::rs_check(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd = rs_calc_syndromes(msg, nsym);
    auto it = std::max_element(synd.begin(), synd.end());
    return (*it == 0);
}

std::vector<int> BMS_RS::rs_find_errata_locator(const std::vector<int>& e_pos) {
    // Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions
    // (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx"
    // with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed
    // since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the
    // erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.

    std::vector<int> e_loc{1};  // just to init because we will multiply, so it must be 1 so that the multiplication starts correctly without nulling any term
    // erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials.
    for (const int i : e_pos) {
        e_loc = bms_gf.gf_poly_mul(e_loc, bms_gf.gf_poly_add({1}, {bms_gf.gf_pow(3, i), 0}));
    }

    return e_loc;
}

std::vector<int> BMS_RS::rs_find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym) {
    // Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)

    std::vector<int> v{1};
    std::vector<int> tmp(nsym + 1, 0);
    v.insert(v.end(), tmp.begin(), tmp.end());

    std::vector<int> remainder = bms_gf.gf_poly_div(bms_gf.gf_poly_mul(synd, err_loc), v);
    return remainder;
}

std::vector<int> BMS_RS::rs_find_errors(std::vector<int>& err_loc, int nmess) {
    int errs = err_loc.size() - 1;
    std::vector<int> err_pos;
    for (int i = 0; i < nmess; ++i) {
        if (bms_gf.gf_poly_eval(err_loc, bms_gf.gf_pow(3, i)) == 0) {
            err_pos.push_back(nmess - 1 - i);
        }
    }
    // Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
    if (err_pos.size() != errs) {
        // couldn't find error locations
        std::cerr << "ERROR: too many errors" << std::endl;
        exit(EXIT_FAILURE);
    }
    return err_pos;
}

std::vector<int> BMS_RS::rs_find_error_locator(const std::vector<int>& synd, int nsym, int erase_count) {
    std::vector<int> err_loc{1};
    std::vector<int> old_loc{1};

    int synd_shift = synd.size() - nsym;
    for (int i = 0; i < nsym - erase_count; i++) {
        int K = i + synd_shift;  //(erase_loc.empty()) ? (i + synd_shift) : (erase_count + i + synd_shift);
        int delta = synd[K];
        for (int j = 1; j < err_loc.size(); j++) {
            delta ^= bms_gf.gf_mul(err_loc.at(err_loc.size() - (j + 1)), synd[K - j]);
        }

        old_loc.push_back(0);

        if (delta != 0) {
            if (old_loc.size() > err_loc.size()) {
                std::vector<int> new_loc = bms_gf.gf_poly_scale(old_loc, delta);
                old_loc = bms_gf.gf_poly_scale(err_loc, bms_gf.gf_inverse(delta));
                err_loc = new_loc;
            }
            err_loc = bms_gf.gf_poly_add(err_loc, bms_gf.gf_poly_scale(old_loc, delta));
        }
    }

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

std::vector<int> BMS_RS::rs_forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess) {
    std::vector<int> erase_pos_reversed;
    for (auto p : pos) {
        erase_pos_reversed.push_back(nmess - 1 - p);
    }

    std::vector<int> fsynd(synd.begin() + 1, synd.end());

    for (int i = 0; i < pos.size(); i++) {
        int x = bms_gf.gf_pow(2, erase_pos_reversed[i]);
        for (int j = 0; j < fsynd.size() - 1; j++) {
            fsynd[j] = bms_gf.gf_mul(fsynd[j], x) ^ fsynd[j + 1];
        }
    }

    return fsynd;
}

std::vector<int> BMS_RS::rs_correct_errata(std::vector<int>& msg_in, std::vector<int>& synd, std::vector<int>& err_pos) {
    std::vector<int> coef_pos;
    for (int p : err_pos) {
        coef_pos.push_back(msg_in.size() - 1 - p);
    }

    std::vector<int> err_loc = rs_find_errata_locator(coef_pos);
    reverse(synd.begin(), synd.end());
    std::vector<int> err_eval = rs_find_error_evaluator(synd, err_loc, err_loc.size() - 1);

    std::vector<int> X;
    for (int i = 0; i < coef_pos.size(); i++) {
        int l = 255 - coef_pos[i];
        X.push_back(bms_gf.gf_pow(3, -l));
    }

    std::vector<int> E(msg_in.size(), 0);
    int Xlength = X.size();
    for (int i = 0; i < X.size(); i++) {
        int Xi_inv = bms_gf.gf_inverse(X[i]);

        std::vector<int> err_loc_prime_tmp;
        for (int j = 0; j < Xlength; j++) {
            if (j != i) {
                err_loc_prime_tmp.push_back(bms_gf.gf_sub(1, bms_gf.gf_mul(Xi_inv, X[j])));
            }
        }
        int err_loc_prime = 1;
        for (auto coef : err_loc_prime_tmp) {
            err_loc_prime = bms_gf.gf_mul(err_loc_prime, coef);
        }

        auto y = bms_gf.gf_poly_eval(err_eval, Xi_inv);
        y = bms_gf.gf_mul(bms_gf.gf_pow(X[i], 1 - 1), y);

        if (err_loc_prime == 0) {
            std::cerr << "ERROR: could not find magnitude" << std::endl;
            exit(EXIT_FAILURE);
        }

        auto magnitude = bms_gf.gf_div(y, err_loc_prime);
        E[err_pos[i]] = magnitude;
    }
    msg_in = bms_gf.gf_poly_add(msg_in, E);
    return msg_in;
}

std::vector<int> BMS_RS::rs_correct_msg(const std::vector<int>& msg_in, int nsym) {
    if (msg_in.size() > 255) {
        std::cerr << "ERROR: message too long" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> erase_pos;
    std::vector<int> msg_out = msg_in;

    if (erase_pos.size() > nsym) {
        std::cerr << "ERROR: too many errors" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto synd = rs_calc_syndromes(msg_out, nsym);

    if (*std::max_element(synd.begin(), synd.end()) == 0) {
        // return {msg_out.begin(), msg_out.end() - nsym}, {msg_out.end() - nsym, msg_out.end()};
        return msg_out;
    }

    auto fsynd = rs_forney_syndromes(synd, erase_pos, msg_out.size());
    auto err_loc = rs_find_error_locator(fsynd, nsym, erase_pos.size());
    reverse(err_loc.begin(), err_loc.end());
    auto err_pos = rs_find_errors(err_loc, msg_out.size());

    if (err_pos.empty()) {
        std::cerr << "ERROR: error location failed" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> tmp = erase_pos;
    tmp.insert(tmp.end(), err_pos.begin(), err_pos.end());
    msg_out = rs_correct_errata(msg_out, synd, tmp);

    if (!rs_check(msg_out, nsym)) {
        std::cerr << "ERROR: message uncorrectable" << std::endl;
        exit(EXIT_FAILURE);
    }

    return msg_out;
}

/*
0110000101100010011000110111100001011111

011000100110110101110011001100000110110001111100

0100100001100101011011000110110001101111001011000010000001110111011011110111001001101100011001000010000110001101000100111111010011111001010000110001000011100101

0100001001100101011110100110010001110010011000010111010001101111011101100110010100100000011000010010000001101101011011110110001001101001011011000110111001101001001000000111001101101001011101000110010101000001100011101110100101101110011010100001110010011110

*/
