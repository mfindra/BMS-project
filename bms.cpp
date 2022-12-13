/*
author: Michal Findra, xfindr00
project: BMS
description:
*/

// other libs
#include <bits/stdc++.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#define BLOCK_LENGTH 8
#define DEBUG 0

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

using namespace std;

template <class T>
T operator+(const T& l, T&& r) {
    T c{};
    c.reserve(l.size() + r.size());
    auto bi = std::back_inserter(c);
    std::copy(l.begin(), l.end(), bi);
    std::move(r.begin(), r.end(), bi);
    return c;
}

// print help message to standard output
void PrintHelp() {
    cout << "ENCRYPT AND DECRYPT MESSAGE USING REED-SOLOMON ALGORITHM - BMS PROJECT 2022" << endl;
    cout << "===========================================================================" << endl
         << endl;
    cout << "Descrition: Encrypt and decrypt message using Reed-Solomon algorithm as described." << endl;
    cout << "Arguments: -r               : file to transfer " << endl;
    cout << "           -s <IP|Hostname> : destination IP address or hostname " << endl;
    cout << "           -l               : runs as server, which listens for incoming ICMP" << endl;
    cout << "                              messages and stores them in current directory" << endl;
    cout << "           -h               : print help" << endl;
    cout << endl;
    cout << "Example usage: " << endl
         << endl;
    cout << "Sending file \"example_file.txt\" to address 192.168.0.1 :" << endl;
    cout << "       server: sudo ./secret -r example_file.txt -s 192.168.0.1" << endl
         << "       reciever: sudo ./secret - l" << endl;
}

void showvector(vector<int> g) {
    vector<int>::iterator it;
    cout << "[";
    for (it = g.begin(); it != g.end(); ++it)
        if (it != g.begin())
            cout << ", " << *it;
        else
            cout << " " << *it;
    cout << "]\n";
}

std::vector<int> split_and_convert(const std::string& s) {
    std::vector<int> res;
    for (size_t i = 0; i < s.size(); i += 8) {
        std::string part = s.substr(i, 8);
        uint8_t val = (uint8_t)std::stoul(part, nullptr, 2);
        res.push_back(val);
    }
    return res;
}

int gf_mul(int x, int y) {
    if (x == 0 || y == 0)
        return 0;
    return GF_field_exp[((GF_field_log[x] + GF_field_log[y]) % 255)];  // should be gf_exp[(gf_log[x]+gf_log[y])%255] if gf_exp wasn't oversized
}

int gf_pow(unsigned int x, int power) {
    if ((GF_field_log[x] * power) < 0) {
        return GF_field_exp[(255 + (GF_field_log[x] * power)) % 255];
    } else
        return GF_field_exp[(GF_field_log[x] * power) % 255];
}

vector<int> gf_poly_mul(vector<int> p, vector<int> q) {
    showvector(p);
    showvector(q);
    // Multiply two polynomials, inside Galois Field
    //  Pre-allocate the result array
    std::vector<int> r(p.size() + q.size() - 1);
    for (int j = 0; j < q.size(); ++j) {
        for (int i = 0; i < p.size(); ++i) {
            cout << "gf_mul: " << gf_mul(p[i], q[j]) << endl;
            r[i + j] ^= gf_mul(p[i], q[j]);  // equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))
                                             // -- you can see it's your usual polynomial multiplication
        }
    }
    return r;
}

int gf_div(int x, int y) {
    if (y == 0) {
        exit;
    }
    if (x == 0) {
        return 0;
    }
    return GF_field_exp[(GF_field_log[x] + 255 - GF_field_log[y]) % 255];
}

int gf_sub(int x, int y) {
    return x ^ y;
}

std::vector<int> gf_poly_add(std::vector<int> p, std::vector<int> q) {
    std::vector<int> r(std::max(p.size(), q.size()), 0);
    for (std::size_t i = 0; i < p.size(); ++i) {
        r[i + r.size() - p.size()] = p[i];
    }
    for (std::size_t i = 0; i < q.size(); ++i) {
        r[i + r.size() - q.size()] ^= q[i];
    }
    return r;
}

int gf_inverse(int x) {
    return GF_field_exp[255 - GF_field_log[x]];
}

std::vector<int> gf_poly_div(const std::vector<int>& dividend, const std::vector<int>& divisor) {
    std::vector<int> msg_out = dividend;  // Copy the dividend
    for (int i = 0; i < dividend.size() - (divisor.size() - 1); i++) {
        int coef = msg_out[i];  // precaching
        if (coef != 0)          // log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization).
        {
            for (int j = 1; j < divisor.size(); j++)  // in synthetic division, we always skip the first coefficient of the divisor,
                                                      // because it's only used to normalize the dividend coefficient
            {
                if (divisor[j] != 0)  // log(0) is undefined
                {
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
    cout << "gf_poly_div: ";
    showvector(b);
    cout << "divident: ";
    showvector(dividend);
    cout << "divisor: ";
    showvector(divisor);
    return b;  // return quotient, remainder
}

vector<int> rs_generator_poly(int nsym) {
    // Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
    vector<int> g;
    g.push_back(1);
    for (int i = 1; i <= nsym; i++) {
        vector<int> tmp;
        tmp.push_back(1);
        tmp.push_back(gf_pow(3, i));
        g = gf_poly_mul(g, tmp);
    }
    if (DEBUG) {
        cout << "generator poly: ";
        showvector(g);
    }
    return g;
}

std::vector<int> rs_encode_msg(const std::vector<int>& msg_in, int nsym) {
    if (msg_in.size() + nsym > 255) {
        exit;
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
                msg_out[i + j] ^= gf_mul(gen[j], coef);
            }
        }
    }
    return msg_out;
}

vector<int> gf_poly_scale(const vector<int>& p, int x) {
    vector<int> r(p.size());
    for (int i = 0; i < p.size(); i++) {
        r[i] = gf_mul(p[i], x);
        cout << ">" << r[i] << " " << p[i] << " " << x << "<" << endl;
    }
    cout << "--" << endl;

    return r;
}

int gf_poly_eval(const std::vector<int>& poly, int x) {
    int y = poly[0];
    for (int i = 1; i < poly.size(); ++i) {
        y = gf_mul(y, x) ^ poly[i];
    }
    return y;
}

std::vector<int> rs_calc_syndromes(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd(nsym);
    for (int i = 0; i < nsym; ++i) {
        synd[i] = gf_poly_eval(msg, gf_pow(3, i + 1));
    }

    // Note the "[0] +" : we add a 0 coefficient for the lowest degree (the constant). This effectively shifts the syndrome, and will shift every computations depending on the syndromes (such as the errors locator polynomial, errors evaluator polynomial, etc. but not the errors positions).
    // This is not necessary, you can adapt subsequent computations to start from 0 instead of skipping the first iteration (ie, the often seen range(1, n-k+1)),
    synd.insert(synd.begin(), 0);

    cout << "rs_calc_syndromes: ";
    showvector(synd);
    return synd;
}

bool rs_check(const std::vector<int>& msg, int nsym) {
    std::vector<int> synd = rs_calc_syndromes(msg, nsym);
    auto it = std::max_element(synd.begin(), synd.end());
    return (*it == 0);
}

std::vector<int> rs_find_errata_locator(const std::vector<int>& e_pos) {
    // Compute the erasures/errors/errata locator polynomial from the erasures/errors/errata positions
    // (the positions must be relative to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to "h_ll_ worldxxxxxxxxx"
    // with xxxxxxxxx being the ecc of length n-k=9, here the string positions are [1, 4], but the coefficients are reversed
    // since the ecc characters are placed as the first coefficients of the polynomial, thus the coefficients of the
    // erased characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified as an argument.

    std::vector<int> e_loc{1};  // just to init because we will multiply, so it must be 1 so that the multiplication starts correctly without nulling any term
    // erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials.
    for (const int i : e_pos) {
        e_loc = gf_poly_mul(e_loc, gf_poly_add({1}, {gf_pow(3, i), 0}));
    }

    return e_loc;
}

std::vector<int> rs_find_error_evaluator(const std::vector<int>& synd, const std::vector<int>& err_loc, int nsym) {
    // Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    cout << "===" << endl;

    std::vector<int> v{1};
    std::vector<int> tmp(nsym + 1, 0);
    v.insert(v.end(), tmp.begin(), tmp.end());

    std::vector<int> remainder = gf_poly_div(gf_poly_mul(synd, err_loc), v);
    cout << "reminder: ";
    showvector(remainder);
    cout << "===" << endl;
    return remainder;
}

std::vector<int> rs_find_errors(std::vector<int>& err_loc, int nmess) {
    int errs = err_loc.size() - 1;
    std::vector<int> err_pos;
    for (int i = 0; i < nmess; ++i) {
        if (gf_poly_eval(err_loc, gf_pow(3, i)) == 0) {
            err_pos.push_back(nmess - 1 - i);
        }
    }
    // Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
    if (err_pos.size() != errs) {
        // couldn't find error locations
        exit;
    }
    return err_pos;
}

std::vector<int> rs_find_error_locator(const std::vector<int>& synd, int nsym, int erase_count) {
    cout << "rs_find_error_args: ";
    showvector(synd);
    cout << nsym << "," << erase_count << endl;
    std::vector<int> err_loc{1};
    std::vector<int> old_loc{1};

    int synd_shift = synd.size() - nsym;
    cout << "synd " << synd_shift << endl;
    for (int i = 0; i < nsym - erase_count; i++) {
        int K = i + synd_shift;  //(erase_loc.empty()) ? (i + synd_shift) : (erase_count + i + synd_shift);
        int delta = synd[K];
        showvector(err_loc);
        for (int j = 1; j < err_loc.size(); j++) {
            delta ^= gf_mul(err_loc.at(err_loc.size() - (j + 1)), synd[K - j]);
        }

        cout << "K: " << delta << endl;
        old_loc.push_back(0);

        if (delta != 0) {
            if (old_loc.size() > err_loc.size()) {
                std::vector<int> new_loc = gf_poly_scale(old_loc, delta);
                old_loc = gf_poly_scale(err_loc, gf_inverse(delta));
                err_loc = new_loc;
            }
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta));
            cout << "old_loc: ";
            showvector(old_loc);
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

    cout << "err_loc: ";
    showvector(err_loc);

    err_loc.erase(err_loc.begin(), err_loc.begin() + count);

    int errs = err_loc.size() - 1;
    if ((errs)*2 > nsym) {
        exit;  // too many errors to correct
    }
    return err_loc;
}

std::vector<int> rs_forney_syndromes(const std::vector<int>& synd, const std::vector<int>& pos, int nmess) {
    std::vector<int> erase_pos_reversed;
    for (auto p : pos) {
        erase_pos_reversed.push_back(nmess - 1 - p);
    }

    std::vector<int> fsynd(synd.begin() + 1, synd.end());

    for (int i = 0; i < pos.size(); i++) {
        int x = gf_pow(2, erase_pos_reversed[i]);
        for (int j = 0; j < fsynd.size() - 1; j++) {
            fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1];
        }
    }

    return fsynd;
}

vector<int> rs_correct_errata(vector<int>& msg_in, vector<int>& synd, vector<int>& err_pos) {
    vector<int> coef_pos;
    for (int p : err_pos) {
        coef_pos.push_back(msg_in.size() - 1 - p);
    }

    vector<int> err_loc = rs_find_errata_locator(coef_pos);
    reverse(synd.begin(), synd.end());
    vector<int> err_eval = rs_find_error_evaluator(synd, err_loc, err_loc.size() - 1);

    showvector(coef_pos);
    vector<int> X;
    for (int i = 0; i < coef_pos.size(); i++) {
        int l = 255 - coef_pos[i];
        X.push_back(gf_pow(3, -l));
    }

    cout << "X:";
    showvector(X);

    vector<int> E(msg_in.size(), 0);
    int Xlength = X.size();
    for (int i = 0; i < X.size(); i++) {
        int Xi_inv = gf_inverse(X[i]);
        cout << "Xi_inv: " << X[i] << endl;

        vector<int> err_loc_prime_tmp;
        for (int j = 0; j < Xlength; j++) {
            if (j != i) {
                err_loc_prime_tmp.push_back(gf_sub(1, gf_mul(Xi_inv, X[j])));
            }
        }
        cout << "err_loc_prime_tmp:";
        showvector(err_loc_prime_tmp);
        int err_loc_prime = 1;
        for (auto coef : err_loc_prime_tmp) {
            err_loc_prime = gf_mul(err_loc_prime, coef);
        }
        cout << "err_loc_prime: " << err_loc_prime << endl;

        cout << "err_eval";
        showvector(err_eval);

        auto y = gf_poly_eval(err_eval, Xi_inv);
        cout << "y: " << y << endl;

        y = gf_mul(gf_pow(X[i], 1 - 1), y);

        cout << "y: " << y << endl;

        if (err_loc_prime == 0) {
            exit;
        }

        auto magnitude = gf_div(y, err_loc_prime);
        E[err_pos[i]] = magnitude;
    }
    msg_in = gf_poly_add(msg_in, E);
    cout << "msg_errata: ";
    showvector(msg_in);
    return msg_in;
}

vector<int> rs_correct_msg(const vector<int>& msg_in, int nsym) {
    if (msg_in.size() > 255) {
        exit(4);
    }

    vector<int> erase_pos;
    vector<int> msg_out = msg_in;

    if (erase_pos.size() > nsym) {
        exit(5);
    }

    auto synd = rs_calc_syndromes(msg_out, nsym);

    if (*std::max_element(synd.begin(), synd.end()) == 0) {
        cout << "message correct" << endl;
        // return {msg_out.begin(), msg_out.end() - nsym}, {msg_out.end() - nsym, msg_out.end()};
        return msg_out;
    }

    auto fsynd = rs_forney_syndromes(synd, erase_pos, msg_out.size());
    cout << "fsynd: ";
    showvector(fsynd);

    auto err_loc = rs_find_error_locator(fsynd, nsym, erase_pos.size());
    cout << "err_loc: ";
    showvector(err_loc);

    reverse(err_loc.begin(), err_loc.end());
    auto err_pos = rs_find_errors(err_loc, msg_out.size());
    cout << "err_pos_: ";
    showvector(err_pos);

    if (err_pos.empty()) {
        exit(6);
    }

    // std::cout << "msg_out: " << msg_out << " synd: " << synd <<
    vector<int> tmp = erase_pos;
    tmp.insert(tmp.end(), err_pos.begin(), err_pos.end());
    msg_out = rs_correct_errata(msg_out, synd, tmp);

    showvector(msg_out);

    if (!rs_check(msg_out, nsym))
        exit(7);

    return msg_out;
}

int main(int argc, char** argv) {
    int n_opt;  // coded message length
    int k_opt;  // message length
    string t_opt;
    string m_opt;  // coded message
    bool encrypt_switch = false;
    bool decrypt_switch = false;
    bool L_opt = false;
    int opt;

    // read and parse program arguments
    while ((opt = getopt(argc, argv, "edn:t:k:m:h")) != -1) {
        switch (opt) {
            case 'e':
                encrypt_switch = true;
                break;
            case 'd':
                decrypt_switch = true;
                break;
            case 'n':
                n_opt = stoi(optarg);
                break;
            case 't':
                t_opt = optarg;
                break;
            case 'k':
                k_opt = stoi(optarg);
                break;
            case 'm':
                m_opt = optarg;
                break;
            case 'h':
                PrintHelp();
                return EXIT_SUCCESS;
                break;
            default:
                cerr << "ERROR - Wrong argument!" << endl;
                return EXIT_FAILURE;
        }
    }

    // check application mode
    if (encrypt_switch) {
        // convert the string to a vector of ordinal values
        std::vector<int> msg_in;
        msg_in.assign(t_opt.begin(), t_opt.end());

        auto a = rs_encode_msg(msg_in, n_opt - msg_in.size());

        if (DEBUG) {
            cout << "message in: ";
            showvector(msg_in);
        }

        a.erase(a.begin(), a.begin() + msg_in.size());
        for (int x : a) msg_in.push_back(x);

        if (DEBUG) {
            cout << "message out: ";
            showvector(msg_in);
        }

        for (int x : msg_in) {
            std::bitset<8> bits(x);         // create a bitset with 8 bits, initialized with the value of x
            std::cout << bits.to_string();  // print the binary representation of the bitset
        }
        cout << endl;

    } else if (decrypt_switch) {
        cout << "Decrypting" << endl;
        std::vector<int> res;

        for (size_t i = 0; i < m_opt.size(); i += 8) {
            std::string part = m_opt.substr(i, 8);
            int val = (int)std::stoul(part, nullptr, 2);
            res.push_back(val);
        }
        showvector(res);
        if (rs_check(res, n_opt - k_opt))
            cout << "yeah" << endl;
        else
            cout << "nay" << endl;

        /* vector<int> tmp{33, 98, 99, 120, 95};
        vector<int> tmo{218, 44, 0};
        vector<int> tmi{0};
        rs_correct_errata(tmp, tmo, tmi); */

        cout << "original message:  ";
        showvector(res);
        vector<int> msg_out;
        msg_out = rs_correct_msg(res, n_opt - k_opt);
        showvector(msg_out);
        for (int i = 0; i < k_opt; i++) {
            cout << char(msg_out[i]);
        }
        cout << endl;
    }

    return EXIT_SUCCESS;
}

/*
0110000101100010011000110111100001011111

011000100110110101110011001100000110110001111100

0100100001100101011011000110110001101111001011000010000001110111011011110111001001101100011001000010000110001101000100111111010011111001010000110001000011100101

0100001001100101011110100110010001110010011000010111010001101111011101100110010100100000011000010010000001101101011011110110001001101001011011000110111001101001001000000111001101101001011101000110010101000001100011101110100101101110011010100001110010011110

*/
