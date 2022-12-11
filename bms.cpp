/*
author: Michal Findra, xfindr00
project: BMS
description:
*/

// other libs
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <bitset>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#define BLOCK_LENGTH 8

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

int gf_mul(int x, int y) {
    if (x == 0 || y == 0)
        return 0;
    // cout << "gf mul: " << GF_field_exp[GF_field_log[x] + GF_field_log[y]] << endl;
    // cout << tmp << endl;
    // cout << "x: " << GF_field_log[x] << "y: " << GF_field_log[y] << "result: " << GF_field_exp[((GF_field_log[x] + GF_field_log[y]) % 255)] << endl;
    return GF_field_exp[((GF_field_log[x] + GF_field_log[y]) % 255)];  // should be gf_exp[(gf_log[x]+gf_log[y])%255] if gf_exp wasn't oversized
}

int gf_pow(int x, int power) {
    cout << ">gf pow: " << GF_field_exp[(GF_field_log[x] * power) % 255] << " with x=" << x << ", power=" << power << endl;
    return GF_field_exp[(GF_field_log[x] * power) % 255];
}

vector<int> gf_poly_mul(vector<int> p, vector<int> q) {
    // Multiply two polynomials, inside Galois Field
    //  Pre-allocate the result array
    vector<int> r;
    vector<int> rr;
    for (int i = 0; i < p.size() + q.size() - 1; i++)
        r.push_back(0);

    // cout << "before bef: ";
    // showvector(r);
    //  Compute the polynomial multiplication (just like the outer product of two vectors,
    //  we multiply each coefficients of p with all coefficients of q)
    for (int j = 0; j < q.size(); j++) {
        for (int i = 0; i < p.size(); i++) {
            auto it = r.begin();
            vector<int>::iterator itp = p.begin();
            vector<int>::iterator itq = q.begin();
            advance(it, i + j);
            advance(itp, i);
            advance(itq, j);
            int itx = *it;
            int itpx = *itp;
            int itqx = *itq;
            // cout << "before===> " << itx << endl;
            //  cout << "item in p at: " << i << " is: " << itpx << endl;
            //  cout << "item in q at: " << j << " is: " << itqx << endl;
            int tmp;
            // cout << "before===> ";
            // showvector(r);

            // auto l = r.erase(it);
            // cout << ">>" << *l << endl;
            tmp = gf_mul(itpx, itqx);
            // cout << "xor(" << itx << ", " << tmp << ") = " << (itx ^ gf_mul(itpx, itqx)) << endl;
            auto k = r.erase(it);
            r.insert(k, itx ^ gf_mul(itpx, itqx));  // equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))
            // cout << "after:  ";
            // cout << "after===> ";
            // showvector(r);
        }
    }
    return r;
}

vector<int> gf_poly_div(const std::vector<int>& dividend, const std::vector<int>& divisor) {
    vector<int> msg_out = dividend;  // copy the dividend

    for (int i = 0; i < dividend.size() - (divisor.size() - 1); ++i) {
        int coef = msg_out[i];  // precaching
        if (coef != 0)          // log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization).
        {
            for (int j = 1; j < divisor.size(); ++j)  // in synthetic division, we always skip the first coefficient of the divisor,
                                                      // because it's only used to normalize the dividend coefficient
            {
                if (divisor[j] != 0)  // log(0) is undefined
                {
                    // equivalent to the more mathematically correct
                    // (but xoring directly is faster): msg_out[i + j] += -divisor[j] * coef
                    msg_out[i + j] ^= gf_mul(divisor[j], coef);
                }
            }
        }
    }

    return msg_out;
}

vector<int> rs_generator_poly(int nsym) {
    // Generate an irreducible generator polynomial (necessary to encode a message into Reed-Solomon)
    vector<int> g;
    g.push_back(1);
    for (int i = 1; i <= nsym; i++) {
        vector<int> tmp;
        tmp.push_back(1);
        tmp.push_back(gf_pow(3, i));
        // cout << "tmp: ";
        // showvector(tmp);
        g = gf_poly_mul(g, tmp);
        // cout << "g after: ";
        // showvector(g);
        // cout << "rs gener: ";
        // showvector(g);
        // cout << "===================" << endl;
    }
    cout << "generator poly: ";
    showvector(g);
    return g;
}

std::vector<int> rs_encode_msg(const std::vector<int>& msg_in, int nsym) {
    if (msg_in.size() + nsym > 255) {
        throw std::invalid_argument("Message is too long");
    }

    cout << nsym << endl;
    std::vector<int> gen = rs_generator_poly(nsym);
    std::vector<int> msg_out(msg_in.size() + gen.size() - 1);

    // initialize the first part of the output message with the input message
    for (int i = 0; i < msg_in.size(); i++) {
        msg_out[i] = msg_in[i];
    }
    cout << "msg_out: ";
    showvector(msg_out);

    cout << "msg_in: ";
    showvector(msg_in);

    for (int i = 0; i < msg_in.size(); i++) {
        int coef = msg_out[i];

        if (coef != 0) {
            for (int j = 1; j < gen.size(); j++) {
                // cout << endl;
                // cout << "msg_out[i,j]: " << msg_out[i + j] << endl;
                // cout << "second: " << gf_mul(gen[j], coef) << endl;
                // cout << "genj: " << gen[j] << "coef: " << coef << endl;
                msg_out[i + j] ^= gf_mul(gen[j], coef);
                // cout << "i: " << i << "j: " << j << endl;
                // showvector(msg_out);
            }
        }
    }

    showvector(msg_out);
    return msg_out;
}

int main(int argc, char** argv) {
    int n_opt;  // coded message length
    string t_opt;
    string k_opt;  // message length
    string m_opt;  // coded message
    bool encrypt_switch = false;
    bool decrypt_switch = false;
    bool L_opt = false;
    int opt;

    bool debug = false;

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
                k_opt = optarg;
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
        showvector(msg_in);

        a.erase(a.begin(), a.begin() + msg_in.size());

        for (int x : a) msg_in.push_back(x);

        showvector(msg_in);
        for (int x : msg_in) {
            std::bitset<8> bits(x);         // create a bitset with 8 bits, initialized with the value of x
            std::cout << bits.to_string();  // print the binary representation of the bitset
        }

    } else if (decrypt_switch) {
    }

    return EXIT_SUCCESS;
}

// 01100001 01100010 011000110 11110000 1011111
// 01100001 01100010 011000110 11110000 1011111

/*
011000100110110101110011001100000110110001111100
011000100110110101110011001100000110110001111100
011000100110110101110011001100000110110001111100

0100100001100101011011000110110001101111001011000010000001110111011011110111001001101100011001000010000110001101000100111111010011111001010000110001000011100101
0100100001100101011011000110110001101111001011000010000001110111011011110111001001101100011001000010000110001101000100111111010011111001010000110001000011100101
0100100001100101011011000110110001101111001011000010000001110111011011110111001001101100011001000010000110001101000100111111010011111001010000110001000011100101

0100001001100101011110100110010001110010011000010111010001101111011101100110010100100000011000010010000001101101011011110110001001101001011011000110111001101001001000000111001101101001011101000110010101000001100011101110100101101110011010100001110010011110
0100001001100101011110100110010001110010011000010111010001101111011101100110010100100000011000010010000001101101011011110110001001101001011011000110111001101001001000000111001101101001011101000110010101000001100011101110100101101110011010100001110010011110
0100001001100101011110100110010001110010011000010111010001101111011101100110010100100000011000010010000001101101011011110110001001101001011011000110111001101001001000000111001101101001011101000110010101000001100011101110100101101110011010100001110010011110

*/
