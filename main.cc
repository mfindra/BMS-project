/**
 * @file main.cc
 * @author Michal Findra, xfindr00
 * @brief main function for Reed-Salomon encoding and decoding
 * @date 12.12.2022
 * Resources used:
 * - Andrew Brown, Stephen Larroque - unireedsolomon
    (https://pypi.org/project/unireedsolomon/)
    MIT licensed (https://mit-license.org)
 * - Wikiversity - Reed–Solomon codes for coders
    (https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)
    licensed under CC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/)
 */

#include <bits/stdc++.h>
#include <getopt.h>

#include <iostream>
#include <string>
#include <vector>

#include "reed_salomon.hh"

// print help message to standard output
void PrintHelp() {
    std::cout << "ENCRYPT AND DECRYPT MESSAGE USING REED-SOLOMON ALGORITHM - BMS PROJECT 2022" << std::endl;
    std::cout << "===========================================================================" << std::endl
              << std::endl;
    std::cout << "Descrition: Encrypt and decrypt message using Reed-Solomon algorithm as described." << std::endl;
    std::cout << "Default generator = 3 and fcr = 1. Galois field tables are precomputed." << std::endl
              << std::endl;
    std::cout << "Arguments: -e                 : encrypt mode " << std::endl;
    std::cout << "           -d                 : decrypt mode " << std::endl;
    std::cout << "           -n                 : code word length " << std::endl;
    std::cout << "           -t <string>        : input text for encoding" << std::endl;
    std::cout << "           -k                 : message length in decoding" << std::endl;
    std::cout << "           -n <binary string> : message in decoding" << std::endl;
    std::cout << "           -h                 : print help" << std::endl;
    std::cout << std::endl;
    std::cout << "Example usage: " << std::endl
              << std::endl;
    std::cout << "./bms -e -n <code word length> -t <input text>" << std::endl;
    std::cout << "./bms -d -n <code word length> -k <message length> -m <message binary>" << std::endl;
}

int main(int argc, char** argv) {
    int n_opt;  // coded message length
    int k_opt;  // message length
    std::string t_opt;
    std::string m_opt;  // coded message
    bool encrypt_switch = false;
    bool decrypt_switch = false;
    bool L_opt = false;
    int opt;
    REED_SOLOMON RS;

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
                n_opt = std::stoi(optarg);
                break;
            case 't':
                t_opt = optarg;
                break;
            case 'k':
                k_opt = std::stoi(optarg);
                break;
            case 'm':
                m_opt = optarg;
                break;
            case 'h':
                PrintHelp();
                return EXIT_SUCCESS;
                break;
            default:
                std::cerr << "ERROR - Wrong argument!" << std::endl;
                return EXIT_FAILURE;
        }
    }

    // check application mode
    if (encrypt_switch) {
        // convert the string to a vector of ordinal values
        std::vector<int> msg_in;
        msg_in.assign(t_opt.begin(), t_opt.end());

        // encode message
        auto a = RS.encode_msg(msg_in, n_opt - msg_in.size());
        a.erase(a.begin(), a.begin() + msg_in.size());
        for (int x : a) msg_in.push_back(x);

        for (int x : msg_in) {
            std::bitset<8> bits(x);         // create a bitset with 8 bits, initialized with the value of x
            std::cout << bits.to_string();  // print the binary representation of the bitset
        }
        std::cout << std::endl;

    } else if (decrypt_switch) {
        std::vector<int> res;
        for (size_t i = 0; i < m_opt.size(); i += 8) {
            std::string part = m_opt.substr(i, 8);
            int val = (int)std::stoul(part, nullptr, 2);
            res.push_back(val);
        }

        std::vector<int> msg_out;
        msg_out = RS.correct_msg(res, n_opt - k_opt);
        for (int i = 0; i < k_opt; i++) {
            std::cout << char(msg_out[i]);
        }
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}