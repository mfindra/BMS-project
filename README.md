# BMS-project

Simple implementation of Reed-Solomon encoder and decoder using CPP11.

Values which are specific for this variation of R-S: 
- Generator: 3
- Galois field Prime: 0x11b
- Galois field Fcr: 1
- Galois field size: 256 (2^8)
- Max length of encoded word: 255

Galois field tables are precomputed. For changes in above values, fields need to be regenerated. 

## Usage

Before usage build app using: 

	make

Output is bms binary.

Print help: 

    ./bms -h

### Encoding message:

    ./bms -e -n <code word length> -t <input text>

Specific example: 

    ./bms -e -n 5 -t "abc"

Output:
 
    0110000101100010011000110111100001011111

### Decoding message:

    ./bms -d -n <code word length> -k <message length> -m <message binary>

Specific example: 

    ./bms -d -n 5 -k 3 -m 0110000101100010011000110111100001011111

Output:

    abc 
    

## Sources 

More info about Reed-Solomon: 

[Reed–Solomon codes for coders](https://en.wikiversity.org/wiki/Reed–Solomon_codes_for_coders) </br>
[Reed–Solomon error correction](https://en.wikipedia.org/wiki/Reed–Solomon_error_correction) </br>
[Pip library for R-S](https://pypi.org/project/unireedsolomon/) </br>

