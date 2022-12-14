/*
author: Michal Findra, xfindr00
project: BMS
*/
#include <vector>
#ifndef RS_H
#define RS_H

#define BLOCK_LENGTH 8

class BMS_RS {
   public:
    BMS_RS();
    ~BMS_RS();
    std::vector<int> rs_correct_msg(const std::vector<int>& msg_in, int nsym);
    std::vector<int> rs_encode_msg(const std::vector<int>& msg_in, int nsym);

   private:
    int GF_field_exp[256];
    int GF_field_log[256];
};

#endif