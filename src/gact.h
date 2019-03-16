/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <queue>
#include "align.h"

#define STRING_BUFFER_LEN 1024
#define PRECOMPILED_BINARY "xl"
#define NUM_DEVICES 1
#define MAX_NUM_DEVICES 16

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1

#include <CL/cl2.hpp>
//#include <ap_int.h>
using std::vector;
#define LENGTH (2)
#define NUM_WORKGROUPS (1)
#define WORKGROUP_SIZE (1)
extern std::vector<cl::Platform> platforms;
extern cl::Platform platform;
extern std::vector<cl::Device> devices;
extern cl::Device device;
extern cl::Context context;
extern cl::CommandQueue q;
extern cl::Program::Binaries bins;
extern cl::Program program;
extern char* xclbinFilename;

extern int forward_kernel_counter;
extern int reverse_kernel_counter;


extern void checkErr(cl_int err, const char * name);

struct Alignment {
    std::string ref_name;
    std::string query_name;
    std::string aligned_ref_str;
    std::string aligned_query_str;
    uint32_t ref_start;
    uint32_t query_start;
    uint32_t aligned_ref_len;
    uint32_t aligned_query_len;
    uint32_t ref_len;
    uint32_t query_len;
    int score;
    int flag;
    char strand;
};

Alignment GACT (char* ref_str, char* query_str, std::string ref_name, std::string query_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold, int mode);

