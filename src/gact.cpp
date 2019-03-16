/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without ion, including without limitation the rights
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

#include "gact.h"
#include <math.h>
#include "xcl2.hpp"
#include <vector>


enum states {Z, D, I, M};
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
extern int forward_kernel_counter;
extern int reverse_kernel_counter;



int dequeue(int *BT_states, int *front, int *queuesize) {
   int data = BT_states[(*front)++];
	
   if(*front == MAX_TILE_SIZE) {
      *front = 0;
   }
	
   (*queuesize)--;
   return data;  
}


std::queue<int> CpuXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states, int *queuesize, char strand, int *rear, int *front);

#ifdef FPGA
std::queue<int> FpgaXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states, int *queuesize, char strand, int *rear, int *front, int mode);
#endif

Alignment GACT (char* ref_str, char* query_str, std::string ref_name, std::string query_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold, int mode) {
    std::queue<int> BT_states_std;

    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";

    Alignment alignment;
    alignment.ref_name = ref_name;
    alignment.query_name = query_name;
    alignment.aligned_ref_str = "";
    alignment.aligned_query_str = "";
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;
    alignment.aligned_ref_len = 0;
    alignment.aligned_query_len = 0;
    alignment.ref_len = ref_length;
    alignment.query_len = query_length;
    alignment.strand = strand;
    alignment.score = 0;
    alignment.flag = 0;
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
   
    int max_ref_pos = 0;
    int max_query_pos = 0;
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;

    while ((ref_pos > 0) && (query_pos > 0) && (((i > 0) && (j > 0)) || first_tile)) {
        //change the tile length if elements less than that of the tile size
        ref_tile_length = (ref_pos > tile_size) ? tile_size : ref_pos;
        query_tile_length = (query_pos > tile_size) ? tile_size : query_pos;

    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode==1) { 
       	       BT_states_std =  FpgaXL (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
        }
#ifdef FPGA
        else if(mode==2) { 
       	       BT_states_std = FpgaXL (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
        }
#endif
        else {
        	BT_states_std = AlignWithBT (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, 
				first_tile, (tile_size - tile_overlap));
	}

        i = 0;
        j = 0;
        int tile_score = BT_states_std.front();
        BT_states_std.pop();

        
        if (first_tile) {
            ref_pos = ref_pos - ref_tile_length + BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos - query_tile_length + BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
            rev_ref_pos = ref_pos;
            rev_query_pos = query_pos;
            first_tile_score = tile_score;
        }

        int num_tb = BT_states_std.size();
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = num_tb-1;
        int query_buf_curr = num_tb-1;

        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr--] = '-';
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = std::string(ref_buf, num_tb) + aligned_ref_str;
            aligned_query_str = std::string(query_buf, num_tb) + aligned_query_str;
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos -= (j);
        query_pos -= (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }
    
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;

    ref_pos = rev_ref_pos;
    query_pos = rev_query_pos;
    
    i =  tile_size;
    j = tile_size;
    
    //starts with the first tile
    while ((ref_pos < ref_length) && (query_pos < query_length) && (((i > 0) && (j > 0)) || first_tile)) {
        ref_tile_length = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
        query_tile_length = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode == 1) { 
                BT_states_std = FpgaXL(ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
        }
#ifdef FPGA
        else if(mode == 2) { 
                BT_states_std = FpgaXL (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
	}
#endif
        else {
        	BT_states_std = AlignWithBT (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, 
				first_tile, (tile_size - tile_overlap));
	}
        i = 0;
        j = 0;
        int tile_score = BT_states_std.front(); BT_states_std.pop(); 

        if (first_tile) {
            ref_pos = ref_pos + ref_tile_length - BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos + query_tile_length - BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
        }

        int num_tb = BT_states_std.size();;
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = 0;
        int query_buf_curr = 0;
        
        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr++] = '-';
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = aligned_ref_str + std::string(ref_buf, num_tb);
            aligned_query_str = aligned_query_str + std::string(query_buf, num_tb);
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos += (j);
        query_pos += (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }

    int total_score = 0;
    bool open = true;
    for (uint32_t j = 0; j < aligned_ref_str.length(); j++) {
        char ref_nt = aligned_ref_str[j];
        char query_nt = aligned_query_str[j];
        if (ref_nt == '-' || query_nt == '-') {
            total_score += (open) ? gap_open : gap_extend;
            open = false;
        }
        else {
            total_score += sub_mat[5*NtChar2Int(query_nt) + NtChar2Int(ref_nt)];
            open = true;
        }
    }
    alignment.aligned_ref_str = aligned_ref_str;
    alignment.aligned_query_str = aligned_query_str;
    alignment.score = total_score;
    return alignment;
}

extern void xl_cpu (const char *  a, const int m,
                                        const char *  b, const int n,
                                        const int gap_open, const int gap_extend,
                                        const int query_pos, const int ref_pos, const int reverse, const int first, const int early_terminate,
                                        int *  sub_mat,  int *  dir_matrix_arg,   int *  max_score_arg,
                                        int *  max_i_arg,  int *  max_j_arg,  int *  pos_score_arg);


std::queue<int> CpuXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate, int* BT_states_param, int *queuesize, char strand, int *rear, int *front) {
    char *a, *b;
    a = ref_str;
    b = query_str;
    int reverse_int = reverse;
    int first_tile_int = first;
    int m = (int)ref_tile_length;
    int n = (int)query_tile_length;

    std::vector<int,aligned_allocator<int>> dir_matrix((MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1));
    std::vector<int,aligned_allocator<int>> max_score(1), max_i(1), max_j(1), pos_score(1);

    for (int i = 0; i < m+1; i++)
           dir_matrix[i*(MAX_TILE_SIZE+1)+0] = ZERO_OP;
    for (int j = 0; j < n + 1; j++)
           dir_matrix[0*(MAX_TILE_SIZE+1)+j] = ZERO_OP;

    std::vector<char,aligned_allocator<char>> a_aligned(m);
    std::vector<char,aligned_allocator<char>> b_aligned(n);
    std::vector<int,aligned_allocator<int>> sub_mat_aligned(25);

    for(int i = 0 ; i < m ; i++)
        a_aligned[i] = a[i];
    for(int i = 0 ; i < n ; i++)
        b_aligned[i] = b[i];
    for(int i = 0 ; i < 25 ; i++)
        sub_mat_aligned[i] = sub_mat[i];

   xl_cpu( a_aligned.data(), m, b_aligned.data(), n, gap_open, gap_extend,  query_pos, ref_pos, 
		reverse_int, first_tile_int, early_terminate, sub_mat_aligned.data(), 
		dir_matrix.data(), max_score.data(), max_i.data(), max_j.data(), pos_score.data());

  if(reverse)
    reverse_kernel_counter++;
  else
    forward_kernel_counter++;

  std::queue<int> BT_states;

  int i_curr=ref_pos, j_curr=query_pos;
  int i_steps = 0, j_steps = 0;


  int open = 0;
  if (first) {
      i_curr = *(max_i.data());
      j_curr = *(max_j.data());
      BT_states.push(*(max_score.data()));
      BT_states.push(i_curr);
      BT_states.push(j_curr);
  }
  else {
      BT_states.push(*(pos_score.data()));
  }

  int state = dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] % 4;

  while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) {
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[(i_curr-1)*(MAX_TILE_SIZE+1)+(j_curr-1)] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
  };

  return BT_states;
}


#ifdef FPGA

std::queue<int> FpgaXL(char* a, long long int m, char* b, long long int n, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first, int early_terminate, int* BT_states_arg, int *queuesize, char strand, int *rear, int *front, int mode) {

    cl_int status;
    int reverse_int = reverse;
    int first_tile_int = first;

    size_t wgSize[3] = {1, 1, 1};
    size_t gSize[3] = {1, 1, 1};

    int d = 0;
    int nbb=1, jj=0;
    int score;
    int m_i = (int) m;
    int n_i = (int)n;

    std::vector<int,aligned_allocator<int>> dir_matrix((MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1));
    std::vector<int,aligned_allocator<int>> max_score(1), max_i(1), max_j(1), pos_score(1);
    std::vector<int,aligned_allocator<int>> sub_mat_aligned(25);

    for (int i = 0; i < m+1; i++)
           dir_matrix[i*(MAX_TILE_SIZE+1)+0] = ZERO_OP;
    for (int j = 0; j < n + 1; j++)
           dir_matrix[0*(MAX_TILE_SIZE+1)+j] = ZERO_OP;

    std::vector<char,aligned_allocator<char>> a_aligned(m);
    std::vector<char,aligned_allocator<char>> b_aligned(n);

    for(int i = 0 ; i < m ; i++)
        a_aligned[i] = a[i];
    for(int i = 0 ; i < n ; i++)
        b_aligned[i] = b[i];
    for(int i = 0 ; i < 25 ; i++)
        sub_mat_aligned[i] = sub_mat[i];
    
   if(mode == 1) {
   xl_cpu( a_aligned.data(), m_i, b_aligned.data(), n_i, gap_open, gap_extend,  query_pos, ref_pos,
                reverse_int, first_tile_int, early_terminate, sub_mat_aligned.data(),
                dir_matrix.data(), max_score.data(), max_i.data(), max_j.data(), pos_score.data());
   } else {

    cl::Buffer cl_a(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, m, a_aligned.data());
    cl::Buffer cl_b(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, n, b_aligned.data());
    cl::Buffer cl_sub_mat(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 25*sizeof(int), sub_mat_aligned.data());
    cl::Buffer cl_dir_matrix(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1), dir_matrix.data());
    cl::Buffer cl_max_score(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int), max_score.data());
    cl::Buffer cl_max_i(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int), max_i.data());
    cl::Buffer cl_max_j(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int), max_j.data());
    cl::Buffer cl_pos_score(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int), pos_score.data());

    std::vector<cl::Memory> inBufVec, outBufVec;
    inBufVec.push_back(cl_a);
    inBufVec.push_back(cl_b);
    inBufVec.push_back(cl_sub_mat);
    outBufVec.push_back(cl_dir_matrix);
    outBufVec.push_back(cl_max_score);
    outBufVec.push_back(cl_max_i);
    outBufVec.push_back(cl_max_j);
    outBufVec.push_back(cl_pos_score);

    cl::Kernel krnl_xl(program, "xl");
    q.enqueueMigrateMemObjects(inBufVec, 0/* 0 means from host*/);

    int narg = 0;
    krnl_xl.setArg(narg++, cl_a);
    krnl_xl.setArg(narg++, m_i);
    krnl_xl.setArg(narg++, cl_b);
    krnl_xl.setArg(narg++, n_i);
    krnl_xl.setArg(narg++, gap_open);
    krnl_xl.setArg(narg++, gap_extend);
    krnl_xl.setArg(narg++, query_pos);
    krnl_xl.setArg(narg++, ref_pos);
    krnl_xl.setArg(narg++, reverse_int);
    krnl_xl.setArg(narg++, first_tile_int);
    krnl_xl.setArg(narg++, early_terminate);
    krnl_xl.setArg(narg++, cl_sub_mat);
    krnl_xl.setArg(narg++, cl_dir_matrix);
    krnl_xl.setArg(narg++, cl_max_score);
    krnl_xl.setArg(narg++, cl_max_i);
    krnl_xl.setArg(narg++, cl_max_j);
    krnl_xl.setArg(narg++, cl_pos_score);

    krnl_xl.setArg(narg++, jj);

    //Launch the Kernel
    q.enqueueTask(krnl_xl);
    q.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
    q.finish();
  }

  if(reverse)
    reverse_kernel_counter++;
  else
    forward_kernel_counter++;

  std::queue<int> BT_states;

  int i_curr=ref_pos, j_curr=query_pos;
  int i_steps = 0, j_steps = 0;


  int open = 0;
  if (first) {
      i_curr = *(max_i.data());
      j_curr = *(max_j.data());
      BT_states.push(*(max_score.data()));
      BT_states.push(i_curr);
      BT_states.push(j_curr);
  }
  else {
      BT_states.push(*(pos_score.data()));
  }

  int state = dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] % 4;

  while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) {
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[(i_curr-1)*(MAX_TILE_SIZE+1)+(j_curr-1)] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr*(MAX_TILE_SIZE+1)+j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
  };

  return BT_states;

}


#endif

