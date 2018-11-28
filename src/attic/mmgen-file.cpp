

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class mm_generator {
public:
  mm_generator() {
    source = fopen("matmat.cpp", "w");
    header = fopen("matmat.hpp", "w");
    table  = fopen("mattable.hpp", "w");
    for (auto f : {header, source, table}) {
      (void)fprintf(f, "#include <cstddef>\n");
      (void)fprintf(f, "#include \"Matrix.hpp\"\n");
    }
    (void)fprintf(source, "#include <future>\n");

    (void)fprintf(table, "#include <string>\n");
    (void)fprintf(table, "#include <tuple>\n");
    (void)fprintf(table, "#include \"matmat.hpp\"\n");
    (void)fprintf(table, "std::tuple<std::string, void (*)(const Matrix&, const Matrix&, Matrix&)> matmat_table[] = {\n");

    for (auto f : {header, source, table}) {
      newline(f);
    }
  }

  ~mm_generator() {
    (void)fprintf(table, "};\n");

    fclose(source);
    fclose(header);
    fclose(table);
  }

private:
  FILE* source;
  FILE* header;
  FILE* table;

  struct loop {
    std::string name;
    std::string b_name;
    std::string begin;
    std::string end;
    std::string last;
    std::string size;
    size_t      the_size;
  };

  std::vector<size_t> inner_loop_perm = {0, 1, 2};
  std::vector<size_t> block_loop_perm = {0, 1, 2};

  loop inner_loop[3] = {{"i", "ib", "i_begin", "i_end", "i_last", "tile_rows", 2},
                        {"j", "jb", "j_begin", "j_end", "j_last", "tile_cols", 2},
                        {"k", "kb", "k_begin", "k_end", "k_last", "tile_inner", 2}};

  loop block_loop[3] = {{"ii", "iib", "ii_begin", "ii_end", "ii_last", "block_rows", 24},
                        {"jj", "jjb", "jj_begin", "jj_end", "jj_last", "block_cols", 32},
                        {"kk", "kkb", "kk_begin", "kk_end", "kk_last", "block_inner", 48}};

  bool transpose_B = true;
  bool hoist = true;
  size_t num_threads = 4;

  void newline(FILE *f) {
    fprintf(f, "\n");
  }
  void newline() {
    newline(source);
  }

  // For loops

  void for_begin_to_last_by_size(loop lp) {
    (void)fprintf(source, "for (");
    (void)fprintf(source, "size_t %s = %s, %s = 0; %s < %s; %s += %lu, %s += %lu",
		  lp.name.c_str(), lp.begin.c_str(), lp.b_name.c_str(), 
		  lp.name.c_str(), lp.last.c_str(),
                  lp.name.c_str(), lp.the_size, lp.b_name.c_str(), lp.the_size);
    (void)fprintf(source, ") {\n");
  }

  void for_begin_to_end_by_size(loop lp) {
    (void)fprintf(source, "for (");
    (void)fprintf(source, "size_t %s = %s, %s = 0; %s < %s; %s += %lu, %s += %lu",
		  lp.name.c_str(), lp.begin.c_str(), lp.b_name.c_str(), 
		  lp.name.c_str(), lp.end.c_str(),
                  lp.name.c_str(), lp.the_size, lp.b_name.c_str(), lp.the_size);
    (void)fprintf(source, ") {\n");
  }

  void for_end_to_last_by_one(loop lp) {
    (void)fprintf(source, "for (");
    (void)fprintf(source, "size_t %s = %s, %s = %s - %s; %s < %s; ++%s, ++%s", 
		  lp.name.c_str(), lp.end.c_str(), lp.b_name.c_str(), lp.end.c_str(), lp.begin.c_str(),
		  lp.name.c_str(), lp.last.c_str(),
                  lp.name.c_str(), lp.b_name.c_str());
    (void)fprintf(source, ") {\n");
  }

  // Loop setup

  void inner_limits(size_t idx) {
    (void)fprintf(source, "size_t %s = %s, %s_step = std::min<size_t>(%lu, %s-%s), %s = %s + %lu*(%s_step/%lu), %s = %s + %s_step;\n",
                  inner_loop[idx].begin.c_str(), block_loop[idx].name.c_str(), inner_loop[idx].name.c_str(), 
		  block_loop[idx].the_size, block_loop[idx].last.c_str(), block_loop[idx].name.c_str(), 
		  inner_loop[idx].end.c_str(), inner_loop[idx].begin.c_str(), inner_loop[idx].the_size, inner_loop[idx].name.c_str(), inner_loop[idx].the_size, 
		  inner_loop[idx].last.c_str(), inner_loop[idx].begin.c_str(), inner_loop[idx].name.c_str());
  }

  // Hoisting

  void hoist_pre() {
    size_t which_to_hoist = inner_loop_perm[2]; // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = (0 == which_to_hoist ? 1 : inner_loop[0].the_size);
    size_t jlimit = (1 == which_to_hoist ? 1 : inner_loop[1].the_size);
    size_t klimit = (2 == which_to_hoist ? 1 : inner_loop[2].the_size);
    
    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {
	  switch(which_to_hoist) {
	  case 2:
	    (void)fprintf(source, "double t%lu%lu = C(%s + %lu, %s + %lu);\n", i, j, inner_loop[0].name.c_str(), i, inner_loop[1].name.c_str(), j);
	    break;
	  case 1:
	    (void)fprintf(source, "double t%lu%lu = A(%s + %lu, %s + %lu);\n", i, k, inner_loop[0].name.c_str(), i, inner_loop[2].name.c_str(), k);
	    break;
	  case 0:
	    if (transpose_B) {
	      (void)fprintf(source, "double t%lu%lu = B_T(%s + %lu, %s + %lu);\n", j, k, inner_loop[1].b_name.c_str(), j, inner_loop[2].b_name.c_str(), k);
	    } else {
	      (void)fprintf(source, "double t%lu%lu = B(%s + %lu, %s + %lu);\n", k, j, inner_loop[2].name.c_str(), k, inner_loop[1].name.c_str(), j);
	    }
	    break;
	  default: throw;
	  }
	}
      }
    }
  }

  void hoist_post() {
    size_t which_to_hoist = inner_loop_perm[2]; // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = (0 == which_to_hoist ? 1 : inner_loop[0].the_size);
    size_t jlimit = (1 == which_to_hoist ? 1 : inner_loop[1].the_size);
    size_t klimit = (2 == which_to_hoist ? 1 : inner_loop[2].the_size);
    
    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {
	  switch(which_to_hoist) {
	  case 2:
	    (void)fprintf(source, "C(%s + %lu, %s + %lu) = t%lu%lu;\n", inner_loop[0].name.c_str(), i, inner_loop[1].name.c_str(), j, i, j);
	    break;
	  case 1:
	  case 0:
	    break;
	  default: throw;
	  }
	}
      }
    }
  }


  // Loop bodies

  void hoisted_loop_body() {
    size_t hoisted = inner_loop_perm[2]; // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = inner_loop[0].the_size;
    size_t jlimit = inner_loop[1].the_size;
    size_t klimit = inner_loop[2].the_size;

    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {
	  
	  switch (hoisted) {
	  case 2:
	    (void)fprintf(source, "t%lu%lu += A(%s + %lu, %s + %lu) * ", i, j, inner_loop[0].name.c_str(), i, inner_loop[2].name.c_str(), k);
	    if (transpose_B) {
	      (void)fprintf(source, "B_T(%s + %lu, %s + %lu);\n", inner_loop[1].b_name.c_str(), j, inner_loop[2].b_name.c_str(), k);
	    } else {
	      (void)fprintf(source, "B(%s + %lu, %s + %lu);\n", inner_loop[2].name.c_str(), k, inner_loop[1].name.c_str(), j);
	    }
	    break;
	  case 1:
	    (void)fprintf(source, "C(%s + %lu, %s + %lu) += t%lu%lu * ", inner_loop[0].name.c_str(), i, inner_loop[1].name.c_str(), j, i, k);
	    if (transpose_B) {
	      (void)fprintf(source, "B_T(%s + %lu, %s + %lu);\n", inner_loop[1].b_name.c_str(), j, inner_loop[2].b_name.c_str(), k);
	    } else {
	      (void)fprintf(source, "B(%s + %lu, %s + %lu);\n", inner_loop[2].name.c_str(), k, inner_loop[1].name.c_str(), j);
	    }
	    break;
	  case 0:
	    if (transpose_B) {
	      (void)fprintf(source, "C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * t%lu%lu;\n", 
			    inner_loop[0].name.c_str(), i, inner_loop[1].name.c_str(), j,
			    inner_loop[0].name.c_str(), i, inner_loop[2].name.c_str(), k, j, k);
	    } else {
	      (void)fprintf(source, "C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * t%lu%lu;\n", 
			    inner_loop[0].name.c_str(), i, inner_loop[1].name.c_str(), j,
			    inner_loop[0].name.c_str(), i, inner_loop[2].name.c_str(), k, k, j);
	    }
	    break;
	  default:
	    throw;
	  }
	}
      }
    }
  }

  void loop_body(bool outer_fringe, bool middle_fringe, bool inner_fringe) {
    size_t ilimit = inner_loop[0].the_size;
    size_t jlimit = inner_loop[1].the_size;
    size_t klimit = inner_loop[2].the_size;

    size_t start = 2;
    for (auto f : {inner_fringe, middle_fringe, outer_fringe}) {
      if (f) {
        switch (inner_loop_perm[start]) {
          case 0:
            ilimit = 1;
            break;
          case 1:
            jlimit = 1;
            break;
          case 2:
            klimit = 1;
            break;
          default:
            throw;
        }
      }
      --start;
    }

    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {

          (void)fprintf(source, "      C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * ", inner_loop[0].name.c_str(), i,
                        inner_loop[1].name.c_str(), j, inner_loop[0].name.c_str(), i, inner_loop[2].name.c_str(), k);

          if (transpose_B) {
	    (void)fprintf(source, "B_T(%s + %lu, %s + %lu);\n", inner_loop[1].b_name.c_str(), j, inner_loop[2].b_name.c_str(), k);
          } else {
            (void)fprintf(source, "B(%s + %lu, %s + %lu);\n", inner_loop[2].name.c_str(), k, inner_loop[1].name.c_str(), j);
          }
        }
      }
    }
  }

  void inner_main_begin(bool outer_fringe, bool middle_fringe) {
    size_t p2 = inner_loop_perm[2];
    for_begin_to_end_by_size(inner_loop[p2]);
  }
  void inner_main_end() { (void)fprintf(source, "}\n"); }

  void inner_fringe_begin(bool outer_fringe, bool middle_fringe) {
    size_t p2 = inner_loop_perm[2];
    for_end_to_last_by_one(inner_loop[p2]);
  }
  void inner_fringe_end() { (void)fprintf(source, "}\n"); }

  void middle_main_begin(bool outer_fringe) {
    size_t p1 = inner_loop_perm[1];
    for_begin_to_end_by_size(inner_loop[p1]);
  }
  void middle_main_end() { (void)fprintf(source, "}\n"); }

  void middle_fringe_begin(bool outer_fringe) {
    size_t p1 = inner_loop_perm[1];
    for_end_to_last_by_one(inner_loop[p1]);
  }
  void middle_fringe_end() { (void)fprintf(source, "}\n"); }

  void outer_main_begin() {
    size_t p0 = inner_loop_perm[0];
    for_begin_to_end_by_size(inner_loop[p0]);
  }
  void outer_main_end() { (void)fprintf(source, "}\n"); }

  void outer_fringe_begin() {
    size_t p0 = inner_loop_perm[0];
    for_end_to_last_by_one(inner_loop[p0]);
  }
  void outer_fringe_end() { (void)fprintf(source, "}\n"); }


  // The inner loop
  void inner(bool outer_fringe, bool middle_fringe) {

    if (hoist && !outer_fringe && !middle_fringe) {
      hoist_pre();
      inner_main_begin(outer_fringe, middle_fringe);
      hoisted_loop_body();
      inner_main_end();
      hoist_post();
    } else {
      inner_main_begin(outer_fringe, middle_fringe);
      loop_body(outer_fringe, middle_fringe, false);
      inner_main_end();
    }

    inner_fringe_begin(outer_fringe, middle_fringe);
    loop_body(outer_fringe, middle_fringe, true);
    inner_fringe_end();
  }


  // The middle loop
  void middle(bool outer_fringe) {
    middle_main_begin(outer_fringe);
    inner(outer_fringe, false);
    middle_main_end();
    middle_fringe_begin(outer_fringe);
    inner(outer_fringe, true);
    middle_fringe_end();
  }

  // The outer loop (register level)
  void outer() {

    (void)fprintf(source, "\n");

    inner_limits(0);
    inner_limits(1);
    inner_limits(2);
    newline();
    

    if (transpose_B) {
      (void)fprintf(source, "Matrix B_T(j_step, k_step);\n");
      (void)fprintf(source, "for (size_t jb = 0, jt = j_begin; jb < j_step; ++jb, ++jt) {\n");
      (void)fprintf(source, "for (size_t kb = 0, kt = k_begin; kb < k_step; ++kb, ++kt) {\n");
      (void)fprintf(source, "B_T(jb, kb) = B(kt, jt);\n");
      (void)fprintf(source, "}\n}\n");
      newline();
    }

    outer_main_begin();
    middle(false);
    outer_main_end();
    outer_fringe_begin();
    middle(true);
    outer_fringe_end();
  }

  void thread_prelude() {
    if (0 != num_threads) {
      (void) fprintf(source, "size_t                         thread_count = 0;\n");
      (void) fprintf(source, "size_t                         num_threads  = %lu;\n", num_threads);
      (void) fprintf(source, "std::vector<std::future<void>> futs(num_threads);\n");
      newline();
    }
  }

  void begin_thread() {
    if (0 != num_threads) {    
      newline();
      (void) fprintf(source, "futs[thread_count++] = std::async(std::launch::async, [&, ii, jj]() -> void {\n");
      newline();
    }
  }

  void end_thread()  {
    if (0 != num_threads) {    
      (void) fprintf(source, "});\n");
      newline();
      (void) fprintf(source, "if (thread_count == num_threads) {\n");
      cleanup_thread();
      (void) fprintf(source, "}\n");
      newline();
    }
  }

  void cleanup_thread() {
    (void) fprintf(source, "for (size_t i = 0; i < thread_count; ++i) { // wait on launched tasks\n");
    (void) fprintf(source, "futs[i].wait();\n");
    (void) fprintf(source, "}\n");
    (void) fprintf(source, "thread_count = 0;\n");
  }


  void end_for(loop lp) {
    (void) fprintf(source, "} // End %s\n", lp.name.c_str());
  }

  // The blocked outer loop
  void block_outer() {

    (void)fprintf(source, "\n");
    (void)fprintf(source, "// block limits \n");
    (void)fprintf(source, "size_t ii_begin = 0, ii_last = C.num_rows();\n");
    (void)fprintf(source, "size_t jj_begin = 0, jj_last = C.num_cols();\n");
    (void)fprintf(source, "size_t kk_begin = 0, kk_last = A.num_cols();\n");
    (void)fprintf(source, "\n");

    thread_prelude();

    for_begin_to_last_by_size(block_loop[block_loop_perm[0]]);
    for_begin_to_last_by_size(block_loop[block_loop_perm[1]]);

    begin_thread();

    for_begin_to_last_by_size(block_loop[block_loop_perm[2]]);

    outer();

    end_for(block_loop[block_loop_perm[2]]);

    end_thread();

    end_for(block_loop[block_loop_perm[1]]);
    end_for(block_loop[block_loop_perm[0]]);
    
    cleanup_thread();
    
    (void)fprintf(source, "\n\n");
  }

public:

  // Generate a matmat function
  void gen(std::vector<size_t> bl_perm, size_t ib, size_t jb, size_t kb, bool tr, std::vector<size_t> il_perm, size_t it,
	   size_t jt, size_t kt, bool ht, size_t nthr = 0) {

    block_loop_perm = bl_perm;
    inner_loop_perm = il_perm;

    block_loop[0].the_size = ib;
    block_loop[1].the_size = jb;
    block_loop[2].the_size = kb;

    inner_loop[0].the_size = it;
    inner_loop[1].the_size = jt;
    inner_loop[2].the_size = kt;

    transpose_B = tr;
    hoist = ht;

    num_threads = nthr;

    char buffer[1024];
    (void)sprintf(buffer, "matmat_%s%s%s_%lux%lux%lu_%s%s%s%s%s_%lux%lux%lu_%lu", block_loop[block_loop_perm[0]].name.c_str(),
                  block_loop[block_loop_perm[1]].name.c_str(), block_loop[block_loop_perm[2]].name.c_str(),
                  block_loop[block_loop_perm[0]].the_size, block_loop[block_loop_perm[1]].the_size,
                  block_loop[block_loop_perm[2]].the_size, (transpose_B ? "BT_" : ""), inner_loop[inner_loop_perm[0]].name.c_str(),
                  inner_loop[inner_loop_perm[1]].name.c_str(), (hoist ? "_H_" : ""), inner_loop[inner_loop_perm[2]].name.c_str(),
                  inner_loop[inner_loop_perm[0]].the_size, inner_loop[inner_loop_perm[1]].the_size,
                  inner_loop[inner_loop_perm[2]].the_size, num_threads);

    (void)fprintf(table, "{ \"%s\", %s },\n", buffer, buffer);
    (void)fprintf(header, "void %s (const Matrix &A, const Matrix& B, Matrix&C);\n", buffer);

    (void)fprintf(source, "\n\nvoid %s (const Matrix &A, const Matrix& B, Matrix&C) {\n", buffer);

    block_outer();

    (void)fprintf(source, "}\n");
  }

};


int main() {

  mm_generator x;

  std::vector<std::vector<size_t>> out = {
    {0, 1, 2}, {1, 0, 2}, {1, 2, 0}, {2, 1, 0}, {2, 0, 1}, {0, 2, 1},
  };

  for (size_t iii = 0; iii < 1; ++iii)
    for (size_t jjj = 0; jjj < 6; ++jjj)

      for (auto ii : {128})
        for (auto jj : {128})
          for (auto kk : {128})
            for (auto i : {2})
              for (auto j : {2})
                for (auto k : {1})
                  for (auto t : {true})
		    for (auto h : {true})		    
		      x.gen(out[iii], ii, jj, kk, t, out[jjj], i, j, k, h);


  return 0;
}
