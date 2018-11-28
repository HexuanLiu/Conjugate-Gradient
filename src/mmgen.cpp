#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class matmat_generator {
public:
  matmat_generator(std::tuple<size_t, size_t, size_t> block_order, std::tuple<size_t, size_t, size_t> block_sizes,
                   std::tuple<size_t, size_t, size_t> inner_order, std::tuple<size_t, size_t, size_t> tile_sizes, bool transpose,
                   bool hoist, bool thread)
      : block_order_(block_order), block_sizes_(block_sizes), inner_order_(block_order), tile_sizes_(block_sizes),
        transpose_(transpose), hoist_(hoist), thread_(thread),
        inner_loop_(std::vector<loop>({
            {"i", "ib", "i_begin", "i_end", "i_last", "tile_rows", std::get<0>(tile_sizes)},
            {"j", "jb", "j_begin", "j_end", "j_last", "tile_cols", std::get<1>(tile_sizes)},
            {"k", "kb", "k_begin", "k_end", "k_last", "tile_inner", std::get<2>(tile_sizes)},
        })),
        block_loop_(std::vector<loop>({
            {"ii", "iib", "ii_begin", "ii_end", "ii_last", "block_rows", std::get<0>(block_sizes)},
            {"jj", "jjb", "jj_begin", "jj_end", "jj_last", "block_cols", std::get<1>(block_sizes)},
            {"kk", "kkb", "kk_begin", "kk_end", "kk_last", "block_inner", std::get<2>(block_sizes)},
        })),
        inner_loop_perm_({std::get<0>(inner_order_), std::get<1>(inner_order_), std::get<2>(inner_order_)}),
        block_loop_perm_({std::get<0>(block_order_), std::get<1>(block_order_), std::get<2>(block_order_)}),
        the_function_name_("matmat") {}

private:
  std::tuple<size_t, size_t, size_t> block_order_, block_sizes_, inner_order_, tile_sizes_;
  bool                               transpose_, hoist_, thread_;

  struct loop {
    std::string name;
    std::string b_name;
    std::string begin;
    std::string end;
    std::string last;
    std::string size;
    size_t      the_size;
  };

  std::vector<loop>   inner_loop_, block_loop_;
  std::vector<size_t> inner_loop_perm_, block_loop_perm_;

  std::vector<size_t> inner_loop_perm = {std::get<0>(inner_order_), std::get<1>(inner_order_), std::get<2>(inner_order_)};
  std::vector<size_t> block_loop_perm = {std::get<0>(block_order_), std::get<1>(block_order_), std::get<2>(block_order_)};

  std::vector<std::string> the_function_;
  std::string              the_function_name_;

  std::string default_function_name() {
    return formatted("matmat_%s%s%s_%lux%lux%lu_%s%s%s%s%s_%lux%lux%lu%s", block_loop_[block_loop_perm[0]].name.c_str(),
                     block_loop_[block_loop_perm[1]].name.c_str(), block_loop_[block_loop_perm[2]].name.c_str(),
                     block_loop_[block_loop_perm[0]].the_size, block_loop_[block_loop_perm[1]].the_size,
                     block_loop_[block_loop_perm[2]].the_size, (transpose_ ? "BT_" : ""),
                     inner_loop_[inner_loop_perm[0]].name.c_str(), inner_loop_[inner_loop_perm[1]].name.c_str(),
                     (hoist_ ? "_H_" : ""), inner_loop_[inner_loop_perm[2]].name.c_str(), inner_loop_[inner_loop_perm[0]].the_size,
                     inner_loop_[inner_loop_perm[1]].the_size, inner_loop_[inner_loop_perm[2]].the_size, (thread_ ? "_MT" : ""));
  }

  // utility
  template<typename... Args>
  void pb(std::vector<std::string>& str_vec, const char* format, Args... args) {
    str_vec.push_back(formatted(format, args...));
  }

  template<typename... Args>
  void pb(const char* format, Args... args) {
    pb(the_function_, format, args...);
  }

  template<typename... Args>
  std::string formatted(const char* format, Args... args) {

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
    size_t                  size = std::snprintf(nullptr, 0, format, args...) + 1;    // Extra space for '\0'
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format, args...);
#pragma clang diagnostic pop

    return std::string(buf.get(), buf.get() + size - 1);    // We don't want the '\0' inside
  }

  void newline(std::vector<std::string>& x) { pb(x, ""); };
  void newline() { pb(the_function_, ""); };

  // For loops
  void for_begin_to_last_by_size(loop lp) {
    pb("for (");
    pb("size_t %s = %s, %s = 0; %s < %s; %s += %lu, %s += %lu", lp.name.c_str(), lp.begin.c_str(), lp.b_name.c_str(),
       lp.name.c_str(), lp.last.c_str(), lp.name.c_str(), lp.the_size, lp.b_name.c_str(), lp.the_size);
    pb(") {");
  }

  void for_begin_to_end_by_size(loop lp) {
    pb("for (");
    pb("size_t %s = %s, %s = 0; %s < %s; %s += %lu, %s += %lu", lp.name.c_str(), lp.begin.c_str(), lp.b_name.c_str(),
       lp.name.c_str(), lp.end.c_str(), lp.name.c_str(), lp.the_size, lp.b_name.c_str(), lp.the_size);
    pb(") {");
  }

  void for_end_to_last_by_one(loop lp) {
    pb("for (");
    pb("size_t %s = %s, %s = %s - %s; %s < %s; ++%s, ++%s", lp.name.c_str(), lp.end.c_str(), lp.b_name.c_str(), lp.end.c_str(),
       lp.begin.c_str(), lp.name.c_str(), lp.last.c_str(), lp.name.c_str(), lp.b_name.c_str());
    pb(") {");
  }

  // Loop setup
  void block_limits() {
    pb("// block limits ");
    pb("size_t ii_begin = 0, ii_last = C.num_rows();");
    pb("size_t jj_begin = 0, jj_last = C.num_cols();");
    pb("size_t kk_begin = 0, kk_last = A.num_cols();");
  }

  void inner_limits() {
    newline();
    inner_limits(0);
    inner_limits(1);
    inner_limits(2);
    newline();
  }

  void inner_limits(size_t idx) {
    pb("size_t %s = %s, %s_step = std::min<size_t>(%lu, %s-%s), %s = %s + %lu*(%s_step/%lu), %s = %s + %s_step;",
       inner_loop_[idx].begin.c_str(), block_loop_[idx].name.c_str(), inner_loop_[idx].name.c_str(), block_loop_[idx].the_size,
       block_loop_[idx].last.c_str(), block_loop_[idx].name.c_str(), inner_loop_[idx].end.c_str(), inner_loop_[idx].begin.c_str(),
       inner_loop_[idx].the_size, inner_loop_[idx].name.c_str(), inner_loop_[idx].the_size, inner_loop_[idx].last.c_str(),
       inner_loop_[idx].begin.c_str(), inner_loop_[idx].name.c_str());
  }

  // Hoisting
  void hoist_pre() {
    size_t which_to_hoist = inner_loop_perm[2];    // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = (0 == which_to_hoist ? 1 : inner_loop_[0].the_size);
    size_t jlimit = (1 == which_to_hoist ? 1 : inner_loop_[1].the_size);
    size_t klimit = (2 == which_to_hoist ? 1 : inner_loop_[2].the_size);

    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {
          switch (which_to_hoist) {
            case 2:
              pb("double t%lu%lu = C(%s + %lu, %s + %lu);", i, j, inner_loop_[0].name.c_str(), i, inner_loop_[1].name.c_str(), j);
              break;
            case 1:
              pb("double t%lu%lu = A(%s + %lu, %s + %lu);", i, k, inner_loop_[0].name.c_str(), i, inner_loop_[2].name.c_str(), k);
              break;
            case 0:
              if (transpose_) {
                pb("double t%lu%lu = B_T(%s + %lu, %s + %lu);", j, k, inner_loop_[1].b_name.c_str(), j,
                   inner_loop_[2].b_name.c_str(), k);
              } else {
                pb("double t%lu%lu = B(%s + %lu, %s + %lu);", k, j, inner_loop_[2].name.c_str(), k, inner_loop_[1].name.c_str(), j);
              }
              break;
            default:
              throw;
          }
        }
      }
    }
  }

  void hoist_post() {
    size_t which_to_hoist = inner_loop_perm[2];    // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = (0 == which_to_hoist ? 1 : inner_loop_[0].the_size);
    size_t jlimit = (1 == which_to_hoist ? 1 : inner_loop_[1].the_size);
    size_t klimit = (2 == which_to_hoist ? 1 : inner_loop_[2].the_size);

    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {
          switch (which_to_hoist) {
            case 2:
              pb("C(%s + %lu, %s + %lu) = t%lu%lu;", inner_loop_[0].name.c_str(), i, inner_loop_[1].name.c_str(), j, i, j);
              break;
            case 1:
            case 0:
              break;
            default:
              throw;
          }
        }
      }
    }
  }

  // Loop bodies

  void hoisted_loop_body() {
    size_t hoisted = inner_loop_perm[2];    // 2 -> C, 1 -> A, 0 -> B

    size_t ilimit = inner_loop_[0].the_size;
    size_t jlimit = inner_loop_[1].the_size;
    size_t klimit = inner_loop_[2].the_size;

    for (size_t i = 0; i < ilimit; ++i) {
      for (size_t j = 0; j < jlimit; ++j) {
        for (size_t k = 0; k < klimit; ++k) {

          switch (hoisted) {
            case 2:
              pb("t%lu%lu += A(%s + %lu, %s + %lu) * ", i, j, inner_loop_[0].name.c_str(), i, inner_loop_[2].name.c_str(), k);
              if (transpose_) {
                pb("B_T(%s + %lu, %s + %lu);", inner_loop_[1].b_name.c_str(), j, inner_loop_[2].b_name.c_str(), k);
              } else {
                pb("B(%s + %lu, %s + %lu);", inner_loop_[2].name.c_str(), k, inner_loop_[1].name.c_str(), j);
              }
              break;
            case 1:
              pb("C(%s + %lu, %s + %lu) += t%lu%lu * ", inner_loop_[0].name.c_str(), i, inner_loop_[1].name.c_str(), j, i, k);
              if (transpose_) {
                pb("B_T(%s + %lu, %s + %lu);", inner_loop_[1].b_name.c_str(), j, inner_loop_[2].b_name.c_str(), k);
              } else {
                pb("B(%s + %lu, %s + %lu);", inner_loop_[2].name.c_str(), k, inner_loop_[1].name.c_str(), j);
              }
              break;
            case 0:
              if (transpose_) {
                pb("C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * t%lu%lu;", inner_loop_[0].name.c_str(), i,
                   inner_loop_[1].name.c_str(), j, inner_loop_[0].name.c_str(), i, inner_loop_[2].name.c_str(), k, j, k);
              } else {
                pb("C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * t%lu%lu;", inner_loop_[0].name.c_str(), i,
                   inner_loop_[1].name.c_str(), j, inner_loop_[0].name.c_str(), i, inner_loop_[2].name.c_str(), k, k, j);
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
    size_t ilimit = inner_loop_[0].the_size;
    size_t jlimit = inner_loop_[1].the_size;
    size_t klimit = inner_loop_[2].the_size;

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

          pb("      C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * ", inner_loop_[0].name.c_str(), i, inner_loop_[1].name.c_str(),
             j, inner_loop_[0].name.c_str(), i, inner_loop_[2].name.c_str(), k);

          if (transpose_) {
            pb("B_T(%s + %lu, %s + %lu);", inner_loop_[1].b_name.c_str(), j, inner_loop_[2].b_name.c_str(), k);
          } else {
            pb("B(%s + %lu, %s + %lu);", inner_loop_[2].name.c_str(), k, inner_loop_[1].name.c_str(), j);
          }
        }
      }
    }
  }

  void inner_main_begin(bool outer_fringe, bool middle_fringe) {
    size_t p2 = inner_loop_perm[2];
    for_begin_to_end_by_size(inner_loop_[p2]);
  }
  void inner_main_end() { pb("}"); }

  void inner_fringe_begin(bool outer_fringe, bool middle_fringe) {
    size_t p2 = inner_loop_perm[2];
    for_end_to_last_by_one(inner_loop_[p2]);
  }
  void inner_fringe_end() { pb("}"); }

  void middle_main_begin(bool outer_fringe) {
    size_t p1 = inner_loop_perm[1];
    for_begin_to_end_by_size(inner_loop_[p1]);
  }
  void middle_main_end() { pb("}"); }

  void middle_fringe_begin(bool outer_fringe) {
    size_t p1 = inner_loop_perm[1];
    for_end_to_last_by_one(inner_loop_[p1]);
  }
  void middle_fringe_end() { pb("}"); }

  void outer_main_begin() {
    size_t p0 = inner_loop_perm[0];
    for_begin_to_end_by_size(inner_loop_[p0]);
  }
  void outer_main_end() { pb("}"); }

  void outer_fringe_begin() {
    size_t p0 = inner_loop_perm[0];
    for_end_to_last_by_one(inner_loop_[p0]);
  }
  void outer_fringe_end() { pb("}"); }

  // The inner loop
  void inner(bool outer_fringe, bool middle_fringe) {

    if (hoist_ && !outer_fringe && !middle_fringe) {
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

  void transpose_block() {
    if (transpose_) {
      newline();
      pb("Matrix B_T(j_step, k_step);");
      pb("for (size_t jb = 0, jt = j_begin; jb < j_step; ++jb, ++jt) {");
      pb("for (size_t kb = 0, kt = k_begin; kb < k_step; ++kb, ++kt) {");
      pb("B_T(jb, kb) = B(kt, jt);");
      pb("}\n}");
      newline();
    }
  }

  // The outer loop (register level)
  void outer() {

    inner_limits();

    transpose_block();

    outer_main_begin();
    middle(false);
    outer_main_end();
    outer_fringe_begin();
    middle(true);
    outer_fringe_end();
  }

  void thread_prelude() {
    if (thread_) {
      newline();
      pb("size_t                         thread_count = 0;");
      pb("std::vector<std::future<void>> futs(num_threads);");

      if (!(hoist_ && (inner_loop_perm_[2] == 2))) {    // Dont need to lock if hoisting C
        pb("std::vector<std::mutex> mtcs((%s + %lu) / %lu);", block_loop_[block_loop_perm[2]].last.c_str(),
           block_loop_[block_loop_perm[2]].the_size - 1, block_loop_[block_loop_perm[2]].the_size);
      }
      newline();
    }
  }

  void lock_guard(loop lp) {
    if (thread_) {
      if (hoist_ && (inner_loop_perm_[2] == 2))    // Dont need to lock if hoisting C
        return;

      newline();
      pb("std::lock_guard<std::mutex> aa(mtcs[%s/%lu]);", lp.b_name.c_str(), lp.the_size);
      newline();
    }
  }

  void begin_thread() {
    if (thread_) {
      newline();
      pb("futs[thread_count++] = std::async(std::launch::async, [&, %s, %s]() -> void {",
         block_loop_[block_loop_perm[0]].name.c_str(), block_loop_[block_loop_perm[1]].name.c_str());

      newline();
    }
  }

  void end_thread() {
    if (thread_) {
      newline();
      pb("});");
      newline();
      pb("if (thread_count == num_threads) {");
      cleanup_thread();
      pb("}");
      newline();
    }
  }

  void cleanup_thread() {
    if (thread_) {
      newline();
      pb("for (size_t i = 0; i < thread_count; ++i) { // wait on launched tasks");
      pb("futs[i].wait();");
      pb("}");
      pb("thread_count = 0;");
      newline();
    }
  }

  void end_for(loop lp) { pb("} // End %s", lp.name.c_str()); }

  // The blocked outer loop
  void block_outer() {

    block_limits();

    thread_prelude();

    for_begin_to_last_by_size(block_loop_[block_loop_perm[0]]);
    for_begin_to_last_by_size(block_loop_[block_loop_perm[1]]);

    begin_thread();

    for_begin_to_last_by_size(block_loop_[block_loop_perm[2]]);

    lock_guard(block_loop_[block_loop_perm[2]]);

    outer();

    end_for(block_loop_[block_loop_perm[2]]);

    end_thread();

    end_for(block_loop_[block_loop_perm[1]]);
    end_for(block_loop_[block_loop_perm[0]]);

    cleanup_thread();
  }

public:
  void gen_function() {
    pb("\n\nvoid %s (const Matrix &A, const Matrix& B, Matrix&C, size_t num_threads) {", the_function_name_.c_str());

    block_outer();

    pb("}");
  }

  std::vector<std::string> get_function() { return the_function_; }

  void gen_name() { the_function_name_ = default_function_name(); }

  std::vector<std::string> get_source_file(const std::string& base_filename) {
    std::vector<std::string> source_file;
    source_file.push_back("#include <cstddef>");
    source_file.push_back("#include <future>");
    source_file.push_back("#include <mutex>");
    source_file.push_back("#include \"Matrix.hpp\"");
    source_file.push_back("#include \"" + base_filename + ".hpp\"");

    for (auto& j : the_function_) {
      source_file.push_back(j);
    }
    return source_file;
  }

  std::vector<std::string> get_header_file(const std::string& base_filename) {
    std::vector<std::string> source_file;
    source_file.push_back("#include <cstddef>");
    source_file.push_back("#include \"Matrix.hpp\"");
    newline(source_file);

    source_file.push_back(
        formatted("void %s (const Matrix &A, const Matrix& B, Matrix&C, size_t num_threads = 1);", the_function_name_.c_str()));

    return source_file;
  }

  std::vector<std::string> get_table_file(const std::string& base_filename) {
    std::vector<std::string> source_file;
    source_file.push_back("#include <cstddef>");
    source_file.push_back("#include <string>");
    source_file.push_back("#include \"Matrix.hpp\"");
    source_file.push_back("#include \"" + base_filename + ".hpp\"");

    return source_file;
  }
};

std::tuple<size_t, size_t, size_t> string_to_tuple(std::string str) {
  for (auto& c : str)
    c = toupper(c);
  if (str == "IJK") return {0, 1, 2};
  if (str == "JIK") return {1, 0, 2};
  if (str == "IKJ") return {0, 2, 1};
  if (str == "KIJ") return {2, 0, 1};
  if (str == "JKI") return {1, 2, 0};
  if (str == "KJI") return {2, 1, 0};
  throw;
};

void usage(char* program_name) {
  std::cout << program_name << " ";
  std::cout << "  [ -l ijk ] // inner loop order, default = ijk" << std::endl;
  std::cout << "  [ -L IJK ] // outer loop order, default = I, J, K" << std::endl;
  std::cout << "  [ -m tile size ] [ -n tile size ] [ -k tile size ], default = 2, 2, 1" << std::endl;
  std::cout << "  [ -M block size ] [ -N block size ] [ -K block size ], default = 64, 64, 64" << std::endl;
  std::cout << "  [ -H ] // don't hoist inner loop invariants, default = true" << std::endl;
  std::cout << "  [ -T ] // don't use block copy transpose, default = true" << std::endl;
  std::cout << "  [ -t ] // multithreaded, default = false " << std::endl;
  std::cout << "  [ -F ] // generate function based on function parameters" << std::endl;
  std::cout << "  [ -f function_name ] // specify function name, default = \"matmat\"" << std::endl;
  std::cout << "  [ -h ] // print this message" << std::endl;
  std::cout << "  [ -o baseneme ] // generate source, header, table files" << std::endl;
  std::cout << "  [ -a baseneme ] // append to source, header, table files" << std::endl;
}

int main(int argc, char* argv[]) {

  std::string                        inner_order   = "ijk";
  std::string                        block_order   = "IJK";
  std::tuple<size_t, size_t, size_t> tile_sizes    = {2, 2, 1};
  std::tuple<size_t, size_t, size_t> block_sizes   = {64, 64, 64};
  bool                               hoist         = true;
  bool                               transpose     = true;
  bool                               thread        = false;
  bool                               generate_name = false;
  std::string                        function_name = "";
  std::string                        base_filename = "";

  try {

    for (size_t arg = 1; arg < argc; ++arg) {

      if (std::string(argv[arg]) == "-o") {
        if (argc == ++arg) usage(argv[0]);
        base_filename = std::string(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-l") {
        if (argc == ++arg) usage(argv[0]);
        inner_order = std::string(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-L") {
        if (argc == ++arg) usage(argv[0]);
        block_order = std::string(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-m") {
        if (argc == ++arg) usage(argv[0]);
        std::get<0>(tile_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-n") {
        if (argc == ++arg) usage(argv[0]);
        std::get<1>(tile_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-k") {
        if (argc == ++arg) usage(argv[0]);
        std::get<2>(tile_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-M") {
        if (argc == ++arg) usage(argv[0]);
        std::get<0>(block_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-N") {
        if (argc == ++arg) usage(argv[0]);
        std::get<1>(block_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-K") {
        if (argc == ++arg) usage(argv[0]);
        std::get<2>(block_sizes) = std::stol(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-f") {
        if (generate_name) {
          usage(argv[0]);
          return -1;
        }
        if (argc == ++arg) {
          usage(argv[0]);
          return -1;
        }
        function_name = std::string(argv[arg]);
      }

      else if (std::string(argv[arg]) == "-t") {
        thread = true;
      }

      else if (std::string(argv[arg]) == "-F") {
        if (function_name != "") {
          usage(argv[0]);
          return -1;
        }
        generate_name = true;
      }

      else if (std::string(argv[arg]) == "-T") {
        transpose = !transpose;
      }

      else if (std::string(argv[arg]) == "-H") {
        hoist = !hoist;
      }

      else {
        usage(argv[0]);
        return -1;
      }
    }
  } catch (int) {
    usage(argv[0]);
    return -1;
  }

  matmat_generator x(string_to_tuple(block_order), block_sizes, string_to_tuple(inner_order), tile_sizes, transpose, hoist, thread);

  x.gen_function();

  if ("" == base_filename) {
    for (auto& y : x.get_function())
      std::cout << y << std::endl;
  } else {
    std::string source_filename = base_filename + ".cpp";
    std::string header_filename = base_filename + ".hpp";
    std::string table_filename  = base_filename + "_table.hpp";

    std::ofstream source_file(source_filename);
    for (auto& y : x.get_source_file(base_filename)) {
      source_file << y << std::endl;
    }

    std::ofstream header_file(header_filename);
    for (auto& y : x.get_header_file(base_filename)) {
      header_file << y << std::endl;
    }

    std::ofstream table_file(table_filename);
    for (auto& y : x.get_table_file(base_filename)) {
      table_file << y << std::endl;
    }
  }

  return 0;
}
