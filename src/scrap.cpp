#if 0

public:
  // Generate a matmat function
  void gen(std::vector<size_t> bl_perm, size_t ib, size_t jb, size_t kb, bool tr, std::vector<size_t> il_perm, size_t it, size_t jt,
           size_t kt, bool ht, size_t nthr = 0) {

    block_loop_perm = bl_perm;
    inner_loop_perm = il_perm;

    block_loop_[0].the_size = ib;
    block_loop_[1].the_size = jb;
    block_loop_[2].the_size = kb;

    inner_loop_[0].the_size = it;
    inner_loop_[1].the_size = jt;
    inner_loop_[2].the_size = kt;

    transpose_ = tr;
    hoist       = ht;

    num_threads = nthr;


    (void)fprintf(table, "{ \"%s\", %s },", buffer, buffer);
    (void)fprintf(header, "void %s (const Matrix &A, const Matrix& B, Matrix&C);", buffer);

    pb("\n\nvoid %s (const Matrix &A, const Matrix& B, Matrix&C) {", buffer);

    block_outer();

    pb("}");
  }

#endif
for (auto j : matmat_table) {

  matmat_generator() {
    source = fopen("matmat.cpp", "w");
    header = fopen("matmat.hpp", "w");
    table  = fopen("mattable.hpp", "w");
    for (auto f : {header, source, table}) {
      (void)fprintf(f, "#include <cstddef>\n");
      (void)fprintf(f, "#include \"Matrix.hpp\"\n");
    }
    pb("#include <future>\n");

    (void)fprintf(table, "#include <string>\n");
    (void)fprintf(table, "#include <tuple>\n");
    (void)fprintf(table, "#include \"matmat.hpp\"\n");
    (void)fprintf(table, "std::tuple<std::string, void (*)(const Matrix&, const Matrix&, Matrix&)> matmat_table[] = {\n");

    for (auto f : {header, source, table}) {
      newline(f);
    }
  }

  ~matma_generator() {
    (void)fprintf(table, "};\n");

    fclose(source);
    fclose(header);
    fclose(table);
  }

  void outer0() {

    transpose_B = false;

    (void)fprintf(source, "\n");
    (void)fprintf(source, "size_t i_begin = 0, i_last = C.num_rows(), i_end = %lu*(i_last/%lu);\n", inner_loop[0].the_size,
                  inner_loop[0].the_size);
    (void)fprintf(source, "size_t j_begin = 0, j_last = C.num_cols(), j_end = %lu*(j_last/%lu);\n", inner_loop[1].the_size,
                  inner_loop[1].the_size);
    (void)fprintf(source, "size_t k_begin = 0, k_last = A.num_cols(), k_end = %lu*(k_last/%lu);\n", inner_loop[2].the_size,
                  inner_loop[2].the_size);
    (void)fprintf(source, "\n");

    outer_main_begin();
    middle(false);
    outer_main_end();
    outer_fringe_begin();
    middle(true);
    outer_fringe_end();
  }
  void gen0(size_t il_perm[3], size_t it, size_t jt, size_t kt) {

    for (size_t i = 0; i < 3; ++i) {
      inner_loop_perm[i] = il_perm[i];
    }

    inner_loop[0].the_size = it;
    inner_loop[1].the_size = jt;
    inner_loop[2].the_size = kt;

    for (auto f : {header, source, table}) {
      if (table != f) (void)fprintf(f, "void ");

      (void)fprintf(f, "matmat_%s%s%s_%lux%lux%lu", inner_loop[inner_loop_perm[0]].name.c_str(),
                    inner_loop[inner_loop_perm[1]].name.c_str(), inner_loop[inner_loop_perm[2]].name.c_str(),
                    inner_loop[inner_loop_perm[0]].the_size, inner_loop[inner_loop_perm[1]].the_size,
                    inner_loop[inner_loop_perm[2]].the_size);

      if (table != f) (void)fprintf(f, "(const Matrix &A, const Matrix& B, Matrix&C)");
      if (table == f) (void)fprintf(f, ",\n");
      if (header == f) (void)fprintf(f, "; // prototype\n\n");
      if (source == f) (void)fprintf(f, " {\n");
    }

    outer0();

    (void)fprintf(source, "}\n");
  }

  if (transpose_B) {
    (void)fprintf(source, "Matrix B_T(j_last-j_begin, k_last-k_begin);\n");
    (void)fprintf(source, "for (size_t j = j_begin, jb = 0; j < j_end; ++j, ++jb) {\n");
    (void)fprintf(source, "for (size_t k = k_begin, kb = 0; k < k_end; ++k, ++kb) {\n");
    (void)fprintf(source, "B_T(jb, kb) = B(k, j);\n");
    (void)fprintf(source, "}\n}\n");
    (void)fprintf(source, "\n");
  }

  size_t i_begin = ii, i_step = std::min<size_t>(32, ii_last - ii), i_end = i_begin + 2 * (i_step / 2), i_last = i_begin + i_step;

  // (void)fprintf(source, "std:: cout << i_begin << " " << i_end << " " << j_begin << " " << j_end << " " << k_begin << " " << k_end << std::endl\n");

  if (transpose_B) {
    (void)fprintf(source, "TransposedMatrixView B_T(B, 0, B.num_rows(), 0, B.num_cols());\n");
  }
  void size_sweep0(mm_generator & x, size_t out[3], size_t in[3]) {
    x.gen0(in, 1, 1, 1);
    x.gen0(in, 1, 2, 2);
    x.gen0(in, 2, 2, 2);
    x.gen0(in, 2, 3, 2);
    x.gen0(in, 2, 3, 5);
    x.gen0(in, 5, 3, 7);
    x.gen0(in, 5, 7, 3);
    x.gen0(in, 3, 7, 3);
    x.gen0(in, 3, 7, 13);
    x.gen0(in, 3, 7, 5);
  }

  void size_sweep2(mm_generator & x, std::vector<size_t> out, std::vector<size_t> in) {

    for (auto ii : {3})
      for (auto jj : {3})
        for (auto kk : {3})
          for (auto i : {1, 2})
            for (auto j : {1, 2})
              for (auto k : {1, 2})
                for (auto t : {true, false})

                  x.gen2(out, ii, jj, kk, t, in, i, j, k);

    for (auto ii : {31, 32, 33})
      for (auto jj : {31, 32, 33})
        for (auto kk : {31, 32, 33})
          for (auto i : {1, 2, 3, 5})
            for (auto j : {1, 2, 3, 5})
              for (auto k : {1, 2, 3, 5})
                for (auto t : {true, false})

                  x.gen2(out, ii, jj, kk, t, in, i, j, k);
  }

  void order_sweep(mm_generator & x, const std::vector<std::vector<size_t>>& out, const std::vector<std::vector<size_t>>& in) {

    for (size_t iii = 0; iii < 5; ++iii)
      for (size_t jjj = 3; jjj < 5; ++jjj)

        for (auto ii : {32, 31})
          for (auto jj : {33, 32})
            for (auto kk : {32, 55})
              for (auto i : {3, 5})
                for (auto j : {1, 3})
                  for (auto k : {2, 3, 4})
                    for (auto t : {true, false}) {

                      x.gen2(out[iii], ii, jj, kk, t, in[jjj], i, j, k);

#if 0
		    for (auto k : out[iii]) 
		      std::cout << k << std::endl;
		    for (auto k : in[jjj]) 
		      std::cout << k << std::endl;

#endif
                    }
  }

  void sanity_sweep(mm_generator & x, const std::vector<std::vector<size_t>>& out, const std::vector<std::vector<size_t>>& in) {

    for (size_t iii = 2; iii < 3; ++iii)
      for (size_t jjj = 4; jjj < 6; ++jjj)

        for (auto ii : {32, 31})
          for (auto jj : {32, 33})
            for (auto kk : {29, 32})
              for (auto i : {1, 2, 3})
                for (auto j : {1, 5})
                  for (auto k : {1, 3})
                    for (auto t : {true, false})
                      x.gen2(out[iii], ii, jj, kk, t, in[jjj], i, j, k);
  }
