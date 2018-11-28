
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>


struct loop {
  std::string name;
  std::string begin;
  std::string end;
  std::string last;
  std::string size;
  size_t the_size;
};


int main() {

  std::string h_path = "matmat.hpp";
  std::string s_path = "matmat.cpp";

  auto header = fopen(h_path.c_str(), "w");
  auto source = fopen(s_path.c_str(), "w");


  for (auto f : { header, source }) {
    (void) fprintf(f, "#include <cstddef>\n");
    (void) fprintf(f, "#include \"Matrix.hpp\"\n\n");
  }

  loop outer[3] = {
    { "ii", "ki_begin", "ii_end", "ii_last", "block_rows", 64 },
    { "jj", "jj_begin", "jj_end", "jj_last", "block_cols", 48 },
    { "kk", "kk_begin", "kk_end", "kk_last", "block_inner", 32 }
  };

  loop inner[3] = {
    { "i", "i_begin", "i_end", "i_last", "tile_rows", 2 },
    { "j", "j_begin", "j_end", "j_last", "tile_cols", 4 },
    { "k", "k_begin", "k_end", "k_last", "tile_inner", 1 }
  };

  size_t inner_perm[3] = { 1, 0, 2 };
  
  for (auto f : { header, source }) {
    (void) fprintf(f, "void matmat_%s%s%s_%lux%lux%lu(const Matrix &A, const Matrix& B, Matrix&C)",
		   inner[inner_perm[0]].name.c_str(), inner[inner_perm[1]].name.c_str(), inner[inner_perm[2]].name.c_str(),
		   inner[inner_perm[0]].the_size, inner[inner_perm[1]].the_size, inner[inner_perm[2]].the_size);
    if (header == f) (void) fprintf(f, "; // prototype\n\n");
    if (source == f) (void) fprintf(f, " {\n");
  }
  
  (void) fprintf(source, "size_t i_begin = 0, i_last = C.num_rows(), i_end = %lu*(i_last/%lu);\n", inner[0].the_size, inner[0].the_size);
  (void) fprintf(source, "size_t j_begin = 0, j_last = C.num_cols(), j_end = %lu*(j_last/%lu);\n", inner[1].the_size, inner[1].the_size);
  (void) fprintf(source, "size_t k_begin = 0, k_last = A.num_cols(), k_end = %lu*(k_last/%lu);\n", inner[2].the_size, inner[2].the_size);
  
  size_t ip = inner_perm[0], jp = inner_perm[1], kp = inner_perm[2];
  (void) fprintf(source, "for (size_t %s = %s; %s < %s; %s += %lu) {\n",
		 inner[ip].name.c_str(), inner[ip].begin.c_str(), inner[ip].name.c_str(), inner[ip].end.c_str(), inner[ip].name.c_str(), inner[ip].the_size);
  (void) fprintf(source, "  for (size_t %s = %s; %s < %s; %s += %lu) {\n",
		 inner[jp].name.c_str(), inner[jp].begin.c_str(), inner[jp].name.c_str(), inner[jp].end.c_str(), inner[jp].name.c_str(), inner[jp].the_size);
  (void) fprintf(source, "    for (size_t %s = %s; %s < %s; %s += %lu) {\n",
		 inner[kp].name.c_str(), inner[kp].begin.c_str(), inner[kp].name.c_str(), inner[kp].end.c_str(), inner[kp].name.c_str(), inner[kp].the_size);
  
  // 0 by 1 by 2
  for (size_t i = 0; i < inner[0].the_size; ++i)  {
    for (size_t j = 0; j < inner[1].the_size; ++j)  {
      for (size_t k = 0; k < inner[2].the_size; ++k) {
	(void) fprintf(source, "      C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * B(%s + %lu, %s + %lu);\n",
		       inner[0].name.c_str(), i, inner[1].name.c_str(), j,
		       inner[0].name.c_str(), i, inner[2].name.c_str(), k,
		       inner[2].name.c_str(), k, inner[1].name.c_str(), j);
      }
    }
  }    
  
  (void) fprintf(source, "    } // %s\n", inner[kp].name.c_str());
  (void) fprintf(source, "    for (size_t %s = %s; %s < %s; ++%s) { // %s cleanup \n", 
		 inner[kp].name.c_str(), inner[kp].end.c_str(), inner[kp].name.c_str(), inner[kp].last.c_str(), inner[kp].name.c_str(), inner[kp].name.c_str());

  // 0 by 1
  for (size_t i = 0; i < inner[ip].the_size; ++i)  {
    for (size_t j = 0; j < inner[jp].the_size; ++j)  {
      for (size_t k = 0; k < 0; ++k)  {
      (void) fprintf(source, "      C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * B(%s + %lu, %s + %lu);\n",
		     inner[0].name.c_str(), i, inner[1].name.c_str(), j,
		     inner[0].name.c_str(), i, inner[2].name.c_str(), k,
		     inner[2].name.c_str(), k, inner[1].name.c_str(), j);
      }
    }
  }

  (void) fprintf(source, "    }\n");
  (void) fprintf(source, "  } // %s\n", inner[jp].name.c_str());


  (void) fprintf(source, "  for (size_t %s = %s; %s < %s; ++%s) { // %s cleanup \n", 
		 inner[jp].name.c_str(), inner[jp].end.c_str(), inner[jp].name.c_str(), inner[jp].last.c_str(), inner[jp].name.c_str(), inner[jp].name.c_str());
  (void) fprintf(source, "    for (size_t %s = %s; %s < %s; ++%s) {\n",
		 inner[kp].name.c_str(), inner[kp].begin.c_str(), inner[kp].name.c_str(), inner[kp].end.c_str(), inner[kp].name.c_str());

  for (size_t i = 0; i < inner[0].the_size; ++i)  {
    size_t k = 0; size_t j = 0;
      (void) fprintf(source, "      C(%s + %lu, %s + %lu) += A(%s + %lu, %s + %lu) * B(%s + %lu, %s + %lu);\n",
		     inner[0].name.c_str(), i, inner[1].name.c_str(), j,
		     inner[0].name.c_str(), i, inner[2].name.c_str(), k,
		     inner[2].name.c_str(), k, inner[1].name.c_str(), j);
  }

  (void) fprintf(source, "    }\n}\n");
  
  (void) fprintf(source, "} // %s\n", inner[ip].name.c_str());
  (void) fprintf(source, "for (size_t %s = %s; %s < %s; ++%s) { // %s cleanup \n", 
		 inner[ip].name.c_str(), inner[ip].end.c_str(), inner[ip].name.c_str(), inner[ip].last.c_str(), inner[ip].name.c_str(), inner[ip].name.c_str());

  (void) fprintf(source, "  for (size_t %s = %s; %s < %s; ++%s) {\n",
		 inner[jp].name.c_str(), inner[jp].begin.c_str(), inner[jp].name.c_str(), inner[jp].end.c_str(), inner[jp].name.c_str());
  (void) fprintf(source, "    for (size_t %s = %s; %s < %s; ++%s) {\n",
		 inner[kp].name.c_str(), inner[kp].begin.c_str(), inner[kp].name.c_str(), inner[kp].end.c_str(), inner[kp].name.c_str());
  (void) fprintf(source, "      C(%s, %s) += A(%s, %s) * B(%s, %s);\n",
		 inner[0].name.c_str(), inner[1].name.c_str(),
		 inner[0].name.c_str(), inner[2].name.c_str(),
		 inner[2].name.c_str(), inner[1].name.c_str());
  (void) fprintf(source, "}\n}\n");

  (void) fprintf(source, "  for (size_t %s = %s; %s < %s; ++%s) {\n",
		 inner[jp].name.c_str(), inner[jp].end.c_str(), inner[jp].name.c_str(), inner[jp].last.c_str(), inner[jp].name.c_str());
  (void) fprintf(source, "    for (size_t %s = %s; %s < %s; ++%s) {\n",
		 inner[kp].name.c_str(), inner[kp].begin.c_str(), inner[kp].name.c_str(), inner[kp].end.c_str(), inner[kp].name.c_str());
  (void) fprintf(source, "      C(%s, %s) += A(%s, %s) * B(%s, %s);\n",
		 inner[0].name.c_str(), inner[1].name.c_str(),
		 inner[0].name.c_str(), inner[2].name.c_str(),
		 inner[2].name.c_str(), inner[1].name.c_str());
  (void) fprintf(source, "}\n}\n");

  (void) fprintf(source, "}\n");
  
  (void) fprintf(source, "}\n");
  
  (void) fclose(header);
  (void) fclose(source);
  
  return 0;
}
