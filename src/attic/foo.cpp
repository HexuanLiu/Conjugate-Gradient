

#include <iostream>
#include <string>

template<typename... Args> 
std::string formatted(const char* format, Args... args) {
    size_t                  size = std::snprintf(nullptr, 0, format, args...) + 1;    // Extra space for '\0'
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format, args...);
    return std::string(buf.get(), buf.get() + size - 1);    // We don't want the '\0' inside
  }




template<typename... Args> 
void f1(const char* format, Args... args) {
  std::cout << formatted(format, args...) << std::endl;
}

template<typename... Args> 
void f2(const char* format, Args... args) {
  f1(format, args...);
}

template<typename... Args> 
void f3(const char* format, Args... args) {
  f2(format, args...);
}


int main() {

  f3("pi = %lf\n", 3.14);
  return 0;
}
