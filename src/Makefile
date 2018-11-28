

#
# This file is part of the course materials for AMATH483/583 at the University of Washington,
# Spring 2018
#
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
# https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Author: Andrew Lumsdaine
# 

LANG		= -std=c++11 -Wc++11-extensions
PICKY		= -Wall

CXX		= c++
#OPTS		= -Ofast -march=native -DNDEBUG # -pthread
DEBUG		= -g
REPORTING   	= -Rpass=.*
INCLUDES	= 
VPATH		= 
LIB		= 

CXXFLAGS	= $(OPTS) $(DEBUG) $(LANG) $(PICKY) $(INCLUDES)

TESTS		= mmtest.cpp mmtime.cpp mmgen.cpp
SOURCES		= matmat.cpp Matrix.cpp

HEADERS		= $(SOURCES:.cpp=.hpp) 
OBJECTS		= $(SOURCES:.cpp=.o) $(TESTS:.cpp=.o)
TARGETS		= $(TESTS:.cpp=.exe)
PCHS		= $(HEADERS:=.gch)

.PHONY		: defreport optreport clean depend all

all		: $(TARGETS)

%.exe        	: %.o
		  $(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

%.o 		: %.cpp
		  $(CXX) -c $(CXXFLAGS) $< -o $@

%.s 		: %.cpp
		  $(CXX) -S $(CXXFLAGS) $<

matmat.cpp	: mmgen.exe
		  ./mmgen.exe -o matmat
		  clang-format -i matmat.cpp
		  clang-format -i matmat.hpp
		  clang-format -i matmat_table.hpp

mmtime.o	: matmat.cpp
mmtest.o	: matmat.cpp

mmtime.exe	: mmtime.o matmat.o Matrix.o
mmtest.exe	: mmtest.o matmat.o Matrix.o

defreport	:
		  $(CXX) -dM -E -x c++ /dev/null

optreport	:
		  echo 'int;' | $(CXX) -xc++ $(CXXFLAGS) - -o /dev/null -\#\#\#

clean		:
		  /bin/rm -f $(TARGETS) $(OBJECTS) $(PCHS) Matrix.s a.out *~

depend: $(TESTS) $(SOURCES)
	@ $(CXX) -MM $(LANG) $(INCLUDES) $(SOURCES) > makedep
	@ echo '/^# DO NOT DELETE THIS LINE/+2,$$d' >eddep
	@ echo '$$r makedep' >>eddep
	@ echo 'w' >>eddep
	@ cp Makefile Makefile.bak
	@ ed - Makefile < eddep
	@ /bin/rm eddep makedep
	@ echo '# DEPENDENCIES MUST END AT END OF FILE' >> Makefile
	@ echo '# IF YOU PUT STUFF HERE IT WILL GO AWAY' >> Makefile
	@ echo '# see make depend above' >> Makefile


# The following 4 (yes 4) lines must be in all the subdirectory makefiles
#-----------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend uses it
# DEPENDENCIES MUST END AT END OF FILE
matmat.o: matmat.cpp Matrix.hpp Vector.hpp
mmgen.o: mmgen.cpp
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
