#!/bin/csh -f

make clean
make mmgen.exe

echo ./mmgen.exe -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -o matmat
./mmgen.exe -m 3 -n 5 -k 7 -M 29 -N 33 -K 31  -o matmat
make mmtest.exe
./mmtest.exe 133 131 313
./mmgen.exe -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -o matmat
make mmtest.exe
./mmtest.exe 133 131 313

./mmgen.exe -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -T -o matmat
make mmtest.exe
./mmtest.exe 133 131 313

./mmgen.exe -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -T -t -o matmat
make mmtest.exe
./mmtest.exe 133 131 313


foreach il ( ijk ikj jik jki kij kji )
  foreach ol ( ijk ikj jik jki kij kji )

    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31  -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    

    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -t -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -t -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    

    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    


    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    

    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -t -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -t -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    


    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -t -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -t -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    

    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -H -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -T -H  -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    


    echo ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -T -t -o matmat
    ./mmgen.exe -l $il -L $ol -m 3 -n 5 -k 7 -M 29 -N 33 -K 31 -H -T -t -o matmat
    make mmtest.exe

foreach i ( 1 27 34 229 256 )
  foreach j ( 1 27 34 229 256 )
    foreach k ( 1 27 34 229 256 )
      ./mmtest.exe $i $j $k
    end
  end
end    


end
end

