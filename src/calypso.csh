#!/bin/csh -f

make clean
make mmgen.exe


foreach ol ( ijk ikj jik jki kij kji )
  foreach il ( ijk ikj jik jki kij kji )
    foreach M (32 64 128)
      foreach N (32 64 128)
         foreach K (32 64 128)
	   foreach m ( 1 2 4 )
  	     foreach n ( 1 2 4 )
  	       foreach k ( 1 2 4 )
        echo ./mmgen.exe -l $il -L $ol -m $m -n $n -k $k -M $M -N $N -K $K -o matmat
	     ./mmgen.exe -l $il -L $ol -m $m -n $n -k $k -M $M -N $N -K $K -o matmat
             make mmtime.exe DEBUG="" OPTS="-march=native -Ofast -DNDEBUG"

            ./mmtime.exe 2048
        echo ./mmgen.exe -l $il -L $ol -m $m -n $n -k $k -M $M -N $N -K $K -t -o matmat
	     ./mmgen.exe -l $il -L $ol -m $m -n $n -k $k -M $M -N $N -K $K -t -o matmat
             make mmtime.exe DEBUG="" OPTS="-march=native -Ofast -DNDEBUG"
            ./mmtime.exe 2048
	       end
            end
          end
        end	    	  
      end	    	  
    end	    	  
  end	    	  
 end	    	  
