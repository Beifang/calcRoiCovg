all:
    ifndef SAMTOOLS_ROOT
	@echo "Please define environment variable SAMDIR to point to your samtools libraries"
    else
	gcc -g -Wall -fopenmp -O2 -I${SAMTOOLS_ROOT} calcRoiCovg.c -o calcRoiCovg -L${SAMTOOLS_ROOT} -lbam -lm -lz -lpthread
    endif
clean:
	rm -f calcRoiCovg

