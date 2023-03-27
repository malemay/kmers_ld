kmers_ld: main.c
	gcc -O2 -I ~/.local/include -o kmers_ld main.c

test: kmers_ld
	./kmers_ld pav_test.txt > kmers_ld_output.txt
