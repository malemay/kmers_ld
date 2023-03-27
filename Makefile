kmers_ld: main.c
	gcc -O2 -I ~/.local/include -o kmers_ld main.c

test: kmers_ld
	./kmers_ld pod_color_blbr_pav_table.txt > kmers_ld_output.txt
