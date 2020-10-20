

rule all:
	input: 'plots/'
	output: touch('.status')

rule PDB_parser: 
        input: 'README.md'
        output: 'data/pdb_ID/pdb_entries_{halide}.txt'
        priority: 80
        shell: 
              'python3 scripts/PDB_parser.py -ian {wildcards.halide} -output {output}'
              
        
rule fetch_structures:
	input: "data/pdb_ID/pdb_entries_{halide}.txt"
	output: directory('data/structures/{halide}_struct')
	priority: 90
	shell: 
		'''
		while read i; do python scripts/fetch_structures.py -pdb_id $i -output {output} &>/dev/null; done < {input}
		i=$(echo {input})
		echo "=============="
		echo "$i: $(cat {input} | wc -l) structures were obtained."
		echo "=============="
		'''

rule filter:
	input: "data/structures/{halide}_struct"
	output: 
		info="data/info/info_{halide}.txt",
		filter="data/filtered_pdb_ID/filtered_{halide}.txt"
	priority: 10
	threads: 4
	shell: "for i in $(ls {input}); do python scripts/filter.py -input {input}/$i -input_type structure -output_info {output.info} -output_filtred {output.filter}; done"


rule get_context_data:
	input: 
		struct='data/structures/{halide}_struct',
		filter='data/filtered_pdb_ID/filtered_{halide}.txt'
	output: 
		one='data/context/{halide}_context.tsv',
		two='data/context/{halide}_context_AA.tsv'
	threads: 4
	shell: "python scripts/get_context.py -input_struct {input.struct} -input_filter {input.filter} -input_type structure -angstrem_radius 5 -output_full {output.one} &>/dev/null -output_AA {output.two} &>/dev/null -ligands y  -C 2 -prot y -sec_struct y -water y"


rule combine_final_data:
	input:	
		one=expand("data/context/{halide}_context.tsv", halide=['BR', 'CL', 'F', 'I']),
		two=expand("data/context/{halide}_context_AA.tsv", halide=['BR', 'CL', 'F', 'I'])
	output: 
                one="data/full_data.tsv",
	        two="data/full_data_AA.tsv"
	shell: 
			'''
			python scripts/parse_context.py -f {input.one} -input data/context/\*_context.tsv -output {output.one}
			cat {input.two} > {output.two}
			'''


rule figures_plotting:
	input: 
		one="data/full_data.tsv",
		two="data/full_data_AA.tsv"
	output: directory('plots/')
	shell: 'Rscript --vanilla scripts/figs.R {input.one} {input.two} {output}  '
