
def metadata_builder(file, SRS_dict = {}, discrepancy = False):
	with open(file) as file:
		for line in file:
			if line[0] == '#':
				continue
			info = line.strip('\n').split('\t')
			if info[0] == 'sample_accession':
				continue
			SRS = info[0]
			if SRS not in SRS_dict:
				SRS_dict[SRS]={'SRR': [info[1]],
					    	  'paired':True if info[2]=='PAIRED' else False, 
					          'organism':info[3].replace(' ', '_'),
		            	      'tech':info[4],
						      'UMI':True if info[5]=='YES' else False,
							  'Study': info[6]}
	return(SRS_dict)

def cellranger_metadata_builder(file):
    SRS_dict={}
    with open(file) as file:
        for line in file:
            info = line.strip('\n').split('\t')
            SRS = info[0]
            if SRS not in SRS_dict:
                SRS_dict[SRS]={'biosample': info[1],'mouse': info[2],'path':info[3]}
    return(SRS_dict)

def lookup_run_from_SRS(SRS, fqp):
	SRR_files=SRS_dict[SRS]['SRR']
	out = []
	for SRR in SRR_files:
		if SRS_dict[SRS]['paired']:
			#PE
			out.append(f'{fqp}/fastq/{SRR}_1.fastq.gz')
			out.append(f'{fqp}/fastq/{SRR}_2.fastq.gz')
		else:
			#SE
			out.append(f'{fqp}/fastq/{SRR}.fastq.gz')
	return(out)

def lookupMapFromFaType(ft):
    if ft == 'cDNA_ONLY':
        return 'references/fa/tr2g.tsv'
    elif ft == 'gencode':
        return 'references/fa/gencode_tx2g.tsv'
    else:
        return 'dummy.yxt'


def get_mouse_fqp(srs, cdict):
    mouse = cdict[srs]['mouse']
    path = cdict[srs]['path']
    return f'bams/{mouse}/{path}/'

bustools_path = '/data/swamyvs/scEiad_DNTX/bustools'
fastq_path = '/data/OGVFB_BG/scEiaD'
whitelist = '10xv2.txt'
SRS_dict = metadata_builder('sampleTable.tsv')
sample_names = list(SRS_dict.keys())
cellranger_dict = cellranger_metadata_builder('cellrangeSampleTable.tsv')
cellranger_samples = list(cellranger_dict.keys())

rule all:
    input: 
        introns= expand('bus_out/cDNA_introns/{sample}/output.{whichCor}.mtx', sample = sample_names, whichCor= ['vanilla', 'corrected'] ),
        other_two = expand('bus_out/{faType}/{sample}/output.{whichCor}.mtx',sample = sample_names, faType=['gencode','cDNA_ONLY' ], whichCor= ['vanilla', 'corrected'] ),
        cellranger = expand('{srs}/', srs= cellranger_samples) ,
        alevin = expand('af_out/{sample}/alevin/quants_mat.mtx', sample  = sample_names)


rule download_annotation:
	output:
		mouse_anno='references/gtf/mm-mus_musculus_anno.gtf.gz',
		mouse_genome = 'references/genome/mm-mus_musculus_genome.fa.gz',
        mouse_anno_uz='references/gtf/mm-mus_musculus_anno.gtf',
		mouse_genome_uz = 'references/genome/mm-mus_musculus_genome.fa'
	shell:
		'''
		mkdir -p references
		wget -O {output.mouse_anno} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz 
        wget -O {output.mouse_genome} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz
        gunzip -c {output.mouse_genome} > {output.mouse_genome_uz}
        gunzip -c {output.mouse_anno} >{output.mouse_anno_uz}
		'''



rule make_gencode_fa:
    output:
        fa= 'references/fa/gencode.fa',
        tx2g = 'references/fa/gencode_tx2g.tsv'
    shell:
        '''
        module load R/4.0.3
        Rscript src/make_gencode_txfa.R
        '''

rule make_cellranger_index:
    input:
        mouse_anno_uz = 'references/gtf/mm-mus_musculus_anno.gtf',
        mouse_genome_uz = 'references/genome/mm-mus_musculus_genome.fa'
    output:
        directory('MusMusculus_genome')
    shell:
        '''
        module load cellranger/5.0.0
        cellranger mkref \
        --genome=MusMusculus_genome \
        --memgb 64 \
        --fasta=references/genome/mm-mus_musculus_genome.fa \
        --genes=references/gtf/mm-mus_musculus_anno.gtf
        '''




rule get_velocity_files:
	input: 
		gtf = 'references/gtf/mm-mus_musculus_anno.gtf.gz'
	output:
		'references/fa/cDNA_introns.fa',
        'references/fa/cDNA_ONLY.fa',
		'references/fa/tr2g.tsv'
	params:
		out_dir = 'references/fa/'
	shell:
		'''
		module load R/3.6
		Rscript src/get_velocity_annotation.R 
		'''


rule kallisto_index:
	input:
		'references/fa/{faType}.fa'
	output:
		'references/kallisto_idx/{faType}.idx'
	shell:
		"""
		module load kallisto/0.46.2
		kallisto index {input} -i {output}
		"""

rule salmon_index:
    input:
        fa= 'references/fa/gencode_gentrom.fa',
        decoys = 'references/decoys.txt'
    output:
        directory('references/salmon_index/gencode/')
    shell:
        '''
        module load salmon/1.4.0
        salmon index -t {input.fa} -d {input.decoys} -p 12 -i {output}
        '''



rule cellrange_count:
    input: 
        genome = 'MusMusculus_genome'
    output:
        results = directory('{SRS}/')
    params:
        fastq_path = lambda wildcards: get_mouse_fqp(wildcards.SRS, cellranger_dict)
    shell:
        '''
        cellranger count \
            --id={wildcards.SRS} \
            --transcriptome={input.genome} \
            --fastqs={params.fastq_path} \
            --localcores 4
        '''


rule rebuild_fastq:
    output:
        fq_l = 'fastq/{SRS}_1.fastq.gz',
        fq_r = 'fastq/{SRS}_2.fastq.gz'
    params:
        srs2fq =lambda wildcards: cellranger_dict[wildcards.SRS]['mouse']
    shell:
        '''
        find bams/{params.srs2fq} -name *_R1_*  | while read p ;
        do 
            gunzip -c $p  >> fastq/{wildcards.SRS}_1.fastq
        done 
        gzip fastq/{wildcards.SRS}_1.fastq

        find bams/{params.srs2fq} -name *_R2_*  | while read p ;
        do 
            gunzip -c $p  >> fastq/{wildcards.SRS}_2.fastq
        done 
        gzip fastq/{wildcards.SRS}_2.fastq
        '''

rule alevin_rad:
    input:
        fastq_l = 'fastq/{sample}_1.fastq.gz' ,
        fastq_r = 'fastq/{sample}_2.fastq.gz',
        idx = 'references/salmon_index/gencode/',
        whitelist = '10xv2.txt', 
        tx2g = 'references/fa/gencode_tx2g.tsv'
    output:
        'alevin_quant/{sample}/map.rad'
    params:
        outdir = lambda wildcards: f'alevin_quant/{wildcards.sample}/'
    shell:
        '''
        module load salmon/1.4.0 
        salmon alevin -l A -1 {input.fastq_l} -2 {input.fastq_r} -i {input.idx} --whitelist {input.whitelist}  --tgMap {input.tx2g} -o {outdir} -p 12 --chromium --rad 

        '''

rule alevin_fry:
    input:
        rad = 'alevin_quant/{sample}/map.rad',
        tx2g = 'references/fa/gencode_tx2g.tsv'
    output:
        'af_out/{sample}/alevin/quants_mat.mtx'
    params:
        idir = lambda wildcards: f'alevin_quant/{wildcards.sample}/',
        odir = lambda wildcards: f'af_out/{wildcards.sample}/'
    shell:
        '''        
        alevin-fry generate-permit-list --input {params.idir} --output-dir {params.idir} --expected-ori fw --knee-distance
        alevin-fry collate --rad-dir {params.idir} --input-dir {params.idir} --max-records 1000000000
        /data/swamyvs/alevin-fry/target/release/alevin-fry quant --input-dir {params.idir} --output-dir {params.odir}  --use-mtx  --resolution full --tg-map {input.tx2g}
        '''




rule kallisto_bus:
	input:
		fastq = [ 'fastq/{sample}_1.fastq.gz','fastq/{sample}_2.fastq.gz' ] ,
		idx = 'references/kallisto_idx/{faType}.idx'
	output:
		sorted_bus = 'quant/{faType}/{sample}/output.vanilla.bus',
        corrected_bus = 'quant/{faType}/{sample}/output.corrected.bus',
		ec = 'quant/{faType}/{sample}/matrix.ec',
		tx_name =  'quant/{faType}/{sample}/transcripts.txt'
	params:
		out_dir = lambda wildcards:  f'quant/{wildcards.faType}/{wildcards.sample}/'
	shell:
		'''
		module load kallisto/0.46.2
		kallisto bus -t 8 -x 10xv2 \
					-i {input.idx} -o {params.out_dir} {input.fastq}
        
        
        {bustools_path}/./bustools sort -t 8 -m 60G -o {output.sorted_bus} {params.out_dir}/output.bus 
        {bustools_path}/./bustools correct -w {whitelist} -o {output.corrected_bus} {output.sorted_bus} 
		'''

rule kallisto_count_noVelo:
    input:
        bus = 'quant/{faType}/{sample}/output.{whichCor}.bus',
        genemap = lambda wildcards: lookupMapFromFaType(wildcards.faType),
        ecmap = 'quant/{faType}/{sample}/matrix.ec',
        txname = 'quant/{faType}/{sample}/transcripts.txt'
    output:
        'bus_out/{faType}/{sample}/output.{whichCor}.mtx'
    params:
        outstem = lambda wildcards: f'bus_out/{wildcards.faType}/{wildcards.sample}/output.{wildcards.whichCor}'
    shell:
        '''
        {bustools_path}/./bustools count -g {input.genemap} -e {input.ecmap} -t {input.txname} --genecounts -o {params.outstem} {input.bus}
        '''
    #^ hot fixed instead of actually running 
rule kallisto_countvelo:
    input:
        bus ='quant/cDNA_introns/{sample}/output.{whichCor}.bus',
        genemap = 'references/fa/tr2g.tsv',
        targets = 'references/fa/introns_tx_to_capture.txt',
        ecmap =  'quant/cDNA_introns/{sample}/matrix.ec',
        txname = 'quant/cDNA_introns/{sample}/transcripts.txt'
    output:
        'bus_out/cDNA_introns/{sample}/output.{whichCor}.mtx'
    params:
        bus_out = lambda wildcards:  f'quant/cDNA_introns/{wildcards.sample}/',
        outstem = lambda wildcards: f'bus_out/cDNA_introns/{wildcards.sample}/output.{wildcards.whichCor}'
    shell:
        '''
        
        {bustools_path}/./bustools capture -s -x -o {params.bus_out}/TMP.spliced.bus -c {input.targets} -e {input.ecmap} -t {input.txname} {input.bus}
        
        {bustools_path}/./bustools count -g {input.genemap} -e {input.ecmap} -t {input.txname} --genecounts -o {params.outstem} {params.bus_out}/TMP.spliced.bus 
        '''



        
        



