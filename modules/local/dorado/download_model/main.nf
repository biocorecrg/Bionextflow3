process DORADO_DOWNLOAD_MODEL {
    tag "$meta.id"
    label 'gpu'

    container "docker.io/ontresearch/dorado:sha00aa724a69ddc5f47d82bd413039f912fdaf4e77"

    input:
    tuple val(meta), path(pod5)
    
    output:
    path("dorado_models")        , type:'dir', emit: modelfolder
    tuple val("${task.process}"), val("dorado"), eval("dorado --version 2>&1 | head -n1 | cut -d '-' -f 2 | cut -d '+' -f 1"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
 
    def args = task.ext.args ?: ''
 
 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def down_pars = args.split(" ").find { it.contains('@') }
    def down_pars2 = args.trim().tokenize()[0]
  
    """
    mkdir dorado_models
    if dorado basecaller ${down_pars2} --max-reads 1 --models-directory \$PWD/dorado_models ./ > test.bam; 
        then
        	echo "Automatic model download succeeded"
    else 
        	echo "Trying the manual download...";
	        dorado download --model ${down_pars} --models-directory \$PWD/dorado_models
	fi
       
   """
}
