#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceGenome = params.referenceGenome
// referencecDNA = params.referencecDNA
referenceGtf = params.referenceGtf
protocol = params.protocol


// REFERENCE_CDNA = Channel.fromPath(referencecDNA,checkIfExists: true).first()
REFERENCE_GTF = Channel.fromPath( referenceGtf,checkIfExists: true ).first()
REFERENCE_GENOME = Channel.fromPath( referenceGenome,checkIfExists: true ).first()

manualDownloadFolder =''
if ( params.containsKey('manualDownloadFolder')){
    manualDownloadFolder = params.manualDownloadFolder
}

fastqProviderConfig = ''
if ( params.containsKey('fastqProviderConfig')){
    fastqProviderConfig = params.fastqProviderConfig
}
Channel
    .fromPath(sdrfFile, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .filter{ row -> (! row.containsKey(params.fields.quality)) || ( row["${params.fields.quality}"].toLowerCase() != 'not ok') }
    .into {
        SDRF_FOR_FASTQS
        SDRF_FOR_STRAND
        SDRF_FOR_TECHREP
        SDRF_FOR_COUNT
    }

//  Read URIs from SDRF, generate target file names, and barcode locations

SDRF_FOR_FASTQS
    .map{ row-> 
      controlled_access='no'
      if (  params.fields.containsKey('controlled_access')){
        controlled_access=row["${params.fields.controlled_access}"]
      }
      tuple(row["${params.fields.run}"], row["${params.fields.cdna_uri}"], row["${params.fields.cell_barcode_uri}"], file(row["${params.fields.cdna_uri}"]).getName(), file(row["${params.fields.cell_barcode_uri}"]).getName(), row["${params.fields.cell_barcode_size}"], row["${params.fields.umi_barcode_size}"], row["${params.fields.end}"], row["${params.fields.cell_count}"], controlled_access) 
    }    
    .set { FASTQ_RUNS }


 


process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 10.hour * task.attempt }
    memory { 20.GB * task.attempt }
    cache 'lenient'

    errorStrategy { task.attempt<=10 & task.exitStatus != 4 ? 'retry' : 'finish' } 
    
    input:
        set runId, cdnaFastqURI, barcodesFastqURI, cdnaFastqFile, barcodesFastqFile, val(barcodeLength), val(umiLength), val(end), val(cellCount), val(controlledAccess) from FASTQ_RUNS

    output:
        set val(runId), file("${cdnaFastqFile}"), file("${barcodesFastqFile}"), val(barcodeLength), val(umiLength), val(end), val(cellCount) into DOWNLOADED_FASTQS

    """
        if [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
           ln -s $manualDownloadFolder/${cdnaFastqFile} ${cdnaFastqFile}
           ln -s $manualDownloadFolder/${barcodesFastqFile} ${barcodesFastqFile}
        elif [ -n "$manualDownloadFolder" ] && [ -e $manualDownloadFolder/${cdnaFastqFile} ] && [ ! -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is available locally, but barcodes file $barcodesFastqFile is not 1>&2
            exit 2    
        elif [ -n "$manualDownloadFolder" ] && [ ! -e $manualDownloadFolder/${cdnaFastqFile} ] && [ -e $manualDownloadFolder/${barcodesFastqFile} ]; then
            echo 'cDNA file $cdnaFastqFile is not available locally, but barcodes file $barcodesFastqFile is 1>&2
            exit 3 
        elif [ "$controlledAccess" = 'yes' ]; then
            echo "One or both of ${cdnaFastqFile}, ${barcodesFastqFile} are not available at $manualDownloadFolder/ for this controlled access experiment" 1>&2
            exit 4   
        else
            confPart=''
            if [ -n "$fastqProviderConfig" ] && [ -e "$fastqProviderConfig" ]; then
                confPart=" -c $fastqProviderConfig"
            fi 
            # Stop fastq downloader from testing different methods -assume the control workflow has done that 
            export NOPROBE=1
        
            fetchFastq.sh -f ${cdnaFastqURI} -t ${cdnaFastqFile} -m ${params.downloadMethod} \$confPart
            
            # Allow for the first download also having produced the second output already
            if [ ! -e ${barcodesFastqFile} ]; then
                fetchFastq.sh -f ${barcodesFastqURI} -t ${barcodesFastqFile} -m ${params.downloadMethod} \$confPart
            fi
        fi
    """
}

// Group read files by run name, or by technical replicate group if specified

if ( params.fields.containsKey('techrep')){

    // If technical replicates are present, create a channel containing that info 

    SDRF_FOR_TECHREP
        .map{ row-> tuple(row["${params.fields.run}"], row["${params.fields.techrep}"]) }
        .groupTuple()
        .map{ row-> tuple( row[0], row[1][0]) }
        .set{ TECHREPS }

    // The target set of results will now be the technical replicate group number

    SDRF_FOR_COUNT
        .map{ row-> tuple(row["${params.fields.techrep}"]) }
        .unique()
        .count()
        .set { TARGET_RESULT_COUNT }
    
    // Now add the tech rep group to the run info, group by it, and create a
    // tuple of files keyed by techrep group

    TECHREPS.join( DOWNLOADED_FASTQS )
        .groupTuple(by: 1)
        .map{ row-> tuple( row[1], row[2].flatten(), row[3].flatten(), row[4][0], row[5][0], row[6][0], row[7][0]) }
        .set{
            FINAL_FASTQS
        }
}else{
    DOWNLOADED_FASTQS.set{ FINAL_FASTQS }
    
    SDRF_FOR_COUNT
      .map{ row-> tuple(row["${params.fields.run}"]) }
      .unique()
      .count()
      .set { TARGET_RESULT_COUNT }
}

FINAL_FASTQS.into{
    FINAL_FASTQS_FOR_CONFIG
    FINAL_FASTQS_FOR_ALEVIN_SPLICI
    FINAL_FASTQS_FOR_ALEVIN_CDNA
    FINAL_FASTQS_FOR_STAR
    FINAL_FASTQS_FOR_KB_TOOLS
    FINAL_FASTQS_FOR_KB_TOOLS_SPLICI
    FINAL_FASTQS_FOR_KB_TOOLS_PRERNA
    FINAL_FASTQS_FOR_ALEVIN_FRY
    FINAL_FASTQS_FOR_ALEVIN_FRY_CDNA
    FINAL_FASTQS_FOR_ALEVIN_FRY_TRANSCRIPTOME
}


process alevin_config {
    cache 'lenient'
    memory { 20.GB * task.attempt }
    cpus 4

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

    output:
        set val(runId), stdout into CONFIG
    
    script:

        def barcodeConfig = ''

        if ( params.containsKey(protocol) ){

            canonicalProtocol = params.get(protocol)
            alevinType = canonicalProtocol.alevinType

            // Non-standard barcode config is supplied as a custom method

            if ( alevinType == 'custom' || "${canonicalProtocol.barcodeLength}" != barcodeLength || "${canonicalProtocol.umiLength}" != umiLength || "${canonicalProtocol.end}" != end ){
                barcodeConfig = "--barcodeLength ${barcodeLength} --umiLength ${umiLength} --end ${end}" 

            }else{
                barcodeConfig = "--$alevinType"
            }
            barcodeConfig = "-l ${canonicalProtocol.libType} $barcodeConfig" 
        }

        """
        if [ -z "$barcodeConfig" ]; then
            echo Input of $protocol results is misconfigured 1>&2
            exit 1
        fi
        # Also check barcode read lengths and return non-0 if they're not what they should be
        targetLen=\$(($umiLength + $barcodeLength))
        barcodesGood=0
        set +e
        while read -r l; do
            checkBarcodeRead.sh -r \$(readlink -f \$l) -b $barcodeLength -u $umiLength -n 1000000 1>&2
            if [ \$? -ne 0 ]; then
                barcodesGood=1
            fi
        done <<< "\$(ls barcodes*.fastq.gz)"
        set -e
        
        echo -n "$barcodeConfig"
        exit \$barcodesGood
        """
}

CONFIG
    .into{
        ALEVIN_CONFIG
        ALEVIN_CONFIG_SPLICI
        ALEVIN_CONFIG_CDNA
        STAR_CONFIG
        KB_CONFIG
        KB_CONFIG_SPLICI
        KB_CONFIG_PRERNA
        ALEVIN_FRY_CONFIG
        ALEVIN_FRY_CONFIG_CDNA
        ALEVIN_FRY_CONFIG_TRANSCRIPTOME
    }



// build index to runSTARSolo
process index_star {
    cache 'lenient'
    memory { 100.GB * task.attempt }
    cpus 4
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/star.yml"

    input:
        path(referenceGenome) from REFERENCE_GENOME
        path(referenceGtf) from REFERENCE_GTF
    output:
        path("STAR_index") into STAR_INDEX
    
    """
    STAR --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${referenceGtf} --genomeSAindexNbases 12 --limitGenomeGenerateRAM 188632691637
                                                                                                                                                                               
    """
}

// run STARSolo 

process run_STARSolo {
   
    cache 'lenient'

    memory { 100.GB * task.attempt }
    cpus 10
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10


    conda "${baseDir}/envs/star.yml"


    input:
    set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_STAR.join(STAR_CONFIG)
    path("STAR_index") from STAR_INDEX

    output:
    set val(runId), path("${runId}_STAR_tmpSolo.out") into STAR_RESULTS
    set val(runId), path(".command.log")  into  MEM_STAR
    // bam into BAM_FILE
    

    script:
    if( barcodeConfig == '10xv3' )
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --soloCBwhitelist '${baseDir}/whitelist/3M-february-2018.txt.gz' --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --readFilesCommand zcat --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out |\
         awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"

        
        
        """
    else if( barcodeConfig == '10xv2' )
  
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --soloCBwhitelist '${baseDir}/whitelist/737K-august-2016.txt' --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --readFilesCommand zcat --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out |\
         awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"

        
        
        """
    else
        """
        STAR --genomeDir STAR_index --readFilesIn \$(ls cdna*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') \$(ls barcodes*.fastq.gz | tr '\\n' ', '| sed  's/,*\$//g') --soloType Droplet --readFilesCommand zcat --soloCBwhitelist None --soloUMIlen ${umiLength} --soloCBlen ${barcodeLength} --soloUMIstart \$(($barcodeLength+1)) --soloCBstart 1 --runThreadN 12 --soloFeatures Gene GeneFull --outFileNamePrefix ${runId}_STAR_tmp --soloBarcodeReadLength 0

        mapping_rate=\$(grep "Uniquely mapped reads %" ${runId}_STAR_tmpLog.final.out | \
        awk '{split(\$0, array, "|"); print array[2]}')
        echo  "\${mapping_rate}"


        
        
        """
}

// make bowtie index


// 
