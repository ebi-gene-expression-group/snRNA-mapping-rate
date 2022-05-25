#!/usr/bin/env nextflow

sdrfFile = params.sdrf
resultsRoot = params.resultsRoot
referenceFasta = params.referenceFasta
referenceGtf = params.referenceGtf
species = params.species
transcriptToGene = params.transcriptToGene
transcriptomeIndex = params.transcriptomeIndex
protocol = params.protocol

manualDownloadFolder =''
if ( params.containsKey('manualDownloadFolder')){
    manualDownloadFolder = params.manualDownloadFolder
}

fastqProviderConfig = ''
if ( params.containsKey('fastqProviderConfig')){
    fastqProviderConfig = params.fastqProviderConfig
}

// Read ENA_RUN column from an SDRF

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

// TRANSCRIPT_TO_GENE = Channel.fromPath( transcriptToGene, checkIfExists: true ).first()
REFERENCE_GENOME = Channel.fromPath(referenceFasta, checkIfExists: true ).first()
REFERENCE_GTF = Channel.fromPath(referenceGtf, checkIfExists: true ).first()


// Read URIs from SDRF, generate target file names, and barcode locations

SDRF_FOR_FASTQS
    .map{ row-> 
      controlled_access='no'
      if (  params.fields.containsKey('controlled_access')){
        controlled_access=row["${params.fields.controlled_access}"]
      }
      tuple(row["${params.fields.run}"], row["${params.fields.cdna_uri}"], row["${params.fields.cell_barcode_uri}"], file(row["${params.fields.cdna_uri}"]).getName(), file(row["${params.fields.cell_barcode_uri}"]).getName(), row["${params.fields.cell_barcode_size}"], row["${params.fields.umi_barcode_size}"], row["${params.fields.end}"], row["${params.fields.cell_count}"], controlled_access) 
    }    
    .set { FASTQ_RUNS }


// generate splici index if it does not exits
process make_splici {
    // publishDir "t2g_alevin_fry/${species}", mode: 'copy', overwrite: true

    cache 'lenient'
   
    conda "${baseDir}/envs/pyroe.yml"

    memory { 20.GB * task.attempt }

    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10

    input:
        path referenceGenome from REFERENCE_GENOME
        path referenceGtf from REFERENCE_GTF

    output:
        path("splici_out/splici_fl*.fa") into SPLICI_FASTA_FOR_FRY
        path("splici_out/splici_fl*.tsv") into T2G_3_FOR_FRY

   
     
    """
    pyroe make-splici ${referenceGenome} ${referenceGtf} 90 splici_out
    """
}

process index_for_alevin_fry {
    // publishDir "index_alevin_fry/${species}", mode: 'copy', overwrite: true
    cache 'deep'
    memory { 20.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"

    input:
        path reference from SPLICI_FASTA_FOR_FRY
        
    output:
        path "alevin_index_splici" into ALEVIN_FRY_INDEX_SPLICI
    
  
    """
    salmon index --transcript ${reference}  -i alevin_index_splici
    """

 }


// Call the download script to retrieve run fastqs

process download_fastqs {
    
    conda "${baseDir}/envs/atlas-fastq-provider.yml"
    
    maxForks params.maxConcurrentDownloads
    time { 1.hour * task.attempt }
    memory { 20.GB * task.attempt }

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
    FINAL_FASTQS_FOR_ALEVIN
}

// Derive Alevin barcodeconfig

process alevin_config {

    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount) from FINAL_FASTQS_FOR_CONFIG

    output:
        set val(runId), stdout into ALEVIN_CONFIG
    
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
// run alevin-fry for quantification with splici index
 process alevin_fry_MR3 {
    publishDir "${resultsRoot}/${params.name}", mode: 'copy', overwrite: true
    cache 'lenient'
    memory { 20.GB * task.attempt }
    errorStrategy { task.exitStatus !=2 && (task.exitStatus == 130 || task.exitStatus == 137 || task.attempt < 3)  ? 'retry' : 'ignore' }
    maxRetries 10
    conda "${baseDir}/envs/alevin-fry_2.yml"
    input:
        set val(runId), file("cdna*.fastq.gz"), file("barcodes*.fastq.gz"), val(barcodeLength), val(umiLength), val(end), val(cellCount), val(barcodeConfig) from FINAL_FASTQS_FOR_ALEVIN.join(ALEVIN_CONFIG)
        path "alevin_index_splici" from ALEVIN_FRY_INDEX_SPLICI
        path "t2g_cDNA.txt" from T2G_3_FOR_FRY
       
    output:
        // publishDir path "${runId}_ALEVIN"
        set val(runId), path("${runId}_ALEVIN_fry_quant") into ALEVIN_FRY_RESULTS_SPLICI
      

    """
    salmon alevin ${barcodeConfig} --sketch -1 \$(ls barcodes*.fastq.gz | tr '\\n' ' ') -2 \$(ls cdna*.fastq.gz | tr '\\n' ' ') \
        -i alevin_index_splici -p ${task.cpus} -o ${runId}_ALEVIN_fry_map t2g_cDNA.txt 

    if (${barcodeConfig} == "--chromium")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/737K-august-2016.txt' --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    elif (${barcodeConfig} == "--chromiumV3")
    then
        alevin-fry generate-permit-list --unfiltered-pl '${baseDir}/whitelist/3M-february-2018.txt.gz' t --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    else
        alevin-fry generate-permit-list --input ${runId}_ALEVIN_fry_map -d fw --output-dir ${runId}_ALEVIN_fry_quant -k --min-reads 1
    fi

    alevin-fry collate -i ${runId}_ALEVIN_fry_quant -r ${runId}_ALEVIN_fry_map -t 16
    alevin-fry quant -i ${runId}_ALEVIN_fry_quant -m t2g_cDNA.txt -t 16 -r cr-like -o ${runId}_ALEVIN_fry_quant --use-mtx

        
    """
}

process mtx_alevin_fry_to_mtx {
    publishDir "${resultsRoot}/${name}", mode: 'copy', overwrite: true
    conda "/nfs/production/irene/ma/users/nnolte/conda/envs/parse_alevin_fry"

    // conda "${baseDir}/envs/parse_alevin_fry.yml"

    memory { 20.GB * task.attempt }
   

    input:
    set val(runId), path("${runId}_ALEVIN_fry_quant") from ALEVIN_FRY_RESULTS_SPLICI

    output:

    set val(runId), path("counts_mtx") into ALEVIN_FRY_MTX


    """
    alevinMtxTo10x.py --cell_prefix ${runId}- ${runId}_ALEVIN_fry_quant counts_mtx
    """      
}





ALEVIN_FRY_RESULTS_SPLICI
    .into{
        ALEVIN_RESULTS_FOR_QC
        ALEVIN_RESULTS_FOR_PROCESSING
        ALEVIN_RESULTS_FOR_OUTPUT
    }

// // Convert Alevin output to MTX. There will be one of these for every run, or
// // technical replicate group of runs


// ALEVIN_FRY_MTX
//     .into{
//         ALEVIN_MTX_FOR_QC
//         ALEVIN_MTX_FOR_EMPTYDROPS
//         ALEVIN_MTX_FOR_OUTPUT
//     }

// Make a diagnostic plot

// ALEVIN_RESULTS_FOR_QC
//     .join(ALEVIN_MTX_FOR_QC)
//     .set{
//         ALEVIN_QC_INPUTS
//     }

// process droplet_qc_plot{
    
//     conda "${baseDir}/envs/alevin.yml"
    
//     memory { 10.GB * task.attempt }
//     errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
//     maxRetries 20

//     input:
//         set val(runId), file(alevinResult), file(rawBarcodeFreq), file(mtx) from ALEVIN_QC_INPUTS

//     output:
//         set val(runId), file("${runId}.png") into ALEVIN_QC_PLOTS

//     """
//     dropletBarcodePlot.R $rawBarcodeFreq $mtx $runId ${runId}.png
//     """ 
// }

// // Remove empty droplets from Alevin results

// process remove_empty_drops {
    
//     conda "${baseDir}/envs/dropletutils.yml"

//     memory { 10.GB * task.attempt }
//     errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'ignore' }
//     maxRetries 20
   
//     input:
//         set val(runId), file(countsMtx) from ALEVIN_MTX_FOR_EMPTYDROPS

//     output:
//         set val(runId), file('nonempty.rds') into NONEMPTY_RDS

//     """
//         dropletutils-read-10x-counts.R -s counts_mtx -c TRUE -o matrix.rds
//         dropletutils-empty-drops.R -i matrix.rds --lower ${params.emptyDrops.lower} --niters ${params.emptyDrops.nIters} --filter-empty ${params.emptyDrops.filterEmpty} \
//             --filter-fdr ${params.emptyDrops.filterFdr} --ignore ${params.minCbFreq} -o nonempty.rds -t nonempty.txt
//     """
// }

// // Convert R matrix object with filtered cells back to .mtx

// process rds_to_mtx{

//     conda "${baseDir}/envs/dropletutils.yml"

//     memory { 10.GB * task.attempt }
//     errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
//     maxRetries 20
   
//     input:
//         set val(runId), file(rds) from NONEMPTY_RDS

//     output:
//         set val(runId), file("counts_mtx_nonempty") into NONEMPTY_MTX

//     """ 
//         #!/usr/bin/env Rscript
        
//         suppressPackageStartupMessages(require(DropletUtils))

//         counts_sce <- readRDS('$rds')
//         write10xCounts(assays(counts_sce)[[1]], path = 'counts_mtx_nonempty', barcodes = colData(counts_sce)\$Barcode, gene.id = rownames(counts_sce))
//     """
// }

// // Compile raw results with raw and emptyDrops-filtered MTX

// ALEVIN_RESULTS_FOR_OUTPUT
//     .join(ALEVIN_MTX_FOR_OUTPUT)
//     .join(NONEMPTY_MTX)
//     .join(ALEVIN_QC_PLOTS)
//     .set{ COMPILED_RESULTS }

// process compile_results{

//     publishDir "$resultsRoot/alevin", mode: 'copy', overwrite: true
    
//     input:
//         set val(runId), file('raw_alevin'), file(rawBarcodeFreq), file(countsMtx), file(countsMtxNonempty), file(qcPlot) from COMPILED_RESULTS

//     output:
//         set val(runId), file("$runId") into RESULTS_FOR_COUNTING

//     """
//         mkdir -p raw_alevin/alevin/mtx
//         cp -P $countsMtx $countsMtxNonempty raw_alevin/alevin/mtx 
//         mkdir -p raw_alevin/alevin/qc
//         cp -P $qcPlot raw_alevin/alevin/qc
//         cp -P raw_alevin $runId
//     """
// }

// // Check the total number of runs we have 

// RESULTS_FOR_COUNTING
//     .count()
//     .set{ ALEVIN_RESULTS_COUNT } 

// process validate_results {
    
//     executor 'local'
    
//     input:
//         val(kallistoResultCount) from ALEVIN_RESULTS_COUNT 
//         val(targetCount) from TARGET_RESULT_COUNT

//     output:
//         stdout DONE

//     """
//     if [ "$kallistoResultCount" -ne "$targetCount" ]; then
//         echo "Alevin results count of $kallistoResultCount does not match expected results number ($targetCount)" 1>&2
//         exit 1
//     else
//         echo "Alevin results count of $kallistoResultCount matches expected results number ($targetCount)"
//     fi
//     """
// }   
