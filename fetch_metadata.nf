



include { GET_GSE_REPORT; GET_CSV_FROM_XML; ASSESS_LIBRARY_STRATEGY                         } from './modules/metadata-tasks.nf'


params.input_csv = "data/ribosome_profiling_superset.csv"


/// Workflow to get metadata. Returns a .csv for each GSE with info on their runs.
workflow metadata_flow {

    take: GSE_inputs

    main:
        gse_report_ch           = GET_GSE_REPORT            ( GSE_inputs )
        csv_ch                  = GET_CSV_FROM_XML          ( gse_report_ch )
        final_ch                = ASSESS_LIBRARY_STRATEGY   ( csv_ch )
}



workflow {

    GSE_inputs = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple("${row.Accession}", "${row.SRA}" )} 

    main:
        metadata_flow(GSE_inputs)

}