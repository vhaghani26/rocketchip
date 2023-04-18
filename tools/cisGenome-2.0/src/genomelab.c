/* ----------------------------------------------------------------------- */
/*  motifsampler.c : Motif discovery through Gibbs sampler algorithm.      */
/*  Author : Ji HongKai ; Time: 2004.07                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "StringLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "MathLib.h"
#include "MotifLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "MicroarrayLib.h"
#include "AffyLib.h"
#include "Alignment.h"
#include "PhysicalNetworkLib.h"
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"
#include "WorkLib.h"

/* -------------- */
/* test functions */
/* -------------- */
int testcoding();
int testloading();
int testnetwork();
int testnormcdf();
int testcel();
int testbar();
int testgeneric();
int testeigen();
int menu_test_chipcluster();
int shhmbtest(int argv, char **argc);

/* -------------- */
/* menu available */
/* -------------- */
int menu_file_cod2bed(int argv, char **argc);
int menu_file_bed2cod(int argv, char **argc);
int menu_file_bed2aln(int argv, char **argc);

/* int menu_codinggenome(); */
int menu_codinggenome(int argv, char **argc);
/* int menu_getseqfromgenome(); */
int menu_getseqfromgenome(int argv, char **argc);
int menu_getseqcsfromgenome(int argv, char **argc);
int menu_getmaskedseqfromgenome(int argv, char **argc);
int menu_getmaskedregionfromgenome(int argv, char **argc);
int menu_getseqfromfasta(int argv, char **argc);
int menu_fasta_soft2hardmask(int argv, char **argc);

int menu_getseqfromgenome_c(int argv, char **argc);
int menu_getseqcsfromgenome_c(int argv, char **argc);
int menu_getmaskedseqfromgenome_c(int argv, char **argc);
int menu_getmaskedregionfromgenome_c(int argv, char **argc);
int menu_seqmask(int argv, char **argc);
int menu_genome_hash(int argv, char **argc);

int menu_codingphastcons(int argv, char **argc);
int menu_codingphastcons_v2(int argv, char **argc);
int menu_maf_conserve_bg(int argv, char **argc);
int menu_maf_conserve_cs(int argv, char **argc);
int menu_maf_footprint(int argv, char **argc);
int menu_cs_getdistn(int argv, char **argc);
int menu_codingCDS(int argv, char **argc);
int menu_csgc_summary(int argv, char **argc);
int menu_regioncs_summary(int argv, char **argc);
int menu_genome_regionextend(int argv, char **argc);

int menu_refgene_sort(int argv, char **argc);
int menu_refgene_gettargettssaround(int argv, char **argc);
int menu_refgene_getortholog(int argv, char **argc);
int menu_refgene_getmultiortholog(int argv, char **argc);
int menu_refgene_getmultiortholog1way(int argv, char **argc);
int menu_refgene_getmultiorthologtssaround(int argv, char **argc);
int menu_refgene_getmultiorthologtssaroundexonmasked(int argv, char **argc);
int menu_refgene_pickspeciesspecific(int argv, char **argc);
int menu_refgene_getnearestgene(int argv, char **argc);
int menu_refgene_getaffy(int argv, char **argc);
int menu_refgene_getmatchedcontrol(int argv, char **argc);
int menu_refgene_getlocationsummary(int argv, char **argc);

int menu_refflat_sort(int argv, char **argc);
int menu_refexon_createmap(int argv, char **argc);
int menu_refmicroarray_createmap(int argv, char **argc);
int menu_reflocus_createmap(int argv, char **argc);
int menu_reflocus_assignvalue(int argv, char **argc);
int menu_reflocus_getneighborgenes(int argv, char **argc);
int menu_reflocus_assignvalue_test(int argv, char **argc);

int menu_refflex_getmultiortholog(int argv, char **argc);
int menu_refflex_getmultiortholog1way(int argv, char **argc);

int menu_affy_getreducedannot();
int menu_affy_getrandomcontrolprobesets();
int menu_affy_loadbar_intensity();
int menu_affy_loadbar_intensity_group();
int menu_affy_loadbpmap_intensity();
int menu_affy_bpmap_filter_gtrans();
int menu_affy_bpmap_filter_ori();
int menu_affy_bpmap_filter_data();
int menu_affy_bar2txt(int argv, char **argc);
int menu_affy_bar2wig(int argv, char **argc);

int menu_network_shortestpath(int argv, char **argc);
int menu_network_createortholognet(int argv, char **argc);
int menu_network_getsubnet(int argv, char **argc);
int menu_network_bind2entrez(int argv, char **argc);
int menu_reflocus_assignedgevalue(int argv, char **argc);

int menu_transloc(int argv, char **argc);

int menu_expression_geneselection(int argv, char **argc);
int menu_expression_rankgenebylocuslink(int argv, char **argc);
int menu_expression_quantilenormalization();
int menu_expression_getspecificprobe(int argv, char **argc);
int menu_expression_getnrprobe(int argv, char **argc);

int menu_tiling_probeselection();
int menu_tiling_ums_fdr(int argv, char **argc);
int menu_tiling_bindingcall_hmm();
int menu_tiling_bindingcall_hmm_explen();
int menu_tiling_bindingcall_hmm_constlen();
int menu_tiling_bindingcall_hmm_baumwelch();

int menu_tilemap_importaffy(int argv, char **argc);
int menu_tilemap_normalization(int argv, char **argc);
int menu_tilemap(int argv, char **argc);
int menu_tilemap_extract(int argv, char **argc);
int menu_tilemap_proberemap(int argv, char **argc);
int menu_tilemapv2_importaffy(int argv, char **argc);
int menu_tilemapv2(int argv, char **argc);
int menu_tilemapv2_regioninfo(int argv, char **argc);
int menu_tilemapv2_regioninfo_integral(int argv, char **argc);
int menu_tilemapv2_txt2bar(int argv, char **argc);
int menu_tilemapv2_quantilenorm(int argv, char **argc);
int menu_tilemapv2_collectprobes(int argv, char **argc);
int menu_tileprobe_buildmodel(int argv, char **argc);
int menu_tileprobe_buildmodel_v2(int argv, char **argc);
int menu_tileprobe_buildmodel_hmmt(int argv, char **argc);
int menu_tileprobe_buildmodel_hmmb(int argv, char **argc);
int menu_tileprobe_buildmodel_hmmm(int argv, char **argc);
int menu_tileprobe_buildmodel_hmmw(int argv, char **argc);
int menu_tileprobe_norm(int argv, char **argc);
int menu_tileprobe_peak(int argv, char **argc);
int menu_tileprobe_mat(int argv, char **argc);
int menu_tileprobe_buildbarmodel(int argv, char **argc);
int menu_tileprobe_buildbarmodel_v2(int argv, char **argc);
int menu_tileprobe_barnorm(int argv, char **argc);
int menu_tileprobe_huberdata(int argv, char **argc);
int menu_tileprobe_resample(int argv, char **argc);

int menu_hts_aln2uniq(int argv, char **argc);
int menu_hts_aln2window(int argv, char **argc);
int menu_hts_aln2diff(int argv, char **argc);
int menu_hts_aln2bar(int argv, char **argc);
int menu_hts_aln2barv2(int argv, char **argc);
int menu_hts_aln2winbar(int argv, char **argc);
int menu_hts_alnshift2bar(int argv, char **argc);
int menu_hts_windowsummary(int argv, char **argc);
int menu_hts_windowsummaryv2(int argv, char **argc);
int menu_hts_onesample_enrich(int argv, char **argc);
int menu_hts_onesamplev2_enrich(int argv, char **argc);
int menu_hts_twosample_windowsummary(int argv, char **argc);
int menu_hts_twosample_windowsummaryv2(int argv, char **argc);
int menu_hts_twosample_enrich(int argv, char **argc);
int menu_hts_twosamplev2_enrich(int argv, char **argc);
int menu_hts_windowsummarypaper(int argv, char **argc);
int menu_hts_createrepeatfilter(int argv, char **argc);
int menu_hts_filterrepeatreads(int argv, char **argc);
int menu_hts_collectreads(int argv, char **argc);
int menu_hts_collectprofile(int argv, char **argc);
int menu_hts_countreads4refgene(int argv, char **argc);
int menu_hts_selectreads(int argv, char **argc);

int menu_seqclust_seg(int argv, char **argc);
int menu_seqclust_count(int argv, char **argc);
int menu_seqclust(int argv, char **argc);
int menu_seqclust_dp(int argv, char **argc);

int menu_seqpeak(int argv, char **argc);

int menu_cnv_aln2window(int argv, char **argc);
int menu_cnv_repeat2window(int argv, char **argc);

int menu_rnaseq_countgeneread(int argv, char **argc);

int menu_flexmodule(int argv, char **argc);
int menu_flexmodule_bgfit(int argv, char **argc);
int menu_flexmodule_extractmotif(int argv, char **argc);

int menu_motifmap_consensusscan(int argv, char **argc);
int menu_motifmap_consensusscan_genome(int argv, char **argc);
int menu_motifmap_matrixscan(int argv, char **argc);
int menu_motifmap_matrixscan_genomebackground(int argv, char **argc);
int menu_motifmap_matrixscan_genome(int argv, char **argc);
int menu_motifmap_matrixintegrate_genome(int argv, char **argc);
int menu_motifmap_multiplematrixscan_genome(int argv, char **argc);
int menu_motifmap_matrixscan_group(int argv, char **argc);
int menu_motifmap_getsitearound(int argv, char **argc);
int menu_motifmap_countkmer(int argv, char **argc);
int menu_motifmap_rescore(int argv, char **argc);
int menu_motifmap_matrixscan_genome_enrich(int argv, char **argc);
int menu_motifmap_consensusscan_genome_enrich(int argv, char **argc);
int menu_motifmap_matrixscan_genome_getcutoff(int argv, char **argc);
int menu_motifmap_matrixscan_genome_summary(int argv, char **argc);
int menu_motifmap_consensusscan_genome_summary(int argv, char **argc);
int menu_motifmap_filter_genome(int argv, char **argc);
int menu_motifmap_getsitearoundcs_genome(int argv, char **argc);
int menu_motifmap_getcluster(int argv, char **argc);

int menu_motif_simuscoredistn_typeI();
int menu_motif_simuscoredistn_typeII();
int menu_motif_simuscoredistn_typeIII();
int menu_motif_getcutoff_typeI();
int menu_motif_getcutoff_typeII();
int menu_motif_getcutoff_typeIII();
int menu_motifsampler_quality();

int menu_malign_blasthit(int argv, char **argc);
int menu_malign_motifmap(int argv, char **argc);
int menu_malign_genome_prepareortholog(int argv, char **argc);
int menu_malign_genome_blasthit(int argv, char **argc);
int menu_malign_modulemap(int argv, char **argc);
int menu_malign_countkmer(int argv, char **argc);
int menu_malign_generatemotifaln(int argv, char **argc);

int menu_mapcMyc();


/* -------------- */
/* batch task     */
/* -------------- */
int pipeline_knownmotifmapping_forgenelist();

/* main function */
/* int main() */
int main(int argv, char **argc)
{
	int nLen;
	int nseed;
	double dTemp;
	
	/* init rand */
	srand( (unsigned)time( NULL ) );

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nseed = (int)(rand()*1000/RAND_MAX);
	}
	else
	{
		nseed = rand()%1000;
	}
	
	
	rand_u_init(nseed);

	
	/* if(argv != 5)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	Expression_GeneSelection_Main(argc[1], argc[2], argc[3], argc[4]); 
	*/

	
	/* ---- */
	/* menu */
	/* ---- */
	/* menu_file_bed2cod(argv, argc); */
	/* menu_file_cod2bed(argv, argc); */
	/* menu_file_bed2aln(argv, argc); */

	/* menu_codinggenome(argv, argc); */
	/* menu_getseqfromgenome(argv, argc); */
	/* menu_codingphastcons(argv, argc); */
	/* menu_codingphastcons_v2(argv, argc); */
	/* menu_maf_conserve_bg(argv, argc); */
	/* menu_maf_conserve_cs(argv, argc); */
	/* menu_maf_footprint(argv, argc); */
	/* menu_cs_getdistn(argv, argc); */
	/* menu_getseqcsfromgenome(argv, argc); */
	/* menu_codingCDS(argv, argc); */
	/* menu_csgc_summary(argv, argc); */
	/* menu_getmaskedseqfromgenome(argv, argc); */
	/* menu_getmaskedregionfromgenome(argv, argc); */
	/* menu_genome_hash(argv, argc); */
	/* menu_regioncs_summary(argv, argc); */
	/* menu_genome_regionextend(argv, argc); */

	/* menu_getseqfromfasta(argv, argc); */

	/* menu_getseqfromgenome_c(argv, argc); */
	/* menu_getseqcsfromgenome_c(argv, argc); */
	/* menu_getmaskedseqfromgenome_c(argv, argc); */
	/* menu_getmaskedregionfromgenome_c(argv, argc); */
	/* menu_seqmask(argv, argc); */
	/* menu_fasta_soft2hardmask(argv, argc); */

	/* menu_refgene_sort(argv, argc); */
	/* menu_refgene_gettargettssaround(argv, argc); */
	/* menu_refgene_getortholog(argv, argc); */
	/* menu_refgene_getmultiortholog(argv, argc); */
	/* menu_refgene_getmultiortholog1way(argv, argc); */
	/* menu_refgene_getmultiorthologtssaround(argv, argc); */
	/* menu_refgene_getmultiorthologtssaroundexonmasked(argv, argc); */
	/* menu_refgene_pickspeciesspecific(argv, argc); */
	/* menu_refgene_getnearestgene(argv, argc); */
	/* menu_refflex_getmultiortholog(argv, argc); */
	/* menu_refflex_getmultiortholog1way(argv, argc); */
	/* menu_refgene_getaffy(argv, argc); */
	/* menu_refgene_getmatchedcontrol(argv, argc); */
	/* menu_reflocus_assignvalue(argv, argc); */
	/* menu_reflocus_getneighborgenes(argv, argc); */
	/* menu_refgene_getlocationsummary(argv, argc); */

	/* menu_refflat_sort(argv, argc); */
	/* menu_reflocus_createmap(argv, argc); */
	/* menu_refexon_createmap(argv, argc); */
	/* menu_refmicroarray_createmap(argv, argc); */

	/* menu_affy_getreducedannot(); */
	/* menu_affy_getrandomcontrolprobesets(); */
	/* menu_affy_loadbar_intensity(); */ 
	/* menu_affy_loadbar_intensity_group(); */
	/* menu_affy_loadbpmap_intensity(); */
	/* menu_affy_bpmap_filter_gtrans(); */
	/* menu_affy_bpmap_filter_ori(); */
	/* menu_affy_bar2txt(argv, argc); */
	/* menu_affy_bar2wig(argv, argc); */

	/* menu_network_shortestpath(argv, argc); */
	/* menu_network_createortholognet(argv, argc); */
	/* menu_network_getsubnet(argv, argc); */
	/* menu_network_bind2entrez(argv, argc); */
	/* menu_reflocus_assignedgevalue(argv, argc); */


	/* menu_transloc(argv, argc); */

	/* menu_shhinsitu(); */
	/* menu_expression_geneselection(argv, argc); */
	/* menu_expression_rankgenebylocuslink(argv, argc); */
	/* menu_expression_quantilenormalization(); */
	/* menu_expression_getspecificprobe(argv, argc); */
	/* menu_expression_getnrprobe(argv, argc); */
	
	/* menu_affychr2122_probematch(); */
	/* menu_tiling_probeselection(); */
	/* menu_tiling_ums_fdr(argv, argc); */
	/* menu_tiling_bindingcall_hmm(); */
	/* menu_affy_bpmap_filter_data(); */
	/* menu_tiling_bindingcall_hmm_explen(); */
	/* menu_tiling_bindingcall_hmm_constlen(); */
	/* menu_tiling_bindingcall_hmm_baumwelch(); */

	/* menu_tilemap_importaffy(argv, argc); */
	/* menu_tilemap_normalization(argv, argc); */
	/* menu_tilemap(argv, argc); */
	/* menu_tilemap_extract(argv, argc); */
	/* menu_tilemap_proberemap(argv, argc); */
	/* menu_tilemapv2_importaffy(argv, argc); */
	/* menu_tilemapv2(argv, argc); */
	/* menu_tilemapv2_regioninfo(argv, argc); */
	/* menu_tilemapv2_regioninfo_integral(argv, argc); */
	/* menu_tilemapv2_txt2bar(argv, argc); */
	/* menu_tilemapv2_quantilenorm(argv, argc); */
	/* menu_tilemapv2_collectprobes(argv, argc); */
	/* menu_tileprobe_buildmodel(argv, argc); */
	/* menu_tileprobe_buildmodel_v2(argv, argc); */
	/* menu_tileprobe_buildmodel_hmmt(argv, argc); */
	/* menu_tileprobe_buildmodel_hmmb(argv, argc); */
	/* menu_tileprobe_buildmodel_hmmm(argv, argc); */
	/* menu_tileprobe_buildmodel_hmmw(argv, argc); */
	/* menu_tileprobe_norm(argv, argc); */
	/* menu_tileprobe_peak(argv, argc); */
	/* menu_tileprobe_mat(argv, argc); */
	/* menu_tileprobe_buildbarmodel(argv, argc); */
	/* menu_tileprobe_buildbarmodel_v2(argv, argc); */
	/* menu_tileprobe_barnorm(argv, argc); */
	/* menu_tileprobe_huberdata(argv, argc); */
	/* menu_tileprobe_resample(argv, argc); */

	/* menu_hts_aln2uniq(argv, argc); */
	/* menu_hts_aln2window(argv, argc); */
	/* menu_hts_aln2diff(argv, argc); */
	/* menu_hts_aln2bar(argv, argc); */
	/* menu_hts_aln2barv2(argv, argc); */
	/* menu_hts_aln2winbar(argv, argc); */ 
	/* menu_hts_alnshift2bar(argv, argc); */
	/* menu_hts_windowsummary(argv, argc); */
	/* menu_hts_windowsummaryv2(argv, argc); */
	/* menu_hts_onesample_enrich(argv, argc); */
	/* menu_hts_onesamplev2_enrich(argv, argc); */
	/* menu_hts_twosample_windowsummary(argv, argc); */
	/* menu_hts_twosample_windowsummaryv2(argv, argc); */
	/* menu_hts_twosample_enrich(argv, argc); */
	/* menu_hts_twosamplev2_enrich(argv, argc); */
	/* menu_hts_windowsummarypaper(argv, argc); */
	/* menu_hts_createrepeatfilter(argv, argc); */
	/* menu_hts_filterrepeatreads(argv, argc); */
	/* menu_hts_collectreads(argv, argc); */
	/* menu_hts_collectprofile(argv, argc); */
	/* menu_hts_countreads4refgene(argv, argc); */

	/* menu_cnv_aln2window(argv, argc); */
	/* menu_cnv_repeat2window(argv, argc); */

	/* menu_rnaseq_countgeneread(argv, argc); */

	/* menu_seqclust_seg(argv, argc); */
	/* menu_seqclust_count(argv, argc); */
	/* menu_seqclust(argv, argc); */
	/* menu_seqclust_dp(argv, argc); */

	/* menu_seqpeak(argv, argc); */

	/* menu_mapcMyc(); */
	/* Count_cMyc_In_TargetRegion_Main(); */
	/* Map_cMyc_Chr2122_Conserved(); */
	/* Count_cMyc_In_ConservedTargetRegion_Main(); */

	/* menu_flexmodule(argv, argc); */
	/* menu_flexmodule_bgfit(argv, argc); */
	/* menu_flexmodule_extractmotif(argv, argc); */

	/* menu_motifmap_consensusscan(argv, argc);*/
	/* menu_motifmap_consensusscan_genome(argv, argc); */
	/* menu_motifmap_matrixscan(argv, argc); */
	/* menu_motifmap_matrixscan_genomebackground(argv, argc); */
	/* menu_motifmap_matrixscan_genome(argv, argc); */
	/* menu_motifmap_matrixintegrate_genome(argv, argc); */
	/* menu_motifmap_multiplematrixscan_genome(argv, argc); */
	/* menu_motifmap_matrixscan_group(argv, argc); */
	/* menu_motifmap_getsitearound(argv, argc); */
	/* menu_motifmap_countkmer(argv, argc); */
	/* menu_motifmap_rescore(argv, argc); */
	/* menu_motifmap_matrixscan_genome_enrich(argv, argc); */
	/* menu_motifmap_consensusscan_genome_enrich(argv, argc); */
	/* menu_motifmap_matrixscan_genome_summary(argv, argc); */
	menu_motifmap_consensusscan_genome_summary(argv, argc);
	/* menu_motifmap_filter_genome(argv, argc); */
	/* menu_motifmap_getsitearoundcs_genome(argv, argc); */
	/* menu_motifmap_matrixscan_genome_getcutoff(argv, argc); */
	/* menu_motifmap_getcluster(argv, argc); */
	/* MotifMap_FilterOverlappingSite("motifs_enrich.txt.1_target_c40.map", "motifs_enrich.txt.1_target_c40.nrmap"); */

	/* menu_motif_simuscoredistn_typeI(); */
	/* menu_motif_simuscoredistn_typeII(); */
	/* menu_motif_simuscoredistn_typeIII(); */
	/* menu_motif_getcutoff_typeI(); */
	/* menu_motif_getcutoff_typeII(); */
	/* menu_motif_getcutoff_typeIII(); */ 
	/* menu_motifsampler_quality(); */

	/* menu_malign_genome_prepareortholog(argv, argc); */
	/* menu_malign_blasthit(argv, argc); */
	/* menu_malign_motifmap(argv, argc); */
	/* menu_malign_modulemap(argv, argc); */
	/* menu_malign_genome_blasthit(argv, argc); */
	/* menu_malign_countkmer(argv, argc); */
	/* menu_malign_generatemotifaln(argv, argc); */

	/* menu_geneset(argv, argc); */

	/* --------- */
	/* pipelines */
	/* --------- */
	/* pipeline_knownmotifmapping_forgenelist(); */

	/* ----------------- */
	/* project workspace */
	/* ----------------- */
	/* menu_MapRPL("C:\\Projects\\research_harvard\\chimp_project\\analysis\\chr4followup\\RPL\\RPL19_blat.txt",
		"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\annotation\\refGene_sorted.txt",
		"C:\\Projects\\research_harvard\\chimp_project\\analysis\\chr4followup\\RPL\\RPL19_map.txt", "human");
	*/

	/* GetEBChIPCod(); */
	/* ConnectSitewithOrtho(argv, argc); */
	/* Oct4GetPromoter(argv, argc); */
	/* GenomeGetConservedSeg_Main(argv, argc); */
	/* HG17GetTSS_Main(argv, argc); */
	/* MM6GetGeneCover_Main(argv, argc); */
	/* EEL_PrepareOrtholog_Main(argv, argc); */
	/* EEL_GetReverseComplement_Main(); */
    /* GliArrayGenerateRand(); */

	/* MapHomoloGene(); */
	/* MapHomoloGene_ByHID(); */
	/* MapGeneID2RefGene(); */
	/* testnetwork(); */
	/* testnormcdf(); */
	/* testcel(); */
	/* testbar(); */
	/* testgeneric(); */
	/* testeigen(); */
	/* menu_test_chipcluster(); */
	/* Feinberg_CollectGenoType(); */
	/* cMyc_AffyMatch(); */
	/* cMyc_AgilentMatch(); */
	/* cMyc_ExonMatch(); */


	/* menu_hts_selectreads(argv, argc); */
	/* dTemp = norminv(0, 1, 0.5); */
	/* IQR_ProbeMatch(); */

	/* book_norm(); */
	/* shhmbtest(argv, argc); */

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_file_bed2cod(int argv, char **argc)
{
	/* define */
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        file_bed2cod             */
	/* -i input                        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          file_bed2cod           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" example: \n");
		printf("    file_bed2cod -i input.bed -o output.cod \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = File_Bed2Cod_Main(strInputPath, strOutputPath);
	}

	return nResult;
}

int menu_file_cod2bed(int argv, char **argc)
{
	/* define */
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        file_cod2bed             */
	/* -i input                        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          file_cod2bed           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" example: \n");
		printf("    file_cod2bed -i input.cod -o output.bed \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = File_Cod2Bed_Main(strInputPath, strOutputPath);
	}

	return nResult;
}

int menu_file_bed2aln(int argv, char **argc)
{
	/* define */
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        file_bed2aln             */
	/* -i input                        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          file_bed2aln           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" example: \n");
		printf("    file_bed2aln -i input.bed -o output.aln \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = File_Bed2Aln_Main(strInputPath, strOutputPath);
	}

	return nResult;
}

int menu_codinggenome(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        codinggenome             */
	/* -d database                     */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          codegenome         \n");
		printf(" -d path of genome database      \n");
		printf(" -o path for saving coded genome (*.sq files) \n");
		printf(" example: \n");
		printf("    codegenome -d /data/genomes/human/b35_hg17/ -o /data/genomes/human/b35_hg17/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded. \n");
		printf("    After coding, a new file chrlen.txt will be created to record the length for each chromosome. \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	/*strcpy(strInPath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");
	strcpy(strOutPath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");*/

	if((dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Fasta_To_Code_4bit_Main(strGenomePath, strOutputPath);
	}

	return nResult;
}

int menu_getseqfromfasta(int argv, char **argc)
{
	/* ------------------------------- */
	/*        fasta_getseq             */
	/* get sequence from FASTA file    */
	/* ------------------------------- */
	char strFASTAPath[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nUp,nDown;
	int nStrandType;
	int ni;
	int dOK,iOK,oOK,rOK,upOK,downOK;
	int nResult;

	/* ------------------------------- */
	/*        fasta_getseq             */
	/* -d FASTA database               */
	/* -i target region, 0-indexed seq */
	/* -o output                       */
	/* -up upstream bp                 */
	/* -down downstream bp             */
	/* -r strand                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        fasta_getseq            \n");
		printf(" -d path of FASTA database      \n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: zero-based seqid; \n");
		printf("		col2: zero-based start position in the assembly \n");
		printf("		col3: zero-based end position in the assembly \n");
		printf("		col4: strand in the assembly, '+' or '-' \n");
		printf(" -o file for saving the sequences\n");
		printf(" -up upstream extension length\n");
		printf(" -down downstream extension length\n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3'; \n");
		printf("		assemblybase [default]: always save the + strand from fasta sequences. \n");
		printf(" example: \n");
		printf("    fasta_getseq -d seq.fa -i coord.txt -o extractedseq.fa -up 2 -down 1 -r genebase \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	/* set default */
	dOK = 0;
	iOK = 0;
	oOK = 0;
	upOK = 0;
	downOK = 0;
	rOK = 0;
	nUp = 0;
	nDown = 0;
	nStrandType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strFASTAPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FastaSequenceExtract_Main(strFASTAPath, strTargetFile, strOutputFile, nUp, nDown, nStrandType); 
	}

	return nResult;
}

int menu_getseqfromgenome(int argv, char **argc)
{
	/* ------------------------------- */
	/*        getseqfromgenome         */
	/* get sequence from genome module */
	/* ------------------------------- */
	char strGenomePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nStrandType;
	int ni;
	int dOK,sOK,iOK,oOK,rOK;
	int nResult;

	/* ------------------------------- */
	/*        getseqfromgenome         */
	/* -d database                     */
	/* -s species                      */
	/* -i target region                */
	/* -o output                       */
	/* -r strand                       */
	/* get sequence from genome module */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        getseqfromgenome         \n");
		printf(" -d path of genome database      \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, D_melanogaster\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in numbers, e.g. use 23 for human chrX, 24 for human chrY \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o file for saving the sequences\n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" example: \n");
		printf("    getseqfromgenome -d /data/genomes/human/b35_hg17/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testseq.fa -r genebase \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;

	/* set default */
	nStrandType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */

	/* nCount = Genome_Code_4bit_GetSeq_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\",
								"human",
								"C:\\Projects\\research_harvard\\genomelab_project\\projects\\xiaoman\\gillData.txt", 
								"C:\\Projects\\research_harvard\\genomelab_project\\projects\\xiaoman\\gillData.fa", 
								0); */

	if( (dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeq_Main(strGenomePath, strSpecies, strTargetFile, strOutputFile, nStrandType); 
	}

	return nResult;
}

int menu_getseqfromgenome_c(int argv, char **argc)
{
	/* ------------------------------- */
	/*        getseqfromgenome         */
	/* get sequence from genome module */
	/* ------------------------------- */
	char strGenomePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nStrandType;
	int ni;
	int dOK,sOK,iOK,oOK,rOK;
	int nResult;

	/* ------------------------------- */
	/*     genome_getseq_c             */
	/* -d database                     */
	/* -s species                      */
	/* -i target region                */
	/* -o output                       */
	/* -r strand                       */
	/* get sequence from genome module */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        genome_getseq_c          \n");
		printf(" -d path of genome database      \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in string, e.g. use chrX, chr1, etc.\n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o file for saving the sequences\n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" example: \n");
		printf("    genome_getseq_c -d /data/genomes/human/b35_hg17/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testseq.fa -r genebase \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;

	/* set default */
	nStrandType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeq_C_Main(strGenomePath, strSpecies, strTargetFile, strOutputFile, nStrandType); 
	}

	return nResult;
}

int menu_getseqcsfromgenome(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	char strCSFormat[LINE_LENGTH];

	int nStrandType;
	int ni;
	int dOK,cOK,sOK,iOK,oOK,aOK,rOK,fOK;
	int nResult;

	/* ------------------------------- */
	/*      getseqcsfromgenome         */
	/* -d database                     */
	/* -c conservation                 */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* -a output sequence file name    */
	/* -r strand                       */
	/* -f conservation format          */
	/* get sequence and conservation   */
	/* from genome                     */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        getseqcsfromgenome       \n");
		printf(" -d path of genome database      \n");
		printf(" -c path of conservation database \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in numbers, e.g. use 23 for human chrX, 24 for human chrY \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" -a file for saving the sequences \n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" -f conservation score format \n");
		printf("        three possible formats: cs, txt, bed\n");
		printf(" example: \n");
		printf("    getseqcsfromgenome -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o ./ -a testseq -r genebase -f txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	cOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	aOK = 0;
	rOK = 0;
	fOK = 0;

	/* set default */
	nStrandType = 0;
	strcpy(strCSFormat, "cs");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strCSFormat, argc[ni]);
			fOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (cOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeqCS_Main(strGenomePath, strConservePath, strSpecies, 
			strTargetFile, strOutputPath, strSeqFile, nStrandType, strCSFormat); 
	}

	return nResult;


	/* nCount = Genome_Code_4bit_GetSeqCS_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\",
								"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\",
								"human",
								"C:\\Projects\\research_harvard\\genomelab_project\\test\\E2F.txt", 
								"C:\\Projects\\research_harvard\\genomelab_project\\test\\", "E2F", 
								1, "cs");
	*/

	/* nCount = Genome_Code_4bit_GetSeqCS_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\",
								"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\",
								"human",
								"C:\\Projects\\research_harvard\\genomelab_project\\test\\chr21.txt", 
								"C:\\Projects\\research_harvard\\genomelab_project\\test\\", "chr21", 
								1, "txt"); */

	/* nCount = Genome_Code_4bit_GetSeqCS_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\",
								"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\",
								"human",
								"C:\\Projects\\research_harvard\\genomelab_project\\projects\\hox\\hoxcord.txt", 
								"C:\\Projects\\research_harvard\\genomelab_project\\projects\\hox\\", "hoxseq", 
								1, "cs");
	*/

}


int menu_getseqcsfromgenome_c(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	char strCSFormat[LINE_LENGTH];

	int nStrandType;
	int ni;
	int dOK,cOK,sOK,iOK,oOK,aOK,rOK,fOK;
	int nResult;

	/* ------------------------------- */
	/*   genome_getseqcs_c             */
	/* -d database                     */
	/* -c conservation                 */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* -a output sequence file name    */
	/* -r strand                       */
	/* -f conservation format          */
	/* get sequence and conservation   */
	/* from genome                     */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_getseqcs_c            \n");
		printf(" -d path of genome database      \n");
		printf(" -c path of conservation database \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in string, e.g. use chrX, chr1, etc. \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" -a file for saving the sequences \n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" -f conservation score format \n");
		printf("        three possible formats: cs, txt, bed\n");
		printf(" example: \n");
		printf("    genome_getseqcs_c -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o ./ -a testseq -r genebase -f txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	cOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	aOK = 0;
	rOK = 0;
	fOK = 0;

	/* set default */
	nStrandType = 0;
	strcpy(strCSFormat, "cs");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strCSFormat, argc[ni]);
			fOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (cOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeqCS_C_Main(strGenomePath, strConservePath, strSpecies, 
			strTargetFile, strOutputPath, strSeqFile, nStrandType, strCSFormat); 
	}

	return nResult;
}

int menu_getmaskedseqfromgenome(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strCdsPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	
	int nStrandType;
	int ni;
	double dC;
	int dOK,cOK,cdOK,cdsOK,sOK,iOK,oOK,rOK;
	int nResult;

	/* ------------------------------- */
	/* genome_getmaskedseq             */
	/* -d database                     */
	/* -c conservation cutoff          */
	/* -cd conservation                */
	/* -cds coding region              */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* -r strand                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_getmaskedseq         \n");
		printf(" -d path of genome database      \n");
		printf(" -c conservation cutoff          \n");
		printf(" -cd path of conservation database \n");
		printf(" -cds path of coding region database\n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in numbers, e.g. use 23 for human chrX, 24 for human chrY \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" example: \n");
		printf("    genome_getmaskedseq -d /data/genomes/human/b35_hg17/ -c 40 -cd /data/genomes/human/b35_hg17/conservation/phastcons/ -cds /data/genomes/human/b35_hg17/conservation/phastcons/cds/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testseq.fa -r genebase\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dC = 0.0;
	dOK = 0;
	cOK = 0;
	cdOK = 0;
	cdsOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	
	/* set default */
	nStrandType = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else if( (cOK == 1) && (cdOK == 0) )
	{
		printf("Error: Input Parameter not correct, no conservation path was specified!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetMaskedSeq_Main(strGenomePath, strSpecies,
			cOK, dC, strConservePath, cdsOK, strCdsPath,
			strTargetFile, strOutputPath, nStrandType); 
	}

	return nResult;
}

int menu_getmaskedseqfromgenome_c(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strCdsPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	
	int nStrandType;
	int ni;
	double dC;
	int dOK,cOK,cdOK,cdsOK,sOK,iOK,oOK,rOK;
	int nResult;

	/* ------------------------------- */
	/* genome_getmaskedseq_c           */
	/* -d database                     */
	/* -c conservation cutoff          */
	/* -cd conservation                */
	/* -cds coding region              */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* -r strand                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_getmaskedseq_c       \n");
		printf(" -d path of genome database      \n");
		printf(" -c conservation cutoff          \n");
		printf(" -cd path of conservation database \n");
		printf(" -cds path of coding region database\n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in string, e.g. use chrX, chr1, etc. \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" example: \n");
		printf("    genome_getmaskedseq_c -d /data/genomes/human/b35_hg17/ -c 40 -cd /data/genomes/human/b35_hg17/conservation/phastcons/ -cds /data/genomes/human/b35_hg17/conservation/phastcons/cds/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testseq.fa -r genebase\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dC = 0.0;
	dOK = 0;
	cOK = 0;
	cdOK = 0;
	cdsOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	
	/* set default */
	nStrandType = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else if( (cOK == 1) && (cdOK == 0) )
	{
		printf("Error: Input Parameter not correct, no conservation path was specified!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetMaskedSeq_C_Main(strGenomePath, strSpecies,
			cOK, dC, strConservePath, cdsOK, strCdsPath,
			strTargetFile, strOutputPath, nStrandType); 
	}

	return nResult;
}

int menu_getmaskedregionfromgenome(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strCdsPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
		
	int ni;
	double dR,dC,dCR,dCDS;
	int dOK,rOK,cOK,crOK,cdOK,cdsOK,cdsdOK,sOK,iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* genome_getmaskedreg             */
	/* -d database                     */
	/* -r repeat cutoff                */
	/* -c conservation cutoff          */
	/* -cd conservation                */
	/* -cds cutoff                     */
	/* -cdsd coding region             */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_getmaskedreg         \n");
		printf(" -d path of genome database      \n");
		printf(" -r repeat cutoff                \n");
		printf(" -c conservation cutoff          \n");
		printf(" -cr minimum percentage of conserved bases \n");
		printf(" -cd path of conservation database \n");
		printf(" -cds CDS cutoff                 \n");
		printf(" -cdsd path of coding region database\n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, cow, chicken, zebrafish\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in numbers, e.g. use 23 for human X, 1, 2, etc. \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    genome_getmaskedreg -d /data/genomes/human/b35_hg17/ -r 0.9 -c 40 -cr 0.9 -cd /data/genomes/human/b35_hg17/conservation/phastcons/ -cds 0.9 -cdsd /data/genomes/human/b35_hg17/conservation/phastcons/cds/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testid_masked.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dR = 0.0;
	dC = 0.0;
	dCR = 0.0;
	dCDS = 0.0;
	dOK = 0;
	rOK = 0;
	cOK = 0;
	crOK = 0;
	cdOK = 0;
	cdsOK = 0;
	cdsdOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cr") == 0)
		{
			ni++;
			dCR = atof(argc[ni]);
			crOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			dCDS = atof(argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-cdsd") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsdOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) || ( (rOK == 0) && (cOK == 0) && (cdsOK == 0) ) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}

    if( (cOK == 1) || (crOK == 1) || (cdOK == 1))
	{
		if( (cOK == 0) || (crOK == 0) || (cdOK == 0))
		{
			printf("Error: Input Parameter not correct, no conservation path was specified!\n");
			exit(EXIT_FAILURE);
		}
	}

	if( (cdsOK == 1) || (cdsdOK == 1))
	{
		if( (cdsOK == 0) || (cdsdOK == 0) )
		{
			printf("Error: Input Parameter not correct, no cds path was specified!\n");
			exit(EXIT_FAILURE);
		}
	}

	nResult = Genome_Code_4bit_GetMaskedReg_Main(strGenomePath, strSpecies,
			rOK, dR, cOK, dC, dCR, strConservePath, cdsOK, dCDS, strCdsPath,
			strTargetFile, strOutputPath); 

	return nResult;
}

int menu_getmaskedregionfromgenome_c(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strCdsPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
		
	int ni;
	double dR,dC,dCR,dCDS;
	int dOK,rOK,cOK,crOK,cdOK,cdsOK,cdsdOK,sOK,iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* genome_getmaskedreg_c           */
	/* -d database                     */
	/* -r repeat cutoff                */
	/* -c conservation cutoff          */
	/* -cd conservation                */
	/* -cds cutoff                     */
	/* -cdsd coding region             */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_getmaskedreg_c       \n");
		printf(" -d path of genome database      \n");
		printf(" -r repeat cutoff                \n");
		printf(" -c conservation cutoff          \n");
		printf(" -cr minimum percentage of conserved bases \n");
		printf(" -cd path of conservation database \n");
		printf(" -cds CDS cutoff                 \n");
		printf(" -cdsd path of coding region database\n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, cow, chicken, zebrafish\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in string, e.g. use chrX, chr1, etc. \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    genome_getmaskedreg_c -d /data/genomes/human/b35_hg17/ -r 0.9 -c 40 -cr 0.9 -cd /data/genomes/human/b35_hg17/conservation/phastcons/ -cds 0.9 -cdsd /data/genomes/human/b35_hg17/conservation/phastcons/cds/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testid_masked.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dR = 0.0;
	dC = 0.0;
	dCR = 0.0;
	dCDS = 0.0;
	dOK = 0;
	rOK = 0;
	cOK = 0;
	crOK = 0;
	cdOK = 0;
	cdsOK = 0;
	cdsdOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cr") == 0)
		{
			ni++;
			dCR = atof(argc[ni]);
			crOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			dCDS = atof(argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-cdsd") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsdOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) || ( (rOK == 0) && (cOK == 0) && (cdsOK == 0) ) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}

    if( (cOK == 1) || (crOK == 1) || (cdOK == 1))
	{
		if( (cOK == 0) || (crOK == 0) || (cdOK == 0))
		{
			printf("Error: Input Parameter not correct, no conservation path was specified!\n");
			exit(EXIT_FAILURE);
		}
	}

	if( (cdsOK == 1) || (cdsdOK == 1))
	{
		if( (cdsOK == 0) || (cdsdOK == 0) )
		{
			printf("Error: Input Parameter not correct, no cds path was specified!\n");
			exit(EXIT_FAILURE);
		}
	}

	nResult = Genome_Code_4bit_GetMaskedReg_C_Main(strGenomePath, strSpecies,
			rOK, dR, cOK, dC, dCR, strConservePath, cdsOK, dCDS, strCdsPath,
			strTargetFile, strOutputPath); 

	return nResult;
}

int menu_genome_hash(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nK = 13;

	int ni;
	int dOK,oOK,kOK;
	int nResult;

	/* ------------------------------- */
	/* genome_hash                     */
	/* -d database                     */
	/* -o output path                  */
	/* -k key length                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_hash                 \n");
		printf(" -d path of genome database      \n");
		printf(" -o output path \n");
		printf(" -k key length \n");
		printf(" example: \n");
		printf("    genome_hash -d /data/genomes/human/b35_hg17/ -o ./ -k 14\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	oOK = 0;
	kOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nK = atoi(argc[ni]);
			kOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (oOK == 0) || (kOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}

	nResult = Genome_CreateHash_Main(strGenomePath, strOutputPath, 
						   nK);

	return nResult;
}

int menu_seqmask(int argv, char **argc)
{
	/* define */
	char strSeqFile[LINE_LENGTH];
	char strMaskFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int nMaskType = 0;
	int iOK,mOK,mtOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* menu_fastaseqmask               */
	/* -i sequence                     */
	/* -m masks                        */
	/* -mt mask type:                  */
	/*     0-soft mask, 1-hard mask    */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_fastaseqmask           \n");
		printf(" -i FASTA sequences                \n");
		printf(" -m masks                          \n");
		printf(" -mt mask type (default 0)         \n");
		printf("     0: soft mask, masking with small letters a, c, g, t.\n");
		printf("     1: hard mask, masking with N \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    genome_fastaseqmask -i test.fa -m test.mask -mt 1 -o testm.fa\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	iOK = 0;
	mOK = 0;
	mtOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMaskFile, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-mt") == 0)
		{
			ni++;
			nMaskType = atoi(argc[ni]);
			mtOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (iOK == 0) ||  (mOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FastaSequenceMask_Main(strSeqFile, strMaskFile,
			nMaskType, strOutFile); 
	}

	return nResult;
}

int menu_fasta_soft2hardmask(int argv, char **argc)
{
	/* define */
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* fasta_soft2hardmask             */
	/* -i sequence                     */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    fasta_soft2hardmask          \n");
		printf(" -i FASTA sequences                \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    fasta_soft2hardmask -i test.fa -o testm.fa\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	iOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FastaSequenceSoft2HardMask_Main(strSeqFile, strOutFile); 
	}

	return nResult;
}

int menu_maf_conserve_bg(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    conserve_bg                  */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_conservebg              \n");
		printf(" example: \n");
		printf("    genome_conservebg conservebg_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = Genome_MafAlign_To_BG_Main(strParamPath);
	
	/* nCount = Genome_MafAlign_To_BG_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\param_bg.txt");
	*/

	return nResult;
}


int menu_maf_conserve_cs(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    conserve_bg                  */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_conservecs              \n");
		printf(" example: \n");
		printf("    genome_conservecs conservecs_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = Genome_MafAlign_To_CS_Main(strParamPath);
	
	/* nCount = Genome_MafAlign_To_CS_Main("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\param_cs.txt");
	*/

	return nResult;
}

int menu_cs_getdistn(int argv, char **argc)
{
	/* define */
	char strDataPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrListFile[LINE_LENGTH];
	char strChrLenFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	
	int ni;
	int dOK,sOK,iOK,oOK,lOK;
	int nResult;

	/* ------------------------------- */
	/*       genome_csgetdistn         */
	/* -d database                     */
	/* -s species                      */
	/* -l chrmosome length file        */
	/* -i target region                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        genome_csgetdistn         \n");
		printf(" -d path of conservation score database \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, D_melanogaster\n");
		printf(" -i target chromosome \n");
		printf(" -o file for saving the output statistics\n");
		printf(" -l file for chromosome length        \n");
		printf(" example: \n");
		printf("    genome_csgetdistn -d /data/genomes/human/b33_hg15/conservation/cs/ -s human -l chrlen.txt -i chrlist.txt -o csstat.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	lOK = 0;

	/* set default */
	ni = 1;
	/* while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strChrListFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	} 

	if( (dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_CS_Get_Distribution(strDataPath, strChrListFile, strChrLenFile, strOutputFile, strSpecies);
	}
	*/
	

	nResult = Genome_CS_Get_Distribution("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\",
			"chrlist.txt", "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\chrlen.txt",
			"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\csstat.txt", "human");


	return nResult;
}

int menu_maf_footprint(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    conserve_bg                  */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_footprint               \n");
		printf(" example: \n");
		printf("    genome_footprint genome_footprint_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = Genome_MafAlign_To_FootPrint_Main(strParamPath);
	
	return nResult;
}

int menu_csgc_summary(int argv, char **argc)
{
	/* define */
	int gdOK,iOK,oOK,cOK,cdOK;
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*   genome_getcsgcsummary         */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   genome_getcsgcsummary    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file    \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" example: \n");
		printf("    genome_getcsgcsummary -gd /data/mm6 -i inputseq.cod -o mm6_summary -c 100 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = 0.0;
	nUseCS = 0;
	strcpy(strCSPath, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((gdOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = Genome_GetGCCS_Summary_Main(strGenomePath, 
					strCodPath,	strOutputPath, 
					cOK, dC, strCSPath);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = Genome_GetGCCS_Summary_Main(strGenomePath, 
					strCodPath,	strOutputPath, 
					cOK, dC, strCSPath);
		}
	}

	/* return */
	return nResult;
}

int menu_regioncs_summary(int argv, char **argc)
{
	/* define */
	int gdOK,iOK,oOK,cOK,cdOK;
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*   genome_regioncssummary        */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   genome_regioncssummary   \n");
		printf(" -gd genome sequence path \n");
		printf(" -i coordinates file    \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" example: \n");
		printf("    genome_regioncssummary -gd /data/mm6 -i inputseq.cod -o inputcs.txt -c 100 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = -0.01;
	nUseCS = 1;
	strcpy(strCSPath, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((gdOK == 0) || (iOK == 0) || (oOK == 0) || (cdOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		cOK = 1;
		nResult = Genome_GetRegionCS_Summary_Main(strGenomePath, 
				strCodPath,	strOutputPath, 
				cOK, dC, strCSPath);
	}

	/* return */
	return nResult;
}


int menu_genome_regionextend(int argv, char **argc)
{
	/* define */
	int iOK,lOK,rOK,dOK,sOK,oOK,cnOK,aOK,tOK;
	char strInputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nL,nR,nCN,nA,nUseStrand;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    genome_regionextend          */
	/* -i input coordinates file       */
	/* -l left extension length        */
	/* -r right extension length       */
	/* -t left/right based on strand info */
	/*    0: ignore strand, assemblybase (default)*/
	/*    1: +/- base                   */
	/* -d genome sequence path         */
	/* -s species                      */
	/* -cn numerical chromosome id     */
	/* -o output file                  */
	/* -a alias type                   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_regionextend  \n");
		printf(" -i input coordinates file   \n");
		printf(" -l left extension length (default = 0) \n");
		printf(" -r right extension length (default = 0) \n");
		printf(" -t left/right based on strand info   \n");
		printf("    0: ignore strand, assemblybase (default)  \n");
		printf("    1: 5' left, 3' right; 5'/3' determined by +/- \n");
		printf(" -d genome sequence path     \n");
		printf(" -s species  (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -cn numerical chromosome id [0=string chr name (default); 1=integer chr id]\n");
		printf(" -o output file \n");
		printf(" -a alias type [0=use original (default); 1=reorder] \n");
		printf(" example: \n");
		printf("    genome_regionextend -i Gli_map.txt -l 100 -r 200 -t 1 -d /data/genomes/mouse/mm6/ -s mouse -cn 0 -a 1 -o Gli_sitearound.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	lOK = 0;
	rOK = 0;
	dOK = 0;
	sOK = 0;
	cnOK = 0;
	oOK = 0;
	aOK = 0;
	tOK = 0;
	
	nL = 0;
	nR = 0;
	nCN = 0;
	nA = 0;
	nUseStrand = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nL = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nR = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nUseStrand = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nA = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-cn") == 0)
		{
			ni++;
			nCN = atoi(argc[ni]);
			cnOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (sOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_RegionExtend_Main(strInputPath, strGenomePath, 
			strOutputPath, strSpecies, nL, nR, nCN, nA, nUseStrand);
	}

	/* return */
	return nResult;
}

int menu_codingphastcons(int argv, char **argc)
{
	/* define */
	char strPhastPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int cOK;
	int ni;
	char strExt[LINE_LENGTH];

	/* ------------------------------- */
	/*        codingphastcons          */
	/* -d genome database              */
	/* -c conservation database        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("       genome_codephastcons        \n");
		printf(" -d path of genome database      \n");
		printf(" -c path of phastCons score database \n");
		printf(" -o path for saving coded score (*.cs files) \n");
		printf(" -e original phastCons file extension (e.g. \"-e .pp\" means chr1.pp is the original score file for chr1). Default = no extension. \n");
		printf(" example: \n");
		printf("    codephastcons -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -o /data/genomes/human/b35_hg17/conservation/phastcons/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded and a chrlen.txt file for chromosome lengths \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	cOK = 0;
	strcpy(strExt, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strPhastPath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			strcpy(strExt, argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	
	if((dOK == 0) || (cOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		AdjustDirectoryPath(strGenomePath);
		AdjustDirectoryPath(strPhastPath);
		AdjustDirectoryPath(strOutputPath);
		nResult = Genome_PhastCons_To_Code_8bit_Main(strGenomePath, strPhastPath, strOutputPath, strExt);
	}

	return nResult;
}

int menu_codingphastcons_v2(int argv, char **argc)
{
	/* define */
	char strPhastPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int cOK;
	int ni;
	char strExt[LINE_LENGTH];

	/* ------------------------------- */
	/*        codingphastcons          */
	/* -d genome database              */
	/* -c conservation database        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_codephastcons_v2        \n");
		printf(" To deal with phastcons format released Nov. 2004 and later. \n"); 
		printf(" -d path of genome database      \n");
		printf(" -c path of phastCons score database \n");
		printf(" -o path for saving coded score (*.cs files) \n");
		printf(" -e original phastCons file extension (e.g. \"-e .pp\" means chr1.pp is the original score file for chr1). Default = no extension. \n");
		printf(" example: \n");
		printf("    genome_codephastcons_v2 -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -o /data/genomes/human/b35_hg17/conservation/phastcons/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded and a chrlen.txt file for chromosome lengths \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	cOK = 0;
	strcpy(strExt, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strPhastPath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			strcpy(strExt, argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	
	if((dOK == 0) || (cOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		AdjustDirectoryPath(strGenomePath);
		AdjustDirectoryPath(strPhastPath);
		AdjustDirectoryPath(strOutputPath);
		nResult = Genome_PhastCons_To_Code_8bit_Main_v2(strGenomePath, strPhastPath, strOutputPath, strExt);
	}

	return nResult;
}

int menu_codingCDS(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int gOK;
	int gtOK;
	int sOK;
	int nGType;
	int ni;

	/* ------------------------------- */
	/*        genome_codingCDS         */
	/* -d genome database              */
	/* -g path of refgene database     */
	/* -gt refgene database type       */
	/* -s species                      */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_codingCDS             \n");
		printf(" -d path of genome database      \n");
		printf(" -g path of refgene database     \n");
		printf(" -gt refgene database type       \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf(" -s species                      \n");
		printf(" -o path for saving coded vector (*.cds files) \n");
		printf(" example: \n");
		printf("    genome_codingCDS -d /data/genomes/human/b35_hg17/ -g /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -gt 1 -s human -o /data/genomes/human/b35_hg17/cds/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded and a chrlen.txt file for chromosome lengths \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	gOK = 0;
	gtOK = 0;
	sOK = 0;
	nGType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-gt") == 0)
		{
			ni++;
			nGType = atoi(argc[ni]);
			gtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	
	if((dOK == 0) || (gOK == 0) || (gtOK == 0) || (sOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_CodeCDS_FromRefGene_Main(strGenomePath, strRefGenePath, nGType, strSpecies, strOutputPath);
	}

	return nResult;
}

int menu_refgene_sort(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;

	/* ------------------------------- */
	/*        refgene_sort             */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          coderefgene            \n");
		printf(" -d path of refgene database      \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" example: \n");
		printf("    coderefgene -d /data/genomes/human/b35_hg17/annotation/refGene.txt -o /data/genomes/human/b35_hg17/annotation/refGene_sorted.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_Convert_UCSC_To_Lab(strRefGenePath, strOutputPath, strSpecies, nChrNum);
	}


	
	/* nCount = RefGene_Convert_UCSC_To_Lab("C:\\Projects\\research_harvard\\genomelab_project\\genomes\\mouse\\b33_mm5\\annotation\\refGene.txt", 
		"C:\\Projects\\research_harvard\\genomelab_project\\genomes\\mouse\\b33_mm5\\annotation\\refGene_sorted.txt", 
		"mouse", 21); */

	return nResult;
}

int menu_refflat_sort(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;

	/* ------------------------------- */
	/*        refflat_encode           */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          refflat_encode          \n");
		printf(" -d path of refgene database      \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" example: \n");
		printf("    refflat_encode -d /data/genomes/human/b35_hg17/annotation/refFlat.txt -o /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefFlat_Convert_UCSC_To_Lab(strRefGenePath, strOutputPath, strSpecies, nChrNum);
	}

	return nResult;
}

int menu_reflocus_createmap(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;
	int nChangeGeneName = 0;

	/* ------------------------------- */
	/*        reflocus_createmap       */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        reflocus_createmap          \n");
		printf(" -d path of sorted refgene&geneid database   \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" -u update gene name (default = 0). If \"-u 1\", the GENEID line should contain 5 fields: GENEID[tab]REFID[tab]ENTREZID[tab]TAXID[tab]GENENAME \n");
		printf(" example: \n");
		printf("    reflocus_createmap -d /data/genomes/human/b35_hg17/annotation/refFlatgeneid.txt -o /data/genomes/human/b35_hg17/annotation/refLocus_sorted.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;
	nChangeGeneName = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nChangeGeneName = atoi(argc[ni]);
			if((nChangeGeneName != 1) && (nChangeGeneName != 0))
				nChangeGeneName = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult =  RefGene_AnnotateWithLocusID(strRefGenePath, strOutputPath, strSpecies, nChrNum, nChangeGeneName);
	}

	return nResult;
}

int menu_refexon_createmap(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;

	/* ------------------------------- */
	/*        refexon_createmap        */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        refexon_createmap          \n");
		printf(" -d path of sorted refgene&exonid database   \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" example: \n");
		printf("    refexon_createmap -d /data/genomes/human/b35_hg17/annotation/refFlatexonid.txt -o /data/genomes/human/b35_hg17/annotation/refExon_map.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult =  RefGene_AnnotateWithExonArrayID(strRefGenePath, strOutputPath, strSpecies, nChrNum);
	}

	return nResult;
}

int menu_refmicroarray_createmap(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nChrNum;
	int nResult;
	int dOK;
	int oOK;
	int sOK;
	int numOK;
	int ni;

	/* ------------------------------- */
	/*    refmicroarray_createmap      */
	/* -d database                     */
	/* -o output                       */
	/* -s species                      */
	/* -n number of chromosome         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     refmicroarray_createmap          \n");
		printf(" -d path of sorted refgene&exonid database   \n");
		printf(" -o path for saving coded and sorted refgene \n");
		printf(" -s species \n");
		printf(" -n number of chromosome \n");
		printf(" example: \n");
		printf("    refmicroarray_createmap -d /data/genomes/human/b35_hg17/annotation/refFlatexonid.txt -o /data/genomes/human/b35_hg17/annotation/refExon_map.txt -s human -n 24\n\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	/* nResult =  RefGene_AnnotateWithMicroArrayID("C:\\Projects\\research_harvard\\hedgehog_project\\microarrays\\Limbud-Steve\\refFlat_Mouse430_rawsorted.txt", 
		"C:\\Projects\\research_harvard\\hedgehog_project\\microarrays\\Limbud-Steve\\refMouse430_map.txt", "mouse", 21);
	*/
	
	dOK = 0;
	oOK = 0;
	sOK = 0;
	numOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nChrNum = atoi(argc[ni]);
			numOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (oOK == 0) || (sOK == 0) || (numOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult =  RefGene_AnnotateWithMicroArrayID(strRefGenePath, strOutputPath, strSpecies, nChrNum);
	}

	return nResult;
}

int menu_refgene_gettargettssaround(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strTargetPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrLen[LINE_LENGTH];
	int nTSSUP;
	int nTSSDOWN;

	int nResult;
	int dOK;
	int tOK;
	int oOK;
	int sOK;
	int upOK;
	int downOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*        refgene_gettssaround     */
	/* -d database                     */
	/* -t target list                  */
	/* -o output                       */
	/* -s species                      */
	/* -up TSS up                      */
	/* -down TSS down                  */
	/* -c chromosome length            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          refgene_gettssaround    \n");
		printf(" -d path of refgene database      \n");
		printf(" -t path of target list           \n");
		printf(" -o path for saving retrieved target coordinates \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -up TSS up \n");
		printf(" -down TSS down \n");
		printf(" -c chromosome length\n");
		printf(" example: \n");
		printf("    refgene_gettssaround -d /data/genomes/human/b35_hg17/annotation/xenoRefGene_sorted.txt -t /data/genomes/human/b35_hg17/annotation/testrefid.txt -o humcod.txt -s human -up 5000 -down 1000 -c /data/genomes/human/b35_hg17/chrlen.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	tOK = 0;
	oOK = 0;
	sOK = 0;
	upOK = 0;
	downOK = 0;
	cOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nTSSUP = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nTSSDOWN = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (tOK == 0) || (oOK == 0) || (sOK == 0) || (upOK == 0) || (downOK == 0) || (cOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetTargetTSSAround(strRefGenePath, strTargetPath,
			nTSSUP, nTSSDOWN, strSpecies, strChrLen, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_refgene_getnearestgene(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nRefType = 0;
	int nUP;
	int nDOWN;

	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int oOK;
	int rOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getnearestgene     */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -i input coordinates            */
	/* -o output                       */
	/* -r reference type               */
	/*    0: TSS-up, TES-down          */
	/*    1: TSS-up, TSS-down          */
	/*    2: TES-up, TES-down          */
	/*    3: CDSS-up, CDSE-down        */
	/*    4: CDSS-up, CDSS-down        */
	/*    5: CDSE-up, CDSE-down        */
	/* -up up distance                 */
	/* -down down distance down        */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        refgene_getnearestgene     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -r reference type \n");
		printf("    0: TSS-up, TES-down (default)    \n");
		printf("    1: TSS-up, TSS-down     \n");
		printf("    2: TES-up, TES-down     \n");
		printf("    3: CDSS-up, CDSE-down   \n");
		printf("    4: CDSS-up, CDSS-down   \n");
		printf("    5: CDSE-up, CDSE-down   \n");
		printf(" -up up distance limit\n");
		printf(" -down down distance limit \n");
		printf(" example: \n");
		printf("    refgene_getnearestgene -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -i target.txt -o target_gene.txt -r 0 -up 5000 -down 1000\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	nUP = 5000;
	nDOWN = 5000;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nRefType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUP = atoi(argc[ni]);
			if(nUP < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDOWN = atoi(argc[ni]);
			if(nDOWN < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) || (upOK == 0) || (downOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetNearestGene_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, strOutputPath,
			nRefType, nUP, nDOWN);
	}

	/* return */
	return nResult;
}

int menu_refgene_getlocationsummary(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nInputType = 0;
	
	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int oOK;
	int rOK;
	int ni;

	/* ------------------------------- */
	/*   refgene_getlocationsummary    */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -i input coordinates            */
	/* -o output                       */
	/* -r input type (0: cod, 1: bed   */
	/*      2: codp, 3: bedp)          */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getlocationsummary    \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -r input type \n");
		printf("    0: cod; 1: bed; 2: codp; 3: bedp \n");
		printf(" example: \n");
		printf("    refgene_getlocationsummary -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -i target.txt -o target_gene.txt -r 0\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetLocationSummary_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, nInputType, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_reflocus_getneighborgenes(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	int nGap = 0;
	int nUP = 1;
	int nDOWN = 1;

	int nResult;
	int dOK;
	int sOK;
	int iOK;
	int oOK;
	int aOK;
	int gOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*    reflocus_getneighborgenes    */
	/* -d database                     */
	/* -s species                      */
	/* -i input coordinates            */
	/* -o output                       */
	/* -a annotation                   */
	/* -g distance upper limit         */
	/* -up no. of upstream genes       */
	/* -down no. of downstream genes   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("      reflocus_getneighborgenes     \n");
		printf(" -d path of refgene database      \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input coordinates              \n");
		printf(" -o path for saving results \n");
		printf(" -a annotation  \n");
		printf(" -g distance upper limit  \n");
		printf(" -up no. of upstream genes \n");
		printf(" -down no. of downstream genes \n");
		printf(" example: \n");
		printf("    reflocus_getneighborgenes -d /data/genomes/human/b35_hg17/annotation/refLocus_sorted.txt -s human -i target.txt -o target_gene.txt -a expression.txt -g 1000000 -up 3 -down 3\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	aOK = 0;
	gOK = 0;
	upOK = 0;
	downOK = 0;
	nUP = 1;
	nDOWN = 1;
	nGap = 1000000000;
	strcpy(strAnnotationPath, "NULL");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUP = atoi(argc[ni]);
			if(nUP < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDOWN = atoi(argc[ni]);
			if(nDOWN < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetNeighborGenes_Main(strDatabasePath, 2,
			strSpecies, strInputPath, strOutputPath,
			strAnnotationPath, nUP, nDOWN, nGap);
	}

	/* nResult = RefGene_GetNeighborGenes_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		2, "mouse", "E:\\Projects\\hedgehog_project\\Test\\Gli3Genome_allt8_reg_mm6.txt", "E:\\Projects\\hedgehog_project\\Test\\peak_annot.txt", "E:\\Projects\\hedgehog_project\\Test\\refLocus_value.txt", 
		5, 5, 10000000); */

	/* return */
	return nResult;
}

int menu_refgene_getmatchedcontrol(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 0;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChrLenPath[LINE_LENGTH];
	int nRepNum = 1;
	int nRegionLen = 2000;
	int nNR = 1;
	
	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int cOK;
	int iOK;
	int oOK;
	int nOK;
	int lOK;
	int nrOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getmatchedcontrol  */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat                  */
	/*     2: refLocus                 */
	/* -c chromosome length            */
	/* -i input coordinates            */
	/* -o output                       */
	/* -n number of replications       */
	/* -l region length                */
	/* -nr remove redundancy           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     refgene_getmatchedcontrol     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format(default) \n");
		printf("     1: UCSC refFlat format          \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -c path of chromosome length \n");
		printf(" -i input coordinates         \n");
		printf(" -o file for saving results   \n");
		printf(" -n number of replications when selecting controls \n");
		printf(" -l control region length     \n");
		printf(" -nr remove redundancy or not \n");
		printf("     0: no                    \n");
		printf("     1: yes (default)         \n");
		printf(" example: \n");
		printf("    refgene_getmatchedcontrol -d /data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt -dt 1 -s human -c /data/genomes/human/b35_hg17/chrlen.txt -i target.txt -o target_ct.txt -n 3 -l 2000 -nr 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	cOK = 0;
	iOK = 0;
	oOK = 0;
	nOK = 0;
	lOK = 0;
	nrOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strChrLenPath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nRepNum = atoi(argc[ni]);
			nOK = 1;

			if(nRepNum <= 0)
			{
				printf("Error: -n must >0! \n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nRegionLen = atoi(argc[ni]);
			lOK = 1;

			if(nRegionLen <= 0)
			{
				printf("Error: -l must >0! \n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-nr") == 0)
		{
			ni++;
			nNR = atoi(argc[ni]);
			nrOK = 1;

			if( (nNR < 0) || (nNR > 1) )
			{
				printf("Error: -nr must be 0 or 1! \n");
				exit(EXIT_FAILURE);
			}
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (dtOK == 0) || (sOK == 0) || (cOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetMatchedControl_Main(strDatabasePath, nDatabaseType,
			strSpecies, strChrLenPath, strInputPath, strOutputPath,
			nRepNum, nRegionLen, nNR);
	}

	/* return */
	return nResult;
}


int menu_refgene_getaffy(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nColumn = 0;

	int nResult;
	int dOK;
	int iOK;
	int oOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*      refgene_getaffy            */
	/* -d affy-refgene map             */
	/* -i input file                   */
	/* -o output file                  */
	/* -c starting column (0-based)    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        refgene_getaffy            \n");
		printf(" -d affy-refgene map               \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -c starting column (0-based)      \n");
		printf(" example: \n");
		printf("    refgene_getaffy -d Mouse430_2_affy2refid.txt -i target.txt -c 1 -o target_affy.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nColumn = atoi(argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetAffy_Main(strDatabasePath, 
			strInputPath, nColumn, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_refgene_getortholog(int argv, char **argc)
{
	/* define */
	char strRefGenePath[LINE_LENGTH];
	char strSourceSpecies[LINE_LENGTH];
	char strSourceRefGenePath[LINE_LENGTH];
	char strDestSpecies[LINE_LENGTH];
	char strMapRefGenePath[LINE_LENGTH];
	char strDestRefGenePath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	int nResult;
	int tOK;
	int ssOK;
	int sdOK;
	int dsOK;
	int mdOK;
	int ddOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        refgene_getortholog      */
	/* -i target list                  */
	/* -ss source species              */
	/* -sd source refgene database     */
	/* -ds dest species                */
	/* -md map refgene database        */
	/* -dd dest refgene database       */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          getrefgeneortholog      \n");
		printf(" -i path of target list           \n");
		printf(" -ss source species               \n");
		printf(" -sd path of source refgene database \n");
		printf(" -ds dest species                 \n");
		printf(" -md path of map refgene database \n");
		printf(" -dd path of dest refgene database \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    getrefgeneortholog -t testrefid.txt -ss mouse -sd /data/genomes/mouse/mm5/annotation/refGene_sorted.txt -ds human -md /data/genomes/human/b35_hg17/annotation/xenoRefGene_sorted.txt -dd /data/genomes/human/b35_hg17/annotation/refGene_sorted.txt -o humortholog\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	tOK = 0;
	ssOK = 0;
	sdOK = 0;
	dsOK = 0;
	mdOK = 0;
	ddOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strRefGenePath, argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-ss") == 0)
		{
			ni++;
			strcpy(strSourceSpecies, argc[ni]);
			ssOK = 1;
		}
		else if(strcmp(argc[ni], "-sd") == 0)
		{
			ni++;
			strcpy(strSourceRefGenePath, argc[ni]);
			sdOK = 1;
		}
		else if(strcmp(argc[ni], "-ds") == 0)
		{
			ni++;
			strcpy(strDestSpecies, argc[ni]);
			dsOK = 1;
		}
		else if(strcmp(argc[ni], "-md") == 0)
		{
			ni++;
			strcpy(strMapRefGenePath, argc[ni]);
			mdOK = 1;
		}
		else if(strcmp(argc[ni], "-dd") == 0)
		{
			ni++;
			strcpy(strDestRefGenePath, argc[ni]);
			ddOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((tOK == 0) || (ssOK == 0) || (sdOK == 0) || (dsOK == 0) || (mdOK == 0) || (ddOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetOrtholog(strRefGenePath, strSourceSpecies, strSourceRefGenePath,
			strDestSpecies, strMapRefGenePath, strDestRefGenePath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_refgene_getmultiortholog(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiortholog     */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getmultiortholog      \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    getrefgenemultiortholog -i testrefid.txt -d orthologsetting.txt -o testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetMultiOrtholog_Main(strTargetPath, strParamPath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_refgene_getmultiortholog1way(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiortholog1way */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getmultiortholog1way   \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    getrefgenemultiortholog1way -i testrefid.msomap -d orthologsetting.txt -o testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_GetMultiOrtholog1way_Main(strTargetPath, strParamPath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_refflex_getmultiortholog(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	
	int nColumn = 0;

	int nResult;
	int iOK;
	int cOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*    refflex_getmultiortholog     */
	/* -i target list                  */
	/* -c starting column (0-based)    */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refflex_getmultiortholog      \n");
		printf(" -i path of target list           \n");
		printf(" -c starting column (0-based)     \n");
		printf("    the column where a UCSC refGene format record starts, e.g. NM_xxxxxx ...\n");
		printf("    default = 0                   \n");
		printf(" -d database infomation           \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    getrefgenemultiortholog -i testrefgene.txt -c 6 -d orthologsetting.txt -o testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	cOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nColumn = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefFlex_GetMultiOrtholog_Main(strTargetPath, nColumn, strParamPath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_refflex_getmultiortholog1way(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*    refflex_getmultiortholog1way */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refflex_getmultiortholog1way   \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" example: \n");
		printf("    refflex_getmultiortholog1way -i testrefid.nsomap -d orthologsetting.txt -o testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefFlex_GetMultiOrtholog1way_Main(strTargetPath, strParamPath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_refgene_getmultiorthologtssaround(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	
	char strSeqFile[LINE_LENGTH];
	int nUp;
	int nDown;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int upOK;
	int downOK;
	int aOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiorthologtss  */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    getmultiorthologtssaround      \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -up TSS up                       \n");
		printf(" -down TSS down                   \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" -a file for saving the sequences \n");
		printf(" example: \n");
		printf("    getmultiorthologtssaround -i testortho.msomap -up -5000 -down 1000 -d genomesetting.txt -o ./ -a testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	upOK = 0;
	downOK = 0;
	aOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0) || (upOK == 0) || (downOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_MultiOrthologGetTSSAround(strTargetPath, nUp, nDown, strParamPath, strOutPath, strSeqFile);
	}

	/* return */
	return nResult;
}


int menu_refgene_getmultiorthologtssaroundexonmasked(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	
	char strSeqFile[LINE_LENGTH];
	int nUp;
	int nDown;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int upOK;
	int downOK;
	int aOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiorthologtss  */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("getmultiorthologtssaroundexonmasked\n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -up TSS up                       \n");
		printf(" -down TSS down                   \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" -a file for saving the sequences \n");
		printf(" example: \n");
		printf("    getmultiorthologtssaroundexonmasked -i testortho.msomap -up -5000 -down 1000 -d genomesetting.txt -o ./ -a testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	upOK = 0;
	downOK = 0;
	aOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0) || (upOK == 0) || (downOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_MultiOrthologGetTSSAroundExonMasked(strTargetPath, nUp, nDown, strParamPath, strOutPath, strSeqFile);
	}

	/* return */
	return nResult;
}


int menu_refgene_pickspeciesspecific(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strDatabasePath[LINE_LENGTH];	
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/* refgene_pickspeciesspecific     */
	/* -i target path                  */
	/* -d database path                */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" refgene_pickspeciesspecific \n");
		printf(" -i path of target list           \n");
		printf(" -d path of database              \n");
		printf(" -o path for output \n");
		printf(" example: \n");
		printf("    refgene_pickspeciesspecific -i mouse_refgene.txt -d human_xenorefgene.txt -o human_mouserefgene.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_PickSpeciesSpecific(strTargetPath, strDatabasePath, strOutPath);
	}

	/* return */
	return nResult;
}

int menu_reflocus_assignvalue(int argv, char **argc)
{
	/* define */
	char strRefLocusDatabasePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];	
	double dTruncateLowerBound; 
	double dTruncateUpperBound;
	char strTransform[LINE_LENGTH];
	int nNormalize = 0;
	double dPostLowerBound; 
	double dPostUpperBound;
	char strPostTransform[LINE_LENGTH]; 
	int nTakeAbsoluteValue = 0;
	char strRefToProbeRowIDPath[LINE_LENGTH];
	char strRefToProbeNamePath[LINE_LENGTH];
	char strNetworkPath[LINE_LENGTH];
	char strNetworkAnnotationPath[LINE_LENGTH];
	int nNetDepth = 1;
		
	int nResult;
	int dOK = 0;
	int sOK = 0;
	int oOK = 0;
	int vOK = 0;
	int vlOK = 0;
	int vuOK = 0;
	int vtOK = 0;
	int nOK = 0;
	int plOK = 0;
	int puOK = 0;
	int ptOK = 0;
	int absOK = 0;
	int d2rOK = 0;
	int d2aOK = 0;
	int netOK = 0;
	int netAOK = 0;
	int netDOK = 0;

	int ni;

	/* ------------------------------- */
	/* menu_reflocus_assignvalue       */
	/* -d reflocus database path       */
	/* -s species                      */
	/* -o output path                  */
	/* -v score path                   */
	/* -vl score truncation lower bound*/
	/* -vu score truncation upper bound*/
	/* -vt score transformation        */
	/* -nOK normalization indicator    */
	/* -pl post normalization truncation lower bound*/
	/* -pu post normalization truncation upper bound*/
	/* -pt post normalization transformation        */
	/* -absOK take absolute value      */
	/* -d2r reflocus to score rowid map*/
	/* -d2a reflocus to probe name map */
	/* -net use physical network to enhance inference*/
	/* -netAOK network annotation      */
	/* -netDOK network search depth    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   reflocus_assignvalue            \n");
		printf(" -d reflocus database path       \n");
		printf(" -s species                      \n");
		printf(" -o output path                  \n");
		printf(" -v score path                   \n");
		printf(" -vl score truncation lower bound\n");
		printf(" -vu score truncation upper bound\n");
		printf(" -vt score transformation        \n");
		printf(" -nOK normalization indicator    \n");
		printf(" -pl post normalization truncation lower bound\n");
		printf(" -pu post normalization truncation upper bound\n");
		printf(" -pt post normalization transformation        \n");
		printf(" -abs take absolute value        \n");
		printf(" -d2r reflocus to score rowid map\n");
		printf(" -d2a reflocus to probe name map \n");
		printf(" -net use physical network to enhance inference\n");
		printf(" -netAOK network annotation      \n");
		printf(" -netDOK network search depth    \n");
		printf(" example: \n");
		printf("    reflocus_assignvalue -d refLocus_sorted.txt -s mouse -o refLocus_POSPOST.txt -v Positive_Post.ori -vl 0.0005 -vu 0.9995 -vt logit -n 0 -pl -1e20 -pu 1e20 -pt -1 -abs 0 -d2r refFlat_exonarray_rowidmap.txt -d2a refFlat_exonarray_transcriptidmap.txt -net HPRD_v6_mouse_enhance.txt -netA geneid2genename_nr.txt -netD 2\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dTruncateLowerBound = -1e20; 
	dTruncateUpperBound = 1e20;
	strcpy(strTransform, "Identity");
	nNormalize = 0;
	dPostLowerBound = -1e20; 
	dPostUpperBound = 1e20;
	strcpy(strPostTransform, "Identity"); 
	nTakeAbsoluteValue = 0;
	strcpy(strNetworkPath, "NULL");
	strcpy(strNetworkAnnotationPath, "NULL");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefLocusDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-v") == 0)
		{
			ni++;
			strcpy(strScorePath, argc[ni]);
			vOK = 1;
		}
		else if(strcmp(argc[ni], "-vl") == 0)
		{
			ni++;
			dTruncateLowerBound = atof(argc[ni]);
			vlOK = 1;
		}
		else if(strcmp(argc[ni], "-vu") == 0)
		{
			ni++;
			dTruncateUpperBound = atof(argc[ni]);
			vuOK = 1;
		}
		else if(strcmp(argc[ni], "-vt") == 0)
		{
			ni++;
			strcpy(strTransform, argc[ni]);
			vtOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nNormalize = atoi(argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-pl") == 0)
		{
			ni++;
			dPostLowerBound = atof(argc[ni]);
			plOK = 1;
		}
		else if(strcmp(argc[ni], "-pu") == 0)
		{
			ni++;
			dPostUpperBound = atof(argc[ni]);
			puOK = 1;
		}
		else if(strcmp(argc[ni], "-pt") == 0)
		{
			ni++;
			strcpy(strPostTransform, argc[ni]);
			ptOK = 1;
		}
		else if(strcmp(argc[ni], "-abs") == 0)
		{
			ni++;
			nTakeAbsoluteValue = atoi(argc[ni]);
			absOK = 1;
		}
		else if(strcmp(argc[ni], "-d2r") == 0)
		{
			ni++;
			strcpy(strRefToProbeRowIDPath, argc[ni]);
			d2rOK = 1;
		}
		else if(strcmp(argc[ni], "-d2a") == 0)
		{
			ni++;
			strcpy(strRefToProbeNamePath, argc[ni]);
			d2aOK = 1;
		}
		else if(strcmp(argc[ni], "-net") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			netOK = 1;
		}
		else if(strcmp(argc[ni], "-netA") == 0)
		{
			ni++;
			strcpy(strNetworkAnnotationPath, argc[ni]);
			netAOK = 1;
		}
		else if(strcmp(argc[ni], "-netD") == 0)
		{
			ni++;
			nNetDepth = atoi(argc[ni]);
			netDOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (sOK == 0) || (oOK == 0) || (vOK == 0) || (vlOK == 0)
		|| (vuOK == 0) || (vtOK == 0) || (nOK == 0) || (plOK == 0) || (puOK == 0)
		|| (ptOK == 0) || (absOK == 0) || (d2rOK == 0) || (d2aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_Database_AssignValues_Main(strRefLocusDatabasePath, 
			strSpecies, strOutPath, strScorePath, nNormalize,
			dTruncateLowerBound, dTruncateUpperBound, strTransform, 
			dPostLowerBound, dPostUpperBound, strPostTransform, nTakeAbsoluteValue,
			strRefToProbeRowIDPath, strRefToProbeNamePath,
			strNetworkPath, strNetworkAnnotationPath, nNetDepth);
	}

	/* return */
	return nResult;
}

int menu_reflocus_assignvalue_test(int argv, char **argc)
{
	RefGene_Database_AssignValues_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		"mouse", "E:\\Projects\\hedgehog_project\\Test\\refLocus_POSPOST.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\Positive_Post.ori", 0,
		0.0005, 0.9995, "logit",
		-1e20, 1e20, "-1", 0,
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_rowidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_transcriptidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\HPRD_v6_mouse_enhance.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\geneid2genename_nr.txt", 2
		);

	/* RefGene_Database_AssignValues_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		"mouse", "E:\\Projects\\hedgehog_project\\Test\\refLocus_M430Neg.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\shh_limb_neg.ori", 0,
		0.0005, 0.9995, "logit",
		-1e20, 1e20, "-1", 0,
		"E:\\Projects\\hedgehog_project\\Test\\refMouse430_rowidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\refMouse430_map.txt",
		"E:\\Projects\\hedgehog_project\\Test\\HPRD_v6_mouse_enhance.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\geneid2genename_nr.txt", 2); */

	/* RefGene_Database_AssignShortest_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		"mouse", "E:\\Projects\\hedgehog_project\\Test\\refLocus_value.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\LimbAlessthanP.ori", 
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_rowidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_transcriptidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\HPRD_v6_mouse_enhance.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\geneid2genename_nr.txt",
		"E:\\Projects\\hedgehog_project\\Test\\exonarraytarget.txt", 10); */

	return 1;
}

int menu_reflocus_assignedgevalue(int argv, char **argc)
{
	/* RefGene_Database_AssignValues_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		"mouse", "E:\\Projects\\hedgehog_project\\Test\\refLocus_GliRepOnly2.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\GliRepOnly_2.ori", 0,
		0.0005, 0.9995, "logit",
		-1e20, 1e20, "-1", 0,
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_rowidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\refFlat_exonarray_transcriptidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\HPRD_v6_mouse_enhance.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\geneid2genename_nr.txt", 2
		); */

	/* RefGene_Database_AssignEdgeValue_Main("E:\\Projects\\hedgehog_project\\Test\\refLocus_sorted.txt", 
		"mouse", "E:\\Projects\\hedgehog_project\\Test\\refLocus_M430Edge_mouse.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\shh_limb_posneg_comb.txt", 0,
		-1e20, 1e20, "identity",
		-1e20, 1e20, "identity", 0,
		"E:\\Projects\\hedgehog_project\\Test\\refMouse430_rowidmap.txt",
		"E:\\Projects\\hedgehog_project\\Test\\refMouse430_map.txt",
		"E:\\Projects\\hedgehog_project\\Test\\mouse_10090_BIND_loc.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\geneid2genename_nr.txt", 0); */

	RefGene_Database_AssignEdgeValue_Main("C:\\Data\\genomes\\human\\b36_hg18\\refLocus_sorted.txt", 
		"human", "C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\Oct-Sox-Nanog\\TTCCCAG\\refLocus_hg18Edge_human.txt", 
		"C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\Oct-Sox-Nanog\\TTCCCAG\\hg17_c40cluster_locnr.txt", 0,
		-1e20, 1e20, "identity",
		-1e20, 1e20, "identity", 0,
		"C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\Oct-Sox-Nanog\\TTCCCAG\\hg18_refFlat_rowid_map.txt",
		"C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\Oct-Sox-Nanog\\TTCCCAG\\hg18_refFlat_locid_map.txt",
		"E:\\Projects\\hedgehog_project\\Test\\human_9606_BIND_loc.txt", 
		"E:\\Projects\\hedgehog_project\\Test\\human_geneid2genename_nr.txt", 0);

	return 1;
}

int menu_affy_bar2txt(int argv, char **argc)
{
	/* ------------------------------- */
	/*        affy_bar2txt             */
	/* convert bar file to txt file    */
	/* ------------------------------- */
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int ni;
	int iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/*        affy_bar2txt             */
	/* -i input file                   */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        affy_bar2txt               \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" example: \n");
		printf("    affy_bar2txt -i input.bar -o output.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Affy_BAR2TXT_Fast(strInputFile, strOutputFile); 
	}

	return nResult;
}

int menu_affy_bar2wig(int argv, char **argc)
{
	/* ------------------------------- */
	/*        affy_bar2wig             */
	/* convert bar file to wig file    */
	/* ------------------------------- */
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int ni;
	int iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/*        affy_bar2wig             */
	/* -i input file                   */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        affy_bar2wig               \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" example: \n");
		printf("    affy_bar2wig -i input.bar -o output.wig\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Affy_BAR2WIG_Main(strInputFile, strOutputFile); 
	}

	return nResult;
}

int menu_affy_getreducedannot()
{
	int nCount;

	nCount = Affy_CSVANNOT_To_Reduced_200408("C:\\Projects\\research_harvard\\affy_project\\data\\NetaffyMm3\\MOE430A_0719\\Mouse430_2_annot.csv", 
		"C:\\Projects\\research_harvard\\genomelab_project\\affymetrix\\mouse\\Mouse430_2_reducedannot_0719.txt", 
		"mouse");

	return nCount;
}

int menu_affy_getrandomcontrolprobesets()
{
	int nCount;
	
	nCount = Affy_PickRandomControlFromCSV("C:\\Projects\\research_harvard\\affy_project\\data\\NetaffyHs\\HG-U133A_annot.csv", 
		"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\egfrm_randomcontrol.txt", 300);

	return nCount;
}

int menu_affy_loadbar_intensity()
{
	int nCount;
	
	/* nCount = Affy_LoadBar_Intensity("C:\\Projects\\research_harvard\\affy_project\\data\\tilingarray\\CD01\\B1_IP_061604.CEL_intensity.bar",
		"C:\\Projects\\research_harvard\\affy_project\\data\\tilingarray\\CD01\\"); */

	nCount = Affy_LoadBar_Intensity("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\IP_8_3A.CEL_intensity.bar",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\");

	return nCount;
}

int menu_affy_loadbar_intensity_group()
{
	int nCount;
	
	/* nCount = Affy_LoadBar_Intensity_Group("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\",
		"cMycChipList.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_norm_PM-MM.txt"); */

	nCount = Affy_LoadBar_Intensity_Group("C:\\Projects\\research_harvard\\estrongen_project\\Raw_supplementary_data\\Carroll_Raw_CEL_files\\",
		"ERChipList.txt",
		"C:\\Projects\\research_harvard\\estrongen_project\\Raw_supplementary_data\\Carroll_Raw_CEL_files\\ER_B_norm_PM.txt");

	return nCount;
}

int menu_affy_loadbpmap_intensity()
{
	int nCount;
	
	nCount = Affy_LoadBPMAP("C:\\Projects\\research_harvard\\affy_project\\data\\tilingarray\\libNCBIv33\\P1_CHIP_A.Anti-Sense.hs.NCBIv33.sary.bpmap",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\ChipA_bpmapv33.txt");

	return nCount;
}

int menu_affy_bpmap_filter_gtrans()
{
	int nCount;
	
	nCount = Affy_BPMAPFilter_GTrans("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_mix_sample2r4_gtrans_3_pvalue.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\ChipA_bpmapv33.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_mix_sample2r4_gtrans_3_pvalue_f.txt");

	return nCount;
}


int menu_affy_bpmap_filter_ori()
{
	int nCount;
	
	nCount = Affy_BPMAPFilter_ORI("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_anova_refvarsh.ori",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\ChipA_bpmapv33.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_anova_refvarsh_f.txt");

	return nCount;
}

int menu_affy_bpmap_filter_data()
{
	int nCount;
	
	nCount = Affy_BPMAPFilter_DATA("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_norm_mat.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\ChipA_bpmapv33.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_norm_mat_f.txt");

	return nCount;
}


int menu_network_shortestpath(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strAnnotPath[LINE_LENGTH];
	char strSourcePath[LINE_LENGTH];	
	char strDestPath[LINE_LENGTH];	
	char strOutPath[LINE_LENGTH];
	int nMaxIteration = 100;
	
	int nResult;
	int nOK;
	int aOK;
	int sOK;
	int dOK;
	int oOK;
	int mOK;
	int ni;

	/* ------------------------------- */
	/*  network_shortestpath           */
	/* -n network path                 */
	/* -a annotation path              */
	/* -s source node path             */
	/* -d dest node path               */
	/* -o output path                  */
	/* -m maximum iteration            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_shortestpath       \n");
		printf(" -n network path            \n");
		printf(" -a annotation path         \n");
		printf(" -s source node path        \n");
		printf(" -d dest node path          \n");
		printf(" -o path for output         \n");
		printf(" -m maximum iteration       \n");
		printf(" example: \n");
		printf("    network_shortestpath -n net.sif -a netannotation.txt -s src.txt -d dst.txt -o out.txt -m 10\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	aOK = 0;
	sOK = 0;
	dOK = 0;
	oOK = 0;
	mOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSourcePath, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDestPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nMaxIteration = atoi(argc[ni]);
			mOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(nMaxIteration <= 0)
	{
		printf("Error: Max iter <= 0!\n");
		exit(EXIT_FAILURE);
	}

	if((nOK == 0) || (aOK == 0) || (sOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_FindShortestPath_Main(strNetworkPath, strAnnotPath,
								  strSourcePath, strDestPath, 
								  strOutPath, nMaxIteration);
	}

	/* return */
	return nResult;
}

int menu_network_createortholognet(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strHomoloPath[LINE_LENGTH];
	int nSrcSpecies,nDestSpecies;
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int nOK;
	int hOK;
	int sOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_createortholognet      */
	/* -n network path                 */
	/* -h homologene path              */
	/* -s source species ID            */
	/* -d dest species ID              */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_createortholognet  \n");
		printf(" -n network path            \n");
		printf(" -h homologene path         \n");
		printf(" -s source species ID       \n");
		printf(" -d dest species ID         \n");
		printf(" -o path for output         \n");
		printf(" example: \n");
		printf("    network_createortholognet -n humannet.sif -h homologene.data -s 9606 -d 10090 -o mousenet.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	hOK = 0;
	sOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-h") == 0)
		{
			ni++;
			strcpy(strHomoloPath, argc[ni]);
			hOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nSrcSpecies = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			nDestSpecies = atoi(argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((nOK == 0) || (hOK == 0) || (sOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_CreateOrthogloNet_Main(strNetworkPath, strHomoloPath,
								  nSrcSpecies, nDestSpecies, strOutPath);
	}
	
	/* return */
	return nResult;
}

int menu_network_getsubnet(int argv, char **argc)
{
	/* define */
	char strNetworkPath[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	char strInPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int nOK;
	int aOK;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_getsubnet              */
	/* -n network path                 */
	/* -a annotation path              */
	/* -i target node path             */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_getsubnet          \n");
		printf(" -n network path            \n");
		printf(" -a annotation path         \n");
		printf(" -i target node path        \n");
		printf(" -o path for output         \n");
		printf(" example: \n");
		printf("    network_getsubnet -n net.sif -h annotation.txt -i targetnode.txt -o subnet.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	nOK = 0;
	aOK = 0;
	iOK = 0;
	oOK = 0;
	strcpy(strAnnotationPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((nOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_GetSubNet_Main(strNetworkPath, strInPath, 
						   aOK, strAnnotationPath, strOutPath);
	}
	
	/* return */
	return nResult;
}

int menu_network_bind2entrez(int argv, char **argc)
{
	/* define */
	char strBindPath[LINE_LENGTH];
	char strEntrezPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	
	int nResult;
	int bOK;
	int eOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*  network_bind2entrez            */
	/* -b BIND path                    */
	/* -e Entrez path                  */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf(" network_bind2entrez        \n");
		printf(" -b BIND path               \n");
		printf(" -e Entrez path             \n");
		printf(" -o output path             \n");
		printf(" example: \n");
		printf("    network_bind2entrez -b bind.dat -e gene2accession.txt -o bind.sif\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	bOK = 0;
	eOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			strcpy(strBindPath, argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			strcpy(strEntrezPath, argc[ni]);
			eOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((bOK == 0) || (eOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Network_LinkBINDEntrez_Main(strBindPath, strEntrezPath, 
						   strOutPath);
	}
	
	/* nResult = Network_LinkBINDEntrez_Main("E:\\Data\\BIND\\GroupbySpecies\\mouse_10090_BIND_nr.txt", 
		"E:\\Data\\BIND\\GroupbySpecies\\mouse_10090_Entrez.txt", 
		"E:\\Data\\BIND\\GroupbySpecies\\mouse_10090_BIND_loc.txt");
	*/

	/* return */
	return nResult;
}

int menu_transloc(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    transloc                     */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    transloc                       \n");
		printf(" example: \n");
		printf("    transloc transloc_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TransLoc_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_expression_rankgenebylocuslink(int argv, char **argc)
{
	int nResult;
	int iOK;
	int ni;
	char strInputPath[LINE_LENGTH];

	/* ------------------------------- */
	/*        powexpressloc            */
	/* -i parameterfile                */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          powexpressloc            \n");
		printf(" -i parameter file      \n");
		printf(" example: \n");
		printf("    powexpressloc -i shhlocuslinkrankinfo.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if(iOK == 0)
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Expression_PowExpressLocusLink_Main(strInputPath);
	}

	/* nCount = Expression_GeneRankByLocusLink_Main("C:\\Projects\\research_harvard\\complexnull_project\\fdrpower\\ebtest\\esloc_locuslinkrankinfo.txt"); */
	/* nCount = Expression_GeneRankByLocusLink_Main("C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_locuslinkrankinfo.txt"); */
	/* nCount = Expression_GeneRankByLocusLink_Main("C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shhlocuslinkrankinfo.txt"); */
	/* nCount = Expression_GeneRankByLocusLink_Main("C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\egfrm_locuslinkrankinfo.txt"); */
	
	/* return */
	return nResult;
}

int menu_expression_getspecificprobe(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nColumn = 0;

	int nResult;
	int dOK;
	int iOK;
	int oOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*  powexpress_getspecificprobe    */
	/* -d raw data                     */
	/* -i input file                   */
	/* -o output file                  */
	/* -c starting column (0-based)    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   powexpress_getspecificprobe     \n");
		printf(" -d raw data                       \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -c starting column (0-based)      \n");
		printf(" example: \n");
		printf("    powexpress_getspecificprobe -d shh.txt -i target.txt -c 1 -o target_affy.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nColumn = atoi(argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Expression_GetSpecificProbe_Main(strDatabasePath, 
			strInputPath, nColumn, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_expression_getnrprobe(int argv, char **argc)
{
	/* define */
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*      expression_getnrprobe      */
	/* -i input file                   */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   powexpress_getnrprobe           \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" example: \n");
		printf("    powexpress_getnrprobe -i target.txt -o target_nr.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Expression_GetNonRedundantProbe_Main(strInputPath, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_expression_geneselection(int argv, char **argc)
{
	int nResult;
	int dOK;
	int aOK;
	int cOK;
	int oOK;
	int ni;
	char strDataPath[LINE_LENGTH];
	char strAnnotationPath[LINE_LENGTH];
	char strCompInfoPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	/* ------------------------------- */
	/*        powexpress               */
	/* -d data file                    */
	/* -a annotation file              */
	/* -c parameter file               */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          powexpress            \n");
		printf(" -d data file      \n");
		printf(" -a annotation file      \n");
		printf(" -c parameter file      \n");
		printf(" -o output file      \n");
		printf(" example: \n");
		printf("    powexpress -d shhdata1.txt -a 3mousechipsinfo.txt -c shhcompinfo1.txt -o shh_8som_pos\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	aOK = 0;
	cOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strAnnotationPath, argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strCompInfoPath, argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (aOK == 0) || (cOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Expression_GeneSelection_Main(strDataPath, strAnnotationPath, strCompInfoPath, strOutPath);
	}

	/* Expression_GeneSelection_Main("C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_data.txt",
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_geneinfo.txt", 
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_compinfo.txt",
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_pos");
	*/

	/* Expression_GeneSelection_Main("C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shhdata1.txt",
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\3mousechipsinfo.txt", 
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shhcompinfo1.txt",
						"C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_8som_pos");
	*/

	/* Expression_GeneSelection_Main("C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\deltaegfrbatchlognorm092004.txt",
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\HG-U133Aannot.txt", 
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\deltaegfrcompinfo.txt",
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\egfrm_vivowd_up");
    */

	/* return */
	return nResult;
}

int menu_expression_quantilenormalization()
{
	/* Expression_Normalization_Quantile_Main(18, 333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_norm.txt", 
		1, 1.0); */

	/* Expression_Normalization_Quantile_Main(2, 372169, 
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\data1\\nimblepilot_fornorm.txt",
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\data1\\nimblepilot_norm.txt", 
		1, 1.0); */

	/* Expression_Normalization_Quantile_Main(14, 372169, 
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\data2\\PairData\\nimble4_fornorm.txt",
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\data2\\PairData\\nimble4_norm.txt", 
		1, 1.0); */

	Expression_Normalization_Quantile_Main(16, 42850, 
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\Agilent\\analysis_step2\\US22502648_RAW_sorted.txt",
		"C:\\Projects\\research_harvard\\hedgehog_project\\analysis\\ChIP-Chip\\Agilent\\analysis_step2\\US22502648_RAW_norm.txt", 
		1, 1.0);
	

	/* Expression_Normalization_Quantile_Main(37, 45025, 
		"C:\\Projects\\research_harvard\\hedgehog_project\\data\\limbbud_steve\\shh_steve_resolver_combined.txt",
		"C:\\Projects\\research_harvard\\hedgehog_project\\data\\limbbud_steve\\shh_steve_resolver_norm.txt", 
		1, 1.0); */

	return(PROC_SUCCESS);
}

int menu_tiling_probeselection()
{
	/* Tiling_ProbeSelection_Main("C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_norm.txt",
		"NULL", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_compinfo.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_sample2naive_Input");
	*/

	Tiling_ProbeSelection_Main("C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_lognorm_anova.txt",
		"NULL", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_mix_compinfo.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMyc_A_anova_refvarsh");


	/* return */
	return PROC_SUCCESS;
}

int menu_tiling_ums_fdr(int argv, char **argc)
{
	int nResult;
	int sOK;
	int dOK;
	int spOK;
	int sqOK;
	int wOK;
	int nOK;
	int oOK;
	int ni;
	char strSelectPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	double dPcut,dQcut;
	int nIntervalNum,nStepSize;

	/* ------------------------------- */
	/*           umsfdr                */
	/* -s selection file               */
	/* -d score file                   */
	/* -sp g0 cutoff                   */
	/* -sq g1 curoff                   */
	/* -w stepsize */
	/* -n interval number              */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          umsfdr            \n");
		printf(" -s selection file   \n");
		printf(" -d score file       \n");
		printf(" -sp g0 cutoff       \n");
		printf(" -sq g1 curoff       \n");
		printf(" -w stepsize         \n");
		printf(" -n interval number  \n");
		printf(" -o output file  \n");
		printf(" example: \n");
		printf("    umsfdr -s selection.txt -d score.txt -sp 0.01 -sq 0.05 -w 1 -n 1000 -o fdr.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	sOK = 0;
	dOK = 0;
	spOK = 0;
	sqOK = 0;
	wOK = 0;
	nOK = 0;
	oOK = 0;
	
	dPcut = 0.01;
	dQcut = 0.05;
	nIntervalNum = 1000;
	nStepSize = 1;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSelectPath, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strScorePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-sp") == 0)
		{
			ni++;
			dPcut = atof(argc[ni]);
			spOK = 1;
		}
		else if(strcmp(argc[ni], "-sq") == 0)
		{
			ni++;
			dQcut = atof(argc[ni]);
			sqOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nStepSize = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nIntervalNum = atoi(argc[ni]);
			nOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((sOK == 0) || (dOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Tiling_UMS_FDR_Main(strSelectPath, strScorePath, 
			dPcut, dQcut, nStepSize, nIntervalNum, strOutPath);
	}

	return nResult;
}

int menu_tiling_bindingcall_hmm()
{
	/* define */
	int nProbeNum;
	double dPrecision;
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strTransitionPath[LINE_LENGTH];
	char strEmissionPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	/* nProbeNum = 330848;
	dPrecision = 0.0001; */

	nProbeNum = 50000;
	dPrecision = 0.001;

	strcpy(strTransform, "identity");
	/* strcpy(strTransform, "invlogit"); */
	
	/* strcpy(strJobName, "cMyc_A_anova_refvarsh_f");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper");
	sprintf(strScorePath, "%s\\%s.txt", strWorkPath, strJobName);
	sprintf(strTransitionPath, "%s\\%s_transitionp.txt", strWorkPath, strJobName);
	sprintf(strEmissionPath, "%s\\%s_emissionp.txt", strWorkPath, strJobName);
	sprintf(strOutPath, "%s\\%s_refbind", strWorkPath, strJobName);
	*/

	strcpy(strJobName, "umstforhmm36");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\ums_simubind_002");
	sprintf(strScorePath, "%s\\%s.txt", strWorkPath, strJobName);
	sprintf(strTransitionPath, "%s\\%s_transitionp.txt", strWorkPath, strJobName);
	sprintf(strEmissionPath, "%s\\%s_emissionp.txt", strWorkPath, strJobName);
	sprintf(strOutPath, "%s\\%s_refbind", strWorkPath, strJobName);

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_probescore.ori",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_bindprob.txt");
	*/

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_ref3naive.ori",
		"logit", 0.0001,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp_ref.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp_ref.txt",
		1000, 0.5,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_refbind");
	*/

	Tiling_BindingRegionSelection_HMM_Main(nProbeNum, 
		strScorePath, strTransform, dPrecision,
		strTransitionPath, strEmissionPath, 
		1000, 0.5,
		strOutPath);

	/* return */
	return PROC_SUCCESS;
}


int menu_tiling_bindingcall_hmm_explen()
{
	/* define */
	int nProbeNum;
	double dPrecision;
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strTransitionPath[LINE_LENGTH];
	char strEmissionPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	nProbeNum = 330848;
	dPrecision = 0.0001;
	strcpy(strTransform, "identity"); 
	/* strcpy(strTransform, "invlogit"); */
	strcpy(strJobName, "cMyc_A_mix_sample3varsh_f");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper");
	sprintf(strScorePath, "%s\\%s.txt", strWorkPath, strJobName);
	sprintf(strTransitionPath, "%s\\%s_translenp.txt", strWorkPath, strJobName);
	sprintf(strEmissionPath, "%s\\%s_emissionp.txt", strWorkPath, strJobName);
	sprintf(strOutPath, "%s\\%s_refbind_explen", strWorkPath, strJobName);

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_probescore.ori",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_bindprob.txt");
	*/

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_ref3naive.ori",
		"logit", 0.0001,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp_ref.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp_ref.txt",
		1000, 0.5,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_refbind");
	*/

	Tiling_BindingRegionSelection_HMM_ExpLen_Main(nProbeNum, 
		strScorePath, strTransform, dPrecision,
		strTransitionPath, strEmissionPath, 
		1000, 0.5,
		strOutPath);

	/* return */
	return PROC_SUCCESS;
}


int menu_tiling_bindingcall_hmm_constlen()
{
	/* define */
	int nProbeNum;
	double dPrecision;
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strTransitionPath[LINE_LENGTH];
	char strEmissionPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	nProbeNum = 330848;
	dPrecision = 0.0001;
	strcpy(strTransform, "identity"); 
	/* strcpy(strTransform, "invlogit"); */
	strcpy(strJobName, "cMyc_A_anova_sample3r4varsh_3_f");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper");
	sprintf(strScorePath, "%s\\%s.txt", strWorkPath, strJobName);
	sprintf(strTransitionPath, "%s\\%s_conslenp.txt", strWorkPath, strJobName);
	sprintf(strEmissionPath, "%s\\%s_emissionp.txt", strWorkPath, strJobName);
	sprintf(strOutPath, "%s\\%s_refbind_conslen", strWorkPath, strJobName);

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_probescore.ori",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp.txt",
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_bindprob.txt");
	*/

	/* Tiling_BindingRegionSelection_HMM_Main(333396, 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_ref3naive.ori",
		"logit", 0.0001,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\transitionp_ref.txt", 
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\emissionp_ref.txt",
		1000, 0.5,
		"C:\\Projects\\research_harvard\\affy_project\\analysis\\tfbs_chr21&22\\chr21_22_refbind");
	*/

	Tiling_BindingRegionSelection_HMM_ConstLen_Main(nProbeNum, 
		strScorePath, strTransform, dPrecision,
		strTransitionPath, strEmissionPath, 
		1000, 0.5,
		strOutPath);

	/* return */
	return PROC_SUCCESS;
}


int menu_tiling_bindingcall_hmm_baumwelch()
{
	/* define */
	int nProbeNum;
	double dPrecision;
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strTransitionPath[LINE_LENGTH];
	char strEmissionPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];

	nProbeNum = 330848;
	dPrecision = 0.0001;
	strcpy(strTransform, "invlogit");
	/* strcpy(strTransform, "identity"); */
	strcpy(strJobName, "cMyc_A_anova_sample2naive_f");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper");
	sprintf(strScorePath, "%s\\%s.txt", strWorkPath, strJobName);
	sprintf(strTransitionPath, "%s\\%s_transitionp.txt", strWorkPath, strJobName);
	sprintf(strEmissionPath, "%s\\%s_emissionp.txt", strWorkPath, strJobName);
	sprintf(strOutPath, "%s\\%s_refbind_EM", strWorkPath, strJobName);

	Tiling_BindingRegionSelection_HMM_BaumWelch_Main(nProbeNum, 
		strScorePath, strTransform, dPrecision,
		strTransitionPath, strEmissionPath,
		1000, 0.9,
		1e-6, 200, 
		strOutPath);

	/* return */
	return PROC_SUCCESS;
}


int menu_tilemap_importaffy(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemap_importaffy           */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemap_importaffy            \n");
		printf(" example: \n");
		printf("    tilemap_importaffy arraylist.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMap_ImportAffy_Normalization_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemap_normalization(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemap_importaffy           */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemap_norm                   \n");
		printf(" example: \n");
		printf("    tilemap_norm tilemap_norm_arg.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMap_Normalization_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemap(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemap_importaffy           */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemap                    \n");
		printf(" example: \n");
		printf("    tilemap tilemap_arg.txt    \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMap_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemap_extract(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemap_importaffy           */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemap                    \n");
		printf(" example: \n");
		printf("    tilemap_extract tilemap_extract_arg.txt    \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMap_Extract_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemapv2_importaffy(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemap_importaffy           */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_importaffy            \n");
		printf(" example: \n");
		printf("    tilemapv2_importaffy tilemapv2_importaffy_arg.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMapv2_ImportAffy_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemapv2(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemapv2                    */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2                      \n");
		printf(" example: \n");
		printf("    tilemapv2 tilemap_arg.txt      \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileMapv2_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemapv2_regioninfo(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	
	char strRegionPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nResult;
	int ni;
	int iOK,oOK,dOK;

	/* ------------------------------- */
	/*    tilemapv2                    */
	/* -i input file                   */
	/* -o output file                  */
	/* -d track info file              */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_regioninfo           \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -d track info file                \n");
		printf(" example: \n");
		printf("    tilemapv2_regioninfo -i region.txt -o output.txt -d tilemap_info_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strRegionPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_RegionInfo_Main(strRegionPath, strParamPath, strOutputPath);
	}
	
	/* return */
	return nResult;
}

int menu_tilemapv2_regioninfo_integral(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	
	char strRegionPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nResult;
	int ni;
	int iOK,oOK,dOK;

	/* ------------------------------- */
	/*    tilemapv2                    */
	/* -i input file                   */
	/* -o output file                  */
	/* -d track info file              */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_regioninfo_integral  \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -d track info file                \n");
		printf(" example: \n");
		printf("    tilemapv2_regioninfo_integral -i region.txt -o output.txt -d tilemap_info_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strRegionPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_RegionInfo_Integral_Main(strRegionPath, strParamPath, strOutputPath);
	}
	
	/* return */
	return nResult;
}

int menu_tilemapv2_collectprobes(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	
	char strRegionPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nResult;
	int ni;
	int iOK,oOK,dOK;

	/* ------------------------------- */
	/*    tilemapv2                    */
	/* -i input file                   */
	/* -o output file                  */
	/* -d track info file              */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_collectprobes        \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -d track info file                \n");
		printf(" example: \n");
		printf("    tilemapv2_collectprobes -i region.txt -o output.txt -d tilemap_info_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strRegionPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_CollectProbes_Main(strRegionPath, strParamPath, strOutputPath);
	}
	
	/* return */
	return nResult;
}


int menu_tilemapv2_txt2bar(int argv, char **argc)
{
	/* ------------------------------- */
	/*     tilemapv2_txt2bar           */
	/* convert bar file to txt file    */
	/* ------------------------------- */
	char strInputFile[LINE_LENGTH];
	char strOutputFolder[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int ni;
	int iOK,oOK,dOK;
	int nResult;

	/* ------------------------------- */
	/*   tilemapv2_txt2bar             */
	/* -i input file                   */
	/* -d output folder                */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_txt2bar               \n");
		printf(" -i input file (full path)         \n");
		printf(" -d output folder (default = current folder) \n");
		printf(" -o output file name               \n");
		printf(" example: \n");
		printf("   tilemapv2_txt2bar -i input.txt -d ./data -o output.cgw\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	strcpy(strOutputFolder, ".");

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_TXT2BAR(strInputFile, strOutputFolder, strOutputFile); 
	}

	return nResult;
}

int menu_tilemapv2_quantilenorm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tilemapv2_quantilenorm        */
	/* quantile normalizatio of a txt  */
	/* ------------------------------- */
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int iOK,oOK,cOK,lOK,uOK,tOK;
	int nResult;
	int nSkipCol = 2;
	double dL = -1000000.0;
	double dU = 1000000.0;
	int nTransform = 0;

	/* ------------------------------- */
	/*  tilemapv2_quantilenorm         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_quantilenorm         \n");
		printf(" -i input file (full path)         \n");
		printf(" -o output file name               \n");
		printf(" -c number of columns to skip (default = 2) \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -1000000) \n");
		printf(" -u truncation upper bound, values bigger than the upper bound will be truncated to the upper bound (default = 1000000) \n");
		printf(" -t transformation (0:identity(default), 1:log2). Transformation will be performed after truncation \n");
		printf(" example: \n");
		printf("   tilemapv2_quantilenorm -i input.txt -o output.txt -c 1 -l 1 -t l\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	cOK = 0;
	lOK = 0;
	uOK = 0;
	tOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nSkipCol = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			dU = atof(argc[ni]);
			uOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_TXT_QuantileNormalization(strInFile, strOutFile,
				nSkipCol, nTransform, dL, dU);

	}

	return nResult;
}

int menu_tileprobe_buildmodel(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	int nShrink = 0;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,sOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel           \n");
		printf(" -i input file that contains the array list for building the background probe model. \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transformation (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -s 1: shrink IQR; 0 (default): not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	sOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_Build_Main(strInFile, strOutFile, nTransform, dL, nShrink, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildmodel_v2(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	int nShrink = 1;
	int nLogAfterNorm = 0;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,sOK,eOK,aOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel_v2        \n");
		printf(" -i input file that contains the array list for building the background probe model. \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transform (0:identity, 1:log2 (default)). Transform will be performed after truncation \n");
		printf(" -a 1: transform after quantile normalization; 0: transform before quantile normalization (default) \n");
		printf(" -s 1: shrink probe variance (default); 0: not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel_v2 -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	aOK = 0;
	sOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nLogAfterNorm = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_Buildv2_Main(strInFile, strOutFile, nTransform, dL, nShrink, nLogAfterNorm, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildmodel_hmmt(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	double dTopPrc = 0.005;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,cOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel_hmmt        \n");
		printf(" -i input file that contains the array list for building the background probe model. (col1=0: control; col1=1: IP) \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transformation (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -c excluding top c*100% probes from the IP arrays (default = 0.005) \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel_hmmt -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	cOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dTopPrc = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BuildHMMT_Main(strInFile, strOutFile, nTransform, dL, dTopPrc, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildmodel_hmmb(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	double dTopPrc = 0.005;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,cOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel_hmmb        \n");
		printf(" -i input file that contains the array list for building the background probe model. (col1=0: control; col1=1: IP) \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transformation (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -c excluding top c*100% probes from the IP arrays (default = 0.005) \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel_hmmb -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	cOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dTopPrc = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BuildHMMB_Main(strInFile, strOutFile, nTransform, dL, dTopPrc, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildmodel_hmmm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	double dTopPrc = 0.005;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,cOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel_hmmm        \n");
		printf(" -i input file that contains the array list for building the background probe model. (col1=0: control; col1=1: IP) \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transformation (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -c excluding top c*100% probes from the IP arrays (default = 0.005) \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel_hmmm -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	cOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dTopPrc = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BuildHMMM_Main(strInFile, strOutFile, nTransform, dL, dTopPrc, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildmodel_hmmw(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = 1.0;
	int nTransform = 1;
	int nShrink = 1;
	double dTopPrc = 0.005;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,sOK,tOK,cOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildmodel_hmmw        \n");
		printf(" -i input file that contains the array list for building the background probe model. (col1=0: control; col1=1: IP) \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transformation (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -c excluding top c*100% probes from the IP arrays (default = 0.005) \n");
		printf(" -s 1: shrink probe variance (default); 0: not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildmodel_hmmw -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	sOK = 0;
	tOK = 0;
	cOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dTopPrc = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BuildHMMW_Main(strInFile, strOutFile, nTransform, dL, dTopPrc, nShrink, nTest);
	}

	return nResult;
}

int menu_tileprobe_norm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_norm                */
	/* normalization based on probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strModelFile[MED_LINE_LENGTH];
	int nInputType = 0;
	double dL = 1.0;
	int nTransform = 1;
	int nLogAfterNorm = 0;
	double dB = 0.0;

	int ni;
	int iOK,oOK,mOK,cOK,lOK,tOK,aOK,bOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_norm                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_norm           \n");
		printf(" -i input cel file or file list \n");
		printf(" -o output folder               \n");
		printf(" -m path and title of the file that contains probe effect models \n");
		printf(" -c 0 (default): input is a single cel file; 1: input is a file that contains a list of array cel files. \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = 1) \n");
		printf(" -t transform (0:identity, 1:log2 (default)). Transformation will be performed after truncation \n");
		printf(" -a 1: transform after quantile normalization; 0: transform before quantile normalization (default).\n");
		printf(" -b shrinkage factor (if <0, no shrinking; 0~1, shrinking; >1 (for future support, automatic shrinking) \n");
		printf(" example: \n");
		printf("   tileprobe_norm -i input.cel -o /home/normalized/ -m MouseProm\n");
		printf("   tileprobe_norm -i input_list.txt -c 1 -o /home/normalized/ -m MouseProm -l -100000 -t 0  -b 0.1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	mOK = 0;
	cOK = 0;
	lOK = 0;
	tOK = 0;
	aOK = 0;
	bOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strModelFile, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nLogAfterNorm = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			dB = atof(argc[ni]);
			bOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (mOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_Norm_Main(strInFile, nInputType, strOutFile, strModelFile, nTransform, dL, dB, nLogAfterNorm);
	}

	return nResult;
}

int menu_tileprobe_peak(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tilemapv2                    */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_peak                 \n");
		printf(" example: \n");
		printf("    tileprobe_peak tileprobe_peak_arg.txt  \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileProbe_Peak_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tileprobe_mat(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tileprobe_mat                */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_mat                  \n");
		printf(" example: \n");
		printf("    tileprobe_mat tileprobe_mat_arg.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileProbe_MATv2_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tileprobe_huberdata(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tileprobe_huberdata          */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_huberdata                 \n");
		printf(" example: \n");
		printf("    tileprobe_huberdata tileprobe_mat_arg.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileProbe_HuberData_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tileprobe_buildbarmodel(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildbarmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = -100000000.0;
	int nTransform = 0;
	int nShrink = 0;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,sOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildbarmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildbarmodel           \n");
		printf(" -i input file that contains the array list for building the background probe model. \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -100000000.0) \n");
		printf(" -t transformation (0=default:identity, 1:log2). Transformation will be performed after truncation \n");
		printf(" -s 1: shrink IQR; 0 (default): not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildbarmodel -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	sOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BARBuild_Main(strInFile, strOutFile, nTransform, dL, nShrink, nTest);
	}

	return nResult;
}

int menu_tileprobe_buildbarmodel_v2(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildbarmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = -100000000.0;
	int nTransform = 0;
	int nShrink = 1;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,sOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildbarmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildbarmodel_v2     \n");
		printf(" -i input file that contains the array list for building the background probe model. \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -100000000.0) \n");
		printf(" -t transformation (0=default:identity, 1:log2). Transformation will be performed after truncation \n");
		printf(" -s 1: shrink variance (default); 0: not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildbarmodel_v2 -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	sOK = 0;
	eOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BARBuildv2_Main(strInFile, strOutFile, nTransform, dL, nShrink, nTest);
	}

	return nResult;
}

int menu_tileprobe_barnorm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_barnorm             */
	/* normalization based on probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strModelFile[MED_LINE_LENGTH];
	int nInputType = 0;
	double dL = -100000000.0;
	int nTransform = 0;
	double dB = 0.0;

	int ni;
	int iOK,oOK,mOK,cOK,lOK,tOK,bOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_norm                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_barnorm           \n");
		printf(" -i input bar file or file list \n");
		printf(" -o output folder               \n");
		printf(" -m path and title of the file that contains probe effect models \n");
		printf(" -c 0 (default): input is a single bar file; 1: input is a file that contains a list of array bar files. \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -100000000.0) \n");
		printf(" -t transformation (0:identity (default), 1:log2). Transformation will be performed after truncation \n");
		printf(" -b shrinkage factor (if <0, no shrinking; 0~1, shrinking; >1 (for future support, automatic shrinking) \n");
		printf(" example: \n");
		printf("   tileprobe_barnorm -i input.bar -o /home/normalized/ -m MouseProm\n");
		printf("   tileprobe_barnorm -i input_list.txt -c 1 -o /home/normalized/ -m MouseProm -l -100000 -t 0  -b 0.1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	mOK = 0;
	cOK = 0;
	lOK = 0;
	tOK = 0;
	bOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strModelFile, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			dB = atof(argc[ni]);
			bOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (mOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BARNorm_Main(strInFile, nInputType, strOutFile, strModelFile, nTransform, dL, dB);
	}

	return nResult;
}

int menu_tileprobe_resample(int argv, char **argc)
{
		/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    tileprobe_huberdata          */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_resample                 \n");
		printf(" example: \n");
		printf("    tileprobe_resample tileprobe_resample_arg.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = TileProbe_Resample_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_tilemap_proberemap(int argv, char **argc)
{
	/* define */
	int nResult = 0;

	TileMapv2_ProbeReMapping_Main("F:\\Data\\genomes\\mouse\\mm8\\Mm35b_P01R_v02-1_NCBIv33.bpmap", "F:\\Data\\genomes\\mouse\\mm8",
						"F:\\Data\\genomes\\mouse\\mm8\\Mm35b_P01R_v02-1_mm8.bpmap", "mouse", 13);
	/* Genome_CreateHash_Chr("C:\\Data\\genomes\\mouse\\mm8\\", "chr19", "C:\\Data\\genomes\\mouse\\mm8\\",
				61321190, 4, 13); */

	/* return */
	return nResult;
}


int menu_mapcMyc()
{
	Map_cMyc_Chr2122();

	return 1;
}

int menu_hts_aln2uniq(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nMax = 1;
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        hts_aln2uniq             */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          hts_aln2uniq           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -m maximal read no. per locus (default = 1) \n");
		printf(" example: \n");
		printf("    hts_aln2uniq -i input.aln -o output.aln -m 5 \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nMax = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}


	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_Aln2Unique(strInputPath, strOutputPath, nMax);
	}

	return nResult;
}

int menu_hts_aln2window(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    hts_aln2window               */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    hts_aln2window                 \n");
		printf(" example: \n");
		printf("    hts_aln2window hts_aln2window_arg.txt  \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = HTS_Aln2Window_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_hts_aln2diff(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    hts_aln2diff                 */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    hts_aln2diff                   \n");
		printf(" example: \n");
		printf("    hts_aln2window hts_aln2window_arg.txt  \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = HTS_Aln2Diff_Main(strParamPath); 
	/* nResult = HTS_Aln2Enrich_Main(strParamPath); */
	
	/* return */
	return nResult;
}

int menu_hts_aln2bar(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        hts_aln2bar              */
	/* -i input                        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          hts_aln2bar           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" example: \n");
		printf("    hts_aln2bar -i input.txt -o output.bar \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_Aln2BAR(strInputPath, strOutputPath);
	}

	return nResult;
}

int menu_hts_aln2barv2(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        hts_aln2barv2            */
	/* -i input                        */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          hts_aln2barv2          \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" example: \n");
		printf("    hts_aln2barv2 -i input.txt -o output.bar \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_Aln2BARv2(strInputPath, strOutputPath);
	}

	return nResult;
}

int menu_hts_aln2winbar(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strBARHeader[LINE_LENGTH];
	int nExtLen = 0;
	int nBinSize = 25;
	char strSpecies[LINE_LENGTH];
	char strChrLenFile[MED_LINE_LENGTH];
	int nStrand = 0;
	int nResult;
	
	int iOK;
	int dOK;
	int oOK;
	int sOK;
	int lOK;

	int ni;

	/* ------------------------------- */
	/*        hts_aln2winbar           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          hts_aln2winbar          \n");
		printf(" -i input aln file \n");
		printf(" -d output folder  \n");
		printf(" -o output bar file header \n");
		printf(" -e extension length \n");
		printf(" -b bin size \n");
		printf(" -strand 1: produce strand specific curves; 0 (default): do not produce strand profiles. \n");
		printf(" -s species \n");
		printf(" -l chromosome length list \n");
		printf(" example: \n");
		printf("    hts_aln2winbar -i input.aln -d . -o output -e 100 -b 25 -s human -l chrlen.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	sOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strBARHeader, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-strand") == 0)
		{
			ni++;
			nStrand = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) || (sOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_Aln2WinBAR(strInputFile, strOutputFolder, strBARHeader, 
				   nExtLen, nBinSize, nStrand,
				   strSpecies, strChrLenFile);
	}

	return nResult;
}

int menu_hts_alnshift2bar(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int sOK;
	int nS = 0;
	int ni;

	/* ------------------------------- */
	/*       hts_alnshift2bar          */
	/* -i input                        */
	/* -o output                       */
	/* -s offset                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("         hts_alnshift2bar          \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -s shift half window (bp)  \n");
		printf(" example: \n");
		printf("    hts_alnshift2bar -i input.bar -o output.bar -s 35\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	sOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(iOK == 0)
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(oOK == 0)
		{
			sprintf(strOutputPath, "%s_C.bar", strInputPath);
		}
		nResult = HTS_AlnShift2BAR(strInputPath, strOutputPath, nS);
	}

	return nResult;
}

int menu_hts_windowsummary(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int ni;
	struct tagBARData *pRepeatData = NULL;

	/* ------------------------------- */
	/*     hts_windowsummary           */
	/* -i input                        */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        hts_windowsummary          \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" example: \n");
		printf("    hts_windowsummary -i input.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_WindowSummary(strInputPath, strChrList, strChrLen, nW, strOutputPath, pRepeatData);
	}

	return nResult;
}

int menu_hts_windowsummaryv2(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nZ = 0;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	char strRepeatFile[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int zOK;
	int mOK;
	int ni;

	/* ------------------------------- */
	/*     hts_windowsummaryv2         */
	/* -i input                        */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -z use combined data after shifting */
	/* -m mask repeats                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        hts_windowsummaryv2        \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -z use combined data after shifting (default = 0) \n");
		printf(" -m mask windows specified in a BAR file which usually represents repeat regions (can be skipped) \n");
		printf(" example: \n");
		printf("    hts_windowsummaryv2 -i input.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt -z 1 -m repeat.bar\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;
	mOK = 0;
	strcpy(strRepeatFile, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strRepeatFile, argc[ni]);
			mOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_WindowSummaryv2(strInputPath, strChrList, strChrLen, nW, strOutputPath, nZ, strRepeatFile);
	}

	return nResult;
}

int menu_hts_windowsummarypaper(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nZ = 0;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int zOK;
	int ni;

	/* ------------------------------- */
	/*    hts_windowsummarypaper       */
	/* -i input                        */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -z use combined data after shifting */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        hts_windowsummarypaper      \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -z use combined data after shifting (default = 0) \n");
		printf(" example: \n");
		printf("    hts_windowsummarypaper -i input.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt -z 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_WindowSummaryPaper(strInputPath, strChrList, strChrLen, nW, strOutputPath, nZ);
	}

	return nResult;
}

int menu_hts_createrepeatfilter(int argv, char **argc)
{
	char strChrListFile[MED_LINE_LENGTH];
	char strChrLenFile[MED_LINE_LENGTH];
	char strGenomeFolder[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	int nW = 100;
	double dR = 0.5;

	int nResult;
	int iOK;
	int lOK;
	int gdOK;
	int dOK;
	int oOK;
	int wOK;
	int rOK;
	int ni;

	/* ------------------------------- */
	/*   hts_createrepeatfilter        */
	/* -i input chrlist file           */
	/* -l input chrlen file			   */
	/* -gd genome sequence folder      */
	/* -d output folder                */
	/* -o output title                 */
	/* -w window size (default = 100)  */
	/* -r repeat cutoff, if repeat percentage>=r, then discard. */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_createrepeatfilter        \n");
		printf(" -i input chrlist file     \n");
		printf(" -l input chrlen file           \n");
		printf(" -gd genome sequence folder           \n");
		printf(" -d output folder  \n");
		printf(" -o output title    \n");
		printf(" -w window size (default = 100)        \n");
		printf(" -r repeat cutoff, if repeat percentage>=r, then discard. \n");
		printf(" example: \n");
		printf("    hts_createrepeatfilter -i chrlist.txt -l chrlen.txt -gd /home/hg17 -d . -o output -w 100 -r 0.5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	lOK = 0;
	gdOK = 0;
	dOK = 0;
	oOK = 0;
	wOK = 0;
	rOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strChrListFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomeFolder, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (lOK == 0) || (gdOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100.\n");
		}
		
		nResult = HTS_CreateRepeatFilter(strGenomeFolder, strChrListFile, strChrLenFile, 
					 nW, dR, strOutputFolder, strOutputTitle);
	}

	return nResult;
}

int menu_hts_filterrepeatreads(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strGenomeFolder[MED_LINE_LENGTH];

	int nResult;
	int iOK;
	int gdOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*   hts_filterrepeatreads         */
	/* -i input file                   */
	/* -o output title                 */
	/* -gd genome sequence folder      */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_filterrepeatreads         \n");
		printf(" -i input file     \n");
		printf(" -o output file    \n");
		printf(" -gd genome sequence folder           \n");
		printf(" example: \n");
		printf("    hts_filterrepeatreads -i input.txt -o output.txt -gd /home/hg17\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	gdOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomeFolder, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (gdOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_FilterRepeatReads(strInputFile, strOutputFile, strGenomeFolder);
	}

	return nResult;
}

int menu_hts_collectreads(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strBARFile[MED_LINE_LENGTH];

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*      hts_collectreads           */
	/* -i input region file            */
	/* -d input read bar file.         */
	/* -o output title                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_collectreads           \n");
		printf(" -i input region file     \n");
		printf(" -d input alignment bar file        \n");
		printf(" -o output file           \n");
		printf(" example: \n");
		printf("    hts_collectreads -i input.txt -o output.txt -d chip.aln.bar \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strBARFile, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_CollectReads(strInputFile, strOutputFile, strBARFile);
	}

	return nResult;
}

int menu_hts_collectprofile(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strBARFile[MED_LINE_LENGTH];
	int nW = 200;
	int nS = 20;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*      hts_collectprofile         */
	/* -i input region file            */
	/* -d input read bar file.         */
	/* -o output title                 */
	/* -w half window size             */
	/* -s step size                    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_collectprofile           \n");
		printf(" -i input region file     \n");
		printf(" -d input alignment bar file        \n");
		printf(" -o output file           \n");
		printf(" -w half window size (default = 200 bp)     \n");
		printf(" -s step size (default = 20 bp)            \n");
		printf(" example: \n");
		printf("    hts_collectprofile -i input.txt -o output.txt -d chip.aln.bar -w 200 -s 50\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strBARFile, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW < nS)
		{
			printf("Warning: nW cannot be smaller than nS. nW is set to nS!\n");
			nW = nS;
		}
		nResult = HTS_CollectProfile_Main(strInputFile, strOutputFile, strBARFile, nW, nS);
	}

	return nResult;
}

int menu_hts_onesample_enrich(int argv, char **argc)
{
	char strInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	int nW = 100;
	int nS = 25;
	int nCut = 10;
	int nMinLen = 1;
	int nMaxGap = 0;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int wOK;
	int sOK;
	int cOK;
	int gOK;
	int lOK;
	int ni;

	/* ------------------------------- */
	/*    hts_onesample_enrich         */
	/* -i input                        */
	/* -d output folder                */
	/* -o output title                 */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -c cutoff (default = 10)        */
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("       hts_peakdetector            \n");
		printf(" -i input      \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -c cutoff (default = 10)        \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" example: \n");
		printf("    hts_peakdetector -i input.bar -d . -o output -w 200 -s 25 -c 10 -g 50 -l 500\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	wOK = 0;
	sOK = 0;
	cOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nCut = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100.\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25.\n");
		}
		if(nCut <= 0)
		{
			nCut = 10;
			printf("Warning: Cutoff<=0! Proceed with default cutoff = 10.\n");
		}

		nResult =  HTS_Enrich_OneSample_Main(strInputPath, nW, nS, nCut,
					  nMinLen, nMaxGap, strOutputFolder, strOutputTitle);
	}

	return nResult;
}

int menu_hts_onesamplev2_enrich(int argv, char **argc)
{
	char strInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	int nW = 100;
	int nS = 25;
	int nCut = 10;
	int nCutF = 0;
	int nCutR = 0;
	int nMinLen = 1;
	int nMaxGap = 0;
	int nBR = 0; 
	int nBRL = 30;
	int nSSF = 0;
	int nZ = 0;
	
	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int wOK;
	int sOK;
	int cOK;
	int cfOK;
	int crOK;
	int gOK;
	int lOK;
	int zOK;
	int ni;

	/* ------------------------------- */
	/*    hts_onesample_enrich         */
	/* -i input                        */
	/* -d output folder                */
	/* -o output title                 */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -c cutoff (default = 10)        */
	/* -cf forward cutoff (default = 10 */
	/* -cr backward cutoff (default = 10 */
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* -z use combined data after shifting (default = 0) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("       hts_peakdetectorv2          \n");
		printf(" -i input      \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -c cutoff (default = 10)        \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" -br apply boundary refinement (0[default]: no; 1: yes) \n");
		printf(" -brl boundary refinement min region length (default=30) \n");
		printf(" -ssf apply single strand filtering (0[default]: no; 1: yes) \n");
		printf(" -cf (currently only for developer's use) forward cutoff (default = c)  \n");
		printf(" -cr (currently only for developer's use) reverse cutoff (default = c)  \n");
		printf(" -z use combined data after shifting (default = 0) \n");
		printf(" example: \n");
		printf("    hts_peakdetectorv2 -i input.bar -d . -o output -w 200 -s 25 -c 10 -g 50 -l 500 -br 1 -brl 30 -ssf 1 -cf 8 -cr 8 -z 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	wOK = 0;
	sOK = 0;
	cOK = 0;
	cfOK = 0;
	crOK = 0;
	gOK = 0;
	lOK = 0;
	zOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nCut = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cf") == 0)
		{
			ni++;
			nCutF = atoi(argc[ni]);
			cfOK = 1;
		}
		else if(strcmp(argc[ni], "-cr") == 0)
		{
			ni++;
			nCutR = atoi(argc[ni]);
			crOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-br") == 0)
		{
			ni++;
			nBR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-brl") == 0)
		{
			ni++;
			nBRL = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ssf") == 0)
		{
			ni++;
			nSSF = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100.\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25.\n");
		}
		if(nCut <= 0)
		{
			nCut = 10;
			printf("Warning: Cutoff<=0! Proceed with default cutoff = 10.\n");
		}
		if(nCutF <= 0)
		{
			nCutF = nCut;
			/* printf("Warning: Forward Cutoff<=0! Proceed with default cutoff = C.\n"); */
		}
		if(nCutR <= 0)
		{
			nCutR = nCut;
			/* printf("Warning: Reverse Cutoff<=0! Proceed with default cutoff = C.\n"); */
		}
		if(nBRL <= 0)
		{
			nBRL = 30;
		}

		nResult =  HTS_Enrich_OneSamplev2_Main(strInputPath, nW, nS, nCut, nCutF, nCutR,
					  nMinLen, nMaxGap, strOutputFolder, strOutputTitle,
					  nBR, nBRL, nSSF, nZ);
	}

	return nResult;
}

int menu_hts_twosample_windowsummary(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strNegInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int nOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int ni;

	/* ------------------------------- */
	/*  hts_windowsummary_twosample    */
	/* -i input1 (positive sample)     */
	/* -n input2 (negative sample)     */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    hts_windowsummary_twosample    \n");
		printf(" -i input (positive)     \n");
		printf(" -n input (negative)     \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" example: \n");
		printf("    hts_windowsummary_twosample -i posinput.bar -n neginput.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_TwoSample_WindowSummary(strInputPath, strNegInputPath, strChrList, strChrLen, nW, strOutputPath);
	}

	return nResult;
}

int menu_hts_twosample_windowsummaryv2(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strNegInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nZ = 0;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int nOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int zOK;
	int ni;

	/* ------------------------------- */
	/*  hts_windowsummary_twosamplev2  */
	/* -i input1 (positive sample)     */
	/* -n input2 (negative sample)     */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -z use combined data after shifting */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    hts_windowsummaryv2_2sample    \n");
		printf(" -i input (positive)     \n");
		printf(" -n input (negative)     \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -z use combined data after shifting (default = 0) \n");
		printf(" example: \n");
		printf("    hts_windowsummaryv2_2sample -i posinput.bar -n neginput.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt -z 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;
	zOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_TwoSample_WindowSummaryv2(strInputPath, strNegInputPath, strChrList, strChrLen, nW, strOutputPath, nZ);
	}

	return nResult;
}

int menu_hts_twosample_enrich(int argv, char **argc)
{
	char strPosInputPath[MED_LINE_LENGTH];
	char strNegInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	char strFDRPath[MED_LINE_LENGTH];
	int nOneSide = 1;
	int nW = 100;
	int nS = 25;
	double dFDRCut = 0.1;
	int nMinLen = 1;
	int nMaxGap = 0;
	int nTCut = 10;
	double dP0 = 0.5;
	
	int nResult;
	int iOK;
	int nOK;
	int dOK;
	int oOK;
	int tOK;
	int wOK;
	int sOK;
	int fOK;
	int cOK;
	int mOK;
	int pOK;
	int gOK;
	int lOK;
	int ni;

	  
	/* ------------------------------- */
	/*    hts_peakdetector_2sample     */
	/* -i input (positive)             */
	/* -n input (negative)             */
	/* -d output folder                */
	/* -o output title                 */
	/* -t one-sided or two-sided test  */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -f FDR file                     */
	/* -c fdr cutoff (default = 0.1)   */
	/* -m min window read count        */
	/* -p p0 = pos/neg background ratio*/
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_peakdetector_2sample      \n");
		printf(" -i input (positive)          \n");
		printf(" -n input (negative)          \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -t 1: one-sided or 0: two-sided test (default=1) \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -f FDR file                     \n");
		printf(" -c fdr cutoff (default = 0.1)   \n");
		printf(" -m min window read_num (default = 10) \n");
		printf(" -p p0 = pos/neg background ratio \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" example: \n");
		printf("    hts_peakdetector_2sample -i pos.bar -n neg.bar -d . -o output -f FDR.txt -c 0.2 -p 0.5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	dOK = 0;
	oOK = 0;
	tOK = 0;
	wOK = 0;
	sOK = 0;
	fOK = 0;
	cOK = 0;
	mOK = 0;
	pOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strPosInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nOneSide = atoi(argc[ni]);
			tOK = 1;

			if(nOneSide != 0)
				nOneSide = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strFDRPath, argc[ni]);
			fOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dFDRCut = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nTCut = atoi(argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-p") == 0)
		{
			ni++;
			dP0 = atof(argc[ni]);
			pOK = 1;

			if( (dP0 <= 0.0) || (dP0 >= 1.0) )
			{
				printf("Error: -p option must be set to a value between 0 and 1 (i.e. 0<p<1)!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (dOK == 0) || (oOK == 0) || (fOK == 0) || (pOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25\n");
		}
		if(nTCut <= 0)
		{
			nTCut = 10;
			printf("Warning: Window read_num cutoff <= 0! Proceed with default window read_num cutoff = 10\n");
		}
		if(dFDRCut <= 0.0)
		{
			dFDRCut = 0.1;
			printf("Warning: FDR cutoff <= 0! Proceed with default fdr cutoff = 0.1\n");
		}
	

		nResult =  HTS_Enrich_TwoSample_Main(strPosInputPath, strNegInputPath, 
					  nOneSide, nW, nS, nTCut, strFDRPath, dFDRCut,
					  nMinLen, nMaxGap, dP0,
					  strOutputFolder, strOutputTitle);
	}

	return nResult;
}

int menu_hts_twosamplev2_enrich(int argv, char **argc)
{
	char strPosInputPath[MED_LINE_LENGTH];
	char strNegInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	char strFDRPath[MED_LINE_LENGTH];
	int nOneSide = 1;
	int nW = 100;
	int nS = 25;
	double dFDRCut = 0.1;
	int nMinLen = 1;
	int nMaxGap = 0;
	int nTCut = 10;
	double dP0 = 0.5;
	int nBR = 0;
	int nBRL = 30; 
	int nSSF = 0;
	int nSSFF = 0;
	int nSSFR = 0;
	int nZ = 0;
	double dFC = 0.0;
	double dTFC = 0.0;
	
	int nResult;
	int iOK;
	int nOK;
	int dOK;
	int oOK;
	int tOK;
	int wOK;
	int sOK;
	int fOK;
	int cOK;
	int mOK;
	int pOK;
	int gOK;
	int lOK;
	int zOK;
	int ni;

	  
	/* ------------------------------- */
	/*    hts_peakdetector_2samplev2   */
	/* -i input (positive)             */
	/* -n input (negative)             */
	/* -d output folder                */
	/* -o output title                 */
	/* -t one-sided or two-sided test  */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -f FDR file                     */
	/* -c fdr cutoff (default = 0.1)   */
	/* -m min window read count        */
	/* -p p0 = pos/neg background ratio*/
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* -z use combined data after shifting (default = 0) */
	/* -fc minimal fold change         */
	/* -tfc minimal region pos_read/neg_read ratio */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_peakdetectorv2_2sample    \n");
		printf(" -i input (positive)          \n");
		printf(" -n input (negative)          \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -t 1: one-sided or 0: two-sided test (default=1) \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -f FDR file                     \n");
		printf(" -c fdr cutoff (default = 0.1)   \n");
		printf(" -m min window read_num (default = 10) \n");
		printf(" -p p0 = pos/neg background ratio \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" -br apply boundary refinement (0[default]: no; 1: yes) \n");
		printf(" -brl boundary refinement min region length (default=30) \n");
		printf(" -ssf apply single strand filtering (0[default]: no; 1: yes) \n");
		printf(" -cf (currently only for developer's use) forward cutoff (default = c)  \n");
		printf(" -cr (currently only for developer's use) reverse cutoff (default = c)  \n");
		printf(" -z use combined data after shifting (default = 0)  \n");
		printf(" -fc minimal fold change  \n");
		printf(" -tfc minimal region pos_read/neg_read ratio  \n");
		printf(" example: \n");
		printf("    hts_peakdetectorv2_2sample -i pos.bar -n neg.bar -d . -o output -f FDR.txt -c 0.2 -p 0.5 -br 1 -brl 30 -ssf 1 -cf 8 -cr 8 -z 1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	dOK = 0;
	oOK = 0;
	tOK = 0;
	wOK = 0;
	sOK = 0;
	fOK = 0;
	cOK = 0;
	mOK = 0;
	pOK = 0;
	gOK = 0;
	lOK = 0;
	zOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strPosInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nOneSide = atoi(argc[ni]);
			tOK = 1;

			if(nOneSide != 0)
				nOneSide = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strFDRPath, argc[ni]);
			fOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dFDRCut = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nTCut = atoi(argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-p") == 0)
		{
			ni++;
			dP0 = atof(argc[ni]);
			pOK = 1;

			if( (dP0 <= 0.0) || (dP0 >= 1.0) )
			{
				printf("Error: -p option must be set to a value between 0 and 1 (i.e. 0<p<1)!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-br") == 0)
		{
			ni++;
			nBR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-brl") == 0)
		{
			ni++;
			nBRL = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ssf") == 0)
		{
			ni++;
			nSSF = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-cf") == 0)
		{
			ni++;
			nSSFF = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-cr") == 0)
		{
			ni++;
			nSSFR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else if(strcmp(argc[ni], "-fc") == 0)
		{
			ni++;
			dFC = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-tfc") == 0)
		{
			ni++;
			dTFC = atof(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (dOK == 0) || (oOK == 0) || (fOK == 0) || (pOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25\n");
		}
		if(nBRL <= 0)
		{
			nBRL = 30;
		}
		if(nTCut <= 0)
		{
			nTCut = 10;
			printf("Warning: Window read_num cutoff <= 0! Proceed with default window read_num cutoff = 10\n");
		}
		if(nSSFF <= 0)
		{
			nSSFF = nTCut;
		}
		if(nSSFR <= 0)
		{
			nSSFR = nTCut;
		}
		if(dFDRCut <= 0.0)
		{
			dFDRCut = 0.1;
			printf("Warning: FDR cutoff <= 0! Proceed with default fdr cutoff = 0.1\n");
		}
	

		nResult =  HTS_Enrich_TwoSamplev2_Main(strPosInputPath, strNegInputPath, 
					  nOneSide, nW, nS, nTCut, strFDRPath, dFDRCut,
					  nMinLen, nMaxGap, dP0,
					  strOutputFolder, strOutputTitle,
					  nBR, nBRL, nSSF, nSSFF, nSSFR, nZ, dFC, dTFC);
	}

	return nResult;
}

int menu_hts_countreads4refgene(int argv, char **argc)
{
	/* define */
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 1;
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nRefType = 0;
	int nUP = 0;
	int nDOWN = 0;
	int nInputType = 0; 
	int nStandardizebyT = 1;
	int nStandardizebyL = 1;

	int nResult;
	int dOK;
	int dtOK;
	int sOK;
	int iOK;
	int itOK;
	int oOK;
	int rOK;
	int upOK;
	int downOK;
	int ni;

	/* ------------------------------- */
	/*      hts_countreads4refgene     */
	/* -d database                     */
	/* -s species                      */
	/* -dt database type 0: refGene    */
	/*     1: refFlat (default)        */
	/*     2: refLocus                 */
	/* -i input bar file               */
	/* -it input type. 0: bar file     */
	/*     1: a txt file listing input bar files */
	/* -o output text file             */
	/* -r reference type               */
	/*    0: TSS-up, TES-down          */
	/*    1: TSS-up, TSS-down          */
	/*    2: TES-up, TES-down          */
	/*    3: CDSS-up, CDSE-down        */
	/*    4: CDSS-up, CDSS-down        */
	/*    5: CDSE-up, CDSE-down        */
	/* -up up distance                 */
	/* -down down distance down        */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        hts_countreads4refgene     \n");
		printf(" -d path of refgene database      \n");
		printf(" -dt database type                \n");
		printf("     0: UCSC refGene format       \n");
		printf("     1: UCSC refFlat format (default) \n");
		printf("     2: refLocus format              \n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -i input bar file              \n");
		printf(" -it input type                 \n");
		printf("	0: bar file (default) \n");
		printf("    1: a txt file listing all input bar files. \n");
		printf(" -o output text file \n");
		printf(" -r reference type \n");
		printf("    0: TSS-up, TES-down (default)    \n");
		printf("    1: TSS-up, TSS-down     \n");
		printf("    2: TES-up, TES-down     \n");
		printf("    3: CDSS-up, CDSE-down   \n");
		printf("    4: CDSS-up, CDSS-down   \n");
		printf("    5: CDSE-up, CDSE-down   \n");
		printf(" -up up distance limit (default = 0)\n");
		printf(" -down down distance limit (default = 0)\n");
		printf(" -normt normalize by total count (1: yes (default), 0: no) \n");
		printf(" -norml normalize by total gene length (1: yes (default), 0: no) \n");
		printf(" example: \n");
		printf("    hts_countreads4refgene -d refFlat_sorted.txt -dt 1 -s human -i K36.bar -o geneK36.txt -r 0 -up 5000 -down 1000 -normt 0\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	dtOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	nUP = 0;
	nDOWN = 0;
	itOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-it") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			itOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nRefType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUP = atoi(argc[ni]);
			if(nUP < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDOWN = atoi(argc[ni]);
			if(nDOWN < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-normt") == 0)
		{
			ni++;
			nStandardizebyT = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-norml") == 0)
		{
			ni++;
			nStandardizebyL = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_CountReads4RefGene_Main(strDatabasePath, nDatabaseType,
			strSpecies, strInputPath, nInputType, strOutputPath,
			nRefType, nUP, nDOWN, nStandardizebyT, nStandardizebyL);
	}

	/* return */
	return nResult;
}

int menu_cnv_aln2window(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nL = 10;
	int nN = 1;
	int nC = 1;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;

	int iOK;
	int oOK;
	int gOK;
	int lOK;
	int wOK;
	int sOK;
	int nOK;
	int cOK;
	int ni;

	/* ------------------------------- */
	/*    cnv_aln2window               */
	/* -i input                        */
	/* -o output                       */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -w window size (default = 100)  */
	/* -s step multiplier (default=10) */
	/* -n step number (default = 3)    */
	/* -c combine forward and reverse reads (default=yes) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    cnv_aln2window               \n");
		printf(" -i input                        \n");
		printf(" -o output                       \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step multiplier (default=10) \n");
		printf(" -n step number (default = 1)    \n");
		printf(" -c combine forward and reverse reads (default=1) or not (0) \n");
		printf(" example: \n");
		printf("    cnv_aln2window -i data.aln -o data_w -g chrlist.txt -l chrlen.txt -w 100 -s 10 -n 3\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	gOK = 0;
	lOK = 0;
	wOK = 0;
	sOK = 0;
	nOK = 0;
	cOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nL = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nN = atoi(argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nC = atoi(argc[ni]);
			cOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nC == 0)
		{
			nResult = CNV_Aln2Window_Main(strInputPath, strOutputPath, strChrList, strChrLen,
						nW, nL, nN);
		}
		else
		{
			nResult = CNV_Aln2WindowC_Main(strInputPath, strOutputPath, strChrList, strChrLen,
						nW, nL, nN);
		}
	}

	return nResult;
}

int menu_cnv_repeat2window(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nL = 10;
	int nN = 1;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;

	int iOK;
	int oOK;
	int gOK;
	int lOK;
	int wOK;
	int sOK;
	int nOK;
	int ni;

	/* ------------------------------- */
	/*    cnv_repeat2window            */
	/* -i genome sequence (*.sq) folder*/
	/* -o output                       */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -w window size (default = 100)  */
	/* -s step multiplier (default=10) */
	/* -n step number (default = 3)    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    cnv_repeat2window               \n");
		printf(" -i genome sequence (*.sq) folder\n");
		printf(" -o output                       \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step multiplier (default=10) \n");
		printf(" -n step number (default = 1)    \n");
		printf(" example: \n");
		printf("    cnv_repeat2window -i /data/hg18 -o hg18_repeat_w -g /data/hg18/chrlist.txt -l /data/hg18/chrlen.txt -w 100 -s 10 -n 3\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	gOK = 0;
	lOK = 0;
	wOK = 0;
	sOK = 0;
	nOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nL = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nN = atoi(argc[ni]);
			nOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = CNV_Repeat2Window_Main(strInputPath, strOutputPath, strChrList, strChrLen,
						nW, nL, nN);
	}

	return nResult;
}

int menu_rnaseq_countgeneread(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	int nInputType = 0;
	char strOutputPath[LINE_LENGTH];
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 1;
	char strSpecies[LINE_LENGTH];
	int nStandardize = 1;
	int nResult;

	int iOK;
	int tOK;
	int oOK;
	int dOK;
	int dtOK;
	int sOK;
	int zOK;
	int ni;

	/* ------------------------------- */
	/*     rnaseq_countgeneread        */
	/* -i input                        */
	/* -t 0(default): input is an alignment file  */
	/*    1: input is a file that contains a list of alignment files */
	/* -o output                       */
	/* -d gene annotatation database   */
	/* -dt annotation type. 0=refGene; 1=refFlat; 2=refLocus */
	/* -s species                      */
	/* -z 1(default): standardize by total read count and length, report # reads/1M reads/1kb; 0: not standardize, report raw counts */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    rnaseq_countgeneread               \n");
		printf(" -i input \n");
		printf(" -t 0(default): input is an alignment file; \n");
		printf("    1: input is a file that contains a list of alignment files. \n");
		printf(" -o output                       \n");
		printf(" -d gene annotatation database             \n");
		printf(" -dt annotation type. 0=refGene; 1=refFlat (default); 2=refLocus \n");
		printf(" -s species \n");
		printf(" -z 1(default): standardize by total read count and length, report # reads/1M reads/1kb; 0: not standardize, report raw counts. \n");
		printf(" example: \n");
		printf("    rnaseq_countgeneread -i input.aln -o output.txt -d /data/hg18/annotation/refFlat_sorted.txt -dt 1 -s human\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	tOK = 0;
	oOK = 0;
	dOK = 0;
	dtOK = 0;
	sOK = 0;
	zOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nStandardize = atoi(argc[ni]);
			zOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0) || (sOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RNASEQ_CountReadPerTranscript_Main(strInputPath, nInputType,
						strOutputPath, strDatabasePath, 
						nDatabaseType, strSpecies, nStandardize);
	}

	return nResult;
}

int menu_seqclust_seg(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nBinSize = 25;
	int nKernelType = 1;
	int nKernelStep = 1;
	int nKernelLen = 300;
	int nKernelBand = 100;
	int nSegType = 0;
	char strSegFile[MED_LINE_LENGTH];
	int nDistType = 0;
	int nUp = 0;
	int nDown = 0;
	int nDatabaseType = 1;
	int nMemBlock = 100000;
	int nCorrBlock = 10;
	int nGridNum = 10000;
	int nCutType = 1;
	double dCutL = 0.99;
	double dCutH = 0.999;
	int nExportBAR = 0;


	int iOK;
	int oOK;
	int dOK;
	int bOK;
	int ksOK;
	int siOK;
	int rOK;
	int upOK;
	int downOK;
	int dtOK;
	int ctOK;
	int clOK;
	int chOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust_seg           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust_seg                \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -b bin size (default = 25 bp) \n");
		printf(" -k smoothing kernel \n");
		printf("    0: no kernel, use raw read count \n");
		printf("    1 (default): one-sided exponential \n");
		printf("    2: two-sided exponential \n"); 
		printf("    3: one-sided gaussian \n");
		printf("    4: two-sided gaussian \n");
		printf(" -ks step size for kernel (a positive integer, default = bin_size) \n");
		printf(" -kw window size for kernel smoothing (a positive integer, default = 300, \n");
		printf("     which means smoothing is only carried out for positions within 300 base \n");
		printf("     pairs if one uses one-sided kernel, or within 300 bp each side if two-sided \n");
		printf("     kernel is used (i.e. 600 bp in total).\n");
		printf(" -kb kernel bandwidth (a positive integer, default = 100 (bp)) \n");
		printf(" -st approach for genome segmentation         \n");
		printf("     0 (default): automatic segmentation (no need to specify -si)      \n");
		printf("     1: user supplied coordinates (*.COD) file for genomic intervals (-si needs to be specified) \n");
		printf("     2: intervals defined based on gene structures (-si is needed; -r, -up, -down, -dt are optional) \n");
		printf(" -si if st=1, si specifies a COD file for genomic intervals \n");
		printf("     if st=2, si specifies a gene annotation file (refFlat_sorted.txt from CisGenome website) \n");
		printf(" -r if st=2, r specifies the reference points to define genomic intervals \n");
		printf("     0: TSS-up, TES-down (default)    \n");
		printf("     1: TSS-up, TSS-down     \n");
		printf("     2: TES-up, TES-down     \n");
		printf("     3: CDSS-up, CDSE-down   \n");
		printf("     4: CDSS-up, CDSS-down   \n");
		printf("     5: CDSE-up, CDSE-down   \n");
		printf(" -up if st=2, up specifies the distance 5' to the first reference point (default = 0) \n");
		printf(" -down if st=2, down specifies the distance 3' to the second reference point (default = 0) \n");
		printf("     Example: -st 2 -si refFlat_sorted.txt -r 0 -up 1000 -down 500 means get intervals \n");
		printf("	   that start from 1000 bp 5' upstream of TSS to 500 bp 3' downstream of TES of    \n");
		printf("       genes specified in the file refFlat_sorted.txt \n");
		printf(" -dt optional if st=2, annotation type: 0=refGene; 1=refFlat (default); 2=refLocus \n");
		printf(" -mb memory block size (default = 100000, which means each processing cycle will handle 100000 bins simultaneously \n");
		printf(" -cb correlation block size (default = 10, which means correlation coef of a bin is determined using its corr with 20 flanking bins (10 from left and 10 from right) \n");
		printf(" -ct cuttype. 0: cut by percentiles of bin statistics; 1 (default): cut by FDR; 2: cut by user supplied cutoffs. \n");
		printf(" -cl lower cutoff, default: 0.99 if cuttype = 0; 0.25 if cuttype = 1. This is the cutoff for defining peak boundaries. \n");
		printf(" -ch higher cutoff, default: 0.999 if cuttype = 0; 0.05 if cuttype = 1. This is the cutoff for initiating peaks. \n");
		printf(" -grid Number of grid points for computing FDR or percentile cutoffs (default = 10000). \n");
		printf(" -bar 1: export bin counts in bar file (may take significantly more time, but bar file can be used for visualization); 0: do not export bar file (default = 0). \n");
		printf(" Example: \n");
		printf("    seqclust_seg -i input_filelist.txt -d /workingfolder -o mycluster\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	bOK = 0;
	ksOK = 0;
	siOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	dtOK = 0;
	ctOK = 0;
	chOK = 0;
	clOK = 0;
	strcpy(strSegFile, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nKernelType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ks") == 0)
		{
			ni++;
			nKernelStep = atoi(argc[ni]);
			ksOK = 1;
		}
		else if(strcmp(argc[ni], "-kw") == 0)
		{
			ni++;
			nKernelLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kb") == 0)
		{
			ni++;
			nKernelBand = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-st") == 0)
		{
			ni++;
			nSegType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-si") == 0)
		{
			ni++;
			strcpy(strSegFile, argc[ni]);
			siOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nDistType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-mb") == 0)
		{
			ni++;
			nMemBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-cb") == 0)
		{
			ni++;
			nCorrBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ct") == 0)
		{
			ni++;
			nCutType = atoi(argc[ni]);
			ctOK = 1;
		}
		else if(strcmp(argc[ni], "-cl") == 0)
		{
			ni++;
			dCutL = atof(argc[ni]);
			clOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			dCutH = atof(argc[ni]);
			chOK = 1;
		}
		else if(strcmp(argc[ni], "-grid") == 0)
		{
			ni++;
			nGridNum = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bar") == 0)
		{
			ni++;
			nExportBAR = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(ksOK == 0)
	{
		nKernelStep = nBinSize;
		if(nKernelStep == 0)
			nKernelStep = 1;
	}
	if(nSegType != 0)
	{
		if(siOK == 0)
		{
			printf("Error: please use -si to specify a genomic interval (if -st 1) or gene annotation file (if -st 2)!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	if(nCutType == 0)
	{
		if(clOK == 0)
			dCutL = 0.99;
		if(chOK == 0)
			dCutH = 0.999;
	}
	else if(nCutType == 1)
	{
		if(clOK == 0)
			dCutL = 0.25;
		if(chOK == 0)
			dCutH = 0.05;
	}
	else
	{
		if((clOK == 0) || (chOK == 0))
		{
			printf("Error: please specify segmentation cutoffs using -cl and/or -ch options!\n");
			exit(EXIT_FAILURE);
		}
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_Seg_Main(strInputPath, strOutputPath, strOutputFile, nBinSize,
			nKernelType, nKernelStep, nKernelLen, nKernelBand, 
			nSegType, strSegFile, nDistType, nUp, nDown, nDatabaseType, 
			nMemBlock, nCorrBlock, nGridNum, nCutType, dCutL, dCutH, nExportBAR);
	}

	return nResult;
}

int menu_seqclust_dp(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nBinSize = 25;
	int nKernelType = 1;
	int nKernelStep = 1;
	int nKernelLen = 300;
	int nKernelBand = 100;
	int nSegType = 0;
	char strSegFile[MED_LINE_LENGTH];
	int nDistType = 0;
	int nUp = 0;
	int nDown = 0;
	int nDatabaseType = 1;
	int nMemBlock = 100000;
	int nCorrBlock = 10;
	int nGridNum = 10000;
	int nCutType = 0;
	double dCutL = 0.99;
	double dCutH = 0.999;
	int nBlockLenCut = 100;
	int nExportBAR = 0;


	int iOK;
	int oOK;
	int dOK;
	int bOK;
	int ksOK;
	int siOK;
	int rOK;
	int upOK;
	int downOK;
	int dtOK;
	int ctOK;
	int clOK;
	int chOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust_dp            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust_dp                \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -b bin size (default = 25 bp) \n");
		printf(" -k smoothing kernel \n");
		printf("    0: no kernel, use raw read count \n");
		printf("    1 (default): one-sided exponential \n");
		printf("    2: two-sided exponential \n"); 
		printf("    3: one-sided gaussian \n");
		printf("    4: two-sided gaussian \n");
		printf(" -ks step size for kernel (a positive integer, default = bin_size) \n");
		printf(" -kw window size for kernel smoothing (a positive integer, default = 300, \n");
		printf("     which means smoothing is only carried out for positions within 300 base \n");
		printf("     pairs if one uses one-sided kernel, or within 300 bp each side if two-sided \n");
		printf("     kernel is used (i.e. 600 bp in total).\n");
		printf(" -kb kernel bandwidth (a positive integer, default = 100 (bp)) \n");
		printf(" -st approach for genome segmentation         \n");
		printf("     0 (default): automatic segmentation (no need to specify -si)      \n");
		printf("     1: user supplied coordinates (*.COD) file for genomic intervals (-si needs to be specified) \n");
		printf("     2: intervals defined based on gene structures (-si is needed; -r, -up, -down, -dt are optional) \n");
		printf(" -si if st=1, si specifies a COD file for genomic intervals \n");
		printf("     if st=2, si specifies a gene annotation file (refFlat_sorted.txt from CisGenome website) \n");
		printf(" -r if st=2, r specifies the reference points to define genomic intervals \n");
		printf("     0: TSS-up, TES-down (default)    \n");
		printf("     1: TSS-up, TSS-down     \n");
		printf("     2: TES-up, TES-down     \n");
		printf("     3: CDSS-up, CDSE-down   \n");
		printf("     4: CDSS-up, CDSS-down   \n");
		printf("     5: CDSE-up, CDSE-down   \n");
		printf(" -up if st=2, up specifies the distance 5' to the first reference point (default = 0) \n");
		printf(" -down if st=2, down specifies the distance 3' to the second reference point (default = 0) \n");
		printf("     Example: -st 2 -si refFlat_sorted.txt -r 0 -up 1000 -down 500 means get intervals \n");
		printf("	   that start from 1000 bp 5' upstream of TSS to 500 bp 3' downstream of TES of    \n");
		printf("       genes specified in the file refFlat_sorted.txt \n");
		printf(" -dt optional if st=2, annotation type: 0=refGene; 1=refFlat (default); 2=refLocus \n");
		printf(" -mb memory block size (default = 100000, which means each processing cycle will handle 100000 bins simultaneously \n");
		printf(" -cb correlation block size (default = 10, which means correlation coef of a bin is determined using its corr with 20 flanking bins (10 from left and 10 from right) \n");
		printf(" -ct cuttype. 0 (default): cut by percentiles of bin statistics; 1: cut by FDR; 2: cut by user supplied cutoffs. \n");
		printf(" -cl lower cutoff, default: 0.99 if cuttype = 0; 0.25 if cuttype = 1. This is the cutoff for defining peak boundaries. \n");
		printf(" -ch higher cutoff, default: 0.999 if cuttype = 0; 0.05 if cuttype = 1. This is the cutoff for initiating peaks. \n");
		printf(" -grid Number of grid points for computing FDR or percentile cutoffs (default = 10000). \n");
		printf(" -bar 1: export bin counts in bar file (may take significantly more time, but bar file can be used for visualization); 0: do not export bar file (default = 0). \n");
		printf(" -minlen minimal region size in base pairs, regions shorter than this size will not be reported\n");
		printf(" Example: \n");
		printf("    seqclust_dp -i input_filelist.txt -d /workingfolder -o mycluster\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	bOK = 0;
	ksOK = 0;
	siOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	dtOK = 0;
	ctOK = 0;
	chOK = 0;
	clOK = 0;
	strcpy(strSegFile, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nKernelType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ks") == 0)
		{
			ni++;
			nKernelStep = atoi(argc[ni]);
			ksOK = 1;
		}
		else if(strcmp(argc[ni], "-kw") == 0)
		{
			ni++;
			nKernelLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kb") == 0)
		{
			ni++;
			nKernelBand = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-st") == 0)
		{
			ni++;
			nSegType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-si") == 0)
		{
			ni++;
			strcpy(strSegFile, argc[ni]);
			siOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nDistType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-mb") == 0)
		{
			ni++;
			nMemBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-cb") == 0)
		{
			ni++;
			nCorrBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ct") == 0)
		{
			ni++;
			nCutType = atoi(argc[ni]);
			ctOK = 1;
		}
		else if(strcmp(argc[ni], "-cl") == 0)
		{
			ni++;
			dCutL = atof(argc[ni]);
			clOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			dCutH = atof(argc[ni]);
			chOK = 1;
		}
		else if(strcmp(argc[ni], "-grid") == 0)
		{
			ni++;
			nGridNum = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bar") == 0)
		{
			ni++;
			nExportBAR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-minlen") == 0)
		{
			ni++;
			nBlockLenCut = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(ksOK == 0)
	{
		nKernelStep = nBinSize;
		if(nKernelStep == 0)
			nKernelStep = 1;
	}
	if(nSegType != 0)
	{
		if(siOK == 0)
		{
			printf("Error: please use -si to specify a genomic interval (if -st 1) or gene annotation file (if -st 2)!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	if(nCutType == 0)
	{
		if(clOK == 0)
			dCutL = 0.99;
		if(chOK == 0)
			dCutH = 0.999;
	}
	else if(nCutType == 1)
	{
		if(clOK == 0)
			dCutL = 0.25;
		if(chOK == 0)
			dCutH = 0.05;
	}
	else
	{
		if((clOK == 0) || (chOK == 0))
		{
			printf("Error: please specify segmentation cutoffs using -cl and/or -ch options!\n");
			exit(EXIT_FAILURE);
		}
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_SegDP_Main(strInputPath, strOutputPath, strOutputFile, nBinSize,
			nKernelType, nKernelStep, nKernelLen, nKernelBand, 
			nSegType, strSegFile, nDistType, nUp, nDown, nDatabaseType, 
			nMemBlock, nCorrBlock, nGridNum, nCutType, dCutL, dCutH, nBlockLenCut, nExportBAR);
	}

	return nResult;
}

int menu_seqpeak(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nPaired = 0;
	int nBinSize = 50;
	int nExtLen = 150;
	int nSegType = 0;
	int nWinSize = 1;
	int nCutType = 0;
	double dCutoff = 3.0;
	int nMaxGap = 50;
	int nMinLen = 100;
	int nExportBAR = 0;
	int nKeepTempFiles = 0;
	int nBoundaryRefine = 1;
	int nBRWin = 5;
	int nCollectRawData = 0;
	int nPoisFilter = 1;
	int nPoisWin = 10000;
	double dPoisCut = 1.0e-5;
	int nTStandardize = 1;
	int nLFCAdj = 0;
	int nLFCWin = 10000;

	int iOK;
	int oOK;
	int dOK;

	int bOK;
	int wOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqpeak                */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqpeak                \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -b bin size, the resolution of peak detection (default = 50 bp) \n");
		printf(" -w half window size, 2*w+1 bins will be used to determine enrichment of the center bin (default = 1) \n");
		printf(" -e read extension length  (default = 150 bp) \n");
		/* printf(" -st approach for genome segmentation         \n");
		printf("     0 (default): automatic segmentation      \n");
		printf(" -ct cuttype.(default = 0) \n");
		*/
		printf(" -ts 1 (default): standardize window statistics. 0: do not standardize. \n");
		printf(" -c cutoff for defining peak boundaries. (default = 3.0) \n");
		printf(" -maxgap (default = 50 bp): two peaks will be merged if the gap (in bp) between them is smaller than this value \n");
		printf(" -minlen minimal region size (default = 100 bp). Regions shorter than this size will not be reported\n");
		printf(" -br 1 (default): refine peak boundaries using +/- strand peak offsets; 0: no boundary refinement \n");
		printf(" -bw boundary refinement resolution (in bp). (default = 5). The smaller the number (i.e. the higher the resolution), the more memory is required. \n");
		printf(" -bar 1: export bin counts in bar file (may take significantly more time, but bar file can be used for visualization); 0 (default): do not export bar file. \n");
		printf(" -dat 1: report normalized read counts for each peak and each sample; 0 (default): not report the normalized raw data. \n");
		printf(" -keeptemp 1: keep intermediate files (may take big disk space); 0 (default): remove intermediate files. \n");
		printf(" -lpois 1 (default): apply local poisson rate filter (need to be careful when only IP sample is available and the signal is long); 0: no local poisson rate filter. \n");
		printf(" -lpwin window size (in bp) used to estimate local poisson rate (default = 10000). \n");
		printf(" -lpcut p-value cutoff for local poisson rate correction (default = 1e-5). \n");
		printf(" -lfc 1: local fold change adjustment (usually do not use it as it may reduce sensitivity); 0 (default): no adjustment. \n");
		printf(" -lfcwin window size (in bp) used to estimate average local coverage for lfc adjustment (default = 10000). \n");
		printf(" Example: \n");
		printf("    seqpeak -i input_filelist.txt -d . -o mypeak -b 25 -w 5 -e 100\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	bOK = 0;
	wOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nWinSize = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-st") == 0)
		{
			ni++;
			nSegType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ct") == 0)
		{
			ni++;
			nCutType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ts") == 0)
		{
			ni++;
			nTStandardize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dCutoff = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-br") == 0)
		{
			ni++;
			nBoundaryRefine = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bw") == 0)
		{
			ni++;
			nBRWin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bar") == 0)
		{
			ni++;
			nExportBAR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-keeptemp") == 0)
		{
			ni++;
			nKeepTempFiles = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-maxgap") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-minlen") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-dat") == 0)
		{
			ni++;
			nCollectRawData = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpois") == 0)
		{
			ni++;
			nPoisFilter = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpwin") == 0)
		{
			ni++;
			nPoisWin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpcut") == 0)
		{
			ni++;
			dPoisCut = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lfc") == 0)
		{
			ni++;
			nLFCAdj = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lfcwin") == 0)
		{
			ni++;
			nLFCWin = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(nBinSize < 0)
		nBinSize = 50;

	if(nWinSize < 0)
		nWinSize = 0;
	
	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqPeak_Main(strInputPath, strOutputPath, strOutputFile, 
			nPaired, nBinSize, nExtLen, 
			nSegType, nWinSize, nCutType, 
			nTStandardize, dCutoff, nMaxGap, nMinLen, 
			nBoundaryRefine, nBRWin,
			nExportBAR, nKeepTempFiles, nCollectRawData, 
			nPoisFilter, nPoisWin, dPoisCut, 
			nLFCAdj, nLFCWin);
	}

	return nResult;
}

int menu_seqclust_count(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strDataPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nExtendLen = 0;
	
	int iOK;
	int dOK;
	int oOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust_count         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust_count                \n");
		printf(" -i input coordinate file \n");
		printf(" -d sample file list  \n");
		printf(" -o output file  \n");
		printf(" -e extension length (default = 0 bp), typically equal to DNA fragment length.\n");
		printf(" Example: \n");
		printf("    seqclust_count -i input_file.cod -d sample_filelist.txt -o mycountdata.txt -e 100\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtendLen = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_Count_Main(strInputPath, strDataPath, strOutputFile, nExtendLen);
	}

	return nResult;
}

int menu_seqclust(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nKmin = 2;
	int nKmax = 20;
	int nBmin = 1;
	int nBmax = 1;
	int nKr = 5;
	int nMethod = 0;
	int nTransform = 0;
	double dTL = 1.0;
	int nRowStandardize = 0;
	double dCut = 0.0;
	int nSeed = 13;
	int nSkipCol = 0;
	int nMaxIter = 1000;
	double dTol = 1e-4;

	int iOK;
	int oOK;
	int dOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust               */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust                    \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -kmin minimal cluster number, default = 2 \n");
		printf(" -kmax maximal cluster number, default = 20 \n");
		printf(" -kr number of trials (using different seeds) for each cluster number, default = 5 \n");
		printf(" -m method for clustering (0: (default) k normals; 1: k dirichlet; \n");
		printf("     2: k multivariate normals (MVN) with user specified group ID. The IDs should be given in \n");
		printf("     the data file as a line starting with GROUPID, samples within a group are modeled as MVN \n");
		printf("     samples from different groups are assumed to be independent. \n");
		printf("     3: k multivariate normals (MVN) with automatically determined group ID.\n");
		printf("     4: k multivariate normals with local factor dimension reduction \n"); 
		printf(" -bmin minimal block/factor level. Required if -m 3 or -m -4 is used. \n");
		printf("     if -m 3, then Block level = 1 means all samples are assumed to be independent. \n");
		printf("     if -m 4, then bmin should be > 0; it is the minimal number of factors used. \n");
		printf("           (note samples are independent only when number of factors = 0. \n");
		printf("     Block level = sample size means all samples are modeled jointly by a MVN. \n");
		printf(" -bmax maximal block level. \n");
		printf(" -t how to transform the data before clustering \n");
		printf("    0: no transform (default)       \n");
		printf("    1: log2 transform               \n");
		printf(" -tl truncation lower bound (default = 1.0), used to avoid log(0) when -t 1 is set \n");
		printf(" -srow 1: standardize rows (x-mean)/sd; 0: no standardization (default) \n");
		printf(" -c cutoff for reporting the clusters, default = 0.0 (i.e. no cut, classify genes to the most likely cluster) \n");
		printf(" -seed seed for generating random numbers, default = 13. \n");
		printf(" -skipcol how many columns to skip before the data start, default = 0 \n");
		printf(" -maxiter maximal iteration for each run, default = 1000 \n");
		printf(" -tol convergence criteria, default tol = 1e-4\n");
		printf(" Example: \n");
		printf("    seqclust -i inputfile.txt -d /workingfolder -o mycluster\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-kmin") == 0)
		{
			ni++;
			nKmin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kmax") == 0)
		{
			ni++;
			nKmax = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kr") == 0)
		{
			ni++;
			nKr = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nMethod = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bmin") == 0)
		{
			ni++;
			nBmin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bmax") == 0)
		{
			ni++;
			nBmax = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-tl") == 0)
		{
			ni++;
			dTL = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-srow") == 0)
		{
			ni++;
			nRowStandardize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dCut = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-seed") == 0)
		{
			ni++;
			nSeed = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-skipcol") == 0)
		{
			ni++;
			nSkipCol = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-maxiter") == 0)
		{
			ni++;
			nMaxIter = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-tol") == 0)
		{
			ni++;
			dTol = atof(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_Main(strInputPath, strOutputPath, strOutputFile, nSkipCol,
			nKmin, nKmax, nKr, nMethod, nBmin, nBmax, nSeed,
			nTransform, dTL, nRowStandardize, dCut, 
			nMaxIter, dTol);
	}

	return nResult;
}

int menu_flexmodule(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    flexmodule                   */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    flexmodule                    \n");
		printf(" example: \n");
		printf("    flexmodule flexmodule_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = FlexModule_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_flexmodule_bgfit(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    flexmodule                   */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    flexmodule_bgfit               \n");
		printf(" example: \n");
		printf("    flexmodule_bgfit flexmodule_bgfit_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = FlexModule_BGFit_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_flexmodule_extractmotif(int argv, char **argc)
{
	/* define */
	char strInPath[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	int iOK,oOK;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*   flexmodule_extractmotif       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    flexmodule_extractmotif                \n");
		printf(" example: \n");
		printf("    flexmodule_extractmotif  -i motif_l.txt -o newmotif_l\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FlexModule_ExtractMotifsFromResult(strInPath, strOutPath);
	}
	
	return nResult;
}


int menu_malign_blasthit(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_blasthit              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_blasthit                \n");
		printf(" example: \n");
		printf("    malign_blasthit malign_blasthit_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_BlastHit_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_malign_motifmap(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_motifmap              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_motifmap                \n");
		printf(" example: \n");
		printf("    malign_motifmap malign_motifmap_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_MotifMap_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_malign_modulemap(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_modulemap              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_modulemap               \n");
		printf(" example: \n");
		printf("    malign_modulemap malign_modulemap_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_ModuleMap_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_malign_generatemotifaln(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	/* ------------------------------- */
	/*    malign_modulemap              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_generatemotifaln        \n");
		printf(" example: \n");
		printf("    malign_generatemotifaln malign_motifaln_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_GenMotifAln_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_malign_countkmer(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_modulemap              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_countkmer               \n");
		printf(" example: \n");
		printf("    malign_countkmer malign_kmerstat_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_KMerStat_Main(strParamPath);
	
	/* return */
	return nResult;
}


int menu_malign_genome_prepareortholog(int argv, char **argc)
{
	/* define */
	char strInPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	int nSkip = 0;
	int nRefType = 0;
	int nUp = 5000;
	int nDown = 5000;
	int nSpeciesNum = 0;

	int iOK,oOK,nOK,sfOK,rOK,upOK,downOK;
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*  malign_genome_prepareortholog  */
	/* ------------------------------- */
	if((argv < 1))
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   malign_genome_prepareortholog   \n");
		printf("    -i input nstmap or nsomap file \n");
		printf("    -o output file                 \n");
		printf("    -n species number              \n");
		printf("    -sf skip the first species (0: No, 1: Yes) \n");
		printf("    -r reference type              \n");
		printf("       0: TSS-up, TES-down         \n");
		printf("       1: TSS-up, TSS-down         \n");
		printf("       2: TES-up, TES-down         \n");
		printf("       3: CDSS-up, CDSE-down       \n");
		printf("       4: CDSS-up, CDSS-down       \n");
		printf("       5: CDSE-up, CDSE-down       \n");
		printf("    -up up distance                \n");
		printf("    -down down distance down       \n");
		printf(" example: \n");
		printf("   malign_genome_prepareortholog -i input.txt -o output.txt -n 4 -sf 1 -r 0 -up 10000 -down 10000 \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	nOK = 0;
	sfOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;


	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nSpeciesNum = atoi(argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-sf") == 0)
		{
			ni++;
			nSkip = atoi(argc[ni]);
			sfOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nRefType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			if(nUp < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			if(nDown < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0) || (nOK == 0) || (upOK == 0) || (downOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MAlign_Genome_PrepareOrtholog (strInPath, strOutPath,
			nSpeciesNum, nSkip, nRefType, nUp, nDown);
	}
	
	/* return */
	return nResult;
}

int menu_malign_genome_blasthit(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_genome_blasthit       */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_genome_blasthit         \n");
		printf(" example: \n");
		printf("    malign_genome_blasthit malign_genome_blasthit_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_Genome_BlastHit_Main(strParamPath);
	
	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_genomebackground(int argv, char **argc)
{
	/* define */
	int dOK,oOK,bOK,sOK,wOK;
	char strGenomePath[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	int nBGOrder;
	int nS;
	int nW;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*  motifmap_matrixscan_genomebg   */
	/* -d genome database path         */
	/* -o exporting directory          */
	/* -b background markov order      */
	/* -s step size                    */
	/* -w window size                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genomebg    \n");
		printf(" -d genome database path \n");
		printf(" -o exporting directory \n");
		printf(" -b background markov order  \n");
		printf(" -s step size (should <= 1000000)\n");
		printf(" -w window size\n");
		printf(" (note: w%s must be 0, i.e. w=a*s, where a is an integer)");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genomebg -d /data/mm6 -o /data/mm6/mcbg/3/ -b 3 -s 100000 -w 1000000\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	dOK = 0;
	oOK = 0;
	bOK = 0;
	sOK = 0;
	wOK = 0;
	
	nBGOrder = 0;
	nS = 100000;
	nW = 1000000;
	
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (oOK == 0) || (bOK == 0) || (sOK == 0) || (wOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_ScanMatrix_GenomeBackground_Main(strGenomePath, 
					strOutPath, nBGOrder, nS, nW);
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_genome(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,oOK,rOK,bOK,btOK,bdOK,bsOK,cOK,cdOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	double dR;
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;

	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_matrixscan_genome    */
	/* -m motif matrix                 */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat (0: no (default); 1: using unmasked genome) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genome    \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file    \n");
		printf(" -o output file (full path) \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genome -m Gli_matrix.txt -gd /data/mm6 -i inputseq.cod -o Gli_map.txt -r 100 -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 100 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dR = 100.0;
	dC = 0.0;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if( (nIncludeRepeat != 0) && (nIncludeRepeat != 1) )
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
					strCodPath,	strOutputPath, dR, 
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
					strCodPath,	strOutputPath, dR, 
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}


int menu_motifmap_matrixintegrate_genome(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,sampOK,oOK,sOK,wOK,rOK,bOK,btOK,bdOK,bsOK,cOK,cdOK,repOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strSampPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	int nW;
	double dR;
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	int nRepeat;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;
	int nResult = 0;
	
	/* ------------------------------- */
	/* motifmap_matrixintegrate_genome */
	/* -m motif matrix                 */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -samp sampling file             */
	/* -o output file                  */
	/* -s species                      */
	/* -w window size                  */
	/* -repeat include repeat          */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixintegrate_genome   \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file    \n");
		printf(" -samp sampling file    \n");
		printf(" -o output file (full path) \n");
		printf(" -s species             \n");
		printf(" -w window size         \n");
		printf(" -repeat include repeat \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" example: \n");
		printf("    motifmap_matrixintegrate_genome -m Gli_matrix.txt -gd /data/mm6 -i inputseq.cod -samp affytile.txt -o Gli_map.txt -s mouse -w 200 -repeat 1 -r -1.0 -b 3 -bt genome -bd /data/mm6/MC/3 -bs 100000 -c -1 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	sampOK = 0;
	sOK = 0;
	wOK = 0;
	repOK = 0;
	
	dR = 100.0;
	dC = 0.0;
	nW = 100;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	nRepeat = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");
	strcpy(strSpecies, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-samp") == 0)
		{
			ni++;
			strcpy(strSampPath, argc[ni]);
			sampOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-repeat") == 0)
		{
			ni++;
			nRepeat = atoi(argc[ni]);
			repOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0)
		|| (sampOK == 0) || (sOK == 0) || (wOK == 0) || (nW < 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Integrate_Genome_Main(strMotifPath, 
					strGenomePath, strCodPath, 
					strSampPath, strOutputPath, 
					strSpecies, nW, nRepeat, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					cOK, dC, strCSPath);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Integrate_Genome_Main(strMotifPath, 
					strGenomePath, strCodPath, 
					strSampPath, strOutputPath, 
					strSpecies, nW, nRepeat, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					cOK, dC, strCSPath);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_genome_enrich(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,nOK,sOK,oOK,rOK,bOK,btOK,bdOK,bsOK,cOK,cdOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strNegCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nTierSize;
	double dR;
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*motifmap_matrixscan_genome_enrich*/
	/* -m motif matrix                 */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -n negative control coordinates */
	/* -s tier size                    */
	/* -o output file                  */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genome_enrich    \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinate file    \n");
		printf(" -n negative control coordinate file \n");
		printf(" -s tier size \n");
		printf(" -o output file (full path) \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genome_enrich -m Gli_matrix.txt -gd /data/mm6 -i inputseq.cod -n control.cod -s 20 -o Gli_tierenrich.txt -r 100 -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 100 -cd /data/mm6/conservation/phastcons/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	nOK = 0;
	sOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	
	nTierSize = 1;
	dR = 100.0;
	dC = 0.0;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strCSPath, "");
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegCodPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nTierSize = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (nOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Genome_Enrich_Main(strMotifPath, strGenomePath, 
					strCodPath,	strNegCodPath, nTierSize, strOutputPath, dR, 
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Genome_Enrich_Main(strMotifPath, strGenomePath, 
					strCodPath,	strNegCodPath, nTierSize, strOutputPath, dR, 
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_consensusscan_genome_enrich(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,nOK,sOK,oOK,mcOK,mdOK,cOK,cdOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strNegCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nTierSize;
	int nMC,nMD;
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*motifmap_consensusscan_genome_enrich*/
	/* -m motif matrix                 */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -n negative control coordinates */
	/* -s tier size                    */
	/* -o output file                  */
	/* -mc mismatch to consensus       */
	/* -md mismatch to degenerate consensus */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_consensusscan_genome_enrich    \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinate file    \n");
		printf(" -n negative control coordinate file \n");
		printf(" -s tier size \n");
		printf(" -o output file (full path) \n");
		printf(" -mc max. mismatch to consensus \n");
		printf(" -md max. mismatch to degenerate consensus  \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_consensusscan_genome_enrich -m Gli_matrix.txt -gd /data/mm6 -i inputseq.cod -n control.cod -s 20 -o Gli_tierenrich.txt -mc 2 -md 1 -c 100 -cd /data/mm6/conservation/phastcons\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	nOK = 0;
	sOK = 0;
	oOK = 0;
	mcOK = 0;
	mdOK = 0;
	cOK = 0;
	cdOK = 0;
	
	nTierSize = 1;
	nMC = 0;
	nMD = 0;
	dC = 0.0;
	nUseCS = 0;
	strcpy(strCSPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegCodPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nTierSize = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-mc") == 0)
		{
			ni++;
			nMC = atoi(argc[ni]);
			mcOK = 1;
		}
		else if(strcmp(argc[ni], "-md") == 0)
		{
			ni++;
			nMD = atoi(argc[ni]);
			mdOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (nOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanConsensus_Genome_Enrich_Main(strMotifPath, strGenomePath, 
					strCodPath,	strNegCodPath, nTierSize, strOutputPath, nMC, nMD, 
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Genome_Enrich_Main(strMotifPath, strGenomePath, 
				strCodPath,	strNegCodPath, nTierSize, strOutputPath, nMC, nMD, 
				cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_genome_getcutoff(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,oOK,bOK,btOK,bdOK,bsOK,cOK,cdOK;
	char strMotifListPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*motifmap_matrixscan_genome_getcutoff*/
	/* -mr motif matrix and required   */
	/*     occurence rate              */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -bs background mc step size     */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genome_getcutoff  \n");
		printf(" -mr motif matrix and required occurence rate (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinate file    \n");
		printf(" -o output file (full path) \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genome_getcutoff -mr motiflist.txt -gd /data/mm6 -i inputseq.cod -o motiflist_cutoff.txt -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 100 -cd /data/mm6/conservation/phastcons\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = 0.0;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-mr") == 0)
		{
			ni++;
			strcpy(strMotifListPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Genome_GetCutoff_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strOutputPath,
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Genome_GetCutoff_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strOutputPath,
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_genome_summary(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,nOK,oOK,bOK,btOK,bdOK,bsOK,cOK,cdOK;
	char strMotifListPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strNegCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*motifmap_matrixscan_genome_summary*/
	/* -mr motif matrix and likelihood */
	/*     ratio list                  */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -n negative control coordinates */
	/* -o output file                  */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -bs background mc step size     */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genome_summary    \n");
		printf(" -mr motif matrix and likelihood ratio list (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinate file    \n");
		printf(" -n negative control coordinate file \n");
		printf(" -o output file (full path) \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genome_summary -mr motiflist.txt -gd /data/mm6 -i inputseq.cod -n control.cod -o motifs_enrich.txt -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 100 -cd /data/mm6/conservation/phastcons/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	nOK = 0;
	oOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = 0.0;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strCSPath, "");
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-mr") == 0)
		{
			ni++;
			strcpy(strMotifListPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegCodPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (nOK == 0) || (oOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Genome_Summary_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strNegCodPath, strOutputPath,
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Genome_Summary_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strNegCodPath, strOutputPath,
					nBGOrder, strBGType, strBGPath, nBGStepSize,
					cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_consensusscan_genome_summary(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,nOK,oOK,cOK,cdOK;
	char strMotifListPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strNegCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*motifmap_consensusscan_genome_summary*/
	/* -mr motif consensus list and    */
	/*     matching criteria           */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -n negative control coordinates */
	/* -o output file                  */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_consensusscan_genome_summary    \n");
		printf(" -mr motif consensus list and matching criteria (full path) \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinate file    \n");
		printf(" -n negative control coordinate file \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_consensusscan_genome_summary -mr motiflist.txt -gd /data/mm6 -i inputseq.cod -n control.cod -o motifs_enrich.txt -c 100 -cd /data/mm6/conservation/phastcons/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	nOK = 0;
	oOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = 0.0;
	nUseCS = 0;
	strcpy(strCSPath, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-mr") == 0)
		{
			ni++;
			strcpy(strMotifListPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegCodPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (nOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanConsensus_Genome_Summary_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strNegCodPath, strOutputPath,
					cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Genome_Summary_Main(strMotifListPath, strGenomePath, 
					strCodPath,	strNegCodPath, strOutputPath,
					cOK, dC, strCSPath, nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_multiplematrixscan_genome(int argv, char **argc)
{
	/* define */
	int sOK,gdOK,dOK,iOK,oOK,bOK,btOK,bdOK,bsOK,cOK,cdOK,wOK;
	char strSpecies[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	char strCSPath[MED_LINE_LENGTH];
	int nW;
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* --------------------------------------- */
	/*   motifmap_multiplematrixscan_genome    */
	/* -s species                              */
	/* -gd genome sequence path                */
	/* -d working directory                    */
	/* -i data files                           */
	/* -o output file                          */
	/* -b background markov order              */
	/* -bt type of background                  */
	/* -bd background mc path                  */
	/* -bs the stepsize used for computing MC  */
	/* -c use conservation score?              */
	/* -cd conservation score data path        */
	/* -w window half width for testing        */
	/*    clustering                           */
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* --------------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ---------------------------------- */\n");
		printf("   motifmap_multiplematrixscan_genome    \n");
		printf(" -s species\n");
		printf(" -gd genome sequence path\n");
		printf(" -d working directory \n");
		printf(" -i data (coordinates&motif) file    \n");
		printf(" -o output file (full path) \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c use conservation? (0: no; 1: yes) \n");
		printf(" -cd conservation score path \n");
		printf(" -w window half width for testing clustering \n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_multiplematrixscan_genome -s human -gd /data/genomes/human/b35_hg17/ -d . -i ChIPmotif.txt -o ChIPmap -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 1 -cd /data/mm6/conservation/phastcons -w 100\n");
		printf("/* ---------------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	sOK = 0;	
	gdOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	wOK = 0;

	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");
	strcpy(strCSPath, "");
	nW = 100;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strWorkPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nUseCS = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}

	if((sOK == 0) || (gdOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (bOK == 0) || (cOK == 0) || (wOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nUseCS == 1)
		{
			if(cdOK == 0)
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}

		nResult = MotifMap_ScanMultipleMatrix_Genome_Main(strSpecies, strGenomePath,  
				strWorkPath, strInputPath, strOutputPath,  
				nBGOrder, strBGType, strBGPath, nBGStepSize,
				nUseCS, strCSPath, nIncludeRepeat, nW);
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan_group(int argv, char **argc)
{
	/* define */
	int dOK,mOK,iOK,oOK,rOK,bOK,cOK,chOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	double dC;
	double dR;
	int nBGOrder;
	int ni;
	int nUseCS;

	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_matrixscan_group     */
	/* -m motif matrix                 */
	/* -d sequence & conservation path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -c conservation cutoff          */
	/* -ch head tag for *.cs file      */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_matrixscan_group      \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file    \n");
		printf(" -o output file (full path) \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -c conservation cutoff \n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_matrixscan_group -m E2F.txt -d . -i inputseq.fa -o outputseq.txt -r 100 -b 3 -c 100 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	cOK = 0;
	chOK = 0;
	
	dC = 0.0;
	dR = 100.0;
	nBGOrder = 0;
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			strcpy(strCSAlias, argc[ni]);
			chOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			nUseCS = 1;
			if(chOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Group_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_ScanMatrix_Group_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, "");
			}
		}
		else
		{
			nUseCS = 0;
			nResult = MotifMap_ScanMatrix_Group_Main(strMotifPath, strInputPath, strSeqFile,
				strOutputPath, dR, nBGOrder, nUseCS, 0.0, "");
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_matrixscan(int argv, char **argc)
{
	/* define */
	int dOK,mOK,iOK,oOK,rOK,bOK,cOK,chOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	double dC;
	double dR;
	int nBGOrder;
	int ni;
	int nUseCS;

	int nResult;
	
	/* ------------------------------- */
	/*    motifmap_matrixscan          */
	/* -m motif matrix                 */
	/* -d sequence & conservation path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -c conservation cutoff          */
	/* -ch head tag for *.cs file      */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_matrixscan            \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file    \n");
		printf(" -o output file (full path) \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -c conservation cutoff \n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_matrixscan -m E2F.txt -d . -i inputseq.fa -o outputseq.txt -r 100 -b 3 -c 100 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	cOK = 0;
	chOK = 0;
	
	dC = 0.0;
	dR = 100.0;
	nBGOrder = 0;
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			strcpy(strCSAlias, argc[ni]);
			chOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			nUseCS = 1;
			if(chOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, "");
			}
		}
		else
		{
			nUseCS = 0;
			nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,
				strOutputPath, dR, nBGOrder, nUseCS, 0.0, "");
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_consensusscan(int argv, char **argc)
{
	/* define */
	int dOK,mOK,iOK,oOK,mcOK,mdOK,cOK,chOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	int nMC,nMD;
	double dC;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    motifmap_consensusscan       */
	/* -m motif consensus              */
	/* -d sequence & conservation score data path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -mc number of stringent         */
	/*    consensus mismatches allowed */
	/* -md number of relaxed mismatches*/
	/*    allowed                      */
	/* -c conservation cutoff\n")      */
	/* -ch head tag for conservation score *.cs file */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_consensusscan         \n");
		printf(" -m motif consensus (full path)    \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file    \n");
		printf(" -o output file (full path) \n");
		printf(" -mc number of stringent consensus mismatches allowed\n");
		printf(" -md number of relaxed consensus mismatches allowed\n");
		printf(" -c conservation cutoff\n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_consensusscan -m E2F.txt -d . -i inputseq.fa -o outputseq.txt -mc 2 -md 0 -c 200 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	mcOK = 0;
	mdOK = 0;
	cOK = 0;
	chOK = 0;
	nMC = 0;
	nMD = 0;
	dC = 0.0;
	strcpy(strCSAlias, "");
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-mc") == 0)
		{
			ni++;
			nMC = atoi(argc[ni]);
			mcOK = 1;
		}
		else if(strcmp(argc[ni], "-md") == 0)
		{
			ni++;
			nMD = atoi(argc[ni]);
			mdOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			strcpy(strCSAlias, argc[ni]);
			chOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (mcOK == 0) || (mdOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(chOK == 1)
			{
				nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
						strOutputPath, nMC, nMD, cOK, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
						strOutputPath, nMC, nMD, cOK, dC, "");
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,
					strOutputPath, nMC, nMD, cOK, 0.0, "");
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_consensusscan_genome(int argv, char **argc)
{
	/* define */
	int gdOK,cdOK,mOK,iOK,oOK,mcOK,mdOK,cOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCSPath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nMC,nMD;
	double dC;
	int nIncludeRepeat = 0;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_consensusscan_genome */
	/* -m motif consensus              */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -mc number of stringent         */
	/*    consensus mismatches allowed */
	/* -md number of relaxed mismatches*/
	/*    allowed                      */
	/* -c conservation cutoff\n")      */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("  motifmap_consensusscan_genome \n");
		printf(" -m motif consensus (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -mc number of stringent consensus mismatches allowed\n");
		printf(" -md number of relaxed consensus mismatches allowed\n");
		printf(" -c conservation cutoff\n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_consensusscan_genome -m E2F.txt -gd mm6 -i inputseq.cod -o outputseq.txt -mc 2 -md 0 -c 200 -cd /mm6/conservation/phastcons\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	mcOK = 0;
	mdOK = 0;
	cOK = 0;
	cdOK = 0;
	nMC = 0;
	nMD = 0;
	dC = 0.0;
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-mc") == 0)
		{
			ni++;
			nMC = atoi(argc[ni]);
			mcOK = 1;
		}
		else if(strcmp(argc[ni], "-md") == 0)
		{
			ni++;
			nMD = atoi(argc[ni]);
			mdOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (mcOK == 0) || (mdOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, strCodPath,  
						strOutputPath, nMC, nMD, cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the path where conservation score is stored need to be specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, strCodPath,
					strOutputPath, nMC, nMD, cOK, 0.0, "", nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_filter_genome(int argv, char **argc)
{
	/* define */
	int iOK,oOK,cOK,cmOK,cdOK,cdsOK,cdsdOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	double dCds;
	char strCdsPath[MED_LINE_LENGTH];
	
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_filter_genome        */
	/* -i input site file              */
	/* -o output file                  */
	/* -c conservation cutoff          */
	/* -cm conservation mask           */
	/* -cd conservation score data path*/
	/* -cds coding sequence cutoff     */
	/* -cdsd coding sequence data path */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("  motifmap_filter_genome \n");
		printf(" -i input site file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff\n");
		printf(" -cm conservation mask\n");
		printf(" -cd conservation score data path\n");
		printf(" -cds coding sequence cutoff     \n");
		printf(" -cdsd coding sequence data path \n");
		printf(" example: \n");
		printf("    motifmap_filter_genome -i Gli.map -o Gli_filtered.map -c 1 -cm Gli_mask.txt -cd mm6/conservation -cds 0.9 -cdsd mm6/cds\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	oOK = 0;
	cOK = 0;
	cmOK = 0;
	cdOK = 0;
	cdsOK = 0;
	cdsdOK = 0;
	dC = 0.0;
	dCds = 0.0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cm") == 0)
		{
			ni++;
			strcpy(strMaskPath, argc[ni]);
			cmOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			dCds = atof(argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-cdsd") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0) || ((cOK == 0) && (cdsOK == 0)) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	

	if(cOK == 1)
	{
		if( (cdOK == 0) || (cmOK == 0) )
		{
			printf("Error: Input Parameter not correct!\n");
			exit(EXIT_FAILURE);
		}
	}

	if(cdsOK == 1)
	{
		if(cdsdOK == 0)
		{
			printf("Error: Input Parameter not correct!\n");
			exit(EXIT_FAILURE);
		}
	}
		
	nResult = MotifMap_Filter_Genome_Main(strInputPath, strOutputPath, 
		cOK, dC, strMaskPath, strCSPath, cdsOK, dCds, strCdsPath);
	
	/* return */
	return nResult;
}

int menu_motifmap_getsitearoundcs_genome(int argv, char **argc)
{
	/* define */
	int iOK,oOK,lOK,wOK,sOK,gdOK,cdOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nW,nMotifLen;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/* motifmap_getsitearoundcs        */
	/* -i input site file              */
	/* -o output file                  */
	/* -l motif length                 */
	/* -w flanking window size         */
	/* -s species                      */
	/* -gd genome sequence data path   */
	/* -cd conservation score data path*/
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("  motifmap_getsitearoundcs \n");
		printf(" -i input site file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -l motif length \n");
		printf(" -w flanking window size\n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -gd genome sequence data path\n");
		printf(" -cd conservation score data path\n");
		printf(" example: \n");
		printf("    motifmap_getsitearoundcs -i Gli.map -o Gli_cscurve.txt -l 12 -w 50 -s mouse -gd mm6/ -cd mm6/conservation\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	oOK = 0;
	lOK = 0;
	wOK = 0;
	sOK = 0;
	gdOK = 0;
	cdOK = 0;
	nW = 0;
	nMotifLen = -1;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMotifLen = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (wOK == 0) || (sOK == 0) || (gdOK == 0) || (cdOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	

	nResult = MotifMap_GetSiteAroundCS_Genome_Main(strInputPath, strOutputPath, 
		nMotifLen, nW, strSpecies, strGenomePath, strCSPath);
	
	/* return */
	return nResult;
}

int menu_motifmap_getsitearound(int argv, char **argc)
{
	/* define */
	int iOK,wOK,dOK,sOK,oOK,cnOK,aOK;
	char strInputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nW,nCN,nA;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*     motifmap_getsitearound      */
	/* -i input coordinates file       */
	/* -w half flanking window size    */ 
	/* -d genome sequence path         */
	/* -s species                      */
	/* -cn numerical chromosome id     */
	/* -o output file                  */
	/* -a alias type                   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_getsitearound  \n");
		printf(" -i input coordinates file   \n");
		printf(" -w half flanking window size\n");
		printf(" -d genome sequence path     \n");
		printf(" -s species                  \n");
		printf(" -cn numerical chromosome id (0: string chr id; 1: integer chr id\n");
		printf(" -o output file \n");
		printf(" -a alias type (0: use original, 1: reorder) \n");
		printf(" example: \n");
		printf("    motifmap_getsitearound -i Gli_map.txt -w 200 -d /data/genomes/mouse/mm6/ -s mouse -cn 1 -a 1 -o Gli_sitearound.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	wOK = 0;
	dOK = 0;
	sOK = 0;
	cnOK = 0;
	oOK = 0;
	aOK = 0;
	
	nW = 0;
	nCN = 0;
	nA = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nA = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-cn") == 0)
		{
			ni++;
			nCN = atoi(argc[ni]);
			cnOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (wOK == 0) || (dOK == 0) || (sOK == 0) || (oOK == 0) || (aOK == 0) || (cnOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_GetSiteAround_Main(strInputPath, strGenomePath, 
			strOutputPath, strSpecies, nW, nCN, nA);
	}

	/* return */
	return nResult;
}

int menu_motifmap_getcluster(int argv, char **argc)
{
	/* define */
	int iOK,wOK,oOK,rOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nW = 500;
	int nInputType = 0;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    menu_motifmap_getcluster     */
	/* -i input coordinates file       */
	/* -w window size                  */ 
	/* -o output file                  */
	/* -r input type, 0:cod, 1:bed,    */
	/*                2:codp,3:bedp    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    menu_motifmap_getcluster  \n");
		printf(" -i input coordinates file   \n");
		printf(" -w window size [default = 500] \n");
		printf(" -o output file \n");
		printf(" -r input type (0:cod, 1:bed, 2:codp,3:bedp) [default = 0] \n");
		printf(" example: \n");
		printf("    menu_motifmap_getcluster -i Gli_map.txt -w 200 -o Gli_clusteredsites.txt -r 0\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	wOK = 0;
	oOK = 0;
	rOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_getcluster_Main(strInputPath, strOutputPath, nW, nInputType);
	}

	/* return */
	return nResult;
}

int menu_motifmap_countkmer(int argv, char **argc)
{
	/* define */
	int dOK,iOK,oOK,kOK,cOK,chOK;
	char strInPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	double dC;
	int nK = 4;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*    menu_motifmap_countkmer      */
	/* -d sequence & conservation score data path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -k kmer length                  */
	/* -c conservation cutoff\n")      */
	/* -ch head tag for conservation score *.cs file */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    menu_motifmap_countkmer        \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file                  \n");
		printf(" -o output file (full path)        \n");
		printf(" -k kmer length                    \n");
		printf(" -c conservation cutoff\n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_countkmer -d . -i inputseq.fa -o kmer.txt -k 5 -c 200 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	dOK = 0;
	iOK = 0;
	oOK = 0;
	kOK = 0;
	cOK = 0;
	chOK = 0;
	nK = 4;
	dC = 0.0;
	strcpy(strCSAlias, "");
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nK = atoi(argc[ni]);
			kOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			strcpy(strCSAlias, argc[ni]);
			chOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (iOK == 0) || (oOK == 0) || (kOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(chOK == 1)
			{
				nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,  
						strOutPath, nK, cOK, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,  
						strOutPath, nK, cOK, dC, "");
			}
		}
		else
		{
			nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,
					strOutPath, nK, cOK, 0.0, "");
		}
	}

	/* return */
	return nResult;
}

int menu_motifmap_rescore(int argv, char **argc)
{
	/* define */
	int mOK,iOK,oOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*       motifmap_rescore          */
	/* -m motif matrix                 */
	/* -i sequence file                */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_rescore            \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -i sequence file  (full path)  \n");
		printf(" -o output file (full path) \n");
		printf(" example: \n");
		printf("    motifmap_rescore -m Gli.txt -i Gli.map -o Gli.remap\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_Rescore_Main(strMotifPath, strInputPath, strOutputPath);
	}

	/* return */
	return nResult;
}

int menu_motif_simuscoredistn_typeI()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	FILE *fpMotifList;
	
	/* background */
	char strWorkPath[LINE_LENGTH];
	char strBGFPath[LINE_LENGTH];
	char strBGRPath[LINE_LENGTH];
	char strBG0Path[LINE_LENGTH];
	struct DOUBLEMATRIX *pBGF;
	struct DOUBLEMATRIX *pBGR;
	struct DOUBLEMATRIX *pBG0;
	int nBGOrder;
	struct DOUBLEMATRIX *pMotifPWM;

	/* parameter */
	int nSimuNum;
	struct DOUBLEMATRIX *pQ;

	/* init */
	nBGOrder = 3;
	nSimuNum = 1000000;
	pQ = NULL;
	pQ = CreateDoubleMatrix(1, 10);
	if(pQ == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot create parameter!\n");
		exit(EXIT_FAILURE);
	}
	pQ->pMatElement[0] = 0.2;
	pQ->pMatElement[1] = 0.1;
	pQ->pMatElement[2] = 0.02;
	pQ->pMatElement[3] = 0.01;
	pQ->pMatElement[4] = 0.002;
	pQ->pMatElement[5] = 0.001;
	pQ->pMatElement[6] = 0.0002;
	pQ->pMatElement[7] = 0.0001;
	pQ->pMatElement[8] = 0.00002;
	pQ->pMatElement[9] = 0.00001;

	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\cutoff\\");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\motifmap\\");
	sprintf(strBGFPath, "%segfr_vivo_mcbg_%d_f.txt", strWorkPath, nBGOrder);
	sprintf(strBGRPath, "%segfr_vivo_mcbg_%d_r.txt", strWorkPath, nBGOrder);
	sprintf(strBG0Path, "%segfr_vivo_mcbg_0_f.txt", strWorkPath);

	/* load background */
	pBG0 = NULL;
	pBG0 = DMLOAD(strBG0Path);
	if(pBG0 == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BG0!\n");
		exit(EXIT_FAILURE);
	}
	pBGF = NULL;
	pBGF = DMLOAD(strBGFPath);
	if(pBGF == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BGF!\n");
		exit(EXIT_FAILURE);
	}
	pBGR = NULL;
	pBGR = DMLOAD(strBGRPath);
	if(pBGR == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BGR!\n");
		exit(EXIT_FAILURE);
	}

	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		printf("%s\n", strLine);

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.txt", strMotifPath, strMotifName);
		pMotifPWM = NULL;
		pMotifPWM = DMLOAD(strLine);
		if(pMotifPWM == NULL)
		{
			printf("Error: menu_motif_simuscoredistn_typeI, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		MotifMap_SimuScoreDistribution_typeI(strMotifName, pMotifPWM, nBGOrder, pBG0, pBGF, pBGR, nSimuNum, pQ, strCutoffPath);

		DestroyDoubleMatrix(pMotifPWM);

	}
	fclose(fpMotifList);

	/* destroy */
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGR);
	DestroyDoubleMatrix(pQ);

	/* return */
	return PROC_SUCCESS;
}

int menu_motif_simuscoredistn_typeII()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	FILE *fpMotifList;
	
	/* background */
	struct DOUBLEMATRIX *pMotifPWM;
	
	/* parameter */
	int nSimuNum;
	struct DOUBLEMATRIX *pQ;

	/* init */
	nSimuNum = 10000;
	pQ = NULL;
	pQ = CreateDoubleMatrix(1, 9);
	if(pQ == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeII, cannot create parameter!\n");
		exit(EXIT_FAILURE);
	}
	pQ->pMatElement[0] = 0.5;
	pQ->pMatElement[1] = 0.4;
	pQ->pMatElement[2] = 0.3;
	pQ->pMatElement[3] = 0.2;
	pQ->pMatElement[4] = 0.1;
	pQ->pMatElement[5] = 0.05;
	pQ->pMatElement[6] = 0.01;
	pQ->pMatElement[7] = 0.005;
	pQ->pMatElement[8] = 0.001;
	
	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\mouse\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\mouse\\cutoff\\");
	
	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeII, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		printf("%s\n", strLine);

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.txt", strMotifPath, strMotifName);
		pMotifPWM = NULL;
		pMotifPWM = DMLOAD(strLine);
		if(pMotifPWM == NULL)
		{
			printf("Error: menu_motif_simuscoredistn_typeII, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		MotifMap_SimuScoreDistribution_typeII(strMotifName, pMotifPWM, nSimuNum, pQ, strCutoffPath);

		DestroyDoubleMatrix(pMotifPWM);

	}
	fclose(fpMotifList);

	/* destroy */
	DestroyDoubleMatrix(pQ);

	/* return */
	return PROC_SUCCESS;
}


int menu_motif_simuscoredistn_typeIII()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	FILE *fpMotifList;
	
	/* background */
	char strWorkPath[LINE_LENGTH];
	char strBGFPath[LINE_LENGTH];
	char strBGRPath[LINE_LENGTH];
	char strBG0Path[LINE_LENGTH];
	struct DOUBLEMATRIX *pBGF;
	struct DOUBLEMATRIX *pBGR;
	struct DOUBLEMATRIX *pBG0;
	int nBGOrder;
	struct DOUBLEMATRIX *pMotifPWM;

	/* parameter */
	int nSimuNum;
	struct DOUBLEMATRIX *pQ;

	/* init */
	nBGOrder = 3;
	nSimuNum = 1000000;
	pQ = NULL;
	pQ = CreateDoubleMatrix(1, 12);
	if(pQ == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeIII, cannot create parameter!\n");
		exit(EXIT_FAILURE);
	}
	pQ->pMatElement[0] = 0.00001;
	pQ->pMatElement[1] = 0.00002;
	pQ->pMatElement[2] = 0.0001;
	pQ->pMatElement[3] = 0.0002;
	pQ->pMatElement[4] = 0.001;
	pQ->pMatElement[5] = 0.002;
	pQ->pMatElement[6] = 0.01;
	pQ->pMatElement[7] = 0.02;
	pQ->pMatElement[8] = 0.1;
	pQ->pMatElement[9] = 0.2;
	pQ->pMatElement[10] = 0.5;
	pQ->pMatElement[11] = 0.9;

	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\cutoff\\");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\motifmap\\");
	sprintf(strBGFPath, "%segfr_vivo_mcbg_%d_f.txt", strWorkPath, nBGOrder);
	sprintf(strBGRPath, "%segfr_vivo_mcbg_%d_r.txt", strWorkPath, nBGOrder);
	sprintf(strBG0Path, "%segfr_vivo_mcbg_0_f.txt", strWorkPath);

	/* load background */
	pBG0 = NULL;
	pBG0 = DMLOAD(strBG0Path);
	if(pBG0 == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BG0!\n");
		exit(EXIT_FAILURE);
	}
	pBGF = NULL;
	pBGF = DMLOAD(strBGFPath);
	if(pBGF == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BGF!\n");
		exit(EXIT_FAILURE);
	}
	pBGR = NULL;
	pBGR = DMLOAD(strBGRPath);
	if(pBGR == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load BGR!\n");
		exit(EXIT_FAILURE);
	}

	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_simuscoredistn_typeI, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		printf("%s\n", strLine);

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.txt", strMotifPath, strMotifName);
		pMotifPWM = NULL;
		pMotifPWM = DMLOAD(strLine);
		if(pMotifPWM == NULL)
		{
			printf("Error: menu_motif_simuscoredistn_typeI, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		MotifMap_SimuScoreDistribution_typeIII(strMotifName, pMotifPWM, nBGOrder, pBG0, pBGF, pBGR, nSimuNum, pQ, strCutoffPath);

		DestroyDoubleMatrix(pMotifPWM);

	}
	fclose(fpMotifList);

	/* destroy */
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGR);
	DestroyDoubleMatrix(pQ);

	/* return */
	return PROC_SUCCESS;
}

int menu_motif_getcutoff_typeI()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	FILE *fpMotifList;
	FILE *fpOut;
	double dR;
	
	/* background */
	struct DOUBLEMATRIX *pCutoff;
	
	/* init */
	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\cutoff\\");
	sprintf(strOutputFile, "%scutoff.txt", strMotifPath);
	
	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot output!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.stat", strCutoffPath, strMotifName);
		pCutoff = NULL;
		pCutoff = DMLOAD(strLine);
		if(pCutoff == NULL)
		{
			printf("Error: menu_expression_geneselection, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		dR = pCutoff->pMatElement[5];
		if(dR < 10.0)
			dR = 10.0;
		fprintf(fpOut, "%f\n", dR);
		DestroyDoubleMatrix(pCutoff);

	}
	fclose(fpMotifList);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

int menu_motif_getcutoff_typeII()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	FILE *fpMotifList;
	FILE *fpOut;
	double dR;
	
	/* background */
	struct DOUBLEMATRIX *pCutoff;
	
	/* init */
	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\mouse\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\mouse\\cutoff\\");
	sprintf(strOutputFile, "%scutoff_typeII.txt", strMotifPath);
	
	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot output!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.true", strCutoffPath, strMotifName);
		pCutoff = NULL;
		pCutoff = DMLOAD(strLine);
		if(pCutoff == NULL)
		{
			printf("Error: menu_expression_geneselection, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		dR = pCutoff->pMatElement[5];
		fprintf(fpOut, "%f\n", dR);
		DestroyDoubleMatrix(pCutoff);
	}
	fclose(fpMotifList);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

int menu_motif_getcutoff_typeIII()
{
	/* define */
	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strCutoffPath[LINE_LENGTH];
	char strMotifName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	FILE *fpMotifList;
	FILE *fpOut;
	double dR;
	
	/* background */
	struct DOUBLEMATRIX *pCutoff;
	
	/* init */
	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\matrix\\");
	strcpy(strMotifList, "MatrixList.txt");
	strcpy(strCutoffPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\cutoff\\");
	sprintf(strOutputFile, "%scutoff_typeIII.txt", strMotifPath);
	
	/* process motif one by one */
	sprintf(strLine, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strLine, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot load motif list!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: menu_motif_getcutoff, cannot output!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LINE_LENGTH, fpMotifList)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		strcpy(strMotifName, strLine);
		sprintf(strLine, "%s%s.pows", strCutoffPath, strMotifName);
		pCutoff = NULL;
		pCutoff = DMLOAD(strLine);
		if(pCutoff == NULL)
		{
			printf("Error: menu_expression_geneselection, cannot load motif!\n");
			exit(EXIT_FAILURE);
		}

		dR = pCutoff->pMatElement[4];
		fprintf(fpOut, "%f\n", dR);
		DestroyDoubleMatrix(pCutoff);

	}
	fclose(fpMotifList);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


int menu_motifsampler_quality()
{
	char strWorkPath[LINE_LENGTH];
	struct DOUBLEMATRIX *pFreqPrior;
	int nMotifNum;
	struct DOUBLEMATRIX *pSeg;
	struct DOUBLEMATRIX *pMatrixPrior[MAX_MOTIF_NUM];
	int nSeedNum;
	struct DOUBLEMATRIX *pSeedPrior[MAX_MOTIF_NUM];
	int ni;
	char strFilePath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];

	/* load data */
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\genomelab_project\\projects\\motifsampler\\");
	sprintf(strFilePath, "%sbgseg.txt", strWorkPath);
	pSeg = DMLOAD(strFilePath);
	sprintf(strFilePath, "%sfreqpriorquality.txt", strWorkPath);
	pFreqPrior = DMLOAD(strFilePath);
	
	
	nMotifNum = 3;
	for(ni=0; ni<MAX_MOTIF_NUM; ni++)
	{
		pMatrixPrior[ni] = NULL;
		pSeedPrior[ni] = NULL;
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		sprintf(strFilePath, "%stfbsprior.txt", strWorkPath);
		pMatrixPrior[ni] = DMLOAD(strFilePath);
	}

	nSeedNum = 1;
	for(ni=0; ni<nSeedNum; ni++)
	{
		sprintf(strFilePath, "%stfbsprior.txt", strWorkPath);
		pSeedPrior[ni] = DMLOAD(strFilePath);
	}

	/* quality */
	sprintf(strSeqFile, "%sFM_3motifs_3.txt", strWorkPath);
	sprintf(strOutFile, "%sFM_3motifs_3_qualityout.txt", strWorkPath);
	MotifSampler_Quality(strSeqFile, strOutFile, 1000, 
						  1, pFreqPrior,
						  3, pMatrixPrior,
						  1, pSeedPrior);

	/* return */
	return PROC_SUCCESS;
}

/* do motif mapping */
int pipeline_knownmotifmapping_forgenelist()
{
	/* define */
	int nClassNum,ni;
	struct tagString **vClassPath;
	char strGenePath[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strAffyCSVAnnotPath[LINE_LENGTH];
	char strReducedAnnotPath[LINE_LENGTH];
	char strRefGenePath[LINE_LENGTH];

	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strHistoryFile[LINE_LENGTH];
	char strProbeRefGeneMapFile[LINE_LENGTH];
	char strUniqRefGeneFile[LINE_LENGTH];
	char strRefForMapFile[LINE_LENGTH];
	char strGenomicTargetFile[LINE_LENGTH];
	char strMarkovBGFile[LINE_LENGTH];
	char strChrLenFile[LINE_LENGTH];

	char strGenomePath[LINE_LENGTH];
	char strCScorePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strVersion[LINE_LENGTH];
	int nChrNum;
	int nBGOrder;

	char strMotifPath[LINE_LENGTH];
	char strMotifList[LINE_LENGTH];
	char strMapOutPath[LINE_LENGTH];
	char strTargetMotifList[LINE_LENGTH];
	char strFilteredOutPath[LINE_LENGTH];
	char strLikehoodRatioCutoffList[LINE_LENGTH];
	char strTypeIICutoffList[LINE_LENGTH];
	char strTypeIICutoffPath[LINE_LENGTH];
	char strScoreOutPath[LINE_LENGTH];

	char strAllScorePath[LINE_LENGTH];
	char strAnnotScorePath[LINE_LENGTH];

	FILE *fpOut;

	double dCSCutoff;
	int nRepeatMask;
	int nTypeIIMask;

	int nTSSUP, nTSSDOWN, nTESUP, nTESDOWN;
	int nIncludeIntron, nIncludeExon;
	

	/* init */
	nBGOrder = 3;
	nChrNum = 24;
	nClassNum = 4;
	vClassPath = (struct tagString **)calloc(nClassNum, sizeof(struct tagString *));
	if(vClassPath == NULL)
	{
		printf("Error: cannot create memory to store class names!\n");
		return PROC_FAILURE;
	}
	for(ni=0; ni<nClassNum; ni++)
	{
		vClassPath[ni] = CreateString(LINE_LENGTH);
	}

	strcpy(strSpecies, "human");
	strcpy(strVersion, "hg16");

	strcpy(strGenePath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\motifmap\\");
	strcpy(vClassPath[0]->m_pString, "egfrm_vivowd_up");
	strcpy(vClassPath[1]->m_pString, "egfrm_vivowd_up_ctr");
	strcpy(vClassPath[2]->m_pString, "egfrm_vivowd_down");
	strcpy(vClassPath[3]->m_pString, "egfrm_vivowd_down_ctr");
	sprintf(strHistoryFile, "%sprobemap_all.txt", strWorkPath);
	sprintf(strProbeRefGeneMapFile, "%segfrm_vivo", strWorkPath);
	sprintf(strUniqRefGeneFile, "%segfrm_vivo_refgene.txt", strWorkPath);
	sprintf(strRefForMapFile, "%segfrm_vivo_uniref.txt", strWorkPath);
	sprintf(strGenomicTargetFile, "%segfrm_vivo_genomictarget.txt", strWorkPath);
	/* sprintf(strGenomicTargetFile, "%stest_genomictarget.txt", strWorkPath); */
	sprintf(strMarkovBGFile, "%segfr_vivo_mcbg", strWorkPath);
	
	strcpy(strAffyCSVAnnotPath, "C:\\Projects\\research_harvard\\affy_project\\data\\NetaffyHs\\HG-U133A_annot.csv");
	strcpy(strReducedAnnotPath, "C:\\Projects\\research_harvard\\genomelab_project\\affymetrix\\human\\HG-U133A_reducedannot.txt");
	strcpy(strRefGenePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\annotation\\refGene_sorted.txt");
	strcpy(strGenomePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");
	strcpy(strCScorePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\");
	sprintf(strChrLenFile, "%schrlen.txt", strGenomePath);
	strcpy(strMotifPath, "C:\\Projects\\research_harvard\\motif_project\\knownmotifs\\human\\matrix\\");
	/*strcpy(strMotifList, "StatMatrixList.txt");
	strcpy(strLikehoodRatioCutoffList, "StatLikeRatioCutoff.txt");
	strcpy(strTypeIICutoffList, "StatTypeIICutoff.txt");
	*/
	strcpy(strMotifList, "MatrixList8.txt");
	strcpy(strLikehoodRatioCutoffList, "cutoff8.txt");
	strcpy(strTypeIICutoffList, "cutoff_typeII_8.txt");
	sprintf(strScoreOutPath, "%segfr_vivo_motifscore_8", strWorkPath);
	sprintf(strAllScorePath, "%segfr_vivo_motifscore_csall_incintron.txt", strWorkPath);
	sprintf(strAnnotScorePath, "%segfr_vivo_motifscore_csall_incintron_annot.txt", strWorkPath);
	
	strcpy(strMapOutPath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\motifmap\\all\\");
	strcpy(strFilteredOutPath, "C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\motifmap\\repfilt\\");

	nRepeatMask = 1;
	nTypeIIMask = 1;
	dCSCutoff = -1.0;
	/* dCSCutoff = 38.0; */
	nTSSUP = 50000;
	nTSSDOWN = 10000;
	nTESUP = 10000;
	nTESDOWN = 10000;
	nIncludeIntron = 1;
	nIncludeExon = 0;

	/* ------------------------------------------------------------------ */
	/* step 0: prepare reduced annotation file from affy's csv annotation */
	/* ------------------------------------------------------------------ */

	/* Affy_CSVANNOT_To_Reduced_200408(strAffyCSVAnnotPath,	strReducedAnnotPath, strSpecies); */

	/* ------------------------------------------------------------------ */
	/* step 1: select genes and random control                            */
	/* ------------------------------------------------------------------ */

	/* Expression_GeneSelection_Main("C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\deltaegfrbatchlognorm092004.txt",
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\HG-U133Aannot.txt", 
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\deltaegfrcompinfo.txt",
						"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\egfrm_vivowd_down");
	*/

	/* ------------------------------- Below no more needed --------------------------------------------- */
	/* Affy_PickRandomControlFromCSV(strAffyCSVAnnotPath, 
		"C:\\Projects\\research_harvard\\gbm_project\\EGFR_Frank\\pooled\\egfrm_randomcontrol.txt", 300); */
	/* -------------------------------------------------------------------------------------------------- */

	/* ------------------------------------------------------------------ */
	/* step 2: link probesets to reduced annotation                       */
	/* ------------------------------------------------------------------ */
	/* fpOut = NULL;
	fpOut = fopen(strHistoryFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: pipeline_knownmotifmapping_forgenelist, cannot archive genemap lists!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nClassNum; ni++)
	{
		sprintf(strInFile, "%s%s.txt", strGenePath, vClassPath[ni]->m_pString);
		sprintf(strOutFile, "%s%s_map.txt", strWorkPath, vClassPath[ni]->m_pString);
		Affy_LinkProbesetToAnnotation(strInFile, strReducedAnnotPath, strOutFile, ni);
		fprintf(fpOut, "%d\t%s\n", ni, strOutFile);
	}
	fclose(fpOut);
	*/

	/* ------------------------------------------------------------------ */
	/* step 3: link reduced annotation with refGene annotation to get     */
	/* regions in the genome to be mapped.                                */
	/* ------------------------------------------------------------------ */

	/* Affy_CreateProbeRefGeneMap_BothDirection(strHistoryFile, strRefGenePath, 
						strProbeRefGeneMapFile, strSpecies); */

	/* RefGene_GetTargetGenomicRegion_TSSTES(strUniqRefGeneFile, 
			strGenomicTargetFile, 50000, 50000, strSpecies, 24, strChrLenFile); */

	/* ------------------------------------------------------------------ */
	/* step 4: map the motifs to the target regions.                      */
	/* ------------------------------------------------------------------ */

	/* get markov backgroud */
	/* Genome_GetMarkovBGForRegion(strGenomePath, strGenomicTargetFile, strSpecies, nBGOrder, strMarkovBGFile); */
	/* Genome_GetMarkovBGForRegion(strGenomePath, strGenomicTargetFile, strSpecies, 0, strMarkovBGFile); */

	/* map motifs one by one */
	/* Genome_MapMotifToRegion(strGenomePath, strGenomicTargetFile, strSpecies, strVersion,
		strMotifPath, strMotifList, nBGOrder, strMarkovBGFile, strLikehoodRatioCutoffList, strMapOutPath);
	*/

	/* filter motifs */
	sprintf(strTargetMotifList, "%s%s", strMotifPath, strMotifList);
	if(nTypeIIMask == 1)
	{
		sprintf(strTypeIICutoffPath, "%s%s", strMotifPath, strTypeIICutoffList);
	}
	/* Genome_FilterMotifSites(strGenomePath, strCScorePath, strSpecies, strVersion,
		strMotifPath, strMotifList, strMapOutPath, strFilteredOutPath, dCSCutoff,
		nRepeatMask, nTypeIIMask, strTypeIICutoffList);
	*/

	/* create *.bed files for display */

	/* ------------------------------------------------------------------ */
	/* step 5: create mapping score for every gene                        */
	/* ------------------------------------------------------------------ */
	/* score refgene */
	/* MotifMap_ScoreRefGene(strSpecies, strVersion, strGenomePath, strCScorePath,
		strTargetMotifList, strFilteredOutPath, dCSCutoff, nRepeatMask,
		strRefForMapFile, nTSSUP, nTSSDOWN, nTESUP, nTESDOWN, 
		nIncludeIntron, nIncludeExon, strScoreOutPath);
	*/

	/* ------------------------------------------------------------------ */
	/* step 6: analyze the output mapping matrix and get significant      */
	/* class-motif relationship.                                          */
	/* ------------------------------------------------------------------ */

	/* ------------------------------------------------------------------ */
	/* step 7: output the class-motif information                         */
	/* ------------------------------------------------------------------ */

	/* link back to gene symbols */
	Affy_LinkRefGeneScoreBackToProbeset(strHistoryFile, strProbeRefGeneMapFile,
		strAllScorePath, strAnnotScorePath);


	/* ------------------------------------------------------------------ */
	/* step 8: verify the result using expression data.                   */
	/* link motif back to the expression of TF.                           */
	/* ------------------------------------------------------------------ */


	/* release memory */
	for(ni=0; ni<nClassNum; ni++)
	{
		DeleteString(vClassPath[ni]);
		vClassPath[ni] = NULL;
	}
	free(vClassPath);

	/* return */
	return PROC_SUCCESS;
}


int testcoding()
{
	int nLen;
	char strFastaFile[LINE_LENGTH];
	char strCodeFile[LINE_LENGTH];

	strcpy(strFastaFile, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b35_hg17\\chr21.fa");
	strcpy(strCodeFile, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b35_hg17\\chr21.sq");

	nLen = Genome_Fasta_To_Code_4bit(strFastaFile, strCodeFile);

	return nLen;
}

int testloading()
{
	struct tagSequence *pSeq;
	char strFastaFile[LINE_LENGTH];
	char strCodeFile[LINE_LENGTH];
	FILE *fpOut;

	strcpy(strFastaFile, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b35_hg17\\test.fa");
	strcpy(strCodeFile, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b35_hg17\\chr21.sq");
	pSeq = NULL;
	pSeq = Genome_Code_4bit_GetSeq(strCodeFile, 10050122, 10051011);
	
	fpOut = NULL;
	fpOut = fopen(strFastaFile, "wt");
	SequenceWriteToFasta(pSeq, fpOut);
	fclose(fpOut);

	SequenceDelete(pSeq);

	return 1;
}

int testnetwork()
{
	struct NETWORKPHYSICAL *pNet = NULL;
	/* pNet = Network_LoadFromSIF("test.sif");
	Network_WriteToSIF(pNet, "test_copy.sif");
	NETWORKPHYSICALDESTROY(&pNet); */

	Network_CreateOrthogloNet_Main("E:\\Data\\HPRD\\HPRD_v6_human.txt", 
		"E:\\Data\\Homologene\\build48\\homologene.data", 9606, 10090,
		"E:\\Data\\HPRD\\HPRD_v6_mouse.txt");

	/* Network_FindShortestPath_Main("test.sif", "NULL",
								  "test_src.txt", "test_dst.txt", 
								  "test.out", 20); */

	return 1;
}

int testnormcdf()
{
	struct DOUBLEMATRIX *pX;
	int ni,nj;
	double dx,dy;

	pX = CreateDoubleMatrix(1001,2);
	for(ni=0; ni<1001; ni++)
	{
		dx = -20.0+ni*0.04;
		dy = normcdf(0.0, 1.0, dx);
		DMSETAT(pX, ni, 0, dx);
		DMSETAT(pX, ni, 1, dy);
	}

	DMSAVE(pX, "E:\\Projects\\hedgehog_project\\Test\\normcdf.txt");
	return 1;
}

int testcel()
{
	struct tagCELData *pCELData = NULL;
	char strFile[255];
 
	/* strcpy(strFile, "D:\\Projects\\estrogen_project\\ER\\hWG_DF1_MCF7_ER_AS_C01_B1_T1.CEL"); */
	strcpy(strFile, "D:\\Data\\test\\Viji_Control_1.CEL"); 
    /* pCELData = Affy_LoadCELv4(strFile); */
	pCELData = Affy_LoadCEL_CmdCslv1(strFile);
	strcpy(strFile, "D:\\Data\\test\\Viji_Control_1_rewrite.CEL");
	Affy_SaveCELv3(strFile, pCELData);
	Affy_CELData_Destroy(&pCELData);

	return 1;
}

int testbar()
{
	struct tagBARData *pBARData = NULL;
	struct tagBARData *pBARData2 = NULL;
	struct tagBARData *pNewBARData = NULL;
	char strFile[255];
	char strTxtFile[255];
	FILE *fpOut;
	int ni,nj;
	struct INTMATRIX *pColIndicator = NULL;

	/*strcpy(strFile, "D:\\Projects\\rest_project\\chip-seq\\chip2002_uniq_hg17.txt.pos.bar");
    pBARData = Affy_LoadBAR(strFile);
	strcpy(strFile, "D:\\Projects\\rest_project\\chip-seq\\mock2002_uniq_hg17.txt.pos.bar");
	pBARData2 = Affy_LoadBAR(strFile); */

	/* strcpy(strFile, "D:\\Projects\\greensolexa_project\\signal\\s_5.eland.30.output.ualn.bar");
    pBARData = Affy_LoadBAR(strFile);
	strcpy(strFile, "D:\\Projects\\greensolexa_project\\signal\\s_6.eland.30.output.ualn.bar");
	pBARData2 = Affy_LoadBAR(strFile); 
	

	pColIndicator = CreateIntMatrix(1, pBARData->nColNum);
	pColIndicator->pMatElement[0] = 2;
	pColIndicator->pMatElement[1] = 1;
	*/

	// Affy_BARData_AbsTS(pNewBARData, pColIndicator);
	/* pNewBARData = Affy_BARData_Clone(pBARData);
	Affy_BARData_Pow(pNewBARData, pBARData2, pColIndicator, 0); */
	/* Affy_BARData_Sub(pBARData, pBARData2, pColIndicator, 0); */

	/* strcpy(strFile, "G:\\Projects\\cMyc_Project\\analysis\\ChIP-chip\\Bcell\\signal\\cMycB_Ab_peak_1.ma.bar"); */
	strcpy(strFile, "G:\\Projects\\cMyc_Project\\ChIP-chip\\GSE11329\\analysis\\bMyc_peak_1.ma.bar");
    pBARData = Affy_LoadBAR(strFile);
	pColIndicator = CreateIntMatrix(1, pBARData->nColNum);
	pColIndicator->pMatElement[1] = 1;
	Affy_BARData_UPowTS(pBARData, 2.0, pColIndicator);
	pColIndicator->pMatElement[0] = 1;
	/* strcpy(strFile, "G:\\Projects\\cMyc_Project\\analysis\\ChIP-chip\\Bcell\\signal\\cMycB_Ab_peak_1.mao.bar"); */
	strcpy(strFile, "G:\\Projects\\cMyc_Project\\ChIP-chip\\GSE11329\\analysis\\bMyc_peak_1.mao.bar");
	Affy_SaveBAR_Columns_Fast(strFile, pBARData, pColIndicator);
	Affy_BARData_Destroy(&pBARData);
	DestroyIntMatrix(pColIndicator);

	/* strcpy(strFile, "D:\\Projects\\rest_project\\chip-seq\\diff2002_uniq_hg17.txt.pos.bar");
	strcpy(strTxtFile, "D:\\Projects\\rest_project\\chip-seq\\diff2002_uniq_hg17.txt");
	Affy_SaveBAR(strFile, pBARData); */

	/* strcpy(strFile, "D:\\Projects\\greensolexa_project\\signal\\diff5-6.bar");
	strcpy(strTxtFile, "D:\\Projects\\greensolexa_project\\signal\\diff5-6.txt");
	Affy_SaveBAR(strFile, pBARData);
	
	
	fpOut = NULL;
	fpOut = fopen(strTxtFile, "w");
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			fprintf(fpOut, "%f\n", 1e6*pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]);
		}
	}
	fclose(fpOut);
	*/

	/* Affy_BARData_Destroy(&pBARData);
	Affy_BARData_Destroy(&pBARData2);
	Affy_BARData_Destroy(&pNewBARData);
	DestroyIntMatrix(pColIndicator); */

	// Affy_BAR2TXT_Fast(strFile, strTxtFile);

	return 1;
}

int testgeneric()
{
	struct tagGenericFileObject *pGenericObj = NULL;
	char strFilePath[1024];

	/* strcpy(strFilePath, "D:\\Projects\\affymetrix_project\\FusionSDK\\SP1_ARR_CEL_AGCC\\4013176_Ab-_2.2.CEL");
	/* strcpy(strFilePath, "D:\\Projects\\estrogen_project\\ER\\hWG_DF1_MCF7_ER_AS_C01_B1_T1.CEL");
	/* strcpy(strFilePath, "D:\\Projects\\affymetrix_project\\FusionSDK\\affy-fusion-release-109\\affy\\sdk\\sample_data\\Test3-1-121502.calvin.CHP");
	/* strcpy(strFilePath, "D:\\Projects\\affymetrix_project\\FusionSDK\\affy-fusion-release-109\\affy\\sdk\\TestDataFiles\\CHP\\StatData-calvin.CHP");
	*/

	/* pGenericObj = Affy_GenericFileObject_Load(strFilePath);

	Affy_GenericFileObject_Delete(&pGenericObj); */

	return 1;
}

int menu_hts_selectreads(int argv, char **argc)
{
	/* define */
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	double dR = 1.0;

	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*      hts_selectreads            */
	/* -i input file                   */
	/* -o output title                 */
	/* -r probability to include a read */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_selectreads           \n");
		printf(" -i input ALN file             \n");
		printf(" -o output ALN file            \n");
		printf(" -r probability to include a read. Default = 1.0 \n");
		printf(" example: \n");
		printf("    hts_selectreads -i input.txt -o output.txt -r 0.5 \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_selectread(strInputFile, strOutputFile, dR);
	}

	return nResult;
}

int menu_test_chipcluster()
{
	struct INTMATRIX *pChrLen = NULL;
	int nChrNum;
	float **vT;
	int ni,nj;

	pChrLen = IMLOAD("D:\\Data\\genomes\\human\\hg18\\chrlen.txt");
	nChrNum = pChrLen->nHeight;

	vT = NULL;
	vT = (float **)calloc(nChrNum, sizeof(float *));
	for(ni=0; ni<nChrNum; ni++)
	{
		nj = pChrLen->pMatElement[ni]/25+1;
		vT[ni] = (float *)calloc(nj, sizeof(float));
		if(vT[ni] == NULL)
		{
			printf("Error!\n");
		}
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		free(vT[ni]);
	}
	free(vT);


	return 1;

}

int testeigen()
{
	struct DOUBLEMATRIX *pX,*pD,*pV,*pI,*pZ,*pW;

	pX = NULL;
	pD = NULL;
	pV = NULL;
	pI = NULL;
	pZ = NULL;
	pW = NULL;

	pX = DMLOAD("G:\\Projects\\SeqClust_project\\Test\\testmatrix.txt");
	DMSYMEIGEN(pX, &pD, &pV);
	DMSYMEIGENQR(pX,&pD,&pV);

	DMSYMINVEIGEN(pD, pV, &pI);
	pZ = DMMUL(pI, pX);
	pW = DMMUL(pX, pI);
	DMSAVE(pD, "G:\\Projects\\SeqClust_project\\Test\\testeigenvalue.txt");
	DMSAVE(pV, "G:\\Projects\\SeqClust_project\\Test\\testeigenvector.txt");
	DMSAVE(pI, "G:\\Projects\\SeqClust_project\\Test\\testinv.txt");
	DMSAVE(pZ, "G:\\Projects\\SeqClust_project\\Test\\testiden1.txt");
	DMSAVE(pW, "G:\\Projects\\SeqClust_project\\Test\\testiden2.txt");
	

	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pD);
	DestroyDoubleMatrix(pV);
	DestroyDoubleMatrix(pI);
	DestroyDoubleMatrix(pZ);
	DestroyDoubleMatrix(pW);




	return 1;
}