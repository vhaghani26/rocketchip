CFLAGS=-O3
LDFLAGS=-lm -O3
BINDIR=../bin
ALL= file_bed2cod file_cod2bed file_eland2aln file_bed2aln \
  genome_encode genome_conservebg genome_conservecs genome_footprint \
  genome_csgetdistn genome_codephastcons genome_codephastcons_v2 \
  genome_getseq genome_getseqcs genome_getmaskedseq \
  genome_getmaskedreg genome_getseq_c genome_getseqcs_c \
  genome_getmaskedseq_c genome_getmaskedreg_c genome_codingCDS \
  genome_getcsgcsummary genome_regioncssummary genome_regionextend \
  genome_fastaseqmask genome_hash fasta_soft2hardmask refgene_encode \
  refgene_pickspeciesspecific refgene_gettssaround refgene_getmultiortholog \
  refgene_getmultiortholog1way refgene_getmultiorthologtssaround \
  refgene_getmultiorthologtssaroundexonmasked refgene_getnearestgene \
  refgene_getaffy refgene_getmatchedcontrol refgene_getlocationsummary \
  refflat_encode reflocus_createmap reflocus_getneighborgenes \
  reflocus_assignvalue refexon_createmap refmicroarray_createmap \
  refflex_getmultiortholog refflex_getmultiortholog1way affy_bar2txt \
  powexpress powexpressloc powexpress_getspecificprobe powexpress_getnrprobe \
  tileprobe_buildmodel tileprobe_buildmodel_v2 tileprobe_buildmodel_hmmt \
  tileprobe_buildmodel_hmmb tileprobe_buildmodel_hmmm \
  tileprobe_buildmodel_hmmw tileprobe_huberdata tileprobe_norm \
  tileprobe_peak tileprobe_mat tileprobe_buildbarmodel tileprobe_barnorm \
  tileprobe_buildbarmodel_v2 tileprobe_resample tilemap tilemap_norm \
  tilemap_importaffy tilemap_extract tilemapv2_importaffy \
  tilemapv2 tilemapv2_regioninfo tilemapv2_regioninfo_integral \
  tilemapv2_txt2bar tilemapv2_quantilenorm hts_aln2window hts_aln2bar \
  hts_alnshift2bar hts_aln2barv2 hts_windowsummary hts_windowsummaryv2 \
  hts_windowsummary_2sample hts_windowsummaryv2_2sample hts_peakdetector \
  hts_peakdetector_2sample hts_peakdetectorv2 hts_peakdetectorv2_2sample \
  hts_createrepeatfilter hts_filterrepeatreads hts_collectreads \
  hts_collectprofile hts_countreads4refgene cnv_aln2window \
  cnv_repeat2window rnaseq_countgeneread seqclust_seg seqclust_count \
  seqclust network_shortestpath network_createortholognet network_getsubnet \
  network_bind2entrez motifmap_consensusscan motifmap_consensusscan_genome \
  motifmap_matrixscan motifmap_matrixscan_genomebg \
  motifmap_matrixscan_genome motifmap_matrixintegrate_genome \
  motifmap_multiplematrixscan_genome motifmap_matrixscan_group \
  motifmap_getsitearound motifmap_getsitearoundcs motifmap_countkmer \
  motifmap_matrixscan_genome_enrich motifmap_consensusscan_genome_enrich \
  motifmap_matrixscan_genome_summary motifmap_consensusscan_genome_summary \
  motifmap_filter_genome motifmap_matrixscan_genome_getcutoff \
  motifmap_getcluster flexmodule_anchordist flexmodule_motif flexmodule_tnum \
  flexmodule_oneeach flexmodule_moreeach flexmodule_tnumrelax \
  flexmodule_oneeachrelax flexmodule_moreeachrelax flexmodule_tworeverse \
  flexmodule flexmodule_bgfit malign_blasthit malign_motifmap \
  malign_genome_prepareortholog malign_genome_blasthit malign_modulemap \
  malign_countkmer malign_generatemotifaln work_esclusgetrefflat \
  work_esgetpromoter work_shhmbtest region_overlap tablesorter tablesorter_str\
  file_eland2aln_v2 file_bowtie2aln hts_aln2uniq hts_aln2winbar \
  hts_selectreads seqclust_dp seqpeak affy_bar2wig
all: $(ALL)
bin: $(ALL)
	cp $(ALL) $(BINDIR)
clean:
	rm -f *.o $(ALL)
file_bed2cod: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o file_bed2cod.o
file_cod2bed: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o file_cod2bed.o
file_eland2aln: file_eland2aln.o
file_bed2aln: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o file_bed2aln.o
genome_encode: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_encode.o
genome_conservebg: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_conservebg.o
genome_conservecs: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_conservecs.o
genome_footprint: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_footprint.o
genome_csgetdistn: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_csgetdistn.o
genome_codephastcons: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_codephastcons.o
genome_codephastcons_v2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_codephastcons_v2.o
genome_getseq: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getseq.o
genome_getseqcs: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getseqcs.o
genome_getmaskedseq: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getmaskedseq.o
genome_getmaskedreg: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getmaskedreg.o
genome_getseq_c: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getseq_c.o
genome_getseqcs_c: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getseqcs_c.o
genome_getmaskedseq_c: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getmaskedseq_c.o
genome_getmaskedreg_c: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getmaskedreg_c.o
genome_codingCDS: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_codingCDS.o
genome_getcsgcsummary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_getcsgcsummary.o
genome_regioncssummary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_regioncssummary.o
genome_regionextend: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_regionextend.o
genome_fastaseqmask: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_fastaseqmask.o
genome_hash: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o genome_hash.o
fasta_soft2hardmask: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o fasta_soft2hardmask.o
refgene_encode: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_encode.o
refgene_pickspeciesspecific: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_pickspeciesspecific.o
refgene_gettssaround: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_gettssaround.o
refgene_getmultiortholog: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getmultiortholog.o
refgene_getmultiortholog1way: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getmultiortholog1way.o
refgene_getmultiorthologtssaround: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getmultiorthologtssaround.o
refgene_getmultiorthologtssaroundexonmasked: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getmultiorthologtssaroundexonmasked.o
refgene_getnearestgene: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getnearestgene.o
refgene_getaffy: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getaffy.o
refgene_getmatchedcontrol: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getmatchedcontrol.o
refgene_getlocationsummary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refgene_getlocationsummary.o
refflat_encode: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refflat_encode.o
reflocus_createmap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o reflocus_createmap.o
reflocus_getneighborgenes: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o reflocus_getneighborgenes.o
reflocus_assignvalue: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o reflocus_assignvalue.o
refexon_createmap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refexon_createmap.o
refmicroarray_createmap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refmicroarray_createmap.o
refflex_getmultiortholog: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refflex_getmultiortholog.o
refflex_getmultiortholog1way: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o refflex_getmultiortholog1way.o
affy_bar2txt: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o affy_bar2txt.o
powexpress: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o powexpress.o
powexpressloc: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o powexpressloc.o
powexpress_getspecificprobe: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o powexpress_getspecificprobe.o
powexpress_getnrprobe: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o powexpress_getnrprobe.o
tileprobe_buildmodel: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel.o
tileprobe_buildmodel_v2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel_v2.o
tileprobe_buildmodel_hmmt: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel_hmmt.o
tileprobe_buildmodel_hmmb: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel_hmmb.o
tileprobe_buildmodel_hmmm: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel_hmmm.o
tileprobe_buildmodel_hmmw: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildmodel_hmmw.o
tileprobe_huberdata: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_huberdata.o
tileprobe_norm: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_norm.o
tileprobe_peak: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_peak.o
tileprobe_mat: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_mat.o
tileprobe_buildbarmodel: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildbarmodel.o
tileprobe_barnorm: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_barnorm.o
tileprobe_buildbarmodel_v2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_buildbarmodel_v2.o
tileprobe_resample: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tileprobe_resample.o
tilemap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o tilemap.o
tilemap_norm: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o tilemap_norm.o
tilemap_importaffy: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o tilemap_importaffy.o
tilemap_extract: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o tilemap_extract.o
tilemapv2_importaffy: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2_importaffy.o
tilemapv2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2.o
tilemapv2_regioninfo: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2_regioninfo.o
tilemapv2_regioninfo_integral: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2_regioninfo_integral.o
tilemapv2_txt2bar: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2_txt2bar.o
tilemapv2_quantilenorm: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o tilemapv2_quantilenorm.o
hts_aln2window: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_aln2window.o
hts_aln2bar: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_aln2bar.o
hts_alnshift2bar: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_alnshift2bar.o
hts_aln2barv2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_aln2barv2.o
hts_windowsummary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_windowsummary.o
hts_windowsummaryv2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_windowsummaryv2.o
hts_windowsummary_2sample: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_windowsummary_2sample.o
hts_windowsummaryv2_2sample: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_windowsummaryv2_2sample.o
hts_peakdetector: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_peakdetector.o
hts_peakdetector_2sample: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_peakdetector_2sample.o
hts_peakdetectorv2: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_peakdetectorv2.o
hts_peakdetectorv2_2sample: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_peakdetectorv2_2sample.o
hts_createrepeatfilter: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_createrepeatfilter.o
hts_filterrepeatreads: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_filterrepeatreads.o
hts_collectreads: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_collectreads.o
hts_collectprofile: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_collectprofile.o
hts_countreads4refgene: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_countreads4refgene.o
cnv_aln2window: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o cnv_aln2window.o
cnv_repeat2window: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o cnv_repeat2window.o
rnaseq_countgeneread: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o rnaseq_countgeneread.o
seqclust_seg: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o seqclust_seg.o
seqclust_count: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o seqclust_count.o
seqclust: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o seqclust.o
network_shortestpath: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o FlexModuleInterface.o network_shortestpath.o
network_createortholognet: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o FlexModuleInterface.o network_createortholognet.o
network_getsubnet: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o FlexModuleInterface.o network_getsubnet.o
network_bind2entrez: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o FlexModuleInterface.o network_bind2entrez.o
motifmap_consensusscan: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_consensusscan.o
motifmap_consensusscan_genome: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_consensusscan_genome.o
motifmap_matrixscan: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan.o
motifmap_matrixscan_genomebg: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_genomebg.o
motifmap_matrixscan_genome: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_genome.o
motifmap_matrixintegrate_genome: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixintegrate_genome.o
motifmap_multiplematrixscan_genome: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_multiplematrixscan_genome.o
motifmap_matrixscan_group: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_group.o
motifmap_getsitearound: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_getsitearound.o
motifmap_getsitearoundcs: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_getsitearoundcs.o
motifmap_countkmer: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_countkmer.o
motifmap_matrixscan_genome_enrich: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_genome_enrich.o
motifmap_consensusscan_genome_enrich: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_consensusscan_genome_enrich.o
motifmap_matrixscan_genome_summary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_genome_summary.o
motifmap_consensusscan_genome_summary: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_consensusscan_genome_summary.o
motifmap_filter_genome: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_filter_genome.o
motifmap_matrixscan_genome_getcutoff: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_matrixscan_genome_getcutoff.o
motifmap_getcluster: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o motifmap_getcluster.o
flexmodule_anchordist: flexmodule
	ln -s $< $@
flexmodule_motif: flexmodule
	ln -s $< $@
flexmodule_tnum: flexmodule
	ln -s $< $@
flexmodule_oneeach: flexmodule
	ln -s $< $@
flexmodule_moreeach: flexmodule
	ln -s $< $@
flexmodule_tnumrelax: flexmodule
	ln -s $< $@
flexmodule_oneeachrelax: flexmodule
	ln -s $< $@
flexmodule_moreeachrelax: flexmodule
	ln -s $< $@
flexmodule_tworeverse: flexmodule
	ln -s $< $@
flexmodule: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o flexmodule.o
flexmodule_bgfit: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o flexmodule_bgfit.o
malign_blasthit: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_blasthit.o
malign_motifmap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_motifmap.o
malign_genome_prepareortholog: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_genome_prepareortholog.o
malign_genome_blasthit: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_genome_blasthit.o
malign_modulemap: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_modulemap.o
malign_countkmer: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_countkmer.o
malign_generatemotifaln: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o malign_generatemotifaln.o
work_esclusgetrefflat: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o work_esclusgetrefflat.o
work_esgetpromoter: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o work_esgetpromoter.o
work_shhmbtest: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o Alignment.o WorkLib.o FlexModuleInterface.o work_shhmbtest.o
tablesorter: tablesorter.cpp
tablesorter_str: tablesorter_str.cpp
region_overlap: region_overlap.cpp
file_eland2aln_v2: file_eland2aln_v2.o
file_bowtie2aln: file_bowtie2aln.o
hts_aln2uniq: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_aln2uniq.o
hts_aln2winbar: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_aln2winbar.o
hts_selectreads: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o hts_selectreads.o
seqclust_dp: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o seqclust_dp.o
seqpeak: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o TilingArrayLib.o HTSequencingLib.o seqpeak.o
affy_bar2wig: AffyLib.o GenomeLib.o MathLib.o MatrixLib.o MicroarrayLib.o MotifLib.o RandomLib.o SequenceLib.o StringLib.o PhysicalNetworkLib.o WorkLib.o FlexModuleInterface.o affy_bar2wig.o

