BASE_LIB_DIR = ../../
BASE_DIR=../../

BASE_EXEC_PATH=../../bin
BASE_EXEC_EXTRAS_PATH=../../bin/extras


TEST_MAP_SOURCE = generateMaps.cxx
TEST_MAP_EXEC = generateMaps

ASTMO_SOURCE = addStrandToMapOfObs.cxx
ASTMO_EXEC = addStrandToMapOfObs

CS_SOURCE = findConsensus-v1.cxx naive_utils.cxx
CS_EXEC = findConsensus

FBP_SOURCE = findBreakPointFromBLAT.cxx
FBP_EXEC = findBreakPointFromBLAT

CTRI_SOURCE = callTrErrorsWithIndelsFromCircSeq.cxx naive_utils.cxx
CTRI_EXEC = callTrErrorsWithIndelsFromCircSeq

RFB_SOURCE = refineBreakpoint.cxx naive_utils.cxx
RFB_EXEC = refineBreakpoint

RFBI_SOURCE = refineBreakpoint-withIndels.cxx naive_utils.cxx
RFBI_EXEC = refineBreakpoint-withIndels

DMT_EXEC = demulti
DMT_SOURCE = demulti.cxx

MOM_SOURCE = mergeObservationMaps.cxx
MOM_EXEC = mergeObservationMaps

MGM_SOURCE = mergeMapsGenome.cxx
MGM_EXEC = mergeMapsGenome


ABO_SOURCE = addBaseToObservations.cxx
ABO_EXEC = addBaseToObservations

ABOS_SOURCE = addBaseAndStrandToObservations.cxx
ABOS_EXEC = addBaseAndStrandToObservations



CXX = g++
CFLAGS = -g



map : $(TEST_MAP_SOURCE)
	$(CXX) $(CFLAGS) $(TEST_MAP_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapGenome -lcsc -lcircMapCovSimple -lcircMapIndels -lutils_fasta -lsamRead -lGFF -o $(BASE_EXEC_PATH)/$(TEST_MAP_EXEC)

cs : $(CS_SOURCE)
	$(CXX) $(CFLAGS) $(CS_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lutils_fasta -lcsc -lsamRead -lGFF -o $(BASE_EXEC_PATH)/$(CS_EXEC)

fbp : $(FBP_SOURCE)
	$(CXX) $(CFLAGS) $(FBP_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lutils_fasta -lparse_blast -o $(BASE_EXEC_PATH)/$(FBP_EXEC)

ctri : $(CTRI_SOURCE)
	$(CXX) $(CFLAGS) $(CTRI_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapGenome -lcircMapCovSimple -lcircMapIndels -lutils_fasta -lparse_blast -lsamRead -lcsc -lGFF -lcircMapObs -o $(BASE_EXEC_PATH)/$(CTRI_EXEC)

astmo : $(ASTMO_SOURCE)
	$(CXX) $(CFLAGS) $(ASTMO_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapGenome -lcircMapCovSimple -lcircMapIndels -lutils_fasta -lparse_blast -lsamRead -lcsc -lGFF -lcircMapObs -o $(BASE_EXEC_PATH)/$(ASTMO_EXEC)

rfb : $(RFB_SOURCE)
	$(CXX) $(CFLAGS) $(RFB_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapCovSimple -lutils_fasta -lparse_blast -lsamRead -lcsc -lGFF -o $(BASE_EXEC_PATH)/$(RFB_EXEC)

rfbi : $(RFBI_SOURCE)
	$(CXX) $(CFLAGS) $(RFBI_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapCovSimple -lutils_fasta -lparse_blast -lsamRead -lcsc -lGFF -o $(BASE_EXEC_PATH)/$(RFBI_EXEC)

dmt : $(DMT_SOURCE)
	$(CXX) $(CFLAGS) $(DMT_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lutils_fasta -o $(BASE_EXEC_PATH)/$(DMT_EXEC)


# Extras:
mom : $(MOM_SOURCE)
	$(CXX) $(CFLAGS) $(MOM_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapGenome -lGFF -lutils_fasta -lcircMapObs -lcircMapIndels -lsamRead -lcircMapCovSimple -o $(BASE_EXEC_EXTRAS_PATH)/$(MOM_EXEC)

mgm : $(MGM_SOURCE)
	$(CXX) $(CFLAGS) $(MGM_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapGenome -lGFF -lutils_fasta -lcircMapObs -lcircMapIndels -lsamRead -lcircMapCovSimple -o $(BASE_EXEC_EXTRAS_PATH)/$(MGM_EXEC)

abo : $(ABO_SOURCE)
	$(CXX) $(CFLAGS) $(ABO_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapObs -lGFF -lutils_fasta -lcircMapGenome -lcircMapIndels -lsamRead -lcircMapCovSimple -o $(BASE_EXEC_EXTRAS_PATH)/$(ABO_EXEC)

abos : $(ABOS_SOURCE)
	$(CXX) $(CFLAGS) $(ABOS_SOURCE) -I$(BASE_DIR)/src/libs -L$(BASE_DIR)/libs -lcircMapObs -lGFF -lutils_fasta -lcircMapGenome -lcircMapIndels -lsamRead -lcircMapCovSimple -o $(BASE_EXEC_EXTRAS_PATH)/$(ABOS_EXEC)

