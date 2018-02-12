#BASE_PATH=/N/u/jgout/Mason/pipeline-for-publication
BASE_PATH := ${CURDIR}

LIBSRC_PATH=$(BASE_PATH)/src/libs

LIBSRCDIRS := $(wildcard $(LIBSRC_PATH)/*/.)

#LIBSRCDIRS = $(LIBSRC_PATH)/CircSeqConsensus/.

PIPELINE_PATH=$(BASE_PATH)/src/pipeline

EXEC_LIB_PATH=$(BASE_PATH)/libs
EXEC_PATH=$(BASE_PATH)/bin
EXEC_EXTRAS_PATH=$(EXEC_PATH)/extras


libs: $(LIBSRCDIRS) force
	if [ ! -d $(EXEC_LIB_PATH) ]; then mkdir $(EXEC_LIB_PATH); fi
	for dir in $(LIBSRCDIRS) ; do \
		echo $$dir ; \
		cd $$dir; \
		make; \
	done

	mv $(LIBSRC_PATH)/*/*.a $(EXEC_LIB_PATH)
#	$(MAKE) $(LIBSRCDIRS)

.PHONY:	force

libsclean: $(LIBSRCDIRS)
	for dir in $(LIBSRCDIRS) ; do \
		echo $$dir ; \
		cd $$dir; \
		make clean; \
	done

findConsensus: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make cs

map: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make map

findBreakPoint: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make fbp

callErrors: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make ctri

refineBreakPoint: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make rfbi

demult: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make dmt

makeBinsDir:	force
	if [ ! -d $(EXEC_PATH) ]; then mkdir $(EXEC_PATH); fi

makeBinsExtraDir:	force
	if [ ! -d $(EXEC_EXTRAS_PATH) ]; then mkdir $(EXEC_EXTRAS_PATH); fi


pipeline: force makeBinsDir findConsensus map findBreakPoint callErrors refineBreakPoint demult


mom: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make mom

abo: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make abo

abos: $(PIPELINE_PATH)
	cd $(PIPELINE_PATH); \
	make abos


extras: force makeBinsExtraDir mom abo abos

all: libs pipeline extras

clean:
	if [ -d $(EXEC_PATH) ]; then find $(EXEC_PATH) -type f -delete; fi
	if [ -d $(EXEC_LIB_PATH) ]; then find $(EXEC_LIB_PATH) -type f -delete; fi
	find $(BASE_PATH) -type f -name "*.o" -delete
	find $(BASE_PATH) -type f -name "*.a" -delete
	find $(BASE_PATH) -type f -name "*~" -delete

