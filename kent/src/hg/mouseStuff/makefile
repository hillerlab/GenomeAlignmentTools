# Build all directories in makeDb 

kentSrc = ../..

include ${kentSrc}/inc/localEnvironment.mk
include ${kentSrc}/inc/common.mk

# List of directories to build for userApps
UTILS_APPLIST = \
	axtChain \
	chainAntiRepeat \
	chainFilter \
	chainMergeSort \
	chainNet \
	chainPreNet \
	chainSort \
	netToAxt \
	axtToMaf \
	mafToAxt \
	axtToPsl \
	netFilter \
	netSyntenic \
	axtSort \
	axtSwap \
	chainBridge \
	chainSplit \
	chainStitchId \
	chainSwap \
	chainToAxt \
	chainToPsl \
	chainToPslBasic \
	lavToAxt \
	lavToPsl \
	mafRanges \
	mafNoAlign \
	mafToPsl \
	mafToFa \
	netSplit \
	netToBed \
	pslMrnaCover \
	chainCleaner \
	chainBridge 

# List of directories to build
USE_DIRS = $(UTILS_APPLIST)

all: ${USE_DIRS:%=%.all}
%.all:
	cd $* && echo $* && ${MAKE}

userApps: $(UTILS_APPLIST:%=%.userApp)
%.userApp:
	cd $* && echo $* && $(MAKE)

clean:: $(USE_DIRS:%=%.cleanAll)
%.cleanAll:
	cd $* && echo $* && $(MAKE) clean
