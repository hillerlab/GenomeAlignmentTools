include inc/localEnvironment.mk
include inc/common.mk

all: libs userApps

topLibs:
	@./checkUmask.sh
	@MACHTYPE=${MACHTYPE} ./machTest.sh
	cd lib && ${MAKE}
	cd jkOwnLib && ${MAKE}
	cd htslib && ${MAKE}

optLib:
	cd optimalLeaf && ${MAKE}

hgLib:
	cd hg/lib && ${MAKE}

libs: topLibs hgLib optLib

userApps: libs
	cd hg/mouseStuff/ && ${MAKE}

clean:
	cd lib && ${MAKE} clean 
	cd jkOwnLib && ${MAKE} clean
	cd htslib && ${MAKE} clean
	cd hg/lib/ && ${MAKE} clean
	cd optimalLeaf/ && ${MAKE} clean
	cd hg/mouseStuff/ && ${MAKE} clean
