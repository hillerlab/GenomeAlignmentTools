SUBDIRS = chainCleaner chainNet scoreChain

.PHONY: all clean $(SUBDIRS)

all:
	@for dir in $(SUBDIRS); do \
		$(MAKE) -C $${dir}; \
	done


clean:
	@for dir in $(SUBDIRS); do \
		$(MAKE) clean -C $${dir}; \
	done
	@rm -f ../bin/chainCleaner ../bin/chainNet ../bin/scoreChain
