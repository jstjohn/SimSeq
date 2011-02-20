CC=gcc
MACHTYPE=x86_64
CFLAGS=-c
COPTS=-I../ThirdParty/kent/inc -I../inc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}
SAMINC=-I../ThirdParty/samtools
LDFLAGS=-lm
SAMLD=-L../ThirdParty/samtools -lbam
FASIZE=faSize.c
SAMTOFQ=samToFastq.c
ERRORPROFILE=getErrorProfile.c
LIBSRC=../ThirdParty/kent/lib/common.c ../ThirdParty/kent/lib/dlist.c ../ThirdParty/kent/lib/bPlusTree.c ../ThirdParty/kent/lib/dystring.c \
	 ../ThirdParty/kent/lib/linefile.c ../ThirdParty/kent/lib/nib.c ../ThirdParty/kent/lib/sqlNum.c ../ThirdParty/kent/lib/base64.c ../ThirdParty/kent/lib/dnaLoad.c \
	 ../ThirdParty/kent/lib/errabort.c ../ThirdParty/kent/lib/hash.c ../ThirdParty/kent/lib/localmem.c ../ThirdParty/kent/lib/obscure.c \
         ../ThirdParty/kent/lib/twoBit.c ../ThirdParty/kent/lib/bits.c ../ThirdParty/kent/lib/dnaseq.c ../ThirdParty/kent/lib/fa.c ../ThirdParty/kent/lib/https.c  ../ThirdParty/kent/lib/mime.c \
         ../ThirdParty/kent/lib/options.c ../ThirdParty/kent/lib/udc.c ../ThirdParty/kent/lib/cheapcgi.c ../ThirdParty/kent/lib/dnautil.c lib/fastq.c \
	 ../ThirdParty/kent/lib/internet.c ../ThirdParty/kent/lib/net.c ../ThirdParty/kent/lib/pipeline.c ../ThirdParty/kent/lib/verbose.c ../ThirdParty/kent/lib/memalloc.c \
	 ../ThirdParty/kent/lib/osunix.c ../ThirdParty/kent/lib/portimpl.c ../ThirdParty/kent/lib/htmshell.c ../ThirdParty/kent/lib/filePath.c \
	 ../ThirdParty/kent/lib/log.c ../ThirdParty/kent/lib/servpws.c ../ThirdParty/kent/lib/servmsII.c ../ThirdParty/kent/lib/servCrunx.c \
	 ../ThirdParty/kent/lib/servBrcMcw.c ../ThirdParty/kent/lib/wildcmp.c ../ThirdParty/kent/lib/servcis.c ../ThirdParty/kent/lib/intExp.c\
	 ../ThirdParty/kent/lib/servcl.c ../ThirdParty/kent/lib/kxTok.c

FASIZEO=$(FASIZE:.c=.o)
SAMTOFQO=$(SAMTOFQ:.c=.o)
ERRORPROFILEO=$(ERRORPROFILE:.c=.o)
LIBOBJECTS=$(LIBSRC:.c=.o)

all: $(SOURCES) $(LIBOBJECTS) getErrorProfile faSize samToFastq

getErrorProfile: $(ERRORPROFILEO) $(LIBOBJECTS)
	$(CC) $(COPTS) $(LDFLAGS) $(ERRORPROFILEO) $(LIBOBJECTS) -o $@

faSize: $(LIBOBJECTS) $(FASIZEO)
	$(CC) $(COPTS) $(LDFLAGS) $(FASIZEO) $(LIBOBJECTS) -o $@

samToFastq: SAMLIB $(LIBOBJECTS) $(SAMTOFQO)
	$(CC) $(COPTS) $(SAMLD) $(LDFLAGS) $(SAMTOFQO) $(LIBOBJECTS) -o $@

SAMLIB: ../ThirdParty/samtools/libbam.a

../ThirdParty/samtools/libbam.a:
	-pushd ../ThirdParty/samtools
	-make
	-popd

../bin:
	-mkdir ../bin

install: all ../bin
	-cp getErrorProfile ../bin
	-cp samToFastq ../bin
	-cp faSize ../bin

.c.o:
	$(CC) $(COPTS) $(CFLAGS) $(SAMINC) $< -o $@

clean:
	-rm $(OBJECTS) $(LIBOBJECTS) $(EXECUTABLE)

check-syntax:
	$(CC) $(CFLAGS) $(SAMINC) -o .nul -S $(CHK_SOURCES)


