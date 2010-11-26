CC=gcc
MACHTYPE=x86_64
CFLAGS=-c
COPTS=-Ikent/inc -Iinc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}
LDFLAGS=-lm
SOURCES=getErrorProfile.c
LIBSRC=kent/lib/common.c kent/lib/dlist.c kent/lib/bPlusTree.c kent/lib/dystring.c \
	 kent/lib/linefile.c kent/lib/nib.c kent/lib/sqlNum.c kent/lib/base64.c kent/lib/dnaLoad.c \
	 kent/lib/errabort.c kent/lib/hash.c kent/lib/localmem.c kent/lib/obscure.c \
         kent/lib/twoBit.c kent/lib/bits.c kent/lib/dnaseq.c kent/lib/fa.c kent/lib/https.c  kent/lib/mime.c \
         kent/lib/options.c kent/lib/udc.c kent/lib/cheapcgi.c kent/lib/dnautil.c lib/fastq.c \
	 kent/lib/internet.c kent/lib/net.c kent/lib/pipeline.c kent/lib/verbose.c kent/lib/memalloc.c \
	 kent/lib/osunix.c kent/lib/portimpl.c kent/lib/htmshell.c kent/lib/filePath.c \
	 kent/lib/log.c kent/lib/servpws.c kent/lib/servmsII.c kent/lib/servCrunx.c \
	 kent/lib/servBrcMcw.c kent/lib/wildcmp.c kent/lib/servcis.c kent/lib/intExp.c\
	 kent/lib/servcl.c kent/lib/kxTok.c


OBJECTS=$(SOURCES:.c=.o)
LIBOBJECTS=$(LIBSRC:.c=.o)
EXECUTABLE=getErrorProfile

all: $(SOURCES) $(LIBOBJECTS) $(EXECUTABLE) $(LIBA)

$(EXECUTABLE): $(OBJECTS) $(LIBOBJECTS) $(LIBA)
	$(CC) $(COPTS) $(LDFLAGS) $(OBJECTS) $(LIBOBJECTS) -o $@

install: all
	-mv $(EXECUTABLE) ../bin

#$(LIBA): $(LIBOBJECTS)
#	ar rcus $(LIBA) $(LIBOBJECTS)

.c.o:
	$(CC) $(COPTS) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS) $(LIBOBJECTS) $(EXECUTABLE)

check-syntax:
	$(CC) ${CFLAGS} -o .nul -S ${CHK_SOURCES}


