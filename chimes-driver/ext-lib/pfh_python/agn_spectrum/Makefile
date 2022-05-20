#OPT    += -DDIRECT_RAY_INTEGRATION

EXEC   = agn_spectrum_py.so

SRCS   = agn_spectrum.c


OBJS   = $(SRCS:.c=.o)
INCL   = 


#CFLAGS =  -shared -m32 $(OPT)
#CFLAGS =  -bundle $(OPT)
CFLAGS =  -shared -fPIC -O3 -g  #  -Wall
#SDK    = /Developer/SDKs/MacOSX10.4u.sdk
#CFLAGS  =  -bundle $(OPT) -isysroot $(SDK) -arch ppc -arch ppc64 -arch i386
#LDFLAGS = -isysroot $(SDK) -syslibroot,$(SDK)
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   = -lm

CC     =  gcc 

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)
#	$(CC) $(CFLAGS) $(OBJS) $(LIBS)  -o $(EXEC)

#$(OBJS): $(SRCS) 
$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)










