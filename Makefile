CC = mpicc
CFLAGS = -lm -Wall -march=native -Ofast
OBJS=utils.o
BINS=onesided_gol sendrecv_gol

all: $(BINS) 

%.o: %.c Makefile
	$(CC) $(CFLAGS) $< -c -o $@

onesided_gol: onesided_gol.c Makefile $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(OBJS)

sendrecv_gol: sendrecv_gol.c Makefile $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(OBJS)
clean:
	$(RM) $(BINS) $(OBJS)
