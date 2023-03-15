CC = /opt/compiler/gcc/11.2.0/bin/gcc
CFLAGS = -O3 -flto -march=native
LDFLAGS = -lm
OBJS = linsys.o
TARGETS = run1 run2 run3 run4 run5

all: $(TARGETS)

clean:
	rm -f *.o $(TARGETS)

run1: main1.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< ${OBJS} $(LDFLAGS)

run2: main2.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< ${OBJS} $(LDFLAGS)

run3: main3.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< ${OBJS} $(LDFLAGS)	

run4: main4.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< ${OBJS} $(LDFLAGS)

run5: main5.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $< ${OBJS} $(LDFLAGS)
