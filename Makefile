CC=g++
DEPS=gumbel_max.o inverse_sampling.o misc.o
LIBS=-lgsl 
experiments: experiments.cc $(DEPS)
	$(CC) -o experiments $(DEPS) experiments.cc $(LIBS)

gumbel_max.o: gumbel_max.cc gumbel_max.h
	$(CC) $(LIBS) -c -o gumbel_max.o gumbel_max.cc 

inverse_sampling.o: inverse_sampling.cc inverse_sampling.h
	$(CC) $(LIBS) -c -o inverse_sampling.o inverse_sampling.cc 

misc.o: misc.cc misc.h
	$(CC) $(LIBS) -c -o misc.o misc.cc 

clean:
	rm -rf *.o *~ experiments
