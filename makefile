# Compilation options
CXXFLAGS += -O3 

# The library is:
LIB=libSAT.a

# and is made from the following object files
LIBOBJS= \
  SAT_DE.o SAT_Filter.o SAT_Force.o SAT_Kepler.o \
  SAT_RefSys.o SAT_Time.o SAT_VecMat.o

# The application files themselves are:
EXERCISES= \
  Exercise_2_1.out Exercise_2_2.out Exercise_2_3.out Exercise_2_4.out \
  Exercise_2_5.out Exercise_2_6.out \
  Exercise_3_1.out Exercise_3_2.out Exercise_3_3.out Exercise_3_4.out \
  Exercise_4_1.out Exercise_4_2.out Exercise_4_3.out \
  Exercise_5_1.out Exercise_5_2.out Exercise_5_3.out \
  Exercise_6_1.out Exercise_6_2.out Exercise_6_3.out Exercise_6_4.out \
  Exercise_7_1.out \
  Exercise_8_1.out Exercise_8_2.out Exercise_8_3.out 

APPLICS= \
  GEODA.out RTOD.out TDRSOD.out

# Our main target is to make all the applications
all: $(APPLICS) $(EXERCISES)

clean:
	@-rm -rf *.o *.out *.a *.log *~ 

# The library is made from the object files by using 'ar'
libSAT.a: $(LIBOBJS)
	ar rc $(LIB) $(LIBOBJS)

# Each of the application files, in addition to being dependant on its
# own source file, needs the library
$(APPLICS): $(LIB)
	$(CXX) $(CXXFLAGS) $(@:.out=.cpp) -o $@ $(LIB)

$(EXERCISES): $(LIB)
	$(CXX) $(CXXFLAGS) $(@:.out=.cpp) -o $@ $(LIB)

test:all
	for exe in $(shell ls Exercise*.out);do ./$$exe;done
	./GEODA.out GEODA_A1.inp
	./RTOD.out RTOD_A.inp
	./TDRSOD.out 

commit:
	@git commit -am "Commit"
	@git push origin master

pull:
	@git reset --hard HEAD
