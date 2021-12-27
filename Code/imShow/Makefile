MPATH=$(realpath $(lastword $(MAKEFILE_LIST)))
CXX = g++
CUDAROOT  = /usr/local/cuda
CUDAINC   = $(CUDAROOT)/lib64

INCLUDES = -Iinclude -I/usr/local/cuda/include -I../shared -I../shared/include -I../shared/FastAC -I../shared/SplineGenerate -L$(CUDAINC) -L/usr/lib/x86_64-linux-gnu -L/usr/local/cuda-10.0/targets/x86_64-linux/lib 
CPPFLAGS = -Wall -lm -Wno-variadic-macros -std=c++11 -g -DGTE_USE_LINUX
CPPFLAGS += $(INCLUDES)
CUDAFLAGS = -lcudart -lcuda -lGL -lGLU -lGLEW -lglut

LFLAG=-L../shared/SplineGenerate/lib/Debug -lgtapplications -lgtgraphics -lgtmathematicsgpu -lX11 -lXext -lGL -lEGL -lpng -lpthread -lm 

NOTHING = 
SPACE = $(NOTHING) $(NOTHING)

CPPFLAGS += -pg -g -O3 -Wall -Wextra -pedantic 
#ifdef GPROF
#	DBG = 1
#	CPPFLAGS += -pg
#endif
#
#ifdef DBG
#	CPPFLAGS += -Wall -pedantic -g 
#else
#	CPPFLAGS += -O3
#endif


#ADD VERBOSE
ifdef VERBOSE
	CPPFLAGS += -DSHOW_DEBUG_OUTPUT=1
endif

OUTDIR = obj/
OUTFILE = show_skeleton

SOURCE = $(wildcard *.cpp)
OBJ = $(OUTDIR)$(subst $(SPACE), $(OUTDIR),$(SOURCE:.cpp=.o))
EXTERNAL = ../shared/CUDASkel2D/CUDASkel2D.a ../shared/FastAC/arithmetic_codec.o ../shared/lodepng/lodepng.o ../shared/SplineGenerate/BSplineCurveFitter/Debug/BSplineCurveFitterWindow3.o

ABSOUT=$(MPATH:Makefile=)$(OUTDIR)

default:all

all: $(OBJ) ../shared/FastAC/arithmetic_codec.o ../shared/SplineGenerate/BSplineCurveFitter/Debug/BSplineCurveFitterWindow3.o
	$(CXX) $(CPPFLAGS) -o $(OUTFILE) $(OBJ) $(EXTERNAL) $(CUDAFLAGS) $(LFLAG)
	
$(OUTDIR)%.o : %.cpp
	@mkdir -p $(OUTDIR)
	$(CXX) $(CPPFLAGS) -c $< -o $@

../shared/CUDASkel2D/CUDASkel2D.a:
	make -C ../shared/CUDASkel2D -f Makefile skeleton

../shared/lodepng/lodepng.o:
	make -C ../shared/lodepng -f Makefile

../shared/SplineGenerate/BSplineCurveFitter/Debug/BSplineCurveFitterWindow3.o: ../shared/SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.cpp ../shared/SplineGenerate/BSplineCurveFitter/BSplineCurveFitterWindow3.h
	make -C ../shared/SplineGenerate/BSplineCurveFitter CFG=Debug APP=BSplineCurveFitter -f ../makesample.gte

../shared/FastAC/arithmetic_codec.o:
	g++ -g ../shared/FastAC/arithmetic_codec.cpp -c -o ../shared/FastAC/arithmetic_codec.o

clean:
	$(RM) *.gch ../shared/SplineGenerate/BSplineCurveFitter/Debug/BSplineCurveFitterWindow3.o
	$(RM) ../shared/FastAC/arithmetic_codec.o
	$(RM) $(OUTDIR)/*
	$(RM) -r $(OUTDIR)
	$(RM) -rf lib
	$(RM) $(OUTFILE)

clean-full: clean
	make -C ../shared/CUDASkel2D -f Makefile clean
	make -C ../shared/lodepng -f Makefile clean
