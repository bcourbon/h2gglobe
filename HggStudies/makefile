# définition des cibles particuliès
.PHONY: clean, mrproper
 
# désactivation des règles implicitesit
.SUFFIXES:
 
# définition des variablesb
CC = g++ 
#CFLAGS = -W -Wall -v -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags`
#CFLAGS = -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags` 
CFLAGS = `root-config --cflags --libs` -L $ROOTSYS/lib -lRooFit -lHtml -lMinuit -lRooFitCore -W

SEDIR = $(shell pwd)
BOOST_PATH=/Applications/boost_1_49_0

USERINCLUDES += -I$(ROOTSYS)/include
USERINCLUDES += -I$(BOOST_PATH)/include/
USERINCLUDES += -I $(CMSSW_BASE)/src/ -I $(CMSSW_RELEASE_BASE)/src/
USERINCLUDES += -I$(ROOFITSYS)/include

INCLUDE +=  -I $(BASEDIR)/interface
INCLUDE += $(USERINCLUDES)

CXXFLAGS = -Wall -g -Wl,-rpath=$(BASEDIR)/lib
CXXFLAGS += $(INCLUDE)

ROOT_DICT=RootDict.cxx

USERLIBS += -L$(BOOST_PATH)/libs/ -lboost_regex -lboost_program_options
USERLIBS += -L$(BOOST_PATH)/lib/ -lboost_regex -lboost_program_options
USERLIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore
USERLIBS += -L$(CMSSW_BASE)/lib
USERLIBS += $(shell root-config --glibs) -lTMVA -lMinuit -lFoam -lRooFit -lRooFitCore -lMathMore
LIBS += $(USERLIBS)
##LIBS += -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lHiggsAnalysisCombinedLimit
LIBS += -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/

LIBNAME = BiasStudy

LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
OBJDIR = $(BASEDIR)/obj





# all
all: TruthChoice.o ZeeFunctionsChoice.o BiasStudy.o BiasStudy_new.o BiasPlotter.o BiasStudy_muSpe.o RooExponentialSum.o RooPowerLawSum.o functions.o makeBkgPlotsGeneric_2.o makeVarPlots.o makeVarPlotsHistos.o massPlot.o BiasPlotter_2.o makeLatexSlides.o makeEffAcc.o ZeePlots.o TestFunctions.o
	$(CC) TruthChoice.o functions.o RooExponentialSum.o RooPowerLawSum.o -o TruthChoice.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) TestFunctions.o functions.o RooExponentialSum.o RooPowerLawSum.o -o TestFunctions.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) makeEffAcc.o functions.o RooExponentialSum.o RooPowerLawSum.o -o makeEffAcc.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)	
	$(CC) makeLatexSlides.o functions.o RooExponentialSum.o RooPowerLawSum.o -o makeLatexSlides.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR) 
	$(CC) ZeeFunctionsChoice.o functions.o RooExponentialSum.o RooPowerLawSum.o -o ZeeFunctionsChoice.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR) 
	$(CC) BiasStudy.o functions.o RooExponentialSum.o RooPowerLawSum.o -o BiasStudy.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR) 
	$(CC) BiasStudy_new.o functions.o RooExponentialSum.o RooPowerLawSum.o -o BiasStudy_new.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) BiasPlotter.o functions.o RooExponentialSum.o RooPowerLawSum.o -o BiasPlotter.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR) 
	$(CC)  BiasStudy_muSpe.o functions.o RooExponentialSum.o RooPowerLawSum.o -o BiasStudy_muSpe.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) BiasPlotter_2.o functions.o RooExponentialSum.o RooPowerLawSum.o -o BiasPlotter_2.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) makeBkgPlotsGeneric_2.o -o makeBkgPlotsGeneric_2.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR) 
	$(CC) makeVarPlots.o -o makeVarPlots.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) makeVarPlotsHistos.o -o makeVarPlotsHistos.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) massPlot.o functions.o RooExponentialSum.o RooPowerLawSum.o -o massPlot.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)
	$(CC) ZeePlots.o functions.o RooExponentialSum.o RooPowerLawSum.o -o ZeePlots.exe $(CXXFLAGS) $(ROOTCLASSES) $(LIBS) -L$(LIBDIR)

TruthChoice.o: TruthChoice.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c TruthChoice.cpp -o TruthChoice.o $(CXXFLAGS)

TestFunctions.o: TestFunctions.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c TestFunctions.cpp -o TestFunctions.o $(CXXFLAGS)

makeEffAcc.o: makeEffAcc.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c makeEffAcc.cpp -o makeEffAcc.o $(CXXFLAGS)

makeLatexSlides.o: makeLatexSlides.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c makeLatexSlides.cpp -o makeLatexSlides.o $(CXXFLAGS)

ZeeFunctionsChoice.o: ZeeFunctionsChoice.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c ZeeFunctionsChoice.cpp -o ZeeFunctionsChoice.o $(CXXFLAGS)

BiasStudy.o: BiasStudy.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c BiasStudy.cpp -o BiasStudy.o $(CXXFLAGS)

BiasStudy_new.o: BiasStudy_new.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c BiasStudy_new.cpp -o BiasStudy_new.o $(CXXFLAGS)

BiasPlotter.o: BiasPlotter.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c BiasPlotter.cpp -o BiasPlotter.o $(CXXFLAGS)

BiasStudy_muSpe.o: BiasStudy_muSpe.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c BiasStudy_muSpe.cpp -o BiasStudy_muSpe.o $(CXXFLAGS)

BiasPlotter_2.o: BiasPlotter_2.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c BiasPlotter_2.cpp -o BiasPlotter_2.o $(CXXFLAGS)

RooExponentialSum.o: RooExponentialSum.cc
	$(CC) -c RooExponentialSum.cc -o RooExponentialSum.o $(CXXFLAGS)

functions.o: functions.cpp RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c functions.cpp -o functions.o $(CXXFLAGS)

RooPowerLawSum.o: RooPowerLawSum.cc
	$(CC) -c RooPowerLawSum.cc -o RooPowerLawSum.o $(CXXFLAGS)

makeBkgPlotsGeneric_2.o: makeBkgPlotsGeneric_2.C
	$(CC) -c makeBkgPlotsGeneric_2.C -o makeBkgPlotsGeneric_2.o $(CXXFLAGS)

makeVarPlots.o: makeVarPlots.cpp
	$(CC) -c makeVarPlots.cpp -o makeVarPlots.o $(CXXFLAGS)

makeVarPlotsHistos.o: makeVarPlotsHistos.cpp
	$(CC) -c makeVarPlotsHistos.cpp -o makeVarPlotsHistos.o $(CXXFLAGS)

massPlot.o: massPlot.C functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c massPlot.C -o massPlot.o $(CXXFLAGS)


ZeePlots.o: ZeePlots.cpp functions.h RooExponentialSum.h RooPowerLawSum.h
	$(CC) -c ZeePlots.cpp -o ZeePlots.o $(CXXFLAGS)

# clean
clean:
	rm -rf *.bak rm -rf *.o

# mrproper
mrproper: clean
	rm -rf Programme

