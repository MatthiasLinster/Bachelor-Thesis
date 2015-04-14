# Name der binären Datei für den Wirkungsquerschnitt
CSBIN = crossSection
# Name der binären Datei für den Zerfall
DECAYBIN = decay

# Name der Library
LIBNAME = libbachelorarbeit

# Compiler & Flags
CC = g++
CFLAGS = -Wall -O3 -std=c++11
LDFLAGS = -lgsl -lgslcblas -lm
LHAPDFFLAGS = `lhapdf-config --cppflags -cxxflags --cflags --ldflags --libs`
ROOTFLAGS = `root-config --cflags --libs`

# Ausgabe-Ordner für die Objektdateien und Library
OBJDIR = obj
LIBDIR = lib

# Objekte
OBJECTS = $(addprefix $(OBJDIR)/, MonteCarlo.o MonteCarloFunction.o ScatteringMonteCarloFunction.o PhaseSpaceMonteCarloFunction.o TwoTwoOneFinalMassiveHadronScattering.o Scattering.o PhaseSpace.o TwoBodyOneFinalMassivePhaseSpace.o ThreeBodyMasslessPhaseSpace.o Decay.o ThreeBodyMasslessDecay.o DecayMonteCarloFunction.o TwoBodyMassivePhaseSpace.o ThreeBodyMassivePhaseSpace.o ThreeBodyMassiveDecay.o)

# Phony
.PHONY: clean library

# Wirkungsquerschnitt & Zerfall erzeugen
all: crosssection decay

# Erzeugt die Ausgabeordner für die Objekte und die Bibliothek
directories:
	mkdir -p $(OBJDIR)
	mkdir -p $(LIBDIR)
	
# Erzeugt das Library-Objekt
library: directories $(OBJECTS)
	ar crf $(LIBDIR)/$(LIBNAME).a $(OBJECTS)
	ranlib $(LIBDIR)/$(LIBNAME).a
	cp -p *.h $(LIBDIR)

# Programm für die Berechnung des Wirkungsquerschnitts
crosssection: $(OBJDIR)/crossSection.o library
	$(CC) $(CFLAGS) $(LHAPDFFLAGS) $(ROOTFLAGS) $(LDFLAGS) $(OBJDIR)/crossSection.o $(LIBDIR)/$(LIBNAME).a -o $(CSBIN)
	
# Programm für die Berechnung der Zerfallsbreite
decay: $(OBJDIR)/decay.o library
	$(CC) $(CFLAGS) $(LHAPDFFLAGS) $(ROOTFLAGS) $(LDFLAGS) $(OBJDIR)/decay.o $(LIBDIR)/$(LIBNAME).a -o $(DECAYBIN)
	
# Einzelne Objektdateien
$(OBJDIR)/crossSection.o: crossSection.cpp library
	$(CC) $(CFLAGS) $(LHAPDFFLAGS) $(ROOTFLAGS) -I$(LIBDIR)/ -o $@ -c $<
	
$(OBJDIR)/decay.o: decay.cpp library
	$(CC) $(CFLAGS) $(ROOTFLAGS) -I$(LIBDIR)/ -o $@ -c $<

$(OBJDIR)/MonteCarlo.o: MonteCarlo.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<

$(OBJDIR)/MonteCarloFunction.o: MonteCarloFunction.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/ScatteringMonteCarloFunction.o: ScatteringMonteCarloFunction.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<
	
$(OBJDIR)/PhaseSpaceMonteCarloFunction.o: PhaseSpaceMonteCarloFunction.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/TwoTwoOneFinalMassiveHadronScattering.o: TwoTwoOneFinalMassiveHadronScattering.cpp
	$(CC) $(CFLAGS) $(LHAPDFFLAGS) $(ROOTFLAGS) -o $@ -c $<

$(OBJDIR)/Scattering.o: Scattering.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<
	
$(OBJDIR)/PhaseSpace.o: PhaseSpace.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/TwoBodyOneFinalMassivePhaseSpace.o: TwoBodyOneFinalMassivePhaseSpace.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/ThreeBodyMasslessPhaseSpace.o: ThreeBodyMasslessPhaseSpace.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/Decay.o: Decay.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<
	
$(OBJDIR)/ThreeBodyMasslessDecay.o: ThreeBodyMasslessDecay.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<
	
$(OBJDIR)/DecayMonteCarloFunction.o: DecayMonteCarloFunction.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<
	
$(OBJDIR)/ThreeBodyMassivePhaseSpace.o: ThreeBodyMassivePhaseSpace.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/TwoBodyMassivePhaseSpace.o: TwoBodyMassivePhaseSpace.cpp
	$(CC) $(CFLAGS) -o $@ -c $<
	
$(OBJDIR)/ThreeBodyMassiveDecay.o: ThreeBodyMassiveDecay.cpp
	$(CC) $(CFLAGS) $(ROOTFLAGS) -o $@ -c $<

# Aufräumen
clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(CSBIN) $(DECAYBIN)

