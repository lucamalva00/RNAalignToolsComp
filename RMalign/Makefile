## set variable

DD = g++ -fPIE -g -O3 -Wall -Wno-sign-compare -fomit-frame-pointer -ffast-math 

OBJECT =  RNAalign.o Kabsch.o function.o

Main:$(OBJECT) main.o TMscore.o needle.o SARATMscore.o
	$(DD) $(OBJECT) main.o -o RNAalign
	$(DD) $(OBJECT) TMscore.o needle.o -o TMscore
	$(DD) $(OBJECT) SARATMscore.o needle.o -o SARATMscore
main.o:main.cpp
	$(DD) -c main.cpp 
RNAalign.o:RNAalign.cpp RNAalign.h global_variable.h
	$(DD) -c RNAalign.cpp
Kabsch.o:Kabsch.cpp Kabsch.h
	$(DD) -c Kabsch.cpp
function.o:function.cpp function.h
	$(DD) -c function.cpp
TMscore.o: TMscore.cpp needle.h
	$(DD) -c TMscore.cpp
needle.o: needle.cpp needle.h
	$(DD) -c needle.cpp
SARATMscore.o:SARATMscore.cpp SARATMscore.h needle.h
	$(DD) -c SARATMscore.cpp
clean:
	rm *.o
