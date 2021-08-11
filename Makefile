all: SolarWind_Distribution.x SolarWind_Profile.x Barometric_Law.x Hartmann.x

SolarWind_Distribution.x: SolarWind_Distribution.cpp
	g++ $^ -o $@

SolarWind_Profile.x: SolarWind_Profile.cpp
	g++ $^ -o $@

Barometric_Law.x: Barometric_Law.cpp
	g++ $^ -o $@

Hartmann.x: Hartmann.cpp
	g++ -o Hartmann.x -std=c++11 Hartmann.cpp

clean:
	-rm -f *~ *# *o *out *x 
