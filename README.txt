Modelo cinético del viento solar (Grupo F)

Instrucciones:

En el cluster se añaden cuatro códigos principales, los cuales corresponden a la ley barométrica, el flujo de Hartmann, y dos códigos del viento solar que calcula el perfil y la distribución de velocidades.

Para hacer un ejecutable de todos los programas haga en la consola
make

Para hacer un ejecutable de un código en específico haga:
make Baromatric_Law.x                   //Si usted desea hacer un ejecutable de la ley barométrica
make Hartmann.x                         //Si usted desea hacer un ejecutable del flujo de Hartmann
make SolarWind_Distribution.x           //Si usted desea hacer un ejecutable de la distribución de velocidades
make SolarWind_Profile.x                //Si usted desea hacer un ejecutable del perfil de velocidades

Es recomendable que haga todos los ejecutables a la vez (escribiendo make) y si usted da en la consola
ls
podrá encontrar 4 archivos que terminan en .x, eso significa que son ejecutables (usualmente en la consola se ponen en color verde).

Para ejecutar los programas haga en la consola:
./nombre_del_ejecutable.x
