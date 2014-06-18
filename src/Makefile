all: snap2web.x

#H5_LIB=/opt/local/lib
#H5_INC=/opt/local/include

snap2lambda: snap2lambda.x

snap2dens: snap2dens.x

snap2fill: snap2fill.x

snap2lambda.x: snap2lambda_histo.c io_arepo.c io_arepo.h  array.h array.c struct.h
	cc -I$(H5_INC) -L$(H5_LIB) -lhdf5 -lm -Wall snap2lambda_histo.c array.c io_arepo.c -o snap2lambda.x  

snap2fill.x: snap2filling.c io_arepo.c io_arepo.h  array.h array.c header.h
	cc -I$(H5_INC) -L$(H5_LIB) -lhdf5 -lm -Wall snap2filling.c array.c io_arepo.c -o snap2fill.x  

snap2dens.x: snap2density_histo.c io_arepo.c io_arepo.h  array.h array.c
	cc -I$(H5_INC) -L$(H5_LIB) -lhdf5 -lm -Wall snap2density_histo.c array.c io_arepo.c -o snap2dens.x  

snap2web.x: snap2web.c io_arepo.c io_arepo.h
	cc -I$(H5_INC) -L$(H5_LIB) -lhdf5 -lm -Wall snap2web.c io_arepo.c -o snap2web.x


clean:
	rm -f *.o
	rm snap2web.x
