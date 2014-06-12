all: snap2web.x

#H5_LIB=/opt/local/lib
#H5_INC=/opt/local/include

snap2web.x: snap2web.c
	cc -I$(H5_INC) -L$(H5_LIB) -lhdf5 -lm -Wall snap2web.c -o snap2web.x

clean:
	rm snap2web.x