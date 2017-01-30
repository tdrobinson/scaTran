## set flags
FFLAGS = -fpp -c -g -O3 $(FMEMMODEL)
LFLAGS = -g -O3 $(FMEMMODEL)

## list of .o files
objs = scatran.o index_refract.o limbdarken.o path_geom.o path_mc.o \
path_rayt.o raytrace.o set_zenith.o transit_init.o transit_interp.o \
transit_optd.o transit_rayt.o transit_scat.o xyinterp.o

scatran: $(objs)
	ifort $(LFLAGS) -o $@ $(objs)

scatran.o: scatran.f param.inc
	ifort $(FFLAGS) scatran.f

index_refract.o: index_refract.f param.inc
	ifort $(FFLAGS) index_refract.f

limbdarken.o: limbdarken.f param.inc
	ifort $(FFLAGS) limbdarken.f

path_geom.o: path_geom.f param.inc
	ifort $(FFLAGS) path_geom.f

path_mc.o: path_mc.f param.inc
	ifort $(FFLAGS) path_mc.f

raytrace.o: raytrace.f param.inc
	ifort $(FFLAGS) raytrace.f

path_rayt.o: path_rayt.f param.inc
	ifort $(FFLAGS) path_rayt.f

set_zenith.o: set_zenith.f param.inc
	ifort $(FFLAGS) set_zenith.f

transit_init.o: transit_init.f param.inc
	ifort $(FFLAGS) transit_init.f

transit_interp.o: transit_interp.f param.inc
	ifort $(FFLAGS) transit_interp.f

transit_optd.o: transit_optd.f param.inc
	ifort $(FFLAGS) transit_optd.f

transit_rayt.o: transit_rayt.f param.inc
	ifort $(FFLAGS) transit_rayt.f

transit_spec.o: transit_spec.f param.inc
	ifort $(FFLAGS) transit_spec.f

transit_scat.o: transit_scat.f param.inc
	ifort $(FFLAGS) transit_scat.f

xyinterp.o: xyinterp.f param.inc
	ifort $(FFLAGS) xyinterp.f


## rules
default:	
	make scatran

install:
	@echo "No installation rules."

clean:
	rm *.o; rm *.i; rm scatran