LIBNAME=MASS
LD=ld

OBJS=	MASS.o

$(LIBNAME).so:  $(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so
