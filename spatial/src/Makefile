LIBNAME=spatial
LD=ld

OBJS=	krc.o \
	pps.o 

spatial:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean:
	@rm -f Makefile *.o *.so
