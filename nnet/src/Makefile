LIBNAME=nnet
LD=ld

OBJS=	nnet.o

$(LIBNAME).so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean: clean
