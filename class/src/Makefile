LIBNAME=class
LD=ld

OBJS=   class.o

$(LIBNAME).so: $(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so
