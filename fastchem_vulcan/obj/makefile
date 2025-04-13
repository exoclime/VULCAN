
# get options which hold for all files of the project
D=..
include $D/make.global_options


# create lists of file names
o-files     = $(shell echo *.o)

#o-files-hydro = $(addprefix ../$(OBJDIR), $(patsubst %,%,$(c-files:.c=.$(OBJEXT))))
#o-files      = $(sort $(o-files-hydro))
#h-files      = $(sort $(shell echo *.h))


# top level targets
all-l:
	@echo "linking program --> $(o-files)"
	@$(CC) $(o-files) $(F77LIBS) $(LDFLAGS) -o ../$(EXECUTABLE_NAME)
	
clean: 
	rm -f $(o-files) 
