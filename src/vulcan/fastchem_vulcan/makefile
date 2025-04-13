
.DEFAULT_GOAL := all

include make.global_options


model_main_: 
	cd $/model_main && $(MAKE) $(MAKEOPTIONS) model_main_
	
chemistry_: 
	cd $/fastchem_src && $(MAKE) $(MAKEOPTIONS) fastchem_
	
	
all-c: model_main_ chemistry_

	@echo "compiling done."

all-l:
	cd $/obj && $(MAKE) $(MAKEOPTIONS) all-l
	
	@echo "linking done."
	
all: all-c all-l

	@echo "everything is done and fine. enjoy your day!"

	
#and here we clean the mess
clean: clean-binary clean-model_main clean-chemistry

clean-binary: 
	rm -f $(EXECUTABLE_NAME)
	
clean-model_main: 
	cd $/model_main && $(MAKE) clean
	
clean-chemistry: 
	cd $/fastchem_src && $(MAKE) clean
	
	@echo "all clean. have a nice day!"
	
