all: clean build_reference build_alignr bin_file

# Finds a `bin` folder inside `samar-lite` and deletes it if it exists
clean:
	@find . -type d -name 'bin' -exec rm -rf {} +

# Builds the reference
build_reference:
	@echo +++ BUILDING REFERENCE...
	@cd reference && cargo build --release && echo +++ DONE BUILDING REFERENCE
	@echo +++

# Builds the alignr 
build_alignr:
	@echo +++ BUILDING ALIGNR...
	@cd alignr && cargo build --release && echo +++ DONE BUILDING ALIGNR
	@echo +++

# Creates the `bin` directory and places the copies of the executable in it
bin_file:
	@echo +++ CREATING bin DIRECTORY...
	@mkdir bin
	@echo +++ COPYING REFERENCE EXECUTBLE TO bin DIRECTORY
	@cp ./reference/target/release/ref-align ./bin
	@echo +++ COPYING ALIGNR EXECUTBLE TO bin DIRECTORY
	@cp ./alignr/target/release/alignr ./bin
	@echo +++ DONE INITIALIZING bin DIRECTORY
	@echo +++
