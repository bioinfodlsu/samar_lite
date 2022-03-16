all: build_reference build_alignr bin_file

build_reference:
	@echo +++ BUILDING REFERENCE...
	@cd reference && cargo build --release && echo +++ DONE BUILDING REFERENCE
	@echo +++

build_alignr:
	@echo +++ BUILDING ALIGNR...
	@cd alignr && cargo build --release && echo +++ DONE BUILDING ALIGNR
	@echo +++

bin_file:
	@echo +++ CREATING BIN FILE...
	@mkdir bin
	@echo +++ COPYING REFERENCE EXECUTBLE TO BIN FILE
	@cp ./reference/target/release/ref-align.exe ./bin
	@echo +++ COPYING ALIGNR EXECUTBLE TO BIN FILE
	@cp ./alignr/target/release/alignr.exe ./bin
	@echo +++ DONE INITIALIZING BIN FILE
