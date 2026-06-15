# Makefile for compiling sim_home_test.c
SRC := sim_test_bckp.c
OUT := sim_test_bckp

CC := gcc
CFLAGS := -fdiagnostics-color=always -g

all: $(OUT)

$(OUT): $(SRC)
	@echo "Compiling $(SRC) -> $(OUT)"
	$(CC) $(CFLAGS) "$(SRC)" -o "$(OUT)"

clean:
	@echo "Cleaning up..."
	rm -f "$(OUT)"
    find . -type f -executable -exec rm {} \;
