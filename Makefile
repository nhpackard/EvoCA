# One shared library for both models (EvoCA + ChemCA) and their shared core.
# The same command, for the shell:
#   gcc -O2 -Wall -fPIC -dynamiclib -o C/libevoca.dylib C/evoca.c C/core.c C/chemca.c
# (python/evoca_build.py states it for Python; keep the three in step.)
#
# zsh note: on N's machine `make` is shadowed by an autoloaded zsh function
# ("make: function definition file not found"), which is a shell setup issue,
# not this Makefile. `command make` or `/usr/bin/make` runs the real make.

CC     = gcc
CFLAGS = -O2 -Wall -fPIC
SRCS   = C/evoca.c C/core.c
HDRS   = C/evoca.h C/core.h

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    SHARED_FLAG = -dynamiclib
    LIB         = C/libevoca.dylib
else
    SHARED_FLAG = -shared
    LIB         = C/libevoca.so
endif

.PHONY: all clean

all: $(LIB)

$(LIB): $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(SHARED_FLAG) -o $@ $(SRCS)

clean:
	rm -f C/libevoca.so C/libevoca.dylib
